#!/usr/bin/env python3
"""
Standardize metadata value spellings in an RNA-seq expression matrix CSV.

This script is a publication-friendly replacement for ad-hoc "hyoki_yure.py"-style cleaners.
It keeps the same basic I/O style used in your current scripts:

  python standardize_rnaseq_metadata.py input.csv output.csv

Key feature:
  - Value replacements for treatment/cultivar/stage/tissue can be provided from outside
    via a TSV mapping file or repeated --replace flags.

#Example of a TSV mapping file
# column	from	to
tissue	Leaves	leaf
tissue	leaves	leaf
tissue	Mature_leaves	leaf
stage	Seedling	seedling_stage
stage	seedling	seedling_stage
treatment	Control-no_stress	control
"""
import argparse
import collections
import re
import sys
from typing import Dict, List, Tuple

import pandas as pd


# -------------------------
# Utilities for reporting
# -------------------------
def norm_key(s: str) -> str:
    """For variant detection: casefold + remove separators/punctuations."""
    s = str(s).strip().casefold()
    s = re.sub(r"[\s_\-\/]+", "", s)
    s = re.sub(r"[^\w]", "", s)
    return s


def report_variants(df: pd.DataFrame, col: str) -> Dict[str, List[str]]:
    """Return groups of values that look equivalent after normalization."""
    vals = df[col].dropna().astype(str)
    groups: Dict[str, set] = collections.defaultdict(set)
    for v in vals:
        groups[norm_key(v)].add(v)
    return {k: sorted(vs) for k, vs in groups.items() if len(vs) > 1}


# -------------------------
# Normalizers (kept close to your current behavior)
# -------------------------
def to_lower_snake(x):
    """
    Convert spaces/slashes to underscores and collapse repeats.
    NOTE: Intentionally does NOT lowercase (to keep backward compatibility with your scripts).
    """
    if pd.isna(x):
        return x
    s = str(x).strip().replace("℃", "")
    s = re.sub(r"[\/\s]+", "_", s)
    s = re.sub(r"[_\-]+", "_", s)
    s = re.sub(r"__+", "_", s)
    return s.strip("_")


def to_lower_snake_then_lower(x):
    """For stage/tissue: snake_case then force lowercase (as in your scripts)."""
    x = to_lower_snake(x)
    if pd.isna(x):
        return x
    return str(x).lower()


def clean_temperature(x):
    """Extract the first numeric token (e.g., '25.0℃' -> '25.0')."""
    if pd.isna(x):
        return x
    m = re.search(r"(\d+(?:\.\d+)?)", str(x))
    return m.group(1) if m else str(x).strip()


# -------------------------
# Replacement mapping I/O
# -------------------------
def _parse_replace_spec(spec: str) -> Tuple[str, str, str]:
    """
    Parse --replace 'column:FROM=TO'
    Example: --replace 'tissue:Leaves=leaf'
    """
    if ":" not in spec or "=" not in spec:
        raise ValueError(f"Invalid --replace format: {spec!r}. Expected 'column:FROM=TO'.")
    col, rest = spec.split(":", 1)
    src, dst = rest.split("=", 1)
    col = col.strip()
    src = src.strip()
    dst = dst.strip()
    if not col:
        raise ValueError(f"Invalid --replace format (empty column): {spec!r}")
    return col, src, dst


def load_map_tsv(path: str) -> Dict[str, Dict[str, str]]:
    """
    TSV format (tab-separated):
        column <TAB> from <TAB> to

    - Lines starting with '#' are ignored.
    - Empty lines are ignored.
    """
    mapping: Dict[str, Dict[str, str]] = collections.defaultdict(dict)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for ln_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"{path}:{ln_no}: expected 3 columns (column, from, to)")
            col, src, dst = parts[0].strip(), parts[1].strip(), parts[2].strip()
            if not col:
                raise ValueError(f"{path}:{ln_no}: empty column name")
            mapping[col][src] = dst
    return dict(mapping)


def merge_mappings(maps: List[Dict[str, Dict[str, str]]]) -> Dict[str, Dict[str, str]]:
    merged: Dict[str, Dict[str, str]] = collections.defaultdict(dict)
    for mp in maps:
        for col, d in mp.items():
            merged[col].update(d)
    return dict(merged)


def apply_mapping(out: pd.DataFrame, mapping: Dict[str, Dict[str, str]]):
    """Apply exact string replacements per column (only if the column exists)."""
    for col, rep in mapping.items():
        if col in out.columns and rep:
            out[col] = out[col].replace(rep)
    return out


# -------------------------
# Column-specific cleanups
# -------------------------
def clean_original_10th_column_series(v):
    """
    Equivalent to the following awk logic applied to the *original CSV* 10th field ($10):

      gsub(/_/, " ", $10)
      gsub(/;[[:space:]]*/, "; ", $10)

    IMPORTANT:
      Because we read with index_col=0, original $10 becomes out.columns[8] inside pandas.
    """
    if pd.isna(v):
        return v
    if not isinstance(v, str):
        return v
    v = v.replace("_", " ")
    v = re.sub(r";\s*", "; ", v)
    return v


def main():
    p = argparse.ArgumentParser(
        description="Standardize metadata spellings in RNA-seq expression matrix CSV (表記ゆれ修正)."
    )
    p.add_argument("input_csv", help="Input CSV matrix (first column is treated as index; like your current scripts).")
    p.add_argument("output_csv", help="Output CSV path.")
    p.add_argument(
        "--map-tsv",
        action="append",
        default=[],
        help="TSV mapping file: column<TAB>from<TAB>to (repeatable).",
    )
    p.add_argument(
        "--replace",
        action="append",
        default=[],
        help="Inline replacement: column:FROM=TO (repeatable). Example: --replace 'tissue:Leaves=leaf'",
    )
    p.add_argument(
        "--meta-cols",
        type=int,
        default=9,
        help="How many columns from the left (after index_col=0) are metadata columns to fill NA for. Default: 9.",
    )
    p.add_argument(
        "--na-token",
        default="NA",
        help="Token to use for missing/blank values in metadata columns. Default: NA",
    )
    p.add_argument(
        "--no-report",
        action="store_true",
        help="Do not print variant report.",
    )
    p.add_argument(
        "--skip-cols4-7-underscore-to-space",
        action="store_true",
        help="Skip '_' -> ' ' conversion for output columns 4-7 (excluding index).",
    )
    p.add_argument(
        "--skip-col10-clean",
        action="store_true",
        help="Skip cleaning of the original CSV 10th column (underscore/semicolon spacing).",
    )
    args = p.parse_args()

    df = pd.read_csv(args.input_csv, index_col=0)

    # 1) Report possible variants (like your current scripts)
    if not args.no_report:
        for c in ["treatment", "tissue", "stage", "cultivar", "temperature"]:
            if c in df.columns:
                v = report_variants(df, c)
                if v:
                    print(f"\n[{c}] variants:")
                    for _, forms in v.items():
                        print(" ", forms)

    out = df.copy()

    # 2) Load user-provided replacement rules
    maps: List[Dict[str, Dict[str, str]]] = []
    for tsv in args.map_tsv:
        maps.append(load_map_tsv(tsv))

    inline_map: Dict[str, Dict[str, str]] = collections.defaultdict(dict)
    for spec in args.replace:
        col, src, dst = _parse_replace_spec(spec)
        inline_map[col][src] = dst
    maps.append(dict(inline_map))

    mapping = merge_mappings(maps)

    # 3) Apply replacements BEFORE normalization
    out = apply_mapping(out, mapping)

    # 4) Normalize (kept close to your current scripts)
    if "treatment" in out.columns:
        out["treatment"] = out["treatment"].map(to_lower_snake)

    if "stage" in out.columns:
        out["stage"] = out["stage"].map(to_lower_snake_then_lower)

    if "tissue" in out.columns:
        out["tissue"] = out["tissue"].map(to_lower_snake_then_lower)

    if "temperature" in out.columns:
        out["temperature"] = out["temperature"].map(clean_temperature)

    # 5) Apply replacements AGAIN (so rules can target normalized values too)
    out = apply_mapping(out, mapping)

    # 6) Fill blanks/NaN -> NA for metadata columns (leftmost N columns)
    meta_cols = list(out.columns[: max(0, args.meta_cols)])
    if meta_cols:
        out[meta_cols] = out[meta_cols].replace(r"^\s*$", pd.NA, regex=True).fillna(args.na_token)

    # 7) Convert '_' -> ' ' for columns 4-7 (excluding index) (as in your current scripts)
    if not args.skip_cols4_7_underscore_to_space and out.shape[1] >= 6:
        cols_4_7 = list(out.columns[2:6])  # 0-index [2,3,4,5] => 4th-7th columns in the original CSV
        out[cols_4_7] = (
            out[cols_4_7]
            .astype("string")
            .apply(lambda s: s.str.replace("_", " ", regex=False))
        )

    # 8) Clean original CSV 10th column (underscore removal + '; ' normalization)
    if not args.skip_col10_clean and out.shape[1] >= 9:
        col10 = out.columns[8]  # because index_col=0 shifts original $10 -> columns[8]
        out[col10] = out[col10].map(clean_original_10th_column_series)

    out.to_csv(args.output_csv)
    print(f"\nSaved: {args.output_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
