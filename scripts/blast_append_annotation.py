#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Append Arabidopsis gene annotation columns to a BLAST outfmt6-like table.

- Joins by AT* ID found in a specified BLAST column (1-based).
- Annotation TSV is assumed to have a key column (1-based) and additional columns.
- Appends ONLY the annotation columns excluding the key (so AT ID is NOT duplicated on the right).
- If exact match not found, tries "gene-level" fallback by stripping isoform suffix:
    AT1G01010.1 -> AT1G01010
  (also strips trailing .<digits> generally)

Usage:
  python blast_append_annotation.py \
    --blast Pg_Tift_proteins_At_proteins_blastp_out.txt \
    --anno  Arab_gene_annotation_list.tsv \
    --out   Pg_Tift_proteins_At_proteins_blastp_out.with_ArabAnno.tsv

Optional:
  --blast-id-col 5
  --anno-key-col 1
  --missing NA
  --blast-sep "\t"  (default: auto detect as tab or whitespace)
  --anno-sep  "\t"  (default: tab)
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def normalize_id(x: str) -> str:
    """Normalize ID for robust matching (strip spaces)."""
    return x.strip()


_iso_re = re.compile(r"\.\d+$")


def gene_level_id(x: str) -> str:
    """AT1G01010.1 -> AT1G01010 (strip trailing .digits)"""
    x = normalize_id(x)
    return _iso_re.sub("", x)


def detect_sep(sample_line: str, preferred: str = "\t") -> str:
    """Detect separator for blast file (tab vs whitespace)."""
    if preferred:
        # If user explicitly passes, use it.
        return preferred
    # Auto: if tabs exist, use tab; else split on any whitespace.
    return "\t" if "\t" in sample_line else None  # None => csv with split() later


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--blast", required=True, help="BLAST out file (TSV/whitespace delimited).")
    p.add_argument("--anno", required=True, help="Annotation TSV file.")
    p.add_argument("--out", required=True, help="Output TSV file.")
    p.add_argument("--blast-id-col", type=int, default=5, help="1-based column index of AT ID in BLAST file.")
    p.add_argument("--anno-key-col", type=int, default=1, help="1-based key column index in annotation TSV.")
    p.add_argument("--missing", default="NA", help="Value to fill when annotation not found.")
    p.add_argument("--blast-sep", default="", help=r"Separator for BLAST file. Use '\t' for tab. Empty = auto.")
    p.add_argument("--anno-sep", default="\t", help=r"Separator for annotation file (default: tab).")
    return p.parse_args()


def read_annotation(anno_path: Path, key_col_1based: int, sep: str):
    """
    Returns:
      anno_map_exact: dict[key -> list[anno_cols_without_key]]
      anno_map_gene : dict[gene_level_key -> list[anno_cols_without_key]] (only if key has isoform)
      n_anno_cols_without_key: int
    """
    key_idx = key_col_1based - 1
    anno_map_exact = {}
    anno_map_gene = {}

    with anno_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter=sep)
        first_row = None

        for row in reader:
            if not row:
                continue
            if first_row is None:
                first_row = row

            if key_idx >= len(row):
                # skip malformed row
                continue

            key = normalize_id(row[key_idx])
            # Build appended columns excluding key col
            appended = [row[i] for i in range(len(row)) if i != key_idx]

            # Store exact
            if key and key not in anno_map_exact:
                anno_map_exact[key] = appended

            # Store gene-level fallback too
            gkey = gene_level_id(key)
            if gkey and gkey not in anno_map_gene:
                anno_map_gene[gkey] = appended

    if first_row is None:
        raise ValueError(f"Annotation file is empty: {anno_path}")

    n_without_key = len(first_row) - 1
    return anno_map_exact, anno_map_gene, n_without_key


def main():
    args = parse_args()

    blast_path = Path(args.blast)
    anno_path = Path(args.anno)
    out_path = Path(args.out)

    if args.blast_id_col < 1:
        raise ValueError("--blast-id-col must be 1-based (>=1)")
    if args.anno_key_col < 1:
        raise ValueError("--anno-key-col must be 1-based (>=1)")

    # Load annotation maps
    anno_map_exact, anno_map_gene, n_anno_cols = read_annotation(
        anno_path, args.anno_key_col, args.anno_sep
    )

    # Prepare BLAST separator (auto if empty)
    blast_sep = args.blast_sep if args.blast_sep != "" else None

    # Peek first non-empty line to decide delimiter if auto
    with blast_path.open("r", encoding="utf-8", errors="replace") as fb:
        for line in fb:
            if line.strip():
                if blast_sep is None:
                    blast_sep = "\t" if "\t" in line else None
                break

    id_idx = args.blast_id_col - 1
    missing_vec = [args.missing] * n_anno_cols

    n_in = 0
    n_found_exact = 0
    n_found_gene = 0
    n_missing = 0

    with blast_path.open("r", encoding="utf-8", errors="replace") as fin, \
         out_path.open("w", encoding="utf-8", newline="") as fout:

        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        for line in fin:
            if not line.strip():
                continue
            n_in += 1

            if blast_sep is None:
                # whitespace split
                row = line.strip().split()
            else:
                row = line.rstrip("\n").split(blast_sep)

            if id_idx >= len(row):
                # Malformed line: still output with missing appended
                writer.writerow(row + missing_vec)
                n_missing += 1
                continue

            at_id = normalize_id(row[id_idx])

            appended = None
            if at_id in anno_map_exact:
                appended = anno_map_exact[at_id]
                n_found_exact += 1
            else:
                g = gene_level_id(at_id)
                if g in anno_map_gene:
                    appended = anno_map_gene[g]
                    n_found_gene += 1

            if appended is None:
                appended = missing_vec
                n_missing += 1

            writer.writerow(row + appended)

    # Report to stderr (so it won't pollute redirected outputs)
    print(
        f"[done] lines_in={n_in}  found_exact={n_found_exact}  "
        f"found_gene_fallback={n_found_gene}  missing={n_missing}  "
        f"anno_cols_appended={n_anno_cols}",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()
