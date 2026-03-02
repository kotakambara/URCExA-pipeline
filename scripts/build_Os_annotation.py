#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import hashlib
import re
from collections import defaultdict

def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def build_oryzabase_map(gene_en_tsv: str):
    """
    MSU ID -> fields
    重要：同一MSU IDに複数行が来る場合、最後の行が勝つ（上書き）= 添付ファイルと一致するルール
    """
    msu_to = {}
    with open_maybe_gz(gene_en_tsv, "rt") as f:
        f.readline()  # header
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            msu_field = parts[11]
            if not msu_field or msu_field == "_":
                continue
            for msu in re.split(r"\s*,\s*", msu_field):
                msu = msu.strip()
                if not msu:
                    continue
                # overwrite (last wins)
                msu_to[msu] = {
                    "CGSNL Gene Symbol": parts[1],
                    "Gene symbol synonym(s)": parts[2],
                    "CGSNL Gene Name": parts[3],
                    "Gene name synonym(s)": parts[4],
                    "Protein Name": parts[5],
                    "Allele": parts[6],
                    "Chromosome No.": parts[7],
                    "Explanation": parts[8],
                    "Trait Class": parts[9],
                    "RAP ID": parts[10],
                }
    return msu_to

def build_go_map(goslim_tsv: str):
    """
    model -> ordered unique GO IDs
    """
    go_map = defaultdict(list)
    seen = defaultdict(set)
    with open_maybe_gz(goslim_tsv, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            model, goid = cols[0], cols[1]
            if goid not in seen[model]:
                seen[model].add(goid)
                go_map[model].append(goid)
    return go_map

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-en", required=True)
    ap.add_argument("--func", required=True)
    ap.add_argument("--goslim", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--check-equal", default=None, help="gold file path (optional)")
    args = ap.parse_args()

    msu_to = build_oryzabase_map(args.gene_en)
    go_map = build_go_map(args.goslim)

    with open_maybe_gz(args.func, "rt") as f, open(args.out, "w", encoding="utf-8", newline="") as out:
        f.readline()  # header: model\tannotation
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            model, ann = line.split("\t", 1)

            # 15 columns, no header
            cols = [""] * 15
            cols[0] = model
            cols[1] = ann

            if model in msu_to:
                d = msu_to[model]
                cols[2]  = d["CGSNL Gene Symbol"]
                cols[3]  = d["Gene symbol synonym(s)"]
                cols[4]  = d["CGSNL Gene Name"]
                cols[5]  = d["Gene name synonym(s)"]
                cols[6]  = d["Gene name synonym(s)"]   # duplicate (matches attached file)
                cols[7]  = d["Protein Name"]
                cols[8]  = d["Allele"]
                cols[9]  = d["Chromosome No."]
                cols[10] = d["Explanation"]
                cols[11] = d["Explanation"]            # duplicate (matches attached file)
                cols[12] = d["Trait Class"]
                cols[13] = d["RAP ID"]

            if model in go_map:
                cols[14] = ", ".join(go_map[model])

            out.write("\t".join(cols) + "\n")

    if args.check_equal:
        got = sha256_file(args.out)
        exp = sha256_file(args.check_equal)
        if got != exp:
            raise SystemExit(f"[ERROR] Not identical.\n  got: {got}\n  exp: {exp}")
        print(f"[OK] Byte-identical (sha256={got})")

if __name__ == "__main__":
    main()