#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#usege
#python make_db_html.py \
#  --in Template.html \
#  --out FoxtailMillet.html \
#  --title-prefix "Foxtail Millet" \
#  --h1-prefix "Foxtail Millet" \
#  --gene-annot "FoxtailMillet_annot.txt" \
#  --ref-label "Yugu1"

import argparse
import re
import sys
from pathlib import Path

def sub_once(pattern: str, repl: str, text: str, label: str) -> str:
    new_text, n = re.subn(pattern, repl, text, count=0, flags=re.DOTALL)
    if n != 1:
        raise RuntimeError(f"[{label}] expected 1 replacement, but got {n}. pattern={pattern}")
    return new_text

def main():
    ap = argparse.ArgumentParser(
        description="Replace key fields in RNA-seq DB HTML template."
    )
    ap.add_argument("--in", dest="infile", required=True, help="Input template HTML (e.g., Template.html)")
    ap.add_argument("--out", dest="outfile", required=True, help="Output HTML path")
    ap.add_argument("--title-prefix", required=True,
                    help='Replaces the prefix in <title>XXX RNA-seq database</title> (e.g., "Foxtail Millet")')
    ap.add_argument("--h1-prefix", required=True,
                    help='Replaces the prefix in <h1>XXX RNA-seq database</h1>')
    ap.add_argument("--gene-annot", required=True,
                    help='Replaces the quoted filename in const PATH_GENE_ANNOT = "..."')
    ap.add_argument("--ref-label", required=True,
                    help='Replaces the label in <option value="ref">XXX</option>')

    args = ap.parse_args()

    inp = Path(args.infile)
    outp = Path(args.outfile)

    html = inp.read_text(encoding="utf-8")

    # 1) <title>...</title>
    html = sub_once(
        r'(<title>)(.*?)(\s*RNA-seq database</title>)',
        r'\1' + args.title_prefix + r'\3',
        html,
        "title"
    )

    # 2) <h1>...</h1>
    html = sub_once(
        r'(<h1>)(.*?)(\s*RNA-seq database</h1>)',
        r'\1' + args.h1_prefix + r'\3',
        html,
        "h1"
    )

    # 3) const PATH_GENE_ANNOT = "..."
    # （中身が Template_annot.txt でも Template.txt でも、とにかくダブルクォート内を差し替え）
    html = sub_once(
        r'(const\s+PATH_GENE_ANNOT\s*=\s*")([^"]*)(")',
        r'\1' + args.gene_annot + r'\3',
        html,
        "PATH_GENE_ANNOT"
    )

    # 4) <option value="ref">...</option>
    html = sub_once(
        r'(<option\s+value="ref"\s*>)(.*?)(</option>)',
        r'\1' + args.ref_label + r'\3',
        html,
        "option_ref_label"
    )

    outp.write_text(html, encoding="utf-8")
    print(f"OK: wrote {outp}", file=sys.stderr)

if __name__ == "__main__":
    main()