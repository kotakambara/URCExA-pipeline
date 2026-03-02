#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="data"
THREADS=4

URL_ORYZABASE="https://shigen.nig.ac.jp/rice/oryzabase/gene/download?classtag=GENE_EN_LIST"
URL_FUNC="https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.functional_annotation.txt.gz"
URL_GOSLIM="https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.GOSlim.txt.gz"
URL_PEP="https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.pep.fa.gz"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) echo "Unknown option: $1" >&2; exit 2;;
  esac
done

command -v python3 >/dev/null
command -v makeblastdb >/dev/null
command -v gzip >/dev/null

if command -v curl >/dev/null; then
  DL() { curl -L --retry 3 --retry-delay 2 -A "Mozilla/5.0" -o "$2" "$1"; }
elif command -v wget >/dev/null; then
  DL() { wget -O "$2" "$1"; }
else
  echo "Need curl or wget" >&2; exit 2
fi

mkdir -p "$DATA_DIR"

echo "[DL] gene_en_os_list.tsv"
DL "$URL_ORYZABASE" "$DATA_DIR/gene_en_os_list.tsv"

echo "[DL] functional_annotation.txt.gz"
DL "$URL_FUNC" "$DATA_DIR/osa1_r7.all_models.functional_annotation.txt.gz"
gzip -dc "$DATA_DIR/osa1_r7.all_models.functional_annotation.txt.gz" > "$DATA_DIR/osa1_r7.all_models.functional_annotation.txt"

echo "[DL] GOSlim.txt.gz"
DL "$URL_GOSLIM" "$DATA_DIR/osa1_r7.all_models.GOSlim.txt.gz"
gzip -dc "$DATA_DIR/osa1_r7.all_models.GOSlim.txt.gz" > "$DATA_DIR/osa1_r7.all_models.GOSlim.txt"

echo "[DL] all_models.pep.fa.gz"
DL "$URL_PEP" "$DATA_DIR/osa1_r7.all_models.pep.fa.gz"
gzip -dc "$DATA_DIR/osa1_r7.all_models.pep.fa.gz" > "$DATA_DIR/Os_proteins.fasta"

echo "[BLASTDB] makeblastdb"
makeblastdb -in "$DATA_DIR/Os_proteins.fasta" -dbtype prot -out "$DATA_DIR/Os_proteins.fasta" >/dev/null

echo "[BUILD] Os_gene_annotation_list.tsv"
python3 build_Os_annotation.py \
  --gene-en "$DATA_DIR/gene_en_os_list.tsv" \
  --func "$DATA_DIR/osa1_r7.all_models.functional_annotation.txt" \
  --goslim "$DATA_DIR/osa1_r7.all_models.GOSlim.txt" \
  --out "$DATA_DIR/Os_gene_annotation_list.tsv"

echo "[DONE] $DATA_DIR/Os_gene_annotation_list.tsv and BLAST DB are ready."