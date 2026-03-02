#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preprocess a huge CSV (meta cols + many gene cols) into:
  - meta.tsv
  - gene_index.tsv
  - manifest.json
  - shards/shard_XXX.bin (Float32, little-endian, row-major by sample)

This is a one-time offline step. You can run it on your local PC (Python3),
then upload the outputs to the FreeBSD web server (static files only).

CSV assumptions:
  - header row exists
  - first N_META columns are metadata (default: 10)
  - remaining columns are genes (TPM numeric; empty -> NaN)

If your CSV contains commas inside quoted fields, Python's csv module handles it.
"""

import csv
import json
import math
import os
import struct
import sys

N_META = 10
SHARD_GENES = 1024   # adjust: 512/1024/2048 depending on your network and gene count

INPUT_CSV = sys.argv[1]
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "out_sharded"

META_TSV = os.path.join(OUT_DIR, "meta.tsv")
IDX_TSV  = os.path.join(OUT_DIR, "gene_index.tsv")
MANIFEST = os.path.join(OUT_DIR, "manifest.json")
SHARD_DIR = os.path.join(OUT_DIR, "shards")

def safe_float(x):
  if x is None:
    return float("nan")
  s = str(x).strip()
  if s == "" or s.lower() == "na":
    return float("nan")
  try:
    return float(s)
  except Exception:
    return float("nan")

def main():
  os.makedirs(OUT_DIR, exist_ok=True)
  os.makedirs(SHARD_DIR, exist_ok=True)

  # --- Read header to get columns ---
  with open(INPUT_CSV, "r", newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    header = next(reader, None)
    if not header:
      raise SystemExit("ERROR: input CSV has no header")

  meta_cols = header[:N_META]
  gene_cols = header[N_META:]
  if not gene_cols:
    raise SystemExit("ERROR: no gene columns detected (check N_META)")

  # --- Plan shards ---
  n_genes = len(gene_cols)
  shard_gene_counts = []
  gene_to_shard = []  # list of (gene, shardId, colInShard)
  shard_id = 0
  for start in range(0, n_genes, SHARD_GENES):
    end = min(start + SHARD_GENES, n_genes)
    cnt = end - start
    shard_gene_counts.append(cnt)
    for j, gene in enumerate(gene_cols[start:end]):
      gene_to_shard.append((gene, shard_id, j))
    shard_id += 1

  n_shards = len(shard_gene_counts)
  print(f"[INFO] meta cols: {len(meta_cols)}, genes: {n_genes}, shards: {n_shards} (SHARD_GENES={SHARD_GENES})")

  # --- Write gene_index.tsv ---
  with open(IDX_TSV, "w", encoding="utf-8") as f:
    for gene, sid, col in gene_to_shard:
      f.write(f"{gene}\t{sid}\t{col}\n")

  # --- First pass: write meta.tsv and count samples (n_samples) ---
  n_samples = 0
  with open(INPUT_CSV, "r", newline="", encoding="utf-8") as fin, open(META_TSV, "w", encoding="utf-8", newline="\n") as mout:
    reader = csv.reader(fin)
    hdr = next(reader)  # skip header
    mout.write("\t".join(meta_cols) + "\n")
    for row in reader:
      if not row:
        continue
      meta = row[:N_META]
      # ensure length
      if len(meta) < N_META:
        meta = meta + [""] * (N_META - len(meta))
      mout.write("\t".join(meta) + "\n")
      n_samples += 1
      if n_samples % 20000 == 0:
        print(f"[INFO] meta pass: {n_samples} samples")

  print(f"[INFO] samples: {n_samples}")

  # --- Second pass: build shards efficiently with buffering ---
  # We buffer per shard as bytearray to reduce syscall overhead.
  shard_files = []
  for sid in range(n_shards):
    path = os.path.join(SHARD_DIR, f"shard_{sid:03d}.bin")
    shard_files.append(open(path, "wb"))

  # Create mapping from global gene index -> (shardId, colInShard)
  # For row-major writing, we split each row's gene values into shard blocks and append.
  # We'll write in chunks of rows to amortize overhead.
  BUFFER_ROWS = 128  # flush every N rows
  buffers = [bytearray() for _ in range(n_shards)]
  written_rows = 0

  with open(INPUT_CSV, "r", newline="", encoding="utf-8") as fin:
    reader = csv.reader(fin)
    _ = next(reader)  # skip header
    for row in reader:
      if not row:
        continue
      genes = row[N_META:]
      # pad if short
      if len(genes) < n_genes:
        genes = genes + [""] * (n_genes - len(genes))

      # Convert once per row into floats (still heavy but only 1 pass)
      # If this is too slow, consider saving as TSV and using faster parsing tools.
      vals = [safe_float(x) for x in genes]

      # Split into shards and pack float32 little-endian
      base = 0
      for sid, cnt in enumerate(shard_gene_counts):
        block = vals[base:base+cnt]
        # pack as little-endian float32
        buffers[sid] += struct.pack("<%sf" % cnt, *block)
        base += cnt

      written_rows += 1
      if written_rows % BUFFER_ROWS == 0:
        for sid, bf in enumerate(buffers):
          shard_files[sid].write(bf)
          buffers[sid].clear()
      if written_rows % 5000 == 0:
        print(f"[INFO] shard pass: {written_rows}/{n_samples} rows")

  # final flush
  for sid, bf in enumerate(buffers):
    if bf:
      shard_files[sid].write(bf)
      bf.clear()

  for f in shard_files:
    f.close()

  # --- Write manifest.json ---
  manifest = {
    "n_samples": n_samples,
    "n_meta": N_META,
    "shard_genes": SHARD_GENES,
    "n_genes": n_genes,
    "n_shards": n_shards,
    "shard_gene_counts": shard_gene_counts,
  }
  with open(MANIFEST, "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2)

  print("[DONE] Outputs written to:", OUT_DIR)
  print("  - meta.tsv")
  print("  - gene_index.tsv")
  print("  - manifest.json")
  print("  - shards/shard_XXX.bin")

if __name__ == "__main__":
  main()
