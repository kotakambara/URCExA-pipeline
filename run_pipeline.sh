#!/usr/bin/env bash
set -euo pipefail

# RNA-seq DB pipeline runner
# - Runs scripts in ./scripts/ in the order described in "スクリプトの流れ.txt"
# - Designed for reproducible, shareable runs (GitHub-friendly)

usage() {
  cat <<'USAGE'
Usage:
  ./run_pipeline.sh \
    --species "Setaria italica" [--species "foxtail millet" ...] \
    --prefix SI_Yugu1 \
    --star-index /path/to/STAR_index_dir \
    --gff /path/to/genes.gff3 \
    --proteins /path/to/species_proteins.fa \
    --title-prefix "Foxtail Millet" \
    --h1-prefix "Foxtail Millet" \
    --ref-label "Yugu1" \
    --outdir out/SI_Yugu1

Required:
  --species      Organism name(s) used for NCBI SRA search (repeatable)
  --prefix       Output prefix (used in filenames)
  --star-index   STAR genomeDir (pre-built index)
  --gff          GFF3 annotation (gene features must have ID=...)
  --proteins     Protein FASTA for the target species (headers must match gene IDs)
  --title-prefix <title> prefix in HTML ("XXX RNA-seq database")
  --h1-prefix    <h1> prefix in HTML ("XXX RNA-seq database")
  --ref-label    Label shown in HTML reference dropdown (e.g., "Yugu1")
  --outdir       Output directory (will be created)

Optional:
  --template     HTML template (default: scripts/Template.html)
  --map-tsv      Replacement rules TSV for standardize_metadata.py (repeatable)
  --replace      Inline replacement "column:FROM=TO" (repeatable)
  --srr-list     Skip Step1 and use this list file (same format as tmp_list)
  --max-samples  Only map the first N lines of the SRR list (for testing)
  --from         Start step: search|map|clean|filter|blast|html|shard (default: search)
  --to           Stop after step: search|map|clean|filter|blast|html|shard (default: shard)
  --force        Re-run steps even if outputs exist
  --keep-work    Do not delete work dir at the end

Notes:
  - This pipeline downloads potentially LARGE SRA datasets. Prefer running on HPC/storage.
  - NCBI eutils has rate limits; consider setting an NCBI API key in your environment.
USAGE
}

log() { echo "[$(date +'%F %T')] $*"; }

die() { echo "ERROR: $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"
}

need_py_mod() {
  local mod="$1"
  python3 - <<PY >/dev/null 2>&1 || die "Missing Python module: ${mod} (install it in your python3 env)"
import ${mod}
PY
}

# --- defaults ---
TEMPLATE=""
OUTDIR=""
PREFIX=""
STAR_INDEX=""
GFF=""
PROTEINS=""
SRR_LIST=""
MAX_SAMPLES=""
FROM="search"
TO="shard"
FORCE=0
KEEP_WORK=0
MAP_TSVS=()
REPLACES=()
SPECIES=()

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --species) SPECIES+=("$2"); shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --star-index) STAR_INDEX="$2"; shift 2 ;;
    --gff) GFF="$2"; shift 2 ;;
    --proteins) PROTEINS="$2"; shift 2 ;;
    --title-prefix) TITLE_PREFIX="$2"; shift 2 ;;
    --h1-prefix) H1_PREFIX="$2"; shift 2 ;;
    --ref-label) REF_LABEL="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --template) TEMPLATE="$2"; shift 2 ;;
    --map-tsv) MAP_TSVS+=("$2"); shift 2 ;;
    --replace) REPLACES+=("$2"); shift 2 ;;
    --srr-list) SRR_LIST="$2"; shift 2 ;;
    --max-samples) MAX_SAMPLES="$2"; shift 2 ;;
    --from) FROM="$2"; shift 2 ;;
    --to) TO="$2"; shift 2 ;;
    --force) FORCE=1; shift 1 ;;
    --keep-work) KEEP_WORK=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$OUTDIR" ]] || { usage; die "--outdir is required"; }
[[ -n "$PREFIX" ]] || { usage; die "--prefix is required"; }
[[ -n "$STAR_INDEX" ]] || { usage; die "--star-index is required"; }
[[ -n "$GFF" ]] || { usage; die "--gff is required"; }
[[ -n "$PROTEINS" ]] || { usage; die "--proteins is required"; }
[[ -n "${TITLE_PREFIX:-}" ]] || { usage; die "--title-prefix is required"; }
[[ -n "${H1_PREFIX:-}" ]] || { usage; die "--h1-prefix is required"; }
[[ -n "${REF_LABEL:-}" ]] || { usage; die "--ref-label is required"; }

if [[ -z "$SRR_LIST" ]]; then
  [[ ${#SPECIES[@]} -gt 0 ]] || { usage; die "Provide at least one --species or give --srr-list"; }
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
[[ -d "$SCRIPTS_DIR" ]] || die "scripts/ not found: $SCRIPTS_DIR"

if [[ -z "$TEMPLATE" ]]; then
  TEMPLATE="${SCRIPTS_DIR}/Template.html"
fi

# --- dependency checks (minimal, add more if you need) ---
need_cmd python3
need_cmd awk
need_cmd sed
need_cmd grep
need_cmd paste
need_cmd cut
need_cmd sort
need_cmd uniq
need_cmd aria2c
need_cmd wget
need_cmd fasterq-dump
need_cmd fastp
need_cmd STAR
need_cmd samtools
need_cmd featureCounts
need_cmd blastp

need_py_mod bs4
need_py_mod pandas

[[ -d "$STAR_INDEX" ]] || die "STAR index dir not found: $STAR_INDEX"
[[ -f "$GFF" ]] || die "GFF3 not found: $GFF"
[[ -f "$PROTEINS" ]] || die "Protein FASTA not found: $PROTEINS"
[[ -f "$TEMPLATE" ]] || die "Template HTML not found: $TEMPLATE"

mkdir -p "$OUTDIR"
RESULTS_DIR="${OUTDIR}/results"
SITE_DIR="${OUTDIR}/site"
WORKDIR="${OUTDIR}/work"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$RESULTS_DIR" "$SITE_DIR" "$WORKDIR" "$LOGDIR"

# stage scripts into workdir (symlinks) so relative calls inside scripts keep working
STAGE_DIR="${WORKDIR}/stage_scripts"
mkdir -p "$STAGE_DIR"
for f in API_RNA-seq.py findSAMN_info.py get_expt_from_srr.py mapping_script.sh Count_to_TPM.py \
         standardize_metadata.py check_sample.sh blastp_annot.sh blast_append_annotation.py \
         merge_annot.sh drop_dot.sh make_html.py preprocess_shard.py; do
  ln -sf "${SCRIPTS_DIR}/${f}" "${STAGE_DIR}/${f}"
done

# Optional: link data/ if users put it at repo root (recommended)
if [[ -d "${SCRIPT_DIR}/data" ]]; then
  ln -sf "${SCRIPT_DIR}/data" "${STAGE_DIR}/data"
fi

step_id() {
  case "$1" in
    search) echo 1 ;;
    map) echo 2 ;;
    clean) echo 3 ;;
    filter) echo 4 ;;
    blast) echo 5 ;;
    html) echo 6 ;;
    shard) echo 7 ;;
    *) die "Unknown step name: $1" ;;
  esac
}

FROM_ID=$(step_id "$FROM")
TO_ID=$(step_id "$TO")
[[ $FROM_ID -le $TO_ID ]] || die "--from must be <= --to"

run_step() {
  local sid="$1"; shift
  local name="$1"; shift
  if [[ $sid -lt $FROM_ID || $sid -gt $TO_ID ]]; then
    log "[SKIP] step ${sid} (${name}) out of range"
    return 0
  fi
  "$@"
}

# -------------------------
# Step 1: Search SRA (API_RNA-seq.py)
# -------------------------
step_search() {
  cd "$STAGE_DIR"
  local list_out="${WORKDIR}/tmp_list"

  if [[ -n "$SRR_LIST" ]]; then
    log "Using provided SRR list: $SRR_LIST"
    cp -f "$SRR_LIST" "$list_out"
  else
    if [[ -f "$list_out" && $FORCE -ne 1 ]]; then
      log "tmp_list exists, skipping search: $list_out"
    else
      log "Running API_RNA-seq.py (this queries NCBI SRA)..."
      python3 API_RNA-seq.py "${SPECIES[@]}" 2>"${LOGDIR}/01_search.stderr.log" 1>"${LOGDIR}/01_search.stdout.log"
      [[ -f tmp_list ]] || die "API_RNA-seq.py did not create tmp_list"
      cp -f tmp_list "$list_out"
    fi
  fi

  # optional truncation for testing
  if [[ -n "${MAX_SAMPLES}" ]]; then
    log "Limiting SRR list to first ${MAX_SAMPLES} lines (for testing)"
    head -n "${MAX_SAMPLES}" "$list_out" > "${WORKDIR}/tmp_list.max${MAX_SAMPLES}"
    mv -f "${WORKDIR}/tmp_list.max${MAX_SAMPLES}" "$list_out"
  fi

  log "SRR list ready: $list_out"
}

# -------------------------
# Step 2-3: Mapping + Count/TPM matrix (mapping_script.sh)
# -------------------------
step_map() {
  cd "$STAGE_DIR"

  local tpm_csv="${STAGE_DIR}/${PREFIX}_TPM_matrix.csv"
  local cnt_csv="${STAGE_DIR}/${PREFIX}_count_matrix.csv"
  local meta_tsv="${STAGE_DIR}/${PREFIX}_metadata.tsv"

  if [[ -f "$tpm_csv" && $FORCE -ne 1 ]]; then
    log "TPM matrix exists, skipping mapping: $tpm_csv"
  else
    # clean previous partial outputs if --force
    if [[ $FORCE -eq 1 ]]; then
      rm -f "$tpm_csv" "$cnt_csv" "$meta_tsv"
    fi

    log "Running mapping_script.sh (downloads SRA, runs fastp/STAR/featureCounts, builds matrices)..."
    # mapping_script.sh uses relative python calls; run in STAGE_DIR.
    # ulimit may fail on some systems; do not crash pipeline.
    (ulimit -n 100000 || true)

    bash mapping_script.sh "${WORKDIR}/tmp_list" "$PREFIX" "$STAR_INDEX" "$GFF" \
      1>"${LOGDIR}/02_map.stdout.log" 2>"${LOGDIR}/02_map.stderr.log" || true

    [[ -f "$tpm_csv" ]] || die "Mapping did not produce TPM matrix: $tpm_csv"
  fi

  cp -f "$tpm_csv" "$RESULTS_DIR/"
  [[ -f "$cnt_csv" ]] && cp -f "$cnt_csv" "$RESULTS_DIR/" || true
  [[ -f "$meta_tsv" ]] && cp -f "$meta_tsv" "$RESULTS_DIR/" || true

  log "Saved matrices to: $RESULTS_DIR"
}

# -------------------------
# Step 4a: Standardize metadata (standardize_metadata.py)
# -------------------------
step_clean() {
  cd "$STAGE_DIR"
  local in_csv="${RESULTS_DIR}/${PREFIX}_TPM_matrix.csv"
  local out_csv="${RESULTS_DIR}/${PREFIX}_TPM_matrix.cleaned.csv"

  [[ -f "$in_csv" ]] || die "Missing input CSV: $in_csv (run map step first)"

  if [[ -f "$out_csv" && $FORCE -ne 1 ]]; then
    log "Cleaned CSV exists, skipping: $out_csv"
    return 0
  fi

  log "Running standardize_metadata.py ..."
  args=( "$in_csv" "$out_csv" )
  for tsv in "${MAP_TSVS[@]}"; do
    args+=( --map-tsv "$tsv" )
  done
  for r in "${REPLACES[@]}"; do
    args+=( --replace "$r" )
  done

  python3 standardize_metadata.py "${args[@]}" \
    1>"${LOGDIR}/03_clean.stdout.log" 2>"${LOGDIR}/03_clean.stderr.log"

  [[ -f "$out_csv" ]] || die "standardize_metadata.py failed to create: $out_csv"
  log "Cleaned matrix: $out_csv"
}

# -------------------------
# Step 4b: Filter obviously bad samples (check_sample.sh)
# -------------------------
step_filter() {
  cd "$STAGE_DIR"
  local in_csv="${RESULTS_DIR}/${PREFIX}_TPM_matrix.cleaned.csv"
  local out_prefix="${RESULTS_DIR}/${PREFIX}_TPM_matrix.cleaned"
  local out_filt="${out_prefix}.filtered.csv"

  [[ -f "$in_csv" ]] || die "Missing input CSV: $in_csv (run clean step first)"

  if [[ -f "$out_filt" && $FORCE -ne 1 ]]; then
    log "Filtered CSV exists, skipping: $out_filt"
    return 0
  fi

  log "Running check_sample.sh ..."
  bash check_sample.sh "$in_csv" "$out_prefix" \
    1>"${LOGDIR}/04_filter.stdout.log" 2>"${LOGDIR}/04_filter.stderr.log"

  [[ -f "$out_filt" ]] || die "check_sample.sh failed to create: $out_filt"
  log "Filtered matrix: $out_filt"
}

# -------------------------
# Step 5: BLASTP annotation vs Arabidopsis/Rice (blastp_annot.sh)
# -------------------------
step_blast() {
  cd "$STAGE_DIR"

  local annot_out="${STAGE_DIR}/${PREFIX}_annot.txt"

  if [[ -f "$annot_out" && $FORCE -ne 1 ]]; then
    log "Annotation exists, skipping BLAST: $annot_out"
  else
    log "Running blastp_annot.sh (requires ./data/* files; see README)..."
    bash blastp_annot.sh "$PROTEINS" "$PREFIX" \
      1>"${LOGDIR}/05_blast.stdout.log" 2>"${LOGDIR}/05_blast.stderr.log"

    [[ -f "$annot_out" ]] || die "blastp_annot.sh did not create: $annot_out"
  fi

  # Copy to site dir (Template.html expects this file next to HTML)
  cp -f "$annot_out" "${SITE_DIR}/${PREFIX}_annot.txt"
  log "Gene annotation TSV: ${SITE_DIR}/${PREFIX}_annot.txt"
}

# -------------------------
# Step 6: Generate HTML (make_html.py)
# -------------------------
step_html() {
  cd "$STAGE_DIR"
  local out_html="${SITE_DIR}/${PREFIX}.html"

  if [[ -f "$out_html" && $FORCE -ne 1 ]]; then
    log "HTML exists, skipping: $out_html"
    return 0
  fi

  log "Generating HTML from template..."
  python3 make_html.py \
    --in "$TEMPLATE" \
    --out "$out_html" \
    --title-prefix "$TITLE_PREFIX" \
    --h1-prefix "$H1_PREFIX" \
    --gene-annot "${PREFIX}_annot.txt" \
    --ref-label "$REF_LABEL" \
    1>"${LOGDIR}/06_html.stdout.log" 2>"${LOGDIR}/06_html.stderr.log"

  [[ -f "$out_html" ]] || die "make_html.py failed to create: $out_html"
  log "HTML: $out_html"
}

# -------------------------
# Step 7: Shard huge matrix for fast web loading (preprocess_shard.py)
# -------------------------
step_shard() {
  cd "$STAGE_DIR"
  local in_csv="${RESULTS_DIR}/${PREFIX}_TPM_matrix.cleaned.filtered.csv"
  [[ -f "$in_csv" ]] || die "Missing filtered matrix: $in_csv (run filter step first)"

  # If manifest already exists and not forcing, skip
  if [[ -f "${SITE_DIR}/manifest.json" && -d "${SITE_DIR}/shards" && $FORCE -ne 1 ]]; then
    log "Sharded outputs exist, skipping preprocess_shard.py"
    return 0
  fi

  log "Running preprocess_shard.py (this can take time for big matrices)..."
  python3 preprocess_shard.py "$in_csv" "$SITE_DIR" \
    1>"${LOGDIR}/07_shard.stdout.log" 2>"${LOGDIR}/07_shard.stderr.log"

  [[ -f "${SITE_DIR}/manifest.json" ]] || die "preprocess_shard.py failed: manifest.json missing"
  log "Sharded outputs written under: $SITE_DIR"
}

# ---------- execute ----------
run_step 1 search step_search
run_step 2 map step_map
run_step 3 clean step_clean
run_step 4 filter step_filter
run_step 5 blast step_blast
run_step 6 html step_html
run_step 7 shard step_shard

log "DONE"
log "Results (CSV/TSV): $RESULTS_DIR"
log "Website folder to upload: $SITE_DIR"
log "Open in browser after upload: ${PREFIX}.html"

if [[ $KEEP_WORK -ne 1 ]]; then
  # keep staged scripts but remove huge temp files, FASTQs etc (workdir can be very large)
  log "Cleaning workdir (use --keep-work to keep)"
  rm -rf "${WORKDIR}/stage_scripts"/*.fastq "${WORKDIR}/stage_scripts"/*.sra 2>/dev/null || true
fi
