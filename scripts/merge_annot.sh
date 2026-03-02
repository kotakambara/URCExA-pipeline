#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  merge_AtOs_annot.sh --arab <Fox_at_annot.txt> --rice <Fox_Os_annot.txt> --output <Fox_annot.txt>

Notes:
  - Input TSV: 1st col = gene_id, 2nd col = annotation (may be empty, duplicates allowed)
  - Output TSV: gene_id <tab> Arabidopsis_annotations <tab> Rice_annotations
  - For backward compatibility: --input <file> is treated as --arab <file>
USAGE
}

arab=""
rice=""
output=""
input=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arab|--at)   arab=${2:-}; shift 2 ;;
    --rice|--os)   rice=${2:-}; shift 2 ;;
    --output|-o)   output=${2:-}; shift 2 ;;
    --input)       input=${2:-}; shift 2 ;;  # compatibility alias
    -h|--help)     usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# Treat --input as alias of --arab if --arab not explicitly provided
if [[ -z "$arab" && -n "$input" ]]; then
  arab="$input"
fi

if [[ -z "$arab" || -z "$rice" || -z "$output" ]]; then
  echo "Error: --arab (or --input), --rice, and --output are required." >&2
  usage >&2
  exit 1
fi

[[ -f "$arab" ]] || { echo "Error: --arab file not found: $arab" >&2; exit 1; }
[[ -f "$rice" ]] || { echo "Error: --rice file not found: $rice" >&2; exit 1; }

awk -F'\t' -v OFS='\t' '
  # add unique v into vals[k], joined by ";"
  function add_val(vals, seen, k, v, key) {
    # trim CR if any
    sub(/\r$/, "", k); sub(/\r$/, "", v)

    if (v == "") return
    key = k SUBSEP v
    if (seen[key]++) return

    if (vals[k] != "") vals[k] = vals[k] ";" v
    else vals[k] = v
  }

  NR==FNR {
    k=$1
    v=(NF>=2 ? $2 : "")

    if (!(k in at_seen)) { at_seen[k]=1; at_order[++at_n]=k }
    add_val(at_vals, at_val_seen, k, v)
    next
  }

  {
    k=$1
    v=(NF>=2 ? $2 : "")

    if (!(k in os_seen)) { os_seen[k]=1; os_order[++os_n]=k }
    add_val(os_vals, os_val_seen, k, v)
  }

  END {
    for (i=1; i<=at_n; i++) {
      k=at_order[i]
      print k, at_vals[k], os_vals[k]
      printed[k]=1
    }
    for (i=1; i<=os_n; i++) {
      k=os_order[i]
      if (!printed[k]) print k, at_vals[k], os_vals[k]
    }
  }
' "$arab" "$rice" > "$output"

echo "Wrote: $output"
