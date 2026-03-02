in=${1:?}
out=${2:?}

awk -F'\t' 'BEGIN{OFS="\t"}
{
  gsub(/\r$/, "", $0)                # 念のためCRLF対策
  sub(/\.[0-9]+$/, "", $1)           # 1列目のみ: 例 dpca0g000640.840.1 -> dpca0g000640.840

  key = $1 OFS $2 OFS $3
  if (!seen[key]++) print $1, $2, $3 # 3列で完全一致した行だけを重複除去
}' "$in" > "$out"