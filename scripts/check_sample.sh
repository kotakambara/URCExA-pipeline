# ---- settings (必要に応じて調整) ----
meta_n=10      # メタデータ列数（元スクリプトと同じ）:contentReference[oaicite:1]{index=1}
thr=0.001     # ← 今回の条件（TPMが0.00001未満）に合わせる

in=${1:?Usage: filter_lowTPM.sh input.csv [out_prefix]}
prefix=${2:-out}

out_list="${prefix}.lowTPM_samples.tsv"   # 条件に該当したサンプルの1,2列目
out_filt="${prefix}.filtered.csv"         # それらのサンプル行を除外したマトリックス

awk -F',' -v meta_n="$meta_n" -v thr="$thr" \
    -v out_list="$out_list" -v out_filt="$out_filt" '
NR==1{
  # ヘッダ
  print $1"\t"$2 > out_list;
  print $0 > out_filt;
  next
}
{
  low=0; total=0;

  # 遺伝子列（meta_n+1 〜 NF）で判定
  for(i=meta_n+1;i<=NF;i++){
    if($i!=""){
      total++;
      if(($i+0) < thr) low++;
    }
  }

  # 90%以上がthr未満なら drop
  drop = (total>0 && low*100 >= total*90);

  if(drop){
    print $1"\t"$2 >> out_list;
  } else {
    print $0 >> out_filt;
  }
}
' "$in"

echo "Wrote: $out_list"
echo "Wrote: $out_filt"
