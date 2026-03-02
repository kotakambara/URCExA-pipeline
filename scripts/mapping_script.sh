#!/usr/bin/env bash
set -euo pipefail

# ----------
# 引数
#   $1 : diff_list / tmp_list (srr,samn,bio,name,lib がタブ区切りで入ったファイル)
#   $2 : 出力ファイルのプレフィックス (例: PM)
#   $3 : STAR genome index
#   $4 : GFF3 アノテーション
# ----------

shopt -s expand_aliases
alias ddbj_sra='aria2c -x 10'
ulimit -n 100000

diff=${1:?}
out=${2:?}
ref=${3:?}
gff=${4:?}


meta_file="${out}_metadata.tsv"

# メタデータTSVが無ければ作成（4列目は sample_name に変更）
if [ ! -f "${meta_file}" ]; then
    echo -e "srr\tsamn\tbio\tsample_name\tlib\tattribute_raw" > "${meta_file}"
fi

# diff を「カンマ区切り」に正規化 (SRR,SAMN,BioProject,tmp_name,LibraryType)
uniq "${diff}" | sed 's/,/_/g' | sed $'s/\t/,/g' > diff_tmp

############################################
# SAMN_info_tmp から Attribute / Sample name を抜き出す関数
############################################

# 特定の attribute_name の値を取得 (例: cultivar, treatment, tissue, dev_stage, temp)
extract_from_samn () {
    local key="$1"
    local samn_xml
    samn_xml=$(sed 's/&lt;/</g; s/&gt;/>/g' SAMN_info_tmp 2>/dev/null | tr -d '\n')
    if [ -z "${samn_xml}" ]; then
        echo "NA"
        return
    fi

    local val
    val=$(printf '%s\n' "${samn_xml}" \
        | grep -o "<Attribute[^>]*attribute_name=\"${key}\"[^>]*>[^<]*" \
        | head -n1 \
        | sed 's/.*>//')

    if [ -z "${val}" ]; then
        echo "NA"
    else
        printf '%s\n' "${val}" \
            | sed 's/,/_/g; s/^[ \t]*//; s/[ \t]*$//'
    fi
}

# すべての Attribute を key:value;key2:value2;... に整形
extract_all_attributes () {
    local samn_xml
    samn_xml=$(sed 's/&lt;/</g; s/&gt;/>/g' SAMN_info_tmp 2>/dev/null | tr -d '\n')
    if [ -z "${samn_xml}" ]; then
        echo "NA"
        return
    fi

    printf '%s\n' "${samn_xml}" \
      | grep -o "<Attribute[^>]*>[^<]*</Attribute>" \
      | sed 's/.*attribute_name=\"\([^\"]*\)\"[^>]*>\([^<]*\).*/\1:\2/' \
      | sed 's/,/_/g; s/ /_/g' \
      | paste -sd';' -
}

# BioSample XML から Sample name を取り出す
extract_sample_name () {
    local samn_xml
    samn_xml=$(sed 's/&lt;/</g; s/&gt;/>/g' SAMN_info_tmp 2>/dev/null | tr -d '\n')
    if [ -z "${samn_xml}" ]; then
        echo ""
        return
    fi
    # <Id db_label="Sample name"></Id> を抜き出す
    printf '%s\n' "${samn_xml}" \
      | sed -n 's/.*<Id db_label="Sample name">\([^<]*\)<\/Id>.*/\1/p' \
      | sed 's/,/_/g; s/^[ \t]*//; s/[ \t]*$//'
}

############################################
# メインループ：1 行 = 1 サンプル
############################################

while IFS=, read -r srr samn bio old_name lib; do
    [ -z "${srr}" ] && continue

    echo "===== Processing ${srr} (lib=${lib}) ====="

    ######################
    # SAMN メタ情報取得
    ######################
    python findSAMN_info.py "${samn}"

    # Sample name を決定（BioSample → LIBRARY_NAME → SRR/ERR）
    if [ "${samn}" != "NA" ] && [ -n "${samn}" ] && [ ! -s SAMN_info_tmp ]; then
    echo "  [SKIP] metadata failed for ${srr} (SAMN=${samn})"

    # Download_failed_list に書き出し
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${srr}" "${samn}" "${bio}" "${old_name}" "${lib}" \
        >> Download_failed_list.txt

    rm -f SAMN_info_tmp "${srr}.sra" *.fastq
    continue
    fi

    sample_name=$(extract_sample_name)

    if [ -z "${sample_name}" ] || [ "${sample_name}" = "NA" ]; then
        sample_name="${old_name}"
    fi

    if [ -z "${sample_name}" ] || [ "${sample_name}" = "NA" ]; then
        sample_name="${srr}"
    fi

    # ファイル名などに使う code を作る
    code=$(printf '%s' "${sample_name}" | sed 's/[[:space:]]\+/_/g')

    # 各 Attribute を抽出
    cultivar=$(extract_from_samn "cultivar")   # line
    treatment=$(extract_from_samn "treatment")
    tissue=$(extract_from_samn "tissue")
    stage=$(extract_from_samn "dev_stage")     # developmental stage
    temperature=$(extract_from_samn "temp")
    attributes=$(extract_all_attributes)

    # 元の display_name ベースの attribute も一応保存しておく
    attribute_raw=$(grep "display_name=" SAMN_info_tmp 2>/dev/null \
        | sed 's/^[ \t]*//; s/[ \t]*$//' \
        | tr '\n' ';' \
        | sed 's/,/_/g')
    [ -z "${attribute_raw}" ] && attribute_raw="${attributes}"

    ############################
    # DDBJ から .sra を取得（リトライ＋wget＋fasterq-dump フォールバック）
    ############################
    expt=$(python get_expt_from_srr.py "${srr}")
    sexpt=$(echo "${expt}" | cut -c1-6)
    ssexpt=$(echo "${expt}" | cut -c1-3)
    sra="${srr}"

    ddbj_url="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/${ssexpt}/${sexpt}/${expt}/${sra}/${sra}.sra"
    echo "Downloading from DDBJ: ${ddbj_url}"

    download_ok=0
    fastq_ready=0

    # set -e を一時的に無効化してエラーを自前でハンドリング
    set +e

    # --- aria2c (alias ddbj_sra) で1回トライ ---
    for n in 1 2 3 ; do
        echo "  [aria2c] attempt ${n}/3 for ${sra}"
        ddbj_sra "https://sra-pub-run-odp.s3.amazonaws.com/sra/${sra}/${sra}"
        mv "${sra}" "${sra}.sra"
        if [ -s "${sra}.sra" ] || [ -s "${sra}" ]; then
            download_ok=1
            echo "  [aria2c] success on attempt ${n}"
            break
        fi
    done

    # --- 2-2. まだダメなら wget ---
    if [ "${download_ok}" -ne 1 ]; then
        echo "  aria2c failed 3 times. Trying wget..."
        rm -f "${sra}.sra"
        wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/${sra}/${sra}"
        mv "${sra}" "${sra}.sra"
        if [ -s "${sra}.sra" ] || [ -s "${sra}" ]; then
            download_ok=1
            echo "  [wget] success."
        else
            echo "  [wget] failed or file empty."
            rm -f "${sra}.sra"
        fi
    fi

    

    # --- それでもダメなら DDBJ (non-mirror) の sralite を試す（主に DRR 向け） ---
    #  ※ submitter 情報が必要になるため、edirect (efetch/xtract) が入っている環境でのみ動きます
    if [ "${download_ok}" -ne 1 ]; then
    if command -v efetch >/dev/null 2>&1 && command -v xtract >/dev/null 2>&1; then
        if [[ "${sra}" =~ ^DRR ]]; then
            echo "  Trying DDBJ sralite (non-mirror) for ${sra}..."
            expt_edirect=$(efetch -db sra -id "${sra}" -format docsum | xtract -pattern DocumentSummarySet -element Experiment@acc 2>/dev/null)
            submitter=$(efetch -db sra -id "${sra}" -format docsum | xtract -pattern DocumentSummarySet -element Submitter@acc 2>/dev/null)
            ssubmitter=$(echo "${submitter}" | cut -c1-6)

            if [ -n "${expt_edirect}" ]; then
                expt="${expt_edirect}"
                sexpt=$(echo "${expt}" | cut -c1-6)
                ssexpt=$(echo "${expt}" | cut -c1-3)

                ddbj_url="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/${ssexpt}/${sexpt}/${expt}/${sra}/${sra}.sra"
                echo "  [DDBJ sralite] ${ddbj_url}"

                if [ ! -f "${sra}.sra" ] && [ ! -f "${sra}.fastq" ] && [ ! -f "${sra}_1.fastq" ] && [ ! -f "${sra}_2.fastq" ]; then
                    ddbj_sra "${ddbj_url}"
                fi

                if [ -s "${sra}.sra" ]; then
                    download_ok=1
                    echo "  [DDBJ sralite] downloaded."
                fi
            else
                echo "  [DDBJ sralite] expt not found by efetch/xtract for ${sra}."
            fi
        fi
    else
        echo "  [SKIP] efetch/xtract not found. (edirect is required for DDBJ non-mirror download)"
    fi
    fi

# --- .sra から fasterq-dump ---
    if [ "${download_ok}" -eq 1 ]; then
        echo "  Converting ${sra}.sra to FASTQ by fasterq-dump..."
        if [ "${lib}" = "PAIRED" ]; then
            fasterq-dump --split-files -p "${sra}.sra" -e 12
            if [ -s "${sra}_1.fastq" ] && [ -s "${sra}_2.fastq" ]; then
                mv "${sra}_1.fastq" "${code}_raw_R1.fastq"
                mv "${sra}_2.fastq" "${code}_raw_R2.fastq"
                fastq_ready=1
            # PAIRED のはずが片方しか出ない場合は SINGLE 扱いへ切替
            elif [ -s "${sra}_1.fastq" ] && [ ! -s "${sra}_2.fastq" ]; then
                echo "  [WARN] lib=PAIRED but only ${sra}_1.fastq was generated. Override lib -> SINGLE and continue as single-end."
                lib="SINGLE"
                mv "${sra}_1.fastq" "${code}_raw_R.fastq"
                rm -f "${sra}_2.fastq"
                fastq_ready=1

            elif [ -s "${sra}.fastq" ]; then
                echo "  [WARN] lib=PAIRED but only ${sra}.fastq was generated. Override lib -> SINGLE and continue as single-end."
                lib="SINGLE"
                mv "${sra}.fastq" "${code}_raw_R.fastq"
                fastq_ready=1
            fi
        elif [ "${lib}" = "SINGLE" ]; then
            fasterq-dump -p "${sra}.sra" -e 12
            if [ -s "${sra}.fastq" ]; then
                mv "${sra}.fastq" "${code}_raw_R.fastq"
                fastq_ready=1
            fi
        else
            echo "Unknown library type: ${lib} (expect PAIRED or SINGLE)."
        fi
    fi

    # --- それでもダメなら accession に対して fasterq-dump を試す ---
    if [ "${fastq_ready}" -ne 1 ]; then
        echo " Failed converting .sra file. Attempting fasterq-dump ..."
        #ファイルを掃除
        rm -f "${sra}.sra" "${sra}.fastq" "${sra}_1.fastq" "${sra}_2.fastq"

        if [ "${lib}" = "PAIRED" ]; then
            prefetch "${sra}"
            fasterq-dump --split-files -p -e 12 "${sra}"
            if [ -s "${sra}_1.fastq" ] && [ -s "${sra}_2.fastq" ]; then
                mv "${sra}_1.fastq" "${code}_raw_R1.fastq"
                mv "${sra}_2.fastq" "${code}_raw_R2.fastq"
                fastq_ready=1

            # accession 直叩きでも片方だけなら SINGLE に切替
            elif [ -s "${sra}_1.fastq" ] && [ ! -s "${sra}_2.fastq" ]; then
                echo "  [WARN] lib=PAIRED but only ${sra}_1.fastq was generated by accession fasterq-dump. Override lib -> SINGLE and continue as single-end."
                lib="SINGLE"
                mv "${sra}_1.fastq" "${code}_raw_R.fastq"
                rm -f "${sra}_2.fastq"
                fastq_ready=1

            elif [ -s "${sra}.fastq" ]; then
                echo "  [WARN] lib=PAIRED but only ${sra}.fastq was generated by accession fasterq-dump. Override lib -> SINGLE and continue as single-end."
                lib="SINGLE"
                mv "${sra}.fastq" "${code}_raw_R.fastq"
                fastq_ready=1
            fi

        elif [ "${lib}" = "SINGLE" ]; then
            prefetch "${sra}"
            fasterq-dump -p -e 12 "${sra}"
            if [ -s "${sra}.fastq" ]; then
                mv "${sra}.fastq" "${code}_raw_R.fastq"
                fastq_ready=1
            fi
        fi
    fi



    # --- DDBJ public FASTQ (.bz2) を試す（submitter 情報が必要／主に DRR 向け） ---
    if [ "${fastq_ready}" -ne 1 ]; then
    if command -v efetch >/dev/null 2>&1 && command -v xtract >/dev/null 2>&1; then
        if [[ "${sra}" =~ ^DRR ]]; then
            echo "  Trying DDBJ public FASTQ (.bz2) for ${sra}..."

            # expt / submitter が未取得なら取得（既に 2-2b で取れていれば再利用）
            if [ -z "${expt:-}" ]; then
                expt=$(efetch -db sra -id "${sra}" -format docsum | xtract -pattern DocumentSummarySet -element Experiment@acc 2>/dev/null)
            fi
            if [ -z "${submitter:-}" ]; then
                submitter=$(efetch -db sra -id "${sra}" -format docsum | xtract -pattern DocumentSummarySet -element Submitter@acc 2>/dev/null)
            fi
            ssubmitter=$(echo "${submitter}" | cut -c1-6)

            if [ -n "${expt:-}" ] && [ -n "${submitter:-}" ]; then
                # FASTQ がまだ無い場合だけダウンロードを試す
                if  [ ! -f "${sra}.fastq" ] && [ ! -f "${sra}_1.fastq" ] && [ ! -f "${sra}_2.fastq" ]; then
                    ddbj_sra "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/${ssubmitter}/${submitter}/${expt}/${sra}.fastq.bz2"
                    ddbj_sra "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/${ssubmitter}/${submitter}/${expt}/${sra}_1.fastq.bz2"
                    ddbj_sra "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/${ssubmitter}/${submitter}/${expt}/${sra}_2.fastq.bz2"
                    ddbj_sra "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/${ssubmitter}/${submitter}/${expt}/${sra}_3.fastq.bz2"

                    if command -v pbzip2 >/dev/null 2>&1; then
                        if ls "${sra}"*.fastq.bz2 >/dev/null 2>&1; then
                            pbzip2 -d -p60 "${sra}"*.fastq.bz2
                        fi
                    else
                        echo "  [WARN] pbzip2 not found. Cannot decompress ${sra}*.fastq.bz2"
                    fi
                fi

                # 展開後に FASTQ を検出して、既存の命名規則へリネーム
                if [ -s "${sra}_1.fastq" ] && [ -s "${sra}_2.fastq" ]; then
                    mv "${sra}_1.fastq" "${code}_raw_R1.fastq"
                    mv "${sra}_2.fastq" "${code}_raw_R2.fastq"
                    fastq_ready=1
                elif [ -s "${sra}.fastq" ]; then
                    lib="SINGLE"
                    mv "${sra}.fastq" "${code}_raw_R.fastq"
                    fastq_ready=1
                elif [ -s "${sra}_1.fastq" ] && [ ! -s "${sra}_2.fastq" ]; then
                    lib="SINGLE"
                    mv "${sra}_1.fastq" "${code}_raw_R.fastq"
                    rm -f "${sra}_2.fastq"
                    fastq_ready=1
                fi
            else
                echo "  [DDBJ public FASTQ] expt/submitter not found for ${sra}."
            fi
        fi
    else
        echo "  [SKIP] efetch/xtract not found. (edirect is required for DDBJ public FASTQ download)"
    fi
    fi

    # set -e を再有効化
    set -e
    # メタデータTSVに追記（lib の上書き後の値を記録できるようにここへ移動）
    echo -e "${srr}\t${samn}\t${bio}\t${sample_name}\t${lib}\t${attribute_raw}" >> "${meta_file}"

    # --- それでも FASTQ が無ければ失敗として記録して次へ ---
    if [ "${fastq_ready}" -ne 1 ]; then
        echo "${sample_name} failed"
        # tmp_list 由来の情報を Download_failed_list.txt に追記
        printf '%s\t%s\t%s\t%s\t%s\n' "${srr}" "${samn}" "${bio}" "${old_name}" "${lib}" >> Download_failed_list.txt
        rm -f SAMN_info_tmp "${sra}.sra" *.fastq
        continue
    fi
    ########################################
    # FASTQ サイズチェック（100 MB 未満ならスキップ）
    ########################################
    # 100 MB = 100 * 1024 * 1024 bytes
    threshold=$((100 * 1024 * 1024))
    total_size=0

    if [ "${lib}" = "PAIRED" ]; then
        raw1="${code}_raw_R1.fastq"
        raw2="${code}_raw_R2.fastq"
        size1=0
        size2=0
        [ -f "${raw1}" ] && size1=$(wc -c < "${raw1}")
        [ -f "${raw2}" ] && size2=$(wc -c < "${raw2}")
        total_size=$(( (size1 + size2) / 2 ))
    else
        raw="${code}_raw_R.fastq"
        [ -f "${raw}" ] && total_size=$(wc -c < "${raw}") || total_size=0
    fi

    echo "  Total raw FASTQ size for ${sample_name}: ${total_size} bytes"

    if [ "${total_size}" -lt "${threshold}" ]; then
        echo "  [SKIP] ${sample_name}: total FASTQ size < 100 MB"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${srr}" "${samn}" "${bio}" "${sample_name}" "${lib}" "${total_size}" \
            >> samples_too_small_size.txt

        # 小さいサンプルはクリーンアップして次へ
        rm -f "${code}_raw_R1.fastq" "${code}_raw_R2.fastq" "${code}_raw_R.fastq"
        rm -f SAMN_info_tmp "${sra}.sra"
        continue
    fi

    ############################
    # fastp （QC後のFASTQ作成）
    ############################
    if [ "${lib}" = "PAIRED" ]; then
        fastp -i "${code}_raw_R1.fastq" -I "${code}_raw_R2.fastq" \
              -o "${code}_R1.fastq" -O "${code}_R2.fastq"
    else
        fastp -i "${code}_raw_R.fastq" -o "${code}_R.fastq"
    fi

    ############################
    # STAR でマッピング
    ############################
    ulimit -n 100000

    #STARが落ちてもループが止まらないように、一次的に-eをoff
    set +e
    if [ "${lib}" = "PAIRED" ]; then
        STAR --runThreadN 24 \
             --genomeDir "${ref}" \
             --readFilesIn "${code}_R1.fastq" "${code}_R2.fastq" \
             --outFileNamePrefix "${code}_" \
             --outSAMtype BAM SortedByCoordinate
    else
        STAR --runThreadN 24 \
             --genomeDir "${ref}" \
             --readFilesIn "${code}_R.fastq" \
             --outFileNamePrefix "${code}_" \
             --outSAMtype BAM SortedByCoordinate
    fi
    star_status=$?
    #ここで-eを復活
    set -e

        # STAR が失敗した場合は、このサンプルをスキップ
    if [ "${star_status}" -ne 0 ]; then
        echo "  [ERROR] STAR mapping failed for ${sample_name} (${srr}). Skip this sample."

        # 失敗サンプルの一覧をダウンロード失敗リストに入れる
        printf '%s\t%s\t%s\t%s\t%s\tSTAR_failed\n' \
            "${srr}" "${samn}" "${bio}" "${old_name}" "${lib}" \
            >> Download_failed_list.txt

        # そのサンプル由来の一時ファイルを掃除（必要に応じて調整）
        rm -f "${code}_R1.fastq" "${code}_R2.fastq" "${code}_R.fastq"
        rm -f "${code}_Aligned.out.bam" "${code}_Aligned.sortedByCoord.out.bam"
        rm -f SAMN_info_tmp
        rm -rf *tmp *out

        # 次のサンプルへ
        continue
    fi

    mv "${code}_Aligned.sortedByCoord.out.bam" "${code}.bam"
    samtools index "${code}.bam"

    ############################
    # featureCounts → カウント表
    ############################
    if [ "${lib}" = "PAIRED" ]; then
        featureCounts -p -T 20 -t gene -g ID -a "${gff}" \
            -o "${code}_counts.txt" "${code}.bam"
    else
        featureCounts -T 20 -t gene -g ID -a "${gff}" \
            -o "${code}_counts.txt" "${code}.bam"
    fi

    # gene length と counts (featureCounts 出力から切り出し)
    cut -f1,6 "${code}_counts.txt" | grep -v "#" > gene_length.tsv
    cut -f7     "${code}_counts.txt" | grep -v "#" > "${code}_counts_only"
    cut -f1     "${code}_counts.txt" | grep -v "#" > gene_id_tmp
    paste gene_id_tmp "${code}_counts_only" > "${code}_counts_for_tpm.txt"

    ############################
    # TPM 計算 (Count_to_TPM.py)
    ############################
    python Count_to_TPM.py "${code}_counts_for_tpm.txt" gene_length.tsv
    tpm_file="${code}_counts_for_tpm.txt_with_TPM"   # Geneid + {code}.bamTPM

    ###########################################
    # TPM マトリクス (行 = サンプル) の構築
    ###########################################
    if [ ! -f "${out}_TPM_matrix.csv" ]; then
        echo "Generate TPM matrix ${out}_TPM_matrix.csv..."
        # Geneid 列から遺伝子名をヘッダに (1行目 "Geneid" は除く)
        genes=$(cut -f1 "${tpm_file}" | sed '1d' | paste -sd',' -)
        echo "BioProject,SRA,Biosample,treatment,tissue,stage,cultivar,code,temperature,attributes,${genes}" \
            > "${out}_TPM_matrix.csv"
    fi

    # このサンプルのTPM値（2列目）を1行のカンマ区切りに
    tpm_values=$(cut -f2 "${tpm_file}" | sed '1d' | paste -sd',' -)

    # メタ情報の空白はアンダースコアに統一（コード列はすでに code で整形済み）
    treatment_s=$(printf '%s' "${treatment}"   | sed 's/ /_/g' | sed 's/,/_/g')
    tissue_s=$(printf '%s'    "${tissue}"      | sed 's/ /_/g' | sed 's/,/_/g')
    stage_s=$(printf '%s'     "${stage}"       | sed 's/ /_/g' | sed 's/,/_/g')
    line_s=$(printf '%s'      "${cultivar}"    | sed 's/ /_/g' | sed 's/,/_/g')
    temp_s=$(printf '%s'      "${temperature}" | sed 's/ /_/g' | sed 's/,/_/g')
    attr_s=$(printf '%s'      "${attributes}"  | sed 's/ /_/g' | sed 's/,/_/g')

    # 1 サンプル分を1行として書き込み
    # 1列目 PRJNA (bio), 8列目 attributes, 9列目以降が遺伝子TPM
    echo "${bio},${srr},${samn},${treatment_s},${tissue_s},${stage_s},${line_s},${code},${temp_s},${attr_s},${tpm_values}" \
        >> "${out}_TPM_matrix.csv"

    ###########################################
    # Count マトリクス (行 = サンプル) の構築（TPMと同構造）
    ###########################################
    if [ ! -f "${out}_count_matrix.csv" ]; then
        echo "Generating Count matrix ${out}_count_matrix.csv..."
        # Geneid をヘッダに（TPMと同じ並び）
        genes_count=$(cut -f1 "${code}_counts_for_tpm.txt" | sed '1d' | paste -sd',' -)
        echo "PRJNA,SRR,SAMN,treatment,tissue,stage,line,code,temperature,attributes,${genes_count}" \
            > "${out}_count_matrix.csv"
    fi

    # counts（2列目）を1行に
    count_values=$(cut -f2 "${code}_counts_for_tpm.txt" | sed '1d' | paste -sd',' -)

    echo "${bio},${srr},${samn},${treatment_s},${tissue_s},${stage_s},${line_s},${code},${temp_s},${attr_s},${count_values}" \
        >> "${out}_count_matrix.csv"


    ############################
    # 一時ファイルの掃除
    ############################
    set +e
    rm *ra*.fastq
    rm -f *.bam *.bam.bai
    rm -f *counts.txt
    rm -f gene_length.tsv gene_id_tmp "${code}_counts_only" "${code}_counts_for_tpm.txt" "${tpm_file}"
    mv "${sra}.sra" sra
    mv ${sra} sra
    rm -f SAMN_info_tmp
    rm -rf *tmp
    rm -rf *out
    set -e

done < diff_tmp

echo "All samples finished."
echo "TPM matrix (sample x gene): ${out}_TPM_matrix.csv"
echo "Metadata file: ${meta_file}"

