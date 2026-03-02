input=${1:?} #input.fa
out=${2:?} #output prefix e.g. Pg/Os/At

blastp -db data/Arab_proteins.fasta -evalue 1e-20 -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send pident bitscore mismatch gapopen evalue" -query ${input} -num_threads 12 -max_target_seqs 5 > ${out}_proteins_At_proteins_blastp_out.txt
blastp -db data/Os_proteins.fasta -evalue 1e-20 -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send pident bitscore mismatch gapopen evalue" -query ${input} -max_target_seqs 5 -num_threads 12 > ${out}_proteins_Os_proteins_blastp_out.txt

python blast_append_annotation.py --blast ${out}_proteins_At_proteins_blastp_out.txt --anno data/Arab_gene_annotation_list.tsv --out ${out}_At_proteins_blastp.txt
python blast_append_annotation.py --blast ${out}_proteins_Os_proteins_blastp_out.txt --anno data/Os_gene_annotation_list.tsv --out ${out}_Os_proteins_blastp.txt

cut -f1,16 ${out}_At_proteins_blastp.txt > ${out}_At_annot.txt; cut -f1,16 ${out}_Os_proteins_blastp.txt > ${out}_Os_annot.txt
./merge_annot.sh --input anything --arab ${out}_At_annot.txt --rice ${out}_Os_annot.txt --output tmp
./drop_dot.sh tmp ${out}_annot.txt
rm tmp
rm *blastp*.txt