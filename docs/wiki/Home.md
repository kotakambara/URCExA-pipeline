# RNA-seq DB pipeline (Wiki draft)

This Wiki is a draft for GitHub Wiki. Copy each page into your GitHub Wiki repository.

## What this pipeline does

This pipeline turns public RNA-seq SRA runs into an interactive web database:

1) Search SRA by organism name(s)
2) Download FASTQ and map reads (STAR)
3) Count reads per gene (featureCounts) and compute TPM
4) Clean metadata spellings and filter obvious low-quality samples
5) Add gene annotations via BLASTP (Arabidopsis + rice)
6) Generate an HTML front-end from a template
7) Shard the large TPM matrix for fast web loading

## Typical use cases

- Create a web portal for a species/cultivar.
- Rebuild the database periodically with new SRA runs.
- Generate standardized matrices for downstream analysis.

## Minimum required inputs

- `--species`: organism/scientific names (repeatable) used for NCBI SRA query
- `--star-index`: STAR genomeDir (pre-built)
- `--gff`: GFF3 annotation (gene features must have `ID=`)
- `--proteins`: protein FASTA (headers must match gene IDs)

## Quick start (copy/paste)

```bash
./run_pipeline.sh \
  --species "Setaria italica" \
  --species "foxtail millet" \
  --prefix SI_Yugu1 \
  --star-index /path/to/star_index \
  --gff /path/to/genes.gff3 \
  --proteins /path/to/target_species_proteins.fa \
  --title-prefix "Foxtail Millet" \
  --h1-prefix "Foxtail Millet" \
  --ref-label "Yugu1" \
  --outdir out/SI_Yugu1
```

After completion, upload:

- `out/SI_Yugu1/site/`

to your web server (static hosting is enough).

## Outputs

### `results/`

- `<prefix>_TPM_matrix.csv` – raw TPM matrix (sample x gene)
- `<prefix>_TPM_matrix.cleaned.csv` – metadata standardized
- `<prefix>_TPM_matrix.cleaned.filtered.csv` – low-quality samples removed
- `<prefix>_count_matrix.csv` – count matrix (optional)
- `<prefix>_metadata.tsv` – metadata log (from BioSample/SRA)

### `site/`

- `<prefix>.html` – interactive UI
- `<prefix>_annot.txt` – gene annotations for UI
- `meta.tsv`, `gene_index.tsv`, `manifest.json`, `shards/` – sharded matrix files

## Common pitfalls

### 1) GFF3 `ID=` mismatch

`featureCounts` is run with `-t gene -g ID`. Your GFF3 gene features must look like:

```
chr1  ...  gene  ...  ID=GENE12345;...
```

If your gene IDs are stored under a different key (e.g., `gene_id`), edit `mapping_script.sh` accordingly.

### 2) BLAST resources missing

`blastp_annot.sh` expects `./data/Arab_proteins.fasta`, `./data/Os_proteins.fasta`, and annotation TSVs.
Create protein BLAST DBs with `makeblastdb`.

### 3) Download/mapping failures

Some SRR runs may fail to download or map.
`mapping_script.sh` records failures to `Download_failed_list.txt`.

## Recommended publication practice

- Keep your metadata normalization rules as a TSV mapping file and version-control it.
- Record software versions and reference genome/annotation versions.
- Avoid committing large reference files unless licenses allow it.
