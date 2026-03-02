# RNA-seq → Interactive Web Database Pipeline

This repository contains a reproducible pipeline that:

1. searches NCBI SRA for RNA-seq/transcriptome runs,
2. downloads FASTQs, maps reads with STAR, counts with featureCounts,
3. builds TPM/count matrices,
4. standardizes metadata and filters obviously bad samples,
5. annotates genes via BLASTP against Arabidopsis and rice,
6. generates an interactive HTML page from a template, and
7. preprocesses a large matrix into sharded binary files for fast web loading.

## What you need

### Required inputs

- **Organism name(s)** for NCBI SRA search (repeatable `--species`).
- **STAR genome index** (built in advance, `--star-index`).
- **GFF3 annotation** (`--gff`) where gene features have `ID=...` (used by featureCounts `-g ID`).
- **Target species protein FASTA** (`--proteins`) whose headers match the gene IDs you want to annotate.
- **HTML template** (`scripts/Template.html`) and parameters to customize it.

### Software dependencies

The mapping step downloads and processes RNA-seq data, so you need common bioinformatics tools:

- `aria2c`, `wget`
- `fasterq-dump` (SRA-Toolkit)
- `fastp`
- `STAR`
- `samtools`
- `featureCounts` (Subread)
- `blastp` (BLAST+)

Python 3 + modules:

- `pandas`
- `bs4` (BeautifulSoup)

> Tip: consider using conda/mamba and recording exact versions for reproducibility.

## Repository layout

```
.
├── run_pipeline.sh              # main entrypoint
├── scripts/                     # scripts
├── data/                        # Arabidopsis dataset
└── docs/wiki/                   # GitHub Wiki draft pages (copy/paste into wiki)
```

## Quick start

### 
0) Prepare reference files

1) Build STAR index (example):

```bash
STAR --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir /path/to/star_index \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile genes.gff3 \
  --sjdbOverhang 100
```

2) Prepare BLAST databases for Arabidopsis/rice.

- Arabidopsis datasets: TAIR10  
  Data were obtained from TAIR (Phoenix Bioinformatics) Public_Data_Releases and are licensed under CC BY 4.0.  
  We have reformatted files for this pipeline.  
  Source: TAIR Public_Data_Releases (download date: 2024-10-01) (Berardini et al., 2015,  https://doi.org/10.1002/dvg.22877).  
  License: CC BY 4.0.  
  
- Rice datasets: RGAP 7, from the Rice Genome Annotation Project (RGAP, Kawahara et al., 2013) (https://doi.org/10.1186/1939-8433-6-4)  
  and Oryzabase (Kurata and Yamazaki, 2006, https://doi.org/10.1104/pp.105.063008).  
  For the rice datasets, please download files from https://rice.uga.edu/download_osa1r7.shtml /https://shigen.nig.ac.jp/rice/oryzabase/download/gene  
  and create the protein BLAST DBs and annotation files using the following command:  

```bash
# Recommended: build rice (RGAP7 + Oryzabase) reference files automatically.
# This will download the required files and create:
#   - data/Os_proteins.fasta (+ BLAST DB files)
#   - data/Os_gene_annotation_list.tsv
bash scripts/build_os_data.sh --data-dir data --threads 8
```

`blastp_annot.sh` expects the following files under `./data/` (repository root):
- `data/Arab_proteins.fasta` and its BLAST DB files (`.pin/.psq/.phr` etc)
- `data/Os_proteins.fasta` and its BLAST DB files
- `data/Arab_gene_annotation_list.tsv`
- `data/Os_gene_annotation_list.tsv`

### 1) Run the pipeline

```bash
./run_pipeline.sh \
  --species "Cenchrus americanus" \
  --species "Pearl millet" \
  --prefix PM_Tift \
  --star-index /path/to/star_index \
  --gff /path/to/genes.gff3 \
  --proteins /path/to/target_species_proteins.fa \
  --title-prefix "Pearl Millet" \
  --h1-prefix "Pearl Millet" \
  --ref-label "Tift" \
  --outdir out/PM_Tift
```

Outputs:

- `out/PM_Tift/results/` – CSV matrices + metadata TSV
- `out/PM_Tift/site/` – **upload this directory to your web server**
  - `PM_Tift.html`
  - `PM_Tift_annot.txt`
  - `meta.tsv`, `gene_index.tsv`, `manifest.json`, `shards/`


## Customizing metadata normalization

To fix spelling/formatting variation, you can pass mapping rules to `standardize_metadata.py`:

- TSV mapping file (repeatable): `--map-tsv replacements.tsv`
  Example `replacements.tsv` (TAB-separated; lines starting with `#` are ignored):

```tab-delimited file 
# column	from	to
tissue	Leaves	Leaf
tissue	ROOTS	Root
treatment	25.0 degrees	25 °C
treatment	42 Degree	42 °C
```
- Inline rules (repeatable): `--replace 'tissue:Leaves=leaf'`

## Partial runs / resume

- Start from a later step:

```bash
./run_pipeline.sh ... --from clean
```

- Stop early:

```bash
./run_pipeline.sh ... --to html
```

- Use your own SRR list (skip NCBI search):

```bash
./run_pipeline.sh ... --srr-list path/to/tmp_list
```

- Test with only first N samples:

```bash
./run_pipeline.sh ... --max-samples 20
```

## Notes / caveats

- **Large downloads:** SRA datasets can be huge. Always check disk quota.
- **NCBI rate limits:** heavy use of eutils may be throttled.
- The mapping and BLAST scripts currently use fixed thread counts internally.
  If you want to tune performance, edit:
  - `scripts/mapping_script.sh` (STAR/featureCounts/fasterq-dump)
  - `scripts/blastp_annot.sh` (blastp)

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
