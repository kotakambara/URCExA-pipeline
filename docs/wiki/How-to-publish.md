# How to publish the web database

## 1) Build the site folder

Run the pipeline. The final site folder is:

- `out/<prefix>/site/`

It contains:

- `<prefix>.html`
- `<prefix>_annot.txt`
- `meta.tsv`, `gene_index.tsv`, `manifest.json`
- `shards/` (binary shard files)

## 2) Upload to your server

Upload the entire `site/` directory to a static hosting location.

Examples:

- GitHub Pages (for small/medium datasets)
- Nginx/Apache static directory
- Any object storage with static hosting

## 3) Confirm in a browser

Open:

- `<prefix>.html`

If you see network errors, confirm that the relative files exist next to the HTML.

## 4) Optional: custom domain + caching

Large shard downloads benefit from caching.

- Enable gzip/brotli for `meta.tsv`, `gene_index.tsv`, `manifest.json`
- Configure long cache headers for `shards/*.bin`
