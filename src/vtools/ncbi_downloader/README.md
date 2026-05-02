## NCBI / GWAS Batch Downloader

This module provides a resilient way to download records from NCBI (e.g. GWAS Catalog, NCBI Entrez) in bulk, with checkpointing and retry logic.

### Key features

- Batch downloads of sequences or annotations by a list of accession IDs or rsIDs
- Automatic retry and backoff when NCBI rate limits or network errors occur
- Checkpoints: if the run is interrupted, the tool resumes from the last processed ID
- Unified CLI and (optional) GUI interface for non-technical users

### Command-line usage

```bash
# Basic usage: download GenBank records by accession IDs
vtools ncbi-download \
  --input ids.txt \
  --db nucleotide \
  --outdir data/ncbi_records

# Resume a previously interrupted run using a checkpoint file
vtools ncbi-download \
  --input ids.txt \
  --db nucleotide \
  --outdir data/ncbi_records \
  --checkpoint checkpoints/ncbi_checkpoint.json \
  --resume
```

Typical arguments:

- `--input` – path to a text file with one ID per line
- `--db` – NCBI database name (e.g. `nucleotide`, `protein`, `snp`)
- `--outdir` – output directory for downloaded records
- `--email` / `--api-key` – optional NCBI contact email and API key
- `--batch-size` – number of IDs per request
- `--checkpoint` – path to a JSON checkpoint file
- `--resume` – resume from the existing checkpoint instead of starting from scratch

### GUI usage

If you prefer a graphical interface, you can launch the GUI:

```bash
vtools ncbi-gui
```

The GUI allows you to:

- select the input file and output directory
- choose the target NCBI database
- configure batch size and checkpoint file
- monitor progress and errors in a simple status panel

This is especially convenient for biologists who do not want to work with the command line.

### Notes and limitations

- Please respect NCBI usage policies and rate limits.
- For very large ID lists, consider running during off-peak hours and using a valid `--email` and `--api-key`.
- The tool is intended for small to medium scale workflows (e.g. thousands–tens of thousands of IDs), not for mirroring entire databases.