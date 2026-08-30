# 🧬 vtools-genomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

A collection of Python utilities for genomic data processing workflows: preparing PGS score files for PLINK, recovering rsIDs and metadata from NCBI, and converting between ISOGG and YFull Y-DNA haplogroup nomenclature.

All three tools are available both as a unified `vtools` command-line interface and as standalone modules.

## 🛠️ The Utilities

### 1. PGS → PLINK Converter (`vtools pgstoplink`)

- Parses PGS Catalog-style score files (`.txt/.tsv/.csv`, tab- or whitespace-separated) and normalizes them into a PLINK-ready `SNP / A1 / BETA` format.
- Recovers missing rsIDs based on `CHR:POS` via the NCBI Variation/Entrez APIs, with resumable checkpoints for long-running batches.
- Converts genome-build coordinates (hg19 ↔ hg38) via `pyliftover`.
- Optionally generates a matching `.bim` file and prints a ready-to-run PLINK scoring command.

### 2. NCBI Batch Downloader (`vtools downloader`)

- Extracts rsIDs from a TSV/CSV/XLSX table and enriches them with NCBI `esummary` metadata (chromosome, position, gene).
- Automatic retry/backoff on NCBI rate limits and transient network errors (via `urllib3.Retry`).
- JSON-based checkpoints to resume interrupted runs without starting from scratch.
- Incremental Excel output, saved every N processed IDs so partial progress is never lost.

### 3. ISOGG ↔ YFull Haplogroup Converter (`vtools converter`)

- SNP index built from the public YFull `current_tree.json`, cached in SQLite so it isn't re-parsed on every run.
- Case-insensitive ISOGG ↔ YFull lookup with automatic fallback to the nearest parent haplogroup when an exact label isn't found.
- Auto-detects input notation: ISOGG (`R1b1a1`) vs. YFull (`R-M269`).
- Three-level search: exact → `LIKE` → fuzzy `difflib`, with haplogroup queries falling back through parent nodes.
- Imports the ISOGG SNP Index CSV into the same SQLite cache for combined lookups.
- Interactive console menu: load the YFull tree, load JSON/CSV, convert a file, search, or clear the database.

## 🏗️ Project Architecture

```text
vtools-genomics/
├── src/
│   └── vtools/
│       ├── cli.py                     # Unified CLI entry point (click)
│       ├── pgs_to_plink/
│       │   └── pgstoplink.py          # PGS -> PLINK converter
│       ├── ncbi_downloader/
│       │   └── downloader.py          # NCBI Entrez batch downloader
│       └── haplo_converter/
│           └── converter.py           # ISOGG <-> YFull haplogroup converter
├── data/                               # Dummy/test datasets
├── tests/                              # Unit tests (pytest)
├── pyproject.toml                      # Packaging metadata, `vtools` console script
├── requirements.txt                    # Dependencies
└── README.md
```

---

## 🚀 Installation

1. Clone the repository:
```bash
git clone https://github.com/ladachr-dot/vtools-genomics.git
cd vtools-genomics
```

2. Install the package (this registers the `vtools` command):
```bash
pip install -e .
```

*(Optional but recommended)* To run the PLINK scoring step after `pgs-convert`, install [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) and make sure it is available in your `PATH`.

## 💻 Quick start

```bash
vtools pgs-convert score.txt --bim --plink-hint
vtools ncbi-download variants.tsv --email you@example.com
vtools haplo-convert
```

Run `vtools --help` or `vtools <command> --help` for the full list of options for each subcommand.

## 🧪 Tests

```bash
pip install -e .[dev]
pytest tests/
```

Current tests cover:

- `pgs_to_plink`: column-name normalization/matching, genome-build header parsing;
- `ncbi_downloader`: rsID extraction, checkpoint save/resume;
- `haplo_converter`: ISOGG ↔ YFull lookup (exact and parent-fallback), notation auto-detection.

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
