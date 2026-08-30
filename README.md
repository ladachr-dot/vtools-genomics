# 🧬 vtools-genomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

A collection of Python utilities designed to streamline and automate routine genomic data processing workflows.

This toolkit bridges the gap between raw genomic data and analytical tools like PLINK, handling missing rsID recovery via NCBI APIs and safe genome-build liftover.

## 🛠️ The Utilities

This repository is structured as a modular package containing several dedicated tools:

### ✅ 1. PGS → PLINK Converter (`vtools.pgs_to_plink`)

- **NCBI integration:** recovers missing rsIDs based on `CHR:POS` using NCBI Entrez/Variation APIs.
- **PGS score file handling:** parses PGS Catalog-style score files (`.txt/.tsv/.csv`, tab- or whitespace-separated) and normalizes them into a PLINK-ready `SNP / A1 / BETA` format.
- **Automated liftover:** converts genome-build coordinates (hg19 ↔ hg38) via `pyliftover`, with resumable checkpoints for long-running rsID recovery batches.
- **Interfaces:** CLI entry point (`argparse`), with both fully-argument-driven and interactive prompt modes.

> ⏳ **Planned, not implemented yet:** direct `.vcf/.vcf.gz` parsing, a `--bfile`/PLINK-execution wrapper that runs the scoring step automatically, and a GUI wrapper. Currently the tool prepares the score/`.bim` files; running PLINK itself is a manual step (see the printed command hint).

### ✅ 2. NCBI Batch Downloader (`vtools.ncbi_downloader`)

A resilient batch downloader for NCBI Entrez `esummary` lookups by rsID, designed to handle large ID lists without breaking on rate limits or network hiccups.

Current features:

- batch enrichment of a TSV/CSV/XLSX table by extracting rsIDs and querying NCBI `esummary` for chromosome/position/gene metadata;
- automatic retry/backoff on NCBI rate limits and transient network errors (via `urllib3.Retry`);
- JSON-based checkpoints to resume interrupted runs without starting from scratch;
- incremental Excel output, saved every N processed IDs so partial progress is never lost;
- core implemented as a reusable Python module with a CLI entry point (interactive prompts when run without arguments).

> ⏳ **Planned, not implemented yet:** downloading raw FASTA/GenBank records by accession ID, and a GUI wrapper. The current implementation covers rsID metadata enrichment only.

### ✅ 3. ISOGG ↔ YFull Haplogroup Converter (`vtools.haplo_converter`)

A bidirectional converter and SNP resolver for Y-DNA haplogroup nomenclature, bridging the ISOGG alphanumeric system and the YFull SNP-based notation.

Current features:

- SNP index built from `current_tree.json` (downloaded from the public YFull tree), cached in SQLite so it isn't re-parsed on every run;
- case-insensitive ISOGG ↔ YFull lookup with automatic fallback to the nearest parent haplogroup when an exact label isn't found;
- auto-detects input notation: ISOGG (`R1b1a1`) vs. YFull (`R-M269`);
- three-level search: exact → `LIKE` → fuzzy `difflib`, with haplogroup queries falling back through parent nodes;
- imports the ISOGG SNP Index CSV into the same SQLite cache for combined lookups;
- console CLI with a menu-driven workflow (load tree, load CSV, convert, search, clear DB).

> ⏳ **Planned, not implemented yet:** batch Excel-to-Excel conversion with auto-sized columns/frozen header, per-row audit status codes (`matched_by_name_snp`, `ambiguous_multiple_yfull_matches`, …), and a desktop GUI with a progress bar. The current tool is CLI-only and operates interactively, one query at a time.

### 🚧 *work in progress*

A unified `vtools` CLI (`cli.py`) that exposes all three tools as subcommands, and GUI wrappers for each tool, are planned but not yet implemented.

## 🏗️ Project Architecture

The project follows the standard Python `src/` layout:

```text
vtools-genomics/
├── src/
│   └── vtools/
│       ├── pgs_to_plink/
│       │   └── pgstoplink.py      # PGS -> PLINK converter (CLI)
│       ├── ncbi_downloader/
│       │   └── downloader.py      # NCBI Entrez batch downloader (CLI)
│       ├── haplo_converter/
│       │   └── converter.py       # ISOGG <-> YFull haplogroup converter (CLI)
│       └── utils/                 # Shared helpers (planned: logging, checkpoint, HTTP session)
├── data/                          # Dummy/test datasets
├── tests/                         # Unit tests (pytest)
├── requirements.txt               # Dependencies
└── README.md
```

A unified `cli.py` entry point (`vtools pgs-convert`, `vtools ncbi-download`, `vtools haplo-convert`) is planned; for now each tool is run as its own module (see Quick start below).

---

## 🚀 Installation

1. Clone the repository:
```bash
git clone https://github.com/ladachr-dot/vtools-genomics.git
cd vtools-genomics
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

*(Optional but recommended)* To use the full PGS → PLINK pipeline, make sure [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) is installed and available in your `PATH`.

## 💻 Quick start

Each tool currently runs as its own script:

```bash
python src/vtools/pgs_to_plink/pgstoplink.py score.txt
python src/vtools/ncbi_downloader/downloader.py variants.tsv
python src/vtools/haplo_converter/converter.py
```

A unified `python -m vtools.cli` entry point is under active development and will be documented once available.

## 🧪 Tests

The project uses `pytest` for unit tests covering the pure logic (no network, no PLINK, no SQLite-tree download required):

```bash
pip install -r requirements.txt
pytest tests/
```

Current tests cover:

- `pgs_to_plink`: column-name normalization/matching, genome-build header parsing;
- `ncbi_downloader`: rsID extraction, checkpoint save/resume;
- `haplo_converter`: ISOGG ↔ YFull lookup (exact and parent-fallback), notation auto-detection.

Integration-level tests (full VCF/PGS conversion pipelines, live NCBI/YFull calls) are not implemented yet — see the "planned" notes under each tool above.

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
