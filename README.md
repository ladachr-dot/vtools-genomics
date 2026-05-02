# 🧬 vtools-genomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

A collection of Python utilities designed to streamline and automate routine genomic data processing workflows. 

This toolkit bridges the gap between raw genomic data and analytical tools like PLINK, handling everything from missing rsID recovery via NCBI APIs, safe genome-build liftover, and synthetic dataset generation for testing.

## 🛠️ The Utilities

This repository is structured as a modular package containing several dedicated tools:

### ✅ 1. PGS → PLINK Converter & Calculator (`pgs_to_plink`)

The repository is structured as a modular package containing several dedicated tools.

- **NCBI integration:** Recovers missing rsIDs based on `CHR:POS` using NCBI APIs.
- **VCF handling:** Converts patient `.vcf/.vcf.gz` files into PLINK binary format (`.bed/.bim/.fam`).
- **Automated liftover:** Safely converts genome builds (e.g. hg38 ↔ hg19) using temporary internal IDs to avoid PLINK exclusion artefacts.
- **Risk calculation:** Runs PLINK under the hood to compute final risk scores for each sample.
- **Interfaces:** core logic in Python, with a CLI entry point and a simple GUI wrapper (work in progress).

### 🚧 2. VCF Parser (`parser`) — *work in progress*

A multiprocessing‑enabled parser for heavy VCF files, aimed at extracting and filtering specific genomic regions or variant types without loading entire files into RAM.

Planned features:

- filtering by gene, region, or variant type;
- per‑chromosome parallel processing;
- export to tidy CSV/TSV or BED‑like intervals.

### ✅ 3. NCBI Batch Downloader (`ncbi_downloader`)

A resilient batch downloader for NCBI Entrez resources, designed to handle large ID lists without breaking on rate limits or network hiccups.

Current features:

- batch downloads by accession IDs or rsIDs (from a plain text file);
- automatic retry/backoff on NCBI rate limits and transient network errors;
- JSON-based checkpoints to resume interrupted runs without starting from scratch;
- logs progress and failures for later inspection;
- core implemented as a reusable Python module, with a simple CLI and GUI wrapper on top.

Typical usage patterns:

- recovering missing rsIDs or annotations for GWAS/PGS workflows;
- downloading FASTA/GenBank records for a list of accessions;
- building small local datasets for downstream analysis.

### 🚧 4. Synthetic Data Generator (`generator`) — *work in progress*

A factory‑pattern module for generating synthetic genomics datasets for development and testing.

Planned features:

- synthetic VCF files with configurable variant density and error rates;
- small `.bed/.bim/.fam` trios for PLINK pipelines;
- pre‑packaged toy datasets for unit tests and examples.

## 🏗️ Project Architecture

The project follows the standard Python `src/` layout best practices for modularity and easy CLI integration:

```text
vtools-genomics/
├── src/
│   └── vtools/
│       ├── cli.py             # Main CLI entry point
│       ├── pgstoplink.py      # PGS Converter module
│       ├── parser.py          # VCF Parser module
│       ├── downloader.py      # Entrez Downloader module
│       ├── generator.py       # Synthetic Data module
│       └── utils/             # Shared helpers (loggers, validators)
├── data/                      # Dummy/test datasets
├── tests/                     # Unit tests
├── requirements.txt           # Dependencies
└── README.md
```

---

## 🚀 Installation

1. Clone the repository:
```bash
git clone https://github.com/YourUsername/vtools-genomics.git
cd vtools-genomics
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

*(Optional but recommended)* To use the full PGS → PLINK pipeline, make sure [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) is installed and available in your `PATH`.

## 💻 Quick start

For now, the PGS → PLINK converter can be launched directly as a module:

```bash
python -m vtools.pgstoplink
```

or

```bash
python src/vtools/pgstoplink.py
```

GUI and unified `vtools` CLI entrypoints via `cli.py` are under active development and will be documented in upcoming releases.


## 🧪 Tests

The project uses `pytest` for basic regression tests.

```bash
pip install -e .[dev]
pytest
```

Current tests cover:

- synthetic VCF generation on small toy datasets;
- basic PGS → PLINK conversion flows.

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
