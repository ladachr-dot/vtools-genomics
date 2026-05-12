# 🧬 vtools-genomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

A collection of Python utilities designed to streamline and automate routine genomic data processing workflows. 

This toolkit bridges the gap between raw genomic data and analytical tools like PLINK, handling everything from missing rsID recovery via NCBI APIs, safe genome-build liftover.

## 🛠️ The Utilities

This repository is structured as a modular package containing several dedicated tools:

### ✅ 1. PGS → PLINK Converter & Calculator (`pgs_to_plink`)

The repository is structured as a modular package containing several dedicated tools.

- **NCBI integration:** Recovers missing rsIDs based on `CHR:POS` using NCBI APIs.
- **VCF handling:** Converts patient `.vcf/.vcf.gz` files into PLINK binary format (`.bed/.bim/.fam`).
- **Automated liftover:** Safely converts genome builds (e.g. hg38 ↔ hg19) using temporary internal IDs to avoid PLINK exclusion artefacts.
- **Risk calculation:** Runs PLINK under the hood to compute final risk scores for each sample.
- **Interfaces:** core logic in Python, with a CLI entry point and a simple GUI wrapper (work in progress).

### ✅ 2. NCBI Batch Downloader (`ncbi_downloader`)

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

### ✅ 3. ISOGG ↔ YFull Haplogroup Converter (`phylo_resolve`)

A bidirectional converter and SNP resolver for Y-DNA haplogroup nomenclature, bridging the ISOGG alphanumeric system and the YFull SNP-based notation.

Current features:

- live SNP index built from `current_tree.json`, cached in SQLite — no re-parsing on subsequent runs;
- auto-detects input notation: ISOGG (`R1b1a1`) vs. YFull (`R-M269`);
- three-level search: exact → `LIKE` → fuzzy difflib, with haplogroup queries falling back through parent nodes;
- converts ISOGG SNP Index workbooks to Excel with added YFull columns, auto-sized columns, and frozen header;
- detailed per-row status codes for auditability (`matched_by_name_snp`, `ambiguous_multiple_yfull_matches`, …);
- desktop GUI with progress bar + CLI entry point.

Typical usage patterns:

- annotating cohorts with standardised YFull labels for population stratification;
- cross-referencing ISOGG SNP Index entries against the live YFull tree to spot reclassified haplogroups;
- batch-converting legacy haplogroup columns from ISOGG to YFull notation (or vice versa).

### 🚧 *work in progress*

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
│       ├── ... there may be other tools ...
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
