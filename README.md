# 🧬 vtools-genomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

A collection of high-performance Python utilities designed to streamline and automate routine genomic data processing workflows. 

This toolkit bridges the gap between raw genomic data and analytical tools like PLINK, handling everything from missing rsID recovery via NCBI APIs to safe genome build liftovers.

## 🛠️ The Utilities

This repository is structured as a modular package containing several dedicated tools:

### ✅ 1. PGS to PLINK Converter & Calculator (`pgstoplink`)
An interactive CLI module for processing Polygenic Risk Score (PGS) files.
- **Multithreaded NCBI Fetching:** Automatically retrieves missing rsIDs based on `chr:pos` coordinates.
- **VCF Integration:** Seamlessly converts patient `.vcf/.vcf.gz` files to PLINK binary format (`.bed/.bim/.fam`).
- **Automated Liftover:** Safely converts genome builds (e.g., hg38 ↔ hg19) using temporary internal IDs to avoid PLINK exclusion bugs.
- **Auto-Calculation:** Computes the final patient risk profile directly via PLINK subprocesses.

### 🚧 2. VCF Parser (`parser`) - *Work in Progress*
A fast, multiprocessing-enabled parser for heavy VCF files, designed to extract, filter, and format specific genomic regions without overloading RAM.

### 🚧 3. Entrez Batch Downloader (`downloader`) - *Work in Progress*
A robust downloading module with built-in resume logic and error handling for fetching massive datasets from NCBI Entrez.

### 🚧 4. Synthetic Data Generator (`generator`) - *Work in Progress*
A factory-pattern-based module to generate synthetic `.bed/.bim/.fam` and `.vcf` datasets for testing and benchmarking bioinformatics pipelines.

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

*(Optional but recommended)* Ensure you have [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) installed and added to your system's `PATH` to utilize the full calculation features of the `pgstoplink` module.

## 💻 Quick Start

To run the interactive PGS to PLINK converter:
```bash
python src/vtools/pgstoplink.py
```
*(Note: Full integrated CLI usage via `cli.py` will be available in future releases).*

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
