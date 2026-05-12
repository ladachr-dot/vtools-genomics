# 🧬 PGS → PLINK Converter & Calculator

This project provides a simple GUI wrapper for `pgs_to_plink.py` so the converter can be run from a window instead of a terminal, while keeping the original workflow and core functionality of the script intact. The underlying converter prepares PGS scoring files, searches for rsIDs, creates `.txt` and `.bim` outputs, and can run additional PLINK-based steps when PLINK is available.

## Files in the folder

The application folder should contain at least these files:

- `pgs_to_plink.py`
- Input files for conversion may already be placed in the same application folder; this is expected and supported, so you can select them directly from there if they were bundled or copied next to the app in advance.

## 🌟 Key Features

- **Smart Parsing:** Automatically identifies required columns (SNP, effect_allele, weight) by recognizing common synonyms in any `.txt`, `.csv`, or `.tsv` PGS score file. [file:1]
- **Configurable Range:** You can restrict processing to a specific SNP index range (start / end) from the GUI to avoid running the full file when you only need a subset. [file:1]
- **Multithreaded rsID Recovery (NCBI API):** Rapidly fetches missing rsIDs from dbSNP using genomic coordinates (`chr:pos`) with checkpointing and resume support. [file:1]
- **Smart Checkpoint System:** Progress is saved every 100 rows. If your internet connection drops, the script resumes exactly where it left off. Checkpoints are automatically cleaned up upon successful completion. [file:1]
- **Interactive & Auto Modes:** Manually select the best rsID from NCBI results, or switch to an automatic mode to let the script pick the newest variants without additional prompts. [file:1]
- **Safe Liftover (hg19 ↔ hg38):** Safely rebuilds patient binary files (`.bed/.bim/.fam`) using temporary IDs to prevent PLINK exclusion bugs when genome builds do not match. [file:1]
- **VCF Support (when PLINK is available):** Can convert raw `.vcf` / `.vcf.gz` files into PLINK binary format on the fly if PLINK is installed and accessible in `PATH`. [file:1]
- **Quality Control (QC):** Computes SNP overlap (coverage percentage) between the PGS risk table and the patient's genotype, and can optionally run a more detailed QC + Liftover step. [file:1]
- **Intelligent File Search:** The original script can search your system for PLINK and Liftover `.chain.gz` files if they are not present in the working folder. With the GUI, you can also point directly to the required inputs. [file:1]
- **Bundled Example Inputs:** The application folder may already contain example PGS score files that you can use immediately for testing and conversion from the GUI, without downloading additional data. [file:1]

## 🛠 Prerequisites

1. **Python 3.7+**.
2. **PLINK 1.9** (Optional but highly recommended for VCF conversion, Liftover of binary files, and automatic PGS calculation).
   - [Download PLINK](https://www.cog-genomics.org/plink/1.9/).
   - Place `plink.exe` somewhere on your PC and add it to `PATH`, or keep it in the same working directory as your genotype data. [file:1][web:35]
3. **Internet Connection** (for installing dependencies, downloading Liftover chains, and NCBI requests). [file:1]

*Note: The tool supports a “file preparation only” mode. If PLINK is not found on your system, it will still prepare all necessary text files, fetch rsIDs, and generate the exact PLINK command for you to run manually later on a PLINK-equipped machine. The GUI exposes this mode via a simple checkbox.* [file:1]

## 🚀 Usage

You can run the script interactively via the terminal:

```powershell
D:/new/python.exe "D:/path/to/pgs_to_plink_gui_always_send.py"
```

After launch:

- Select the `PGS file`.
- Enable `File preparation only` if PLINK is not installed.
- Choose the `rsID search` mode.
- Optionally select a `Genotype file`.
- Click `Run`.

---
*Developed as part of the `vtools-genomics` toolkit.*
