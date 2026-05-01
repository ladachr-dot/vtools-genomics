# 🧬 PGS → PLINK Converter & Calculator

An interactive, high-performance Python utility designed to automate the processing of Polygenic Risk Score (PGS) files and calculate the genetic risk of a patient based on their genotype.

This tool completely eliminates the most tedious bioinformatics tasks: recovering missing rsIDs from NCBI, aligning mismatched genome builds (Liftover), converting raw VCF files, and managing PLINK subprocesses.

## 🌟 Key Features

- **Smart Parsing:** Automatically identifies required columns (SNP, effect_allele, weight) by recognizing common synonyms in any `.txt`, `.csv`, or `.tsv` file.
- **Multithreaded rsID Recovery (NCBI API):** Rapidly fetches missing rsIDs from dbSNP using genomic coordinates (`chr:pos`).
- **Smart Checkpoint System:** Progress is saved every 100 rows. If your internet connection drops, the script resumes exactly where it left off. Checkpoints are automatically cleaned up upon successful completion.
- **Interactive & Auto Modes:** Manually select the best rsID from NCBI results, or instantly switch to "Auto-mode" mid-process to let the script pick the newest variants automatically.
- **Safe Liftover (hg19 ↔ hg38):** Safely rebuilds patient binary files (`.bed/.bim/.fam`) using temporary IDs to prevent PLINK's exclusion flag bugs when genome builds mismatch.
- **VCF Support:** Automatic, on-the-fly conversion of raw `.vcf`/`.vcf.gz` files into PLINK binary format.
- **Quality Control (QC):** Automatically calculates SNP overlap (coverage percentage) between the PGS risk table and the patient's genotype.
- **Intelligent File Search:** No need to configure system variables. If `plink.exe` or Liftover `.chain.gz` files are missing from the folder, the script will scan your entire computer to find and use them.

## 🛠 Prerequisites

1. **Python 3.7+**.
2. **PLINK 1.9** (Optional but highly recommended).
   - [Download PLINK](https://www.cog-genomics.org/plink/1.9/).
   - Just place `plink.exe` anywhere on your PC; the script will find it automatically.
3. **Internet Connection** (for installing dependencies, downloading Liftover chains, and NCBI requests).

*Note: The script features a "Standalone Mode". If PLINK is not found on your system, it will still prepare all necessary text files, fetch rsIDs, and generate the exact terminal command for you to run manually on a PLINK-equipped machine.*

## 🚀 Usage

You can run the script interactively via the terminal:

```bash
python pgstoplink.py
```

### Interactive Workflow:
1. **Provide the PGS file:** Enter the name of your score file (e.g., `my_score` or `my_score.txt`). The script supports extensionless input and will find the file in the current or subdirectories.
2. **rsID Fetching:** If rsIDs are missing, specify the row range to process. Choose between interactive selection or multithreaded auto-fetching.
3. **Patient Genotype:** Provide the path to the patient's `.vcf` or `.bim` file.
4. **QC & Liftover:** Review the SNP overlap report. If coverage is critically low (<10%), the script will offer to safely Liftover the genome build.
5. **Final Calculation:** Confirm the final step to let the script invoke PLINK and generate the `_PGS_results.profile` file containing the patient's calculated risk.

---
*Developed as part of the `vtools-genomics` toolkit.*
