
"""
PGS → PLINK Converter & Calculator v5.0
=======================================

Core module for preparing Polygenic Risk Score (PGS) files and patient 
genotypes for analysis using PLINK.

Workflow:
1. Score File Parsing:
   - Automatically detects SNP, A1 (Effect Allele), and BETA (Weight) columns.
   - If rsIDs are missing (only chr:pos available), triggers multithreaded fetch.
2. rsID Fetching (NCBI dbSNP):
   - Asynchronous API requests to NCBI based on genomic coordinates.
   - Checkpoint system: progress is saved every 100 SNPs to prevent data loss.
3. Patient Genotype Integration:
   - Converts .vcf / .vcf.gz to PLINK binary format (--make-bed).
   - Performs QC: checks the SNP overlap (coverage) percentage.
4. Liftover (Genome Build Conversion):
   - If coverage is below 10%, safely lifts over patient coordinates (e.g., hg38 -> hg19)
     using temporary IDs (TEMP_SNP_X) to avoid PLINK --exclude flag errors.
5. PGS Calculation:
   - Automatically runs `plink --score` and generates the final *.profile file.
6. "No PLINK" Fallback Mode:
   - If PLINK is not installed, the script switches to preparation mode (parsing, 
     rsID fetching, generating .txt and .bim files), outputting the final 
     terminal command for manual execution on a PLINK-equipped machine.
"""

import sys
import os
import subprocess
import requests
import pandas as pd
import gzip
import time
import json
import urllib.request
import urllib.error
import ssl
import re
import shutil
import concurrent.futures
from pathlib import Path
from io import BytesIO

# Global Variables
LIFTOVER_AVAILABLE = False
PLINK_AVAILABLE = False
MAX_ROWS = 10
PLINK_CMD = "plink"

##########################################################
# AUTO-INSTALLATION BLOCK
##########################################################

def install_package(package_name):
    """Installs a Python package via pip."""
    print(f"📦 Installing {package_name}...")
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", package_name,
            "--quiet", "--disable-pip-version-check"
        ])
        print(f"   ✅ {package_name} installed")
        return True
    except subprocess.CalledProcessError:
        print(f"   ❌ Failed to install {package_name}")
        return False


def check_and_install_dependencies():
    """Checks and installs all required packages."""
    print("=" * 60)
    print("🔧 CHECKING DEPENDENCIES")
    print("=" * 60)

    dependencies = {
        'requests': 'requests',
        'pandas': 'pandas',
        'pyliftover': 'pyliftover',
    }

    all_ok = True

    for import_name, package_name in dependencies.items():
        try:
            __import__(import_name)
            print(f"✅ {package_name} is already installed")
        except ImportError:
            print(f"⚠️  {package_name} not found")
            if not install_package(package_name):
                all_ok = False

    if all_ok:
        try:
            global LiftOver, LIFTOVER_AVAILABLE
            from pyliftover import LiftOver
            LIFTOVER_AVAILABLE = True
            print("✅ pyliftover is ready")
        except ImportError:
            LIFTOVER_AVAILABLE = False
            print("❌ Failed to load pyliftover")

    print()
    return all_ok

def download_liftover_chain(chain_name):
    """Downloads the liftover chain from UCSC to the script's directory."""

    script_dir = os.path.dirname(os.path.abspath(__file__))
    chain_path = os.path.join(script_dir, chain_name)
    
    if 'hg19ToHg38' in chain_name:
        url = f"http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/{chain_name}"
    elif 'hg38ToHg19' in chain_name:
        url = f"http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/{chain_name}"
    else:
        print(f"❌ Unknown chain: {chain_name}")
        return False

    if os.path.exists(chain_path):
        size = os.path.getsize(chain_path)
        if size > 1000:
            print(f"✅ {chain_name} is already downloaded ({size:,} bytes)")
            return True

    print(f"📥 Downloading {chain_name}...")
    print(f"   {url}")

    try:
        resp = requests.get(url, timeout=120, stream=True)
        resp.raise_for_status()

        total_size = int(resp.headers.get('content-length', 0))
        downloaded = 0

        with open(chain_path, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)

                if total_size > 0:
                    percent = (downloaded / total_size) * 100
                    print(f"\r   Progress: {percent:.1f}% ({downloaded:,}/{total_size:,} bytes)", end='')

        print(f"\n   ✅ Downloaded: {os.path.getsize(chain_path):,} bytes")
        return True

    except Exception as e:
        print(f"\n   ❌ Download error: {e}")
        if os.path.exists(chain_path):
            os.remove(chain_path)
        return False

def check_and_download_chains():
    """Checks for chains locally, searches the PC, or downloads if missing."""
    print("=" * 60)
    print("🔗 CHECKING LIFTOVER CHAINS")
    print("=" * 60)

    chains = [
        'hg19ToHg38.over.chain.gz',
        'hg38ToHg19.over.chain.gz'
    ]

    root_path = 'C:\\' if os.name == 'nt' else '/'
    all_ok = True

    script_dir = os.path.dirname(os.path.abspath(__file__))

    for chain in chains:
        chain_path = os.path.join(script_dir, chain)
        if os.path.exists(chain) and os.path.getsize(chain) > 1000:
            print(f"   ✅ {chain} is already in the current directory.")
            continue
            
        print(f"   ⚠️ {chain} not found locally. Initiating search...")
        found_path = find_file(chain, root_path)
        
        if found_path:
            print(f"   📋 Copying {chain} to current directory...")
            shutil.copy(found_path, chain)
        else:
            print(f"   ❌ Failed to copy: {e}")
            if not download_liftover_chain(chain):
                all_ok = False

    if not all_ok:
        print("\n❌ Not all chains were downloaded or found!")
        return False

    print()
    return True


def setup_environment():
    """Full environment setup with a silent pre-check to avoid unnecessary logs."""
    global LiftOver, LIFTOVER_AVAILABLE
    
    # --- 1. SILENT PRE-CHECK ---
    env_ready = True
    try:
        __import__('requests')
        __import__('pandas')
        __import__('pyliftover')
    except ImportError:
        env_ready = False
        
    chains = [
        'hg19ToHg38.over.chain.gz',
        'hg38ToHg19.over.chain.gz'
    ]
    for chain in chains:
        if not os.path.exists(chain) or os.path.getsize(chain) < 1000:
            env_ready = False
            
    if env_ready:
        from pyliftover import LiftOver
        LIFTOVER_AVAILABLE = True
        return True
        
    # --- 2. VERBOSE SETUP ---
    print("\n" + "=" * 60)
    print("🚀 SETTING UP ENVIRONMENT")
    print("=" * 60 + "\n")

    print("📋 Checking pip...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "--version"],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("   ✅ pip is working\n")
    except:
        print("   ❌ pip is not working!\n")
        return False

    if not check_and_install_dependencies():
        print("❌ Not all packages are installed.\n")
        return False

    if not check_and_download_chains():
        print("❌ Liftover chains are not ready.\n")
        return False

    print("✅ Environment is ready to go!\n")
    return True

###################################################################################
# MULTITHREADING 
###################################################################################

def process_single_snp(args):
    i, chrom, pos_hg38, pos_hg19, auto_mode = args
    if pos_hg38 == '.' or pd.isna(pos_hg38):
        return i, '.'
    
    # Slight delay to prevent overloading the NCBI API during multithreading
    time.sleep(0.1)
    rsid = find_rsid_interactive(chrom, pos_hg38, auto_mode=auto_mode)
    return i, rsid if rsid else '.'


###################################################################################
# ФУНКЦИИ ДЛЯ РАБОТЫ С КОЛОНКАМИ (из ноутбука)
###################################################################################

COLUMN_SYNONYMS = {
    "SNP": [
        "rsid", "rs_id", "snp", "variant_id", "marker", "markername", "id",
        "hm_rsid"
    ],
    "A1": [
        "effect_allele", "ea", "a1", "tested_allele", "allele",
        "risk_allele", "hm_effect_allele"
    ],
    "BETA": [
        "beta", "effect_weight", "weight", "score", "log_odds", "logor",
        "hm_beta"
    ]
}

BIM_COLUMN_SYNONYMS = {
    "chr_name": ["chr_name", "chromosome", "chr", "chrom"],
    "chr_position": ["chr_position", "position", "pos", "hm_pos", "bp"],
    "A2": ["other_allele", "a2", "reference_allele", "ref", "hm_other_allele"]
}


def normalize_column_name(name: str) -> str:
    """Normalizes column names to a standard format."""
    return (
        str(name)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )

def find_column(columns, synonyms):
    """Finds the first column that matches one of the synonyms."""
    normalized_columns = {
        normalize_column_name(col): col
        for col in columns
    }

    for synonym in synonyms:
        normalized_synonym = normalize_column_name(synonym)
        if normalized_synonym in normalized_columns:
            return normalized_columns[normalized_synonym]

    return None

def extract_genome_build(file_path: str) -> str:
    """Extracts genome_build from the PGS header."""
    file_path = Path(file_path)

    if not file_path.exists():
        return "unknown"

    with open(file_path, "r", encoding="utf-8", errors="replace") as file:
        for line in file:
            line = line.strip()

            if line.startswith("#genome_build"):
                return line.split("=", 1)[1].strip()

            if line and not line.startswith("#"):
                break

    return "unknown"

###################################################################################
# COORDINATE CONVERSION FUNCTIONS
###################################################################################

def liftover_coordinates(df, from_build, to_build, pos_col='chr_position'):
    """Converts coordinates between genome builds."""

    def normalize(build):
        b = str(build).lower().replace('grch', '').replace('hg', '')
        if '37' in b or '19' in b:
            return '37'
        elif '38' in b:
            return '38'
        return b

    from_b = normalize(from_build)
    to_b = normalize(to_build)

    if from_b == to_b:
        print(f"   ✅ Builds are identical (hg{from_b}). Conversion not required.")
        return df, True

    if from_b == '37' and to_b == '38':
        chain_file = 'hg19ToHg38.over.chain.gz'
        direction = 'hg19 → hg38'
    elif from_b == '38' and to_b == '37':
        chain_file = 'hg38ToHg19.over.chain.gz'
        direction = 'hg38 → hg19'
    else:
        print(f"   ❌ Unsupported conversion: {from_build} → {to_build}")
        return df, False

    if not os.path.exists(chain_file):
        print(f"   ❌ Chain file not found: {chain_file}")
        return df, False

    print(f"   🔄 Converting {direction}...")
    lo = LiftOver(chain_file)
    new_coords = []

    for idx, row in df.iterrows():
        chrom = str(row.get('chr_name', '')).replace('chr', '').replace('CHR', '')
        pos_val = row.get(pos_col, '')

        try:
            pos = int(float(pos_val))
            result = lo.convert_coordinate(f'chr{chrom}', pos - 1)
            if result:
                new_coords.append(str(result[0][1] + 1))
            else:
                new_coords.append('.')
        except (ValueError, TypeError):
            new_coords.append('.')

    df = df.copy()

    if to_b == '38':
        new_col = 'chr_position_hg38'
    else:
        new_col = 'chr_position_hg19'

    df[new_col] = new_coords

    converted = sum(1 for x in new_coords if x != '.')
    print(f"   ✅ Converted: {converted}/{len(new_coords)}")

    return df, True


###################################################################################
# BIM, BED, FAM CONVERSION
###################################################################################

def create_bim_from_score(df, output_path, genome_build="hg19"):
    """
    Creates a BIM file from the score file data.
    BIM format: chr, SNP, cm, pos, A1, A2
    """
    print(f"\n📝 Creating BIM file...")

    # Look for required columns
    chr_col = find_column(df.columns, BIM_COLUMN_SYNONYMS["chr_name"])
    pos_col = find_column(df.columns, BIM_COLUMN_SYNONYMS["chr_position"])
    snp_col = find_column(df.columns, COLUMN_SYNONYMS["SNP"])
    a1_col = find_column(df.columns, COLUMN_SYNONYMS["A1"])
    a2_col = find_column(df.columns, BIM_COLUMN_SYNONYMS["A2"])

    # If A2 is missing, try to find other_allele
    if a2_col is None:
        for col in df.columns:
            if 'other_allele' in col.lower() or col.lower() == 'a2':
                a2_col = col
                break

    # Verify everything was found
    missing = []
    if chr_col is None:
        missing.append("chr_name")
    if pos_col is None:
        missing.append("chr_position")
    if snp_col is None:
        missing.append("SNP/rsID")
    if a1_col is None:
        missing.append("effect_allele/A1")

    if missing:
        print(f"   ⚠️  Required columns for BIM not found: {missing}")
        print(f"   Available columns: {list(df.columns)}")
        return None

    # Create BIM
    bim_df = pd.DataFrame()
    bim_df['chr'] = df[chr_col].astype(str).str.replace('chr', '', case=False)
    bim_df['SNP'] = df[snp_col].astype(str)
    bim_df['cm'] = 0  # Genetic distance (usually 0)
    bim_df['pos'] = pd.to_numeric(df[pos_col], errors='coerce')
    bim_df['A1'] = df[a1_col].astype(str).str.upper()

    # A2: if other_allele exists use it, otherwise use '.'
    if a2_col:
        bim_df['A2'] = df[a2_col].astype(str).str.upper()
    else:
        bim_df['A2'] = '.'

    # Cleaning
    before = len(bim_df)

    # Remove rows without coordinates
    bim_df = bim_df.dropna(subset=['pos'])

    # Remove rows with '.' in chr, SNP, A1
    bim_df = bim_df[bim_df['chr'] != '.']
    bim_df = bim_df[bim_df['SNP'] != '.']
    bim_df = bim_df[bim_df['A1'] != '.']

    # Remove duplicates by SNP
    bim_df = bim_df.drop_duplicates(subset=['SNP'], keep='first')

    after = len(bim_df)

    # Save to file
    bim_df.to_csv(output_path, sep='\t', index=False, header=False)

    print(f"   ✅ BIM file created: {output_path}")
    print(f"   Rows: {after} (removed: {before - after})")
    print(f"   Build: {genome_build}")
    print(f"   Format: chr, SNP, cm, pos, A1, A2")

    # Statistics
    print(f"\n   📊 BIM Statistics:")
    print(f"   Chromosomes: {sorted(bim_df['chr'].unique())}")
    print(f"   Position range: {bim_df['pos'].min():,} - {bim_df['pos'].max():,}")
    print(f"   A1 alleles: {sorted(bim_df['A1'].unique())}")

    return bim_df


def read_bim_sample(bim_path: str, nrows: int = 1000) -> pd.DataFrame:
    """Reads a small part of a .bim file without fully loading it."""
    bim_path = Path(bim_path)
    if not bim_path.exists():
        raise FileNotFoundError(f"BIM file not found: {bim_path}")

    bim = pd.read_csv(bim_path, sep=r"\s+", header=None, nrows=nrows)
    bim.columns = ["chr", "SNP", "cm", "pos", "A1", "A2"]

    bim["SNP"] = bim["SNP"].astype(str)
    bim["A1"] = bim["A1"].astype(str).str.upper()
    bim["A2"] = bim["A2"].astype(str).str.upper()
    bim["pos"] = pd.to_numeric(bim["pos"], errors="coerce")

    return bim

def light_bim_report(bim_path: str, nrows: int = 1000):
    """Quick report for a .bim file."""
    bim_sample = read_bim_sample(bim_path, nrows=nrows)

    rsid_rate = bim_sample["SNP"].str.match(r"^rs\d+$", na=False).mean()

    print("\n===== LIGHT BIM REPORT =====")
    print(f"BIM file: {bim_path}")
    print(f"Rows checked: {len(bim_sample)}")
    print(f"rsID-like IDs: {rsid_rate:.2%}")
    print(f"Chromosomes: {sorted(bim_sample['chr'].astype(str).unique())}")
    print(f"Position range: {bim_sample['pos'].min():,} - {bim_sample['pos'].max():,}")

    return bim_sample


def run_full_qc_report(score_path: str, bim_path: str, genome_build: str):
    """Full QC: checking SNP overlap and allele compatibility."""
    score = pd.read_csv(score_path, sep="\t")
    bim = pd.read_csv(bim_path, sep=r"\s+", header=None)
    bim.columns = ["chr", "SNP", "cm", "pos", "A1", "A2"]

    score["SNP"] = score["SNP"].astype(str)
    score["A1"] = score["A1"].astype(str).str.upper()
    bim["SNP"] = bim["SNP"].astype(str)
    bim["A1"] = bim["A1"].astype(str).str.upper()

    score_snps = set(score["SNP"])
    bim_snps = set(bim["SNP"])
    common_snps = score_snps & bim_snps

    coverage = len(common_snps) / len(score_snps) if score_snps else 0

    print("\n===== FULL QC REPORT =====")
    print(f"BIM file used: {bim_path}")
    print(f"SNP in score: {len(score_snps):,}")
    print(f"SNP in genotype: {len(bim_snps):,}")
    print(f"Common SNP: {len(common_snps):,}")
    print(f"Coverage: {coverage:.2%}")

    if coverage < 0.10:
        print("⚠️  WARNING: Very low SNP coverage! Possible genome build mismatch (hg19/hg38).")
        do_liftover = input("🔄 Do you want to try converting your PLINK files to the hg19 build? (y/n): ").strip().lower()
        
        if do_liftover == 'y':
            bfile_prefix = bim_path.rsplit('.bim', 1)[0]
            target = genome_build if genome_build and genome_build != "unknown" else "hg19"
            source = "hg38" if target == "hg19" else "hg19"
            
            print(f"   Attempting Liftover from presumed build {source} to {target}...")
            new_prefix = convert_user_bfile_liftover(bfile_prefix, source_build=source, target_build=target)
            
            if new_prefix:
                print(f"\n✅ Recalculating QC for the new files {new_prefix}...")
                run_full_qc_report(score_path, f"{new_prefix}.bim", genome_build)
                return f"{new_prefix}.bim"  
            else:
                print("❌ Liftover failed.")
                return bim_path
        return bim_path
                
    elif coverage < 0.50:
        print("⚠️  WARNING: Moderate SNP coverage.")
        return bim_path
    else:
        print("✅ SNP coverage: OK")
        return bim_path


def convert_user_bfile_liftover(bfile_prefix, source_build='hg38', target_build='hg19'):
    """
    Safe conversion of user .bed/.bim/.fam files via temporary IDs 
    to avoid PLINK exclusion bugs with missing names.
    """
    global PLINK_AVAILABLE
    if not PLINK_AVAILABLE:
        print("\n❌ ERROR: Rebuilding binary files after Liftover requires PLINK to be installed!")
        return None

    print(f"\n{'=' * 60}")
    print(f"🧬 GENOTYPE LIFTOVER ({source_build} -> {target_build})")
    print(f"   Files: {bfile_prefix}.*")
    print(f"{'=' * 60}")

    bim_file = f"{bfile_prefix}.bim"
    if not os.path.exists(bim_file):
        print(f"❌ File {bim_file} not found!")
        return None

    # Read BIM
    bim_df = pd.read_csv(bim_file, sep=r'\s+', header=None, 
                         names=['chr', 'rsid_orig', 'cm', 'pos', 'a1', 'a2'], dtype=str)

    print(f"📊 {len(bim_df):,} SNPs loaded from BIM file.")

    # STEP 1: Assign a unique temporary ID to ALL SNPs
    bim_df['temp_id'] = [f"TEMP_SNP_{i}" for i in range(len(bim_df))]
    
    temp_bim_file = f"{bfile_prefix}_temp.bim"
    bim_df[['chr', 'temp_id', 'cm', 'pos', 'a1', 'a2']].to_csv(
        temp_bim_file, sep='\t', index=False, header=False
    )

    # STEP 2: Rename columns for liftover function
    lift_df = bim_df.rename(columns={'chr': 'chr_name', 'pos': 'chr_position'})

    source_num = '38' if source_build == 'hg38' else '37'
    target_num = '37' if target_build == 'hg19' else '38'

    print(f"⏳ Converting coordinates...")
    lifted_df, ok = liftover_coordinates(lift_df, source_num, target_num, 'chr_position')
    
    if not ok:
        print("❌ Error converting coordinates.")
        if os.path.exists(temp_bim_file): os.remove(temp_bim_file)
        return None

    target_col = f'chr_position_{target_build}'
    
    # Split into mapped and unmapped
    mapped_mask = (lifted_df[target_col] != '.') & (lifted_df[target_col].notna())
    mapped_snps = lifted_df[mapped_mask].copy()
    unmapped_snps = lifted_df[~mapped_mask].copy()

    print(f"✅ Successfully converted: {len(mapped_snps):,} SNPs")
    print(f"❌ Failed to convert: {len(unmapped_snps):,} SNPs (these will be removed)")

    out_prefix = f"{bfile_prefix}_{target_build}"
    update_map_file = f"{bfile_prefix}_update_map.txt"
    exclude_file = f"{bfile_prefix}_exclude.txt"
    update_name_file = f"{bfile_prefix}_update_name.txt"

    # STEP 3: Files for PLINK (using temp_id)
    mapped_snps[['temp_id', target_col]].to_csv(update_map_file, sep='\t', index=False, header=False)
    unmapped_snps['temp_id'].to_csv(exclude_file, index=False, header=False)
    mapped_snps[['temp_id', 'rsid_orig']].to_csv(update_name_file, sep='\t', index=False, header=False)

    # STEP 4: Run PLINK
    print("\n🚀 Running PLINK to rebuild binary files...")
    plink_cmd = [
        "plink",
        "--bfile", bfile_prefix,
        "--bim", temp_bim_file,      
        "--exclude", exclude_file,   
        "--update-map", update_map_file, 
        "--update-name", update_name_file, 
        "--make-bed",
        "--out", out_prefix
    ]

    try:
        subprocess.run(plink_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"✅ Conversion complete! New files: {out_prefix}.bed/bim/fam")
    except subprocess.CalledProcessError:
        print("❌ Error running PLINK! Ensure PLINK is installed.")
        return None

    # Cleanup temp files
    for tmp_f in [update_map_file, exclude_file, update_name_file, temp_bim_file]:
        if os.path.exists(tmp_f): os.remove(tmp_f)

    return out_prefix

###################################################################################
# rsID SEARCH FUNCTIONS
###################################################################################

def get_rsid_info(snp_uid):
    """Fetches detailed SNP info from dbSNP."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    fetch_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                 f"db=snp&id={snp_uid}&retmode=json")

    try:
        with urllib.request.urlopen(fetch_url, timeout=10, context=ctx) as resp:
            data = json.loads(resp.read().decode())

        snp_info = data.get('result', {}).get(str(snp_uid), {})

        rsid = (snp_info.get('snp_id') or
                snp_info.get('refsnp_id') or
                snp_info.get('uid'))

        if rsid:
            rsid = str(rsid)
            if rsid.isdigit():
                rsid = 'rs' + rsid

        info = {
            'rsid': rsid if rsid and rsid.startswith('rs') else None,
            'uid': snp_uid,
            'clinical_significance': snp_info.get('clinical_significance', 'not specified'),
            'snp_class': snp_info.get('snp_class', 'not specified'),
            'gene_name': snp_info.get('gene_name', 'not specified'),
            'chr': snp_info.get('chr', 'not specified'),
            'chr_pos': snp_info.get('chrpos', 'not specified'),
            'alleles': snp_info.get('alleles', 'not specified'),
            'validated': snp_info.get('validated', 'not specified'),
        }

        return info

    except Exception as e:
        return None

def find_all_rsids_with_info(chrom, pos):
    """Finds ALL rsIDs for a position with detailed info."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    search_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                  f"db=snp&term={chrom}[Chromosome]+AND+{pos}[Base+Position]"
                  f"&retmax=10&retmode=json")

    try:
        with urllib.request.urlopen(search_url, timeout=10, context=ctx) as resp:
            data = json.loads(resp.read().decode())

        id_list = data.get('esearchresult', {}).get('idlist', [])

        if not id_list:
            return []

        all_snp_info = []
        for uid in id_list:
            info = get_rsid_info(uid)
            if info and info['rsid']:
                all_snp_info.append(info)
            time.sleep(0.1)

        seen = set()
        unique_info = []
        for info in all_snp_info:
            if info['rsid'] not in seen:
                seen.add(info['rsid'])
                unique_info.append(info)

        return unique_info

    except Exception as e:
        return []


def display_snp_options(snp_list, chrom, pos):
    """Displays a list of SNPs with detailed info to the user."""
    print(f"\n   🔬 Position chr{chrom}:{pos}")
    print(f"   {'─' * 50}")
    print(f"   SNPs found: {len(snp_list)}")
    print(f"   {'─' * 50}")

    for i, snp in enumerate(snp_list, 1):
        rsid = snp['rsid']
        clin_sig = snp['clinical_significance']
        snp_class = snp['snp_class']
        gene = snp['gene_name']

        if 'pathogenic' in str(clin_sig).lower():
            sig_icon = "🔴"
        elif 'benign' in str(clin_sig).lower():
            sig_icon = "🟢"
        elif 'uncertain' in str(clin_sig).lower() or 'conflicting' in str(clin_sig).lower():
            sig_icon = "🟡"
        elif 'not provided' in str(clin_sig).lower():
            sig_icon = "⚪"
        else:
            sig_icon = "🔵"

        print(f"   [{i}] {sig_icon} {rsid}")
        print(f"       Class: {snp_class}")
        print(f"       Gene: {gene}")
        print()

    print(f"   {'─' * 50}")
    print(f"   [0] Skip (do not assign rsID)")
    print(f"   [A] Auto-select this one (newest)")
    print(f"   [ALL] SWITCH TO AUTO MODE FOR EVERYTHING")
    print(f"   {'─' * 50}")

    while True:
        choice = input(f"   Your choice: ").strip().upper()

        if choice == 'ALL':
            return 'AUTO'

        if choice == 'A':
            def rs_number(snp):
                try:
                    return int(snp['rsid'].replace('rs', ''))
                except:
                    return 0

            best = max(snp_list, key=rs_number)
            print(f"   ✅ Auto-selected: {best['rsid']} (newest)\n")
            return best['rsid']

        if choice == '0':
            print(f"   ⏭️  Skipped\n")
            return None

        if choice.isdigit():
            idx = int(choice) - 1
            if 0 <= idx < len(snp_list):
                print(f"   ✅ Selected: {snp_list[idx]['rsid']}\n")
                return snp_list[idx]['rsid']

        print(f"   ❌ Invalid choice. Try again.")

def find_rsid_interactive(chrom, pos, auto_mode=False):
    """Searches for rsID with interactive selection."""
    all_snp_info = find_all_rsids_with_info(chrom, pos)

    if not all_snp_info:
        return None

    if len(all_snp_info) == 1:
        print(f"✅ {all_snp_info[0]['rsid']}")
        return all_snp_info[0]['rsid']

    if auto_mode:
        def rs_number(snp):
            try:
                return int(snp['rsid'].replace('rs', ''))
            except:
                return 0

        best = max(all_snp_info, key=rs_number)
        print(f"⚠️  Multiple rsIDs found, auto-selecting: {best['rsid']}")
        return best['rsid']
    else:
        return display_snp_options(all_snp_info, chrom, pos)

###################################################################################
# SCORE FILE CONVERSION
###################################################################################

def convert_score_file(input_path, output_path):
    """Converts score file to PLINK format: SNP / A1 / BETA."""
    input_path = Path(input_path)
    output_path = Path(output_path)

    # Read file
    genome_build = extract_genome_build(input_path)

    if genome_build != "unknown":
        df = pd.read_csv(input_path, sep=r"\s+", comment="#", engine="python")
    else:
        df = pd.read_csv(input_path, sep=None, engine="python")

    print(f"\n📊 Original columns: {list(df.columns)}")

    # Find target columns
    found_columns = {}
    for target_name, synonyms in COLUMN_SYNONYMS.items():
        original_column = find_column(df.columns, synonyms)
        if original_column is None:
            raise ValueError(
                f"Could not find a column for {target_name}. "
                f"Check column names in the file."
            )
        found_columns[target_name] = original_column

    print("Found columns:")
    for target_name, original_column in found_columns.items():
        print(f"  {target_name}: {original_column}")

    # Form the result
    result = df[[found_columns["SNP"], found_columns["A1"], found_columns["BETA"]]].copy()
    result.columns = ["SNP", "A1", "BETA"]

    # Cleaning
    rows_before = len(result)
    result = result.dropna(subset=["SNP", "A1", "BETA"])

    duplicates_before = result["SNP"].duplicated().sum()
    result = result.drop_duplicates(subset=["SNP"], keep="first")
    duplicates_after = result["SNP"].duplicated().sum()

    rows_after = len(result)

    print(f"\n📊 Cleaning statistics:")
    print(f"   Rows before: {rows_before}")
    print(f"   Rows after: {rows_after}")
    print(f"   Duplicate SNPs removed: {duplicates_before}")

    # Save score file
    result.to_csv(output_path, sep="\t", index=False)
    print(f"\n✅ Score file saved: {output_path}")

    return result, df, genome_build

###################################################################################
# ALTERNATIVE METHOD (rsID SEARCH)
###################################################################################

def alternative(filepath, genome_build, skip_lines):
    """Alternative method with rsID search, checkpoints, multithreading, and ranges."""
    print(f"\n{'=' * 60}")
    print(f"🔄 RSID SEARCH:")
    print(f"   File: {os.path.basename(filepath)}")
    print(f"{'=' * 60}")

    if not LIFTOVER_AVAILABLE:
        print("❌ pyliftover is not installed!")
        return

    # Ask for range
    start_idx_input = input("Enter starting SNP index (Enter = 0): ").strip()
    end_idx_input = input("Enter ending SNP index (Enter = to end): ").strip()
    
    start_idx = int(start_idx_input) if start_idx_input.isdigit() else 0
    end_idx = int(end_idx_input) if end_idx_input.isdigit() else None

    # Load data
    try:
        df_full = pd.read_csv(filepath, sep='\t', skiprows=skip_lines, dtype=str)
        print(f"\n📊 Total rows in file: {len(df_full)}")
        
        if end_idx is None:
            end_idx = len(df_full)
        df = df_full.iloc[start_idx:end_idx].copy()
        
        df.reset_index(drop=True, inplace=True)
        print(f"📊 Selected row range [{start_idx}:{end_idx}]. To process: {len(df)}")
    except Exception as e:
        print(f"❌ Read error: {e}")
        return

    if 'chr_name' not in df.columns or 'chr_position' not in df.columns:
        print("❌ chr_name/chr_position columns are missing!")
        return

    # INIT CHECKPOINT
    dir_name = os.path.dirname(filepath)
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    checkpoint_file = os.path.join(dir_name, f"{base_name}_checkpoint_{start_idx}_{end_idx}.csv")

    if os.path.exists(checkpoint_file):
        print(f"♻️ Checkpoint file found: {checkpoint_file}. Loading saved progress...")
        df = pd.read_csv(checkpoint_file, sep='\t', dtype=str)
    else:
        if 'rsID' not in df.columns:
            df['rsID'] = '.'
        df['chr_position_hg19_original'] = df['chr_position']
        
        print(f"\n🔹 Step 1: Converting hg19 → hg38...")
        df, ok = liftover_coordinates(df, '37', '38', 'chr_position')
        if not ok:
            print("❌ Conversion failed!")
            return

    unprocessed_mask = df['rsID'] == '.'
    unprocessed_indices = df[unprocessed_mask].index.tolist()

    print(f"\n🔹 Step 2: Searching for rsID in dbSNP")
    if not unprocessed_indices:
        print("✅ All SNPs in this range have already been processed!")
    else:
        print(f"   [1] Automatic (multithreaded, newest rsID)")
        print(f"   [2] Interactive (SINGLE-THREADED, manual selection)")
        mode_choice = input(f"   Mode (1/2, Enter=1): ").strip()
        auto_mode = mode_choice != '2'

        print(f"\n🔹 RSID Search (Remaining to find: {len(unprocessed_indices)}):")
        print("-" * 60)

        tasks = []
        for i in unprocessed_indices:
            chrom = str(df.at[i, 'chr_name'])
            pos_hg38 = str(df.at[i, 'chr_position_hg38'])
            pos_hg19 = str(df.at[i, 'chr_position_hg19_original'])
            tasks.append((i, chrom, pos_hg38, pos_hg19, auto_mode))

        processed_in_session = 0
        checkpoint_interval = 100

        max_workers = 5 if auto_mode else 1

        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            for i in unprocessed_indices:
                chrom = str(df.at[i, 'chr_name'])
                pos_hg38 = str(df.at[i, 'chr_position_hg38'])
                pos_hg19 = str(df.at[i, 'chr_position_hg19_original'])
                
                # Если авто-режим, то 5 потоков, если интерактивный - ждем ответа
                if auto_mode:
                    rsid = find_rsid_interactive(chrom, pos_hg38, auto_mode=True)
                else:
                    rsid = find_rsid_interactive(chrom, pos_hg38, auto_mode=False)
                    
                    # Если пользователь вернул 'AUTO', переключаем режим для всех остальных
                    if rsid == 'AUTO':
                        print("\n⚡ Switching to AUTOMATIC mode for the rest of the SNPs...")
                        auto_mode = True
                        # Вызываем поиск еще раз уже в авто-режиме
                        rsid = find_rsid_interactive(chrom, pos_hg38, auto_mode=True)
                
                df.at[i, 'rsID'] = rsid
                processed_in_session += 1
                
                if rsid and rsid != '.':
                    print(f"✅ Success: row {i} -> {rsid}")

                if processed_in_session % checkpoint_interval == 0:
                    df.to_csv(checkpoint_file, sep='\t', index=False)
                    print(f"💾 [Checkpoint] Progress saved ({processed_in_session} items)")

        df.to_csv(checkpoint_file, sep='\t', index=False)
        print("-" * 60)
        print(f"✅ Search complete. Checkpoint updated.")

    # Create PLINK and BIM
    print(f"\n🔹 Step 3: Creating output files...")
    
    weight_col = next((c for c in df.columns if 'weight' in c.lower()), None)
    allele_col = next((c for c in df.columns if 'effect_allele' in c.lower()), None)

    if weight_col and allele_col:
        plink_df = df[df['rsID'] != '.'].copy()
        plink_df = plink_df[['rsID', allele_col, weight_col]].copy()
        plink_df[weight_col] = pd.to_numeric(plink_df[weight_col], errors='coerce')
        plink_df = plink_df.dropna()
        plink_df = plink_df.drop_duplicates(subset=['rsID', allele_col])

        score_out = os.path.join(dir_name, f"{base_name}_plink_hg19_{start_idx}_{end_idx}.txt")
        plink_df.to_csv(score_out, sep='\t', index=False, header=False)
        print(f"   ✅ Score file: {score_out}")

        bim_out = os.path.join(dir_name, f"{base_name}_plink_hg19_{start_idx}_{end_idx}.bim")
        create_bim_from_score(df, bim_out, genome_build or "hg19")

        if os.path.exists(checkpoint_file):
            try:
                os.remove(checkpoint_file)
                print(f"\n   🧹 Cleaned up temporary checkpoint: {os.path.basename(checkpoint_file)}")
            except Exception as e:
                pass

        print(f"\n🚀 PLINK command:")
        print(f"   plink --bfile your_data --score {score_out} 1 2 3 header")
    else:
        print("⚠️ Allele/weight columns not found")

###################################################################################
# MAIN PROCESSING
###################################################################################

def find_file(filename, search_path):
    """Searches for a file on the computer. Supports extensionless input."""
    print(f"🔍 Searching for file '{filename}' in '{search_path}'...")
    
    # Расширения, которые мы будем проверять, если пользователь их не ввел
    valid_extensions = ['.txt', '.csv', '.tsv']
    search_names = [filename]
    
    if not any(filename.endswith(ext) for ext in valid_extensions):
        search_names.extend([f"{filename}{ext}" for ext in valid_extensions])

    for root, dirs, files in os.walk(search_path):
        for name in search_names:
            if name in files:
                found_path = os.path.join(root, name)
                print(f"   ✅ Found: {found_path}")
                return found_path
    return None

def process_table(filepath):
    """Main file processing logic."""
    print(f"\n📂 Processing file: {filepath}")

    skip_lines = 0
    genome_build = None
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#'):
                    skip_lines += 1
                    if 'genome_build=' in line.lower():
                        genome_build = line.strip().split('=')[1].strip()
                        print(f"📊 Original build: {genome_build}")
                else:
                    break
    except Exception as e:
        print(f"❌ Header read error: {e}")
        return

    try:
        df = pd.read_csv(filepath, sep='\t', skiprows=skip_lines, low_memory=False)
        print(f"📊 Total rows: {len(df)}")
    except Exception as e:
        print(f"❌ Data read error: {e}")
        return

    top_10_rows = df.head(10)
    rs_col = None

    for col in df.columns:
        if top_10_rows[col].astype(str).str.contains(r'^rs\d+', case=False, na=False, regex=True).any():
            rs_col = col
            break

    dir_name = os.path.dirname(filepath)
    base_name = os.path.splitext(os.path.basename(filepath))[0]

    if rs_col:
        print(f"\n✅ rsID already exists in column: '{rs_col}'")

        score_out = os.path.join(dir_name, f"{base_name}_plink.txt")
        converted_df, full_df, build = convert_score_file(filepath, score_out)

        user_data_check = input("\n🔍 Do you have a user genotype file to calculate PGS? (y/n): ").strip().lower()
        if user_data_check == 'y':
            
            bim_out = os.path.join(dir_name, f"{base_name}_plink.bim")
            create_bim_from_score(full_df, bim_out, genome_build or "hg19")
            
            user_data_check = input("\n🔍 Do you have a user genotype file to calculate PGS? (y/n): ").strip().lower()
            if user_data_check == 'y':
                user_path = input("   File path (.bim, .bed, or .vcf/.vcf.gz): ").strip()
                
                if not os.path.exists(user_path):
                    print("❌ File not found!")
                    return

                if user_path.endswith('.vcf') or user_path.endswith('.vcf.gz'):
                    if not PLINK_AVAILABLE:
                        print("❌ PLINK is not installed. VCF conversion is impossible!")
                        return
                    
                    print("\n🔄 VCF file detected. Converting to PLINK binary format...")
                    bfile_prefix = user_path.replace('.vcf.gz', '').replace('.vcf', '')
                    vcf_cmd = [
                        "plink", "--vcf", user_path, 
                        "--make-bed", "--out", bfile_prefix,
                        "--allow-extra-chr"
                    ]
                    try:
                        subprocess.run(vcf_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        print(f"   ✅ VCF converted to: {bfile_prefix}.bed/bim/fam")
                        user_bim_path = f"{bfile_prefix}.bim"
                    except subprocess.CalledProcessError:
                        print("   ❌ VCF conversion error. Check PLINK installation.")
                        return
                else:
                    user_bim_path = user_path.rsplit('.', 1)[0] + '.bim'
                
                if os.path.exists(user_bim_path):
                    light_bim_report(user_bim_path)

                    full_qc = input("   Run QC check (and possible genome Liftover)? (y/n): ").strip().lower()
                    final_bim_path = user_bim_path
                    
                    if full_qc == 'y':
                        final_bim_path = run_full_qc_report(score_out, user_bim_path, genome_build or "hg19")
                        
                    calc_pgs = input(f"\n🚀 Calculate Polygenic Risk Score (PGS) based on {final_bim_path}? (y/n): ").strip().lower()
                    if calc_pgs == 'y':
                        final_bfile = final_bim_path.rsplit('.bim', 1)[0]
                        pgs_out = os.path.join(dir_name, f"{base_name}_PGS_results")
                        
                        pgs_cmd = [
                            "plink",
                            "--bfile", final_bfile,
                            "--score", score_out, "1", "2", "3", "header",
                            "--out", pgs_out
                        ]
                        
                        print("\n🚀 FINAL PLINK COMMAND:")
                        print(" ".join(pgs_cmd))
                        
                        if PLINK_AVAILABLE:
                            print("\n⏳ Calculating PGS...")
                            try:
                                subprocess.run(pgs_cmd, check=True)
                                print(f"✅ Done! Results saved to: {pgs_out}.profile")
                            except subprocess.CalledProcessError:
                                print("❌ Error calculating PGS. Check PLINK logs.")
                        else:
                            print("\n⚠️ PLINK is not installed! The script generated the command above.")
                            print("   Copy and run it manually in the terminal/console once PLINK is installed.")
        else:
            print("✅ Skipping PGS calculation. Only score file was created.")

    else:
        print(f"\n⚠️  rsID not found in the first 10 rows")
        print(f"   Starting alternative method...")
        alternative(filepath, genome_build, skip_lines)

###################################################################################
# MAIN EXECUTION
###################################################################################

def check_plink_installed():
    """Checks if PLINK is available in the command line."""
    try:
        subprocess.run(["plink", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def main():
    global PLINK_AVAILABLE, MAX_ROWS
    print("\n" + "=" * 60)
    print("🧬 PGS → PLINK + BIM CONVERTER v5.0")
    print("   Score file: SNP/A1/BETA")
    print("   BIM file: chr/SNP/cm/pos/A1/A2")
    print("=" * 60)

    PLINK_AVAILABLE = check_plink_installed()
    
    if not PLINK_AVAILABLE:
        print("\n⚠️ WARNING: PLINK not found in the system (not installed or not in PATH).")
        print("   Without it, the script CAN:")
        print("     ✅ Convert tables")
        print("     ✅ Search for rsID via NCBI dbSNP")
        print("     ✅ Create ready-to-use text .txt and .bim files")
        print("   But CANNOT:")
        print("     ❌ Convert VCFs")
        print("     ❌ Perform Liftover on user binary files")
        print("     ❌ Automatically calculate final risk (PGS)")
        
        choice = input("\nContinue in 'File preparation only' mode? (y/n): ").strip().lower()
        if choice != 'y':
            print("Exiting.")
            return

    if not setup_environment():
        print("\n❌ Failed to set up the environment.")
        return

    target_filename = input("\n📁 Enter the file name (PGS score file): ").strip()

    root_path = 'C:\\' if os.name == 'nt' else '/'
    found_filepath = find_file(target_filename, root_path)

    if found_filepath:
        try:
            df_temp = pd.read_csv(found_filepath, sep='\t', comment='#')
            total_rows = len(df_temp)
            print(f"📊 Total rows: {total_rows:,}")
        except:
            pass

        process_table(found_filepath)

if __name__ == "__main__":

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    try:
        main()
    except Exception as e:
        print(f"\n❌ A system error occurred during execution: {e}")
    finally:
        print("\n" + "-" * 60)
        input("Press Enter to close the window...")
