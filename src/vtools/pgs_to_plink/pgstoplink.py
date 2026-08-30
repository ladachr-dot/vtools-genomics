#!/usr/bin/env python3
"""
PGS -> PLINK Converter (CLI)
Terminal-only tool for preparing PGS score files with optional rsID lookup.

Usage examples:
    python pgs_to_plink.py score.txt
    python pgs_to_plink.py score.txt --bim --rsid-mode interactive
    python pgs_to_plink.py score.txt --start-idx 0 --end-idx 5000 --no-keep-checkpoint
    python pgs_to_plink.py score.txt --plink-hint
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import requests

# ============================================================================
# Globals
# ============================================================================
PLINK_AVAILABLE = False
LIFTOVER_AVAILABLE = False
LiftOver = None  # populated once pyliftover is imported

# ============================================================================
# Config
# ============================================================================
@dataclass
class RunConfig:
    score_file: str
    output_bim: bool = False
    rsid_mode: str = "auto"
    start_idx: int = 0
    end_idx: int | None = None
    keep_checkpoint: bool = True
    show_plink_hint: bool = False

# ============================================================================
# Column mapping
# ============================================================================
COLUMN_SYNONYMS = {
    "SNP": ["rsid", "rs_id", "snp", "variant_id", "marker", "markername", "id", "hm_rsid"],
    "A1": ["effect_allele", "ea", "a1", "tested_allele", "allele", "risk_allele", "hm_effect_allele"],
    "BETA": ["beta", "effect_weight", "weight", "score", "log_odds", "logor", "hm_beta"],
}

BIM_COLUMN_SYNONYMS = {
    "chr_name": ["chr_name", "chromosome", "chr", "chrom"],
    "chr_position": ["chr_position", "position", "pos", "hm_pos", "bp"],
    "A2": ["other_allele", "a2", "reference_allele", "ref", "hm_other_allele"],
}

# ============================================================================
# Helpers
# ============================================================================
def normalize_column_name(name: str) -> str:
    return str(name).strip().lower().replace('-', '_').replace(' ', '_')


def find_column(columns, synonyms):
    normalized_columns = {normalize_column_name(col): col for col in columns}
    for synonym in synonyms:
        s = normalize_column_name(synonym)
        if s in normalized_columns:
            return normalized_columns[s]
    return None


def extract_genome_build(file_path: str) -> str:
    path = Path(file_path)
    if not path.exists():
        return "unknown"
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#genome_build'):
                return line.split('=', 1)[1].strip()
            if line and not line.startswith('#'):
                break
    return 'unknown'


def check_plink_installed() -> bool:
    try:
        subprocess.run(["plink", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception:
        return False


def install_package(package_name: str) -> bool:
    print(f"📦 Installing {package_name}...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package_name, "--quiet", "--disable-pip-version-check"])
        print(f"   ✅ {package_name} installed")
        return True
    except subprocess.CalledProcessError:
        print(f"   ❌ Failed to install {package_name}")
        return False


def check_and_install_dependencies() -> bool:
    print("=" * 60)
    print("🔧 CHECKING DEPENDENCIES")
    print("=" * 60)
    dependencies = {
        'requests': 'requests',
        'pandas': 'pandas',
        'pyliftover': 'pyliftover',
    }
    ok = True
    for import_name, package_name in dependencies.items():
        try:
            __import__(import_name)
            print(f"✅ {package_name} is already installed")
        except ImportError:
            print(f"⚠️ {package_name} not found")
            if not install_package(package_name):
                ok = False
    if ok:
        try:
            global LIFTOVER_AVAILABLE, LiftOver
            from pyliftover import LiftOver as _LiftOver
            LiftOver = _LiftOver
            LIFTOVER_AVAILABLE = True
            print("✅ pyliftover is ready")
        except ImportError:
            print("❌ Failed to load pyliftover")
            ok = False
    print()
    return ok


def chain_path(chain_name: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), chain_name)


def download_liftover_chain(chain_name: str) -> bool:
    path = chain_path(chain_name)
    if 'hg19ToHg38' in chain_name:
        url = f"http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/{chain_name}"
    elif 'hg38ToHg19' in chain_name:
        url = f"http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/{chain_name}"
    else:
        print(f"❌ Unknown chain: {chain_name}")
        return False
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"✅ {chain_name} is already available")
        return True
    print(f"📥 Downloading {chain_name}...")
    try:
        r = requests.get(url, timeout=120, stream=True)
        r.raise_for_status()
        with open(path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"   ✅ Downloaded: {os.path.getsize(path):,} bytes")
        return True
    except Exception as e:
        print(f"   ❌ Download error: {e}")
        if os.path.exists(path):
            os.remove(path)
        return False


def setup_environment() -> bool:
    global LIFTOVER_AVAILABLE, LiftOver
    env_ready = True
    try:
        __import__('requests')
        __import__('pandas')
        __import__('pyliftover')
    except ImportError:
        env_ready = False
    for chain in ['hg19ToHg38.over.chain.gz', 'hg38ToHg19.over.chain.gz']:
        p = chain_path(chain)
        if not os.path.exists(p) or os.path.getsize(p) < 1000:
            env_ready = False
    if env_ready:
        from pyliftover import LiftOver as _LiftOver
        LiftOver = _LiftOver
        LIFTOVER_AVAILABLE = True
        return True

    print("\n" + "=" * 60)
    print("🚀 SETTING UP ENVIRONMENT")
    print("=" * 60 + "\n")
    try:
        subprocess.check_call([sys.executable, '-m', 'pip', '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("📋 pip is working")
    except Exception:
        print("❌ pip is not working")
        return False
    if not check_and_install_dependencies():
        return False
    for chain in ['hg19ToHg38.over.chain.gz', 'hg38ToHg19.over.chain.gz']:
        if not download_liftover_chain(chain):
            return False
    print("✅ Environment is ready\n")
    return True

# ============================================================================
# Core conversion helpers
# ============================================================================
def read_score_table(input_path: str) -> pd.DataFrame:
    path = Path(str(input_path).strip().strip('"'))
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    try:
        df = pd.read_csv(path, sep='\t', comment='#', engine='python', dtype=str)
        if len(df.columns) <= 1:
            raise ValueError('TSV parse produced one column')
        return df
    except Exception:
        try:
            df = pd.read_csv(path, sep=r'\s+', comment='#', engine='python', dtype=str)
            if len(df.columns) <= 1:
                raise ValueError('Whitespace parse produced one column')
            return df
        except Exception:
            return pd.read_csv(path, sep=None, comment='#', engine='python', dtype=str)


def convert_score_file(input_path: str, output_path: str):
    genome_build = extract_genome_build(input_path)
    df = read_score_table(input_path)
    print(f"\n📄 Input file: {input_path}")
    print(f"📊 Original columns: {list(df.columns)}")
    print(f"📊 Parsed column count: {len(df.columns)}")

    snp_col = find_column(df.columns, COLUMN_SYNONYMS['SNP'])
    a1_col = find_column(df.columns, COLUMN_SYNONYMS['A1'])
    beta_col = find_column(df.columns, COLUMN_SYNONYMS['BETA'])

    missing = []
    if snp_col is None:
        missing.append('SNP')
    if a1_col is None:
        missing.append('A1')
    if beta_col is None:
        missing.append('BETA')
    if missing:
        raise ValueError(f"Required columns not found: {missing}; available: {list(df.columns)}")

    out_df = pd.DataFrame({
        'SNP': df[snp_col].astype(str),
        'A1': df[a1_col].astype(str).str.upper(),
        'BETA': pd.to_numeric(df[beta_col], errors='coerce')
    })
    out_df = out_df.dropna(subset=['SNP', 'A1', 'BETA'])
    out_df = out_df[out_df['SNP'] != '.']
    out_df = out_df.drop_duplicates(subset=['SNP', 'A1'])
    out_df.to_csv(output_path, sep='\t', index=False)
    print(f"✅ Score file created: {output_path}")
    return out_df, df, genome_build


def create_bim_from_score(df: pd.DataFrame, output_path: str, genome_build: str = 'unknown'):
    print(f"\n📝 Creating BIM file...")
    chr_col = find_column(df.columns, BIM_COLUMN_SYNONYMS['chr_name'])
    pos_col = find_column(df.columns, BIM_COLUMN_SYNONYMS['chr_position'])
    snp_col = find_column(df.columns, COLUMN_SYNONYMS['SNP'])
    a1_col = find_column(df.columns, COLUMN_SYNONYMS['A1'])
    a2_col = find_column(df.columns, BIM_COLUMN_SYNONYMS['A2'])
    if a2_col is None:
        for col in df.columns:
            low = col.lower()
            if 'other_allele' in low or low == 'a2':
                a2_col = col
                break

    missing = []
    if chr_col is None:
        missing.append('chr_name')
    if pos_col is None:
        missing.append('chr_position')
    if snp_col is None:
        missing.append('SNP/rsID')
    if a1_col is None:
        missing.append('effect_allele/A1')
    if missing:
        print(f"⚠️ BIM was not created; missing columns: {missing}")
        return None

    bim_df = pd.DataFrame()
    bim_df['chr'] = df[chr_col].astype(str).str.replace('chr', '', case=False, regex=False)
    bim_df['SNP'] = df[snp_col].astype(str)
    bim_df['cm'] = 0
    bim_df['pos'] = pd.to_numeric(df[pos_col], errors='coerce')
    bim_df['A1'] = df[a1_col].astype(str).str.upper()
    bim_df['A2'] = df[a2_col].astype(str).str.upper() if a2_col else '.'
    before = len(bim_df)
    bim_df = bim_df.dropna(subset=['pos'])
    bim_df = bim_df[bim_df['chr'] != '.']
    bim_df = bim_df[bim_df['SNP'] != '.']
    bim_df = bim_df[bim_df['A1'] != '.']
    bim_df = bim_df.drop_duplicates(subset=['SNP'], keep='first')
    after = len(bim_df)
    bim_df.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"✅ BIM file created: {output_path}")
    print(f"Rows: {after} (removed: {before - after})")
    print(f"Build: {genome_build}")
    return output_path

# ============================================================================
# rsID / liftover
# ============================================================================
def liftover_coordinates(df: pd.DataFrame, from_build: str, to_build: str, pos_col: str = 'chr_position'):
    def normalize(build):
        b = str(build).lower().replace('grch', '').replace('hg', '')
        if '37' in b or '19' in b:
            return '37'
        if '38' in b:
            return '38'
        return b

    from_b = normalize(from_build)
    to_b = normalize(to_build)
    if from_b == to_b:
        print(f"✅ Builds are identical (hg{from_b}). Conversion not required.")
        return df, True
    if from_b == '37' and to_b == '38':
        chain_file = 'hg19ToHg38.over.chain.gz'
        direction = 'hg19 → hg38'
    elif from_b == '38' and to_b == '37':
        chain_file = 'hg38ToHg19.over.chain.gz'
        direction = 'hg38 → hg19'
    else:
        print(f"❌ Unsupported conversion: {from_build} → {to_build}")
        return df, False

    path = chain_path(chain_file)
    if not os.path.exists(path):
        print(f"❌ Chain file not found: {path}")
        return df, False

    print(f"🔄 Converting {direction}...")
    lo = LiftOver(path)
    new_coords = []
    for _, row in df.iterrows():
        chrom = str(row.get('chr_name', '')).replace('chr', '').replace('CHR', '')
        pos_val = row.get(pos_col, '')
        try:
            pos = int(float(pos_val))
            result = lo.convert_coordinate(f'chr{chrom}', pos - 1)
            new_coords.append(str(result[0][1] + 1) if result else '.')
        except Exception:
            new_coords.append('.')
    out = df.copy()
    new_col = 'chr_position_hg38' if to_b == '38' else 'chr_position_hg19'
    out[new_col] = new_coords
    converted = sum(1 for x in new_coords if x != '.')
    print(f"✅ Coordinate conversion: {converted}/{len(new_coords)} rows converted")
    return out, True


def search_rsid_candidates(chrom: str, pos: str):
    if not chrom or not pos or pos == '.':
        return []

    headers = {'User-Agent': 'Mozilla/5.0'}
    candidates = []
    seen = set()

    urls = [
        f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/?chromosome={chrom}&position={pos}",
        f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/?chromosome={chrom}&position={pos}",
        f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{chrom}/{pos}",
    ]

    for url in urls:
        try:
            r = requests.get(url, headers=headers, timeout=20)
            if r.status_code != 200:
                continue
            data = r.json()
            refs = []
            if isinstance(data, dict):
                if isinstance(data.get('refsnp_set'), list):
                    refs.extend(data['refsnp_set'])
                elif isinstance(data.get('data', {}).get('refsnp_set'), list):
                    refs.extend(data['data']['refsnp_set'])
                elif 'primary_snapshot_data' in data and data.get('refsnp_id'):
                    refs.append(data)
            for item in refs:
                rsid = item.get('refsnp_id') or item.get('id')
                if rsid:
                    rs = f"rs{rsid}"
                    if rs not in seen:
                        seen.add(rs)
                        candidates.append(rs)
        except Exception:
            continue

    if candidates:
        return candidates

    try:
        term = f"{chrom}[CHR] AND {pos}[CPOS]"
        r = requests.get(
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
            params={'db': 'snp', 'retmode': 'json', 'retmax': 200, 'term': term},
            headers=headers,
            timeout=20,
        )
        if r.status_code == 200:
            data = r.json()
            ids = data.get('esearchresult', {}).get('idlist', [])
            for snp_id in ids:
                rs = f"rs{snp_id}"
                if rs not in seen:
                    seen.add(rs)
                    candidates.append(rs)
    except Exception:
        pass

    return candidates


def choose_rsid(chrom: str, pos_hg38: str, auto_mode: bool = True):
    candidates = search_rsid_candidates(chrom, pos_hg38)

    if auto_mode:
        if candidates:
            print(f"✅ chr{chrom}:{pos_hg38} -> {candidates[0]} (auto, found {len(candidates)} candidate(s))")
            return candidates[0]
        print(f"⚠️ chr{chrom}:{pos_hg38} -> no rsID candidates found")
        return '.'

    print(f"\n📍 chr{chrom}:{pos_hg38}")
    if candidates:
        if len(candidates) == 1:
            print(f"   Single candidate found: {candidates[0]} (selected automatically in interactive mode)")
            return candidates[0]
        print(f"   Found {len(candidates)} candidate rsIDs:")
        for i, rs in enumerate(candidates, 1):
            print(f"   [{i}] {rs}")
        print("   [a] automatic mode for the rest")
        print("   [0] skip")
        while True:
            choice = input("Your choice: ").strip().lower()
            if choice == 'a':
                return 'AUTO'
            if choice == '0' or choice == '':
                return '.'
            if choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(candidates):
                    return candidates[idx]
            if choice.startswith('rs'):
                return choice
            print("Invalid choice. Enter a number from the list, 0, a, or type an rsID.")
    else:
        print("   No rsID candidates found for this coordinate.")
        print("   [Enter/0] skip   [a] automatic mode for the rest   [m] manual rsID")
        while True:
            choice = input("Your choice: ").strip().lower()
            if choice in {'', '0'}:
                return '.'
            if choice == 'a':
                return 'AUTO'
            if choice == 'm':
                manual = input("Enter rsID manually (example: rs12345): ").strip()
                if manual:
                    return manual
                return '.'
            if choice.startswith('rs'):
                return choice
            print("Invalid choice. Use Enter, 0, a, m, or type an rsID.")

# ============================================================================
# Main pipeline
# ============================================================================
def score_has_rsid(filepath: str, max_rows: int = 10) -> bool:
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#', nrows=max_rows, dtype=str)
        for col in df.columns:
            series = df[col].astype(str).str.lower()
            if series.str.startswith('rs').any():
                return True
    except Exception:
        return False
    return False


def prepare_with_rsid(filepath: str, config: RunConfig):
    dir_name = os.path.dirname(filepath)
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    score_out = os.path.join(dir_name, f"{base_name}_plink.txt")
    converted_df, full_df, build = convert_score_file(filepath, score_out)
    print(f"\n✅ Preparation finished.")
    if config.output_bim:
        bim_out = os.path.join(dir_name, f"{base_name}.bim")
        create_bim_from_score(full_df, bim_out, build)
    if config.show_plink_hint:
        print("\n💡 PLINK command example:")
        print(f"   plink --bfile your_data --score {score_out} 1 2 3 header")


def prepare_without_rsid(filepath: str, config: RunConfig):
    genome_build = extract_genome_build(filepath)
    print("\n⚠️ rsID not found in the first rows")
    print("   Starting rsID recovery workflow...")
    df_full = pd.read_csv(filepath, sep='\t', comment='#', dtype=str)
    total_rows = len(df_full)
    start_idx = max(0, int(config.start_idx or 0))
    end_idx = config.end_idx if config.end_idx is not None else total_rows
    end_idx = min(end_idx, total_rows)
    if start_idx >= end_idx:
        raise ValueError(f"Invalid range: start={start_idx}, end={end_idx}, total={total_rows}")

    print(f"\n📊 Total rows in file: {total_rows}")
    print(f"📊 Selected row range [{start_idx}:{end_idx}]. To process: {end_idx - start_idx}")

    df = df_full.iloc[start_idx:end_idx].copy().reset_index(drop=True)
    if 'chr_name' not in df.columns or 'chr_position' not in df.columns:
        raise ValueError("chr_name/chr_position columns are missing")

    checkpoint_file = os.path.join(os.path.dirname(filepath), f"{Path(filepath).stem}_checkpoint_{start_idx}_{end_idx}.csv")
    if os.path.exists(checkpoint_file):
        print(f"♻️ Found checkpoint: {checkpoint_file}")
        print("   Resuming saved progress...")
        df = pd.read_csv(checkpoint_file, sep='\t', dtype=str)
    else:
        if 'rsID' not in df.columns:
            df['rsID'] = '.'
        df['chr_position_hg19_original'] = df['chr_position']
        print("\n🔹 Step 1/3: Coordinate conversion")
        df, ok = liftover_coordinates(df, '37', '38', 'chr_position')
        if not ok:
            raise RuntimeError("Coordinate conversion failed")
        df.to_csv(checkpoint_file, sep='\t', index=False)
        print(f"💾 Checkpoint saved: {checkpoint_file}")

    unprocessed = df[df['rsID'] == '.'].index.tolist()
    print("\n🔹 Step 2/3: rsID search")
    print(f"   Remaining rows: {len(unprocessed)}")
    auto_mode = config.rsid_mode != 'interactive'
    print(f"   Mode: {'automatic' if auto_mode else 'interactive'}")

    processed_in_session = 0
    for i in unprocessed:
        chrom = str(df.at[i, 'chr_name'])
        pos_hg38 = str(df.at[i, 'chr_position_hg38'])
        rsid = choose_rsid(chrom, pos_hg38, auto_mode=auto_mode)
        if rsid == 'AUTO':
            print("⚡ Switching to automatic mode for remaining rows")
            auto_mode = True
            rsid = choose_rsid(chrom, pos_hg38, auto_mode=True)
        df.at[i, 'rsID'] = rsid
        processed_in_session += 1
        if rsid and rsid != '.':
            print(f"✅ Row {i}: {rsid}")
        if processed_in_session % 25 == 0:
            df.to_csv(checkpoint_file, sep='\t', index=False)
            print(f"💾 Checkpoint updated ({processed_in_session} rows processed this session)")

    df.to_csv(checkpoint_file, sep='\t', index=False)
    print("✅ rsID search complete")

    print("\n🔹 Step 3/3: Output files")
    weight_col = next((c for c in df.columns if 'weight' in c.lower()), None)
    allele_col = next((c for c in df.columns if 'effect_allele' in c.lower()), None)
    if not weight_col or not allele_col:
        raise ValueError("Allele/weight columns not found after rsID recovery")

    dir_name = os.path.dirname(filepath)
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    score_out = os.path.join(dir_name, f"{base_name}_plink_hg19_{start_idx}_{end_idx}.txt")
    plink_df = df[df['rsID'] != '.'].copy()
    plink_df = plink_df[['rsID', allele_col, weight_col]].copy()
    plink_df[weight_col] = pd.to_numeric(plink_df[weight_col], errors='coerce')
    plink_df = plink_df.dropna()
    plink_df = plink_df.drop_duplicates(subset=['rsID', allele_col])
    plink_df.to_csv(score_out, sep='\t', index=False, header=False)
    print(f"✅ Score file created: {score_out}")
    print(f"📊 Output rows with rsID: {len(plink_df)}")

    if config.output_bim:
        bim_out = os.path.join(dir_name, f"{base_name}_plink_hg19_{start_idx}_{end_idx}.bim")
        create_bim_from_score(df, bim_out, genome_build or 'unknown')

    if config.keep_checkpoint:
        print(f"💾 Checkpoint kept: {checkpoint_file}")
    else:
        try:
            os.remove(checkpoint_file)
            print(f"🧹 Checkpoint removed: {checkpoint_file}")
        except Exception as e:
            print(f"⚠️ Could not remove checkpoint: {e}")

    if config.show_plink_hint:
        print("\n💡 PLINK command example:")
        print(f"   plink --bfile your_data --score {score_out} 1 2 3 header")


def process_score_file(filepath: str, config: RunConfig):
    print("\n" + "=" * 60)
    print("🧬 PGS → PLINK Converter")
    print("=" * 60)
    print(f"📂 File: {filepath}")
    print(f"⚙️  Mode: {'interactive rsID recovery' if config.rsid_mode == 'interactive' else 'automatic rsID recovery'}")
    print(f"⚙️  Range: start={config.start_idx}, end={config.end_idx if config.end_idx is not None else 'ALL'}")
    print(f"⚙️  Create BIM: {'yes' if config.output_bim else 'no'}")
    print(f"⚙️  Keep checkpoint: {'yes' if config.keep_checkpoint else 'no'}")
    print(f"⚙️  Show PLINK hint: {'yes' if config.show_plink_hint else 'no'}")

    if score_has_rsid(filepath):
        print("\nℹ️ rsID detected in input file")
        prepare_with_rsid(filepath, config)
    else:
        print("\nℹ️ rsID not detected in input file")
        prepare_without_rsid(filepath, config)

# ============================================================================
# CLI entry point
# ============================================================================
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a PGS score file into PLINK-ready format, with optional rsID lookup and BIM creation.",
    )
    parser.add_argument("score_file", nargs="?", default=None,
                         help="Path to the input PGS score file (.txt/.tsv/.csv, optionally .gz). "
                              "If omitted, the script will ask for it interactively.")
    parser.add_argument("--rsid-mode", choices=["auto", "interactive"], default="auto",
                         help="How to resolve ambiguous rsID candidates (default: auto)")
    parser.add_argument("--start-idx", type=int, default=0, help="First row index to process (default: 0)")
    parser.add_argument("--end-idx", type=int, default=None, help="Last row index (exclusive) to process (default: all rows)")
    parser.add_argument("--bim", dest="output_bim", action="store_true", help="Also create a .bim file")
    parser.add_argument("--no-keep-checkpoint", dest="keep_checkpoint", action="store_false",
                         help="Delete the checkpoint file after a successful run")
    parser.add_argument("--plink-hint", dest="show_plink_hint", action="store_true",
                         help="Print an example PLINK command after finishing")
    parser.set_defaults(keep_checkpoint=True)
    return parser


def prompt_for_score_file() -> str:
    """Interactively ask the user for a score file when none was given on the command line."""
    while True:
        raw = input("Enter path to the PGS score file: ").strip().strip('"')
        if not raw:
            print("Please enter a path.")
            continue
        candidate = os.path.abspath(raw)
        if os.path.isfile(candidate):
            return candidate
        print(f"❌ File not found: {candidate}")


def prompt_yes_no(question: str, default: bool = False) -> bool:
    suffix = " [Y/n]: " if default else " [y/N]: "
    while True:
        raw = input(question + suffix).strip().lower()
        if not raw:
            return default
        if raw in ("y", "yes"):
            return True
        if raw in ("n", "no"):
            return False
        print("Please answer y or n.")


def prompt_for_config_interactively(score_path: str) -> RunConfig:
    """Ask the user for run options one by one, with sensible defaults on Enter."""
    print("\nNo command-line options were given — let's configure the run interactively.")
    print("(Press Enter at any prompt to accept the default shown in brackets.)\n")

    output_bim = prompt_yes_no("Create a .bim file too?", default=False)

    mode_raw = input("rsID resolution mode - 'auto' or 'interactive' [auto]: ").strip().lower()
    rsid_mode = mode_raw if mode_raw in ("auto", "interactive") else "auto"

    start_raw = input("Start row index [0]: ").strip()
    start_idx = int(start_raw) if start_raw.isdigit() else 0

    end_raw = input("End row index, exclusive (Enter = all rows): ").strip()
    end_idx = int(end_raw) if end_raw.isdigit() else None

    keep_checkpoint = prompt_yes_no("Keep the checkpoint file after a successful run?", default=True)
    show_plink_hint = prompt_yes_no("Show an example PLINK command at the end?", default=False)

    return RunConfig(
        score_file=score_path,
        output_bim=output_bim,
        rsid_mode=rsid_mode,
        start_idx=start_idx,
        end_idx=end_idx,
        keep_checkpoint=keep_checkpoint,
        show_plink_hint=show_plink_hint,
    )


def pause_before_exit():
    """Keep the console window open when the script was double-clicked instead of run from a shell."""
    try:
        input("\nPress Enter to close this window...")
    except (EOFError, KeyboardInterrupt):
        pass


def main(argv=None):
    parser = build_arg_parser()
    # score_file is optional now so the script can be double-clicked and ask interactively
    args, unknown = parser.parse_known_args(argv)

    if args.score_file:
        score_path = os.path.abspath(args.score_file.strip().strip('"'))
        if not os.path.isfile(score_path):
            print(f"❌ File not found: {score_path}")
            pause_before_exit()
            return 1
        # Options were (or could have been) given on the command line
        cfg = RunConfig(
            score_file=score_path,
            output_bim=args.output_bim,
            rsid_mode=args.rsid_mode,
            start_idx=args.start_idx,
            end_idx=args.end_idx,
            keep_checkpoint=args.keep_checkpoint,
            show_plink_hint=args.show_plink_hint,
        )
    else:
        # No arguments at all (e.g. the file was double-clicked) -> go fully interactive
        score_path = prompt_for_score_file()
        cfg = prompt_for_config_interactively(score_path)

    global PLINK_AVAILABLE
    PLINK_AVAILABLE = check_plink_installed()
    if PLINK_AVAILABLE:
        print("✅ plink is installed and available in PATH")
    else:
        print("ℹ️ plink was not found in PATH (only needed to actually run the scoring step)")

    exit_code = 0
    try:
        if not setup_environment():
            print("\n❌ Failed to set up the environment.")
            exit_code = 1
        else:
            process_score_file(score_path, cfg)
            print("\n✅ Finished.")
    except KeyboardInterrupt:
        print("\n⏹️ Interrupted by user.")
        exit_code = 1
    except Exception:
        import traceback
        traceback.print_exc()
        exit_code = 1

    pause_before_exit()
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
