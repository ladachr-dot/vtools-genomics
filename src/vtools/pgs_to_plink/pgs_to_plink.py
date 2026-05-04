# -*- coding: utf-8 -*-
"""
PGS → PLINK Converter (all-in-one)
Single-file GUI application for preparing PGS score files with optional rsID lookup.

Goals:
- one file only
- clean GUI-driven behavior
- no mixed CLI/GUI logic
- predictable outputs
- preserve core functionality where practical
"""

from __future__ import annotations

import concurrent.futures
import io
import os
import queue
import subprocess
import sys
import threading
import time
import traceback
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import requests
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext

# ============================================================================
# Globals
# ============================================================================
PLINK_AVAILABLE = False
LIFTOVER_AVAILABLE = False
APP = None

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
# Stream bridge
# ============================================================================
class GuiOutput(io.TextIOBase):
    def __init__(self, widget: scrolledtext.ScrolledText):
        self.widget = widget
    def write(self, text):
        if text:
            self.widget.after(0, self._append, str(text))
        return len(text)
    def _append(self, text):
        self.widget.configure(state='normal')
        self.widget.insert(tk.END, text)
        self.widget.see(tk.END)
        self.widget.configure(state='disabled')
    def flush(self):
        pass

class GuiInput(io.TextIOBase):
    def __init__(self, app: "App"):
        self.app = app
    def readline(self):
        self.app.after(0, lambda: self.app.set_waiting(True))
        value = self.app.input_queue.get()
        self.app.console.after(0, self.app.echo_input, value)
        self.app.after(0, lambda: self.app.set_waiting(False))
        return value + "\n"
    def flush(self):
        pass

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
            from pyliftover import LiftOver
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
        from pyliftover import LiftOver
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
# GUI
# ============================================================================
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('PGS → PLINK Converter')
        self.geometry('1020x820')
        self.minsize(860, 620)
        self.configure(bg='#1e1e2e')
        self.input_queue = queue.Queue()
        self.worker = None
        self.waiting_for_input = False

        self.folder_var = tk.StringVar()
        self.score_var = tk.StringVar()
        self.rsid_mode_var = tk.StringVar(value='auto')
        self.output_bim_var = tk.BooleanVar(value=False)
        self.keep_checkpoint_var = tk.BooleanVar(value=True)
        self.show_plink_hint_var = tk.BooleanVar(value=False)
        self.start_idx_var = tk.StringVar(value='0')
        self.end_idx_var = tk.StringVar(value='')
        self.status = tk.StringVar(value='Ready.')

        self._build()

    def _row(self, parent):
        f = tk.Frame(parent, bg='#1e1e2e')
        f.pack(fill='x', pady=(0, 6))
        return f
    def _label(self, parent, text, width=18):
        return tk.Label(parent, text=text, bg='#1e1e2e', fg='#cdd6f4', width=width, anchor='w')
    def _btn(self, parent, text, cmd, color='#313244', fg='#cdd6f4'):
        return tk.Button(parent, text=text, command=cmd, bg=color, fg=fg, relief='flat')
    def _entry(self, parent, var, expand=True, width=None):
        kwargs = dict(textvariable=var, bg='#313244', fg='#cdd6f4', insertbackground='white', relief='flat')
        if width:
            kwargs['width'] = width
        e = tk.Entry(parent, **kwargs)
        e.pack(side='left', fill='x' if expand else None, expand=expand)
        return e

    def _build(self):
        top = tk.Frame(self, bg='#181825', padx=10, pady=8)
        top.pack(fill='x')
        tk.Label(top, text='PGS → PLINK Converter', bg='#181825', fg='#cdd6f4', font=('Segoe UI', 13, 'bold')).pack(side='left')
        self._btn(top, 'Restart', self.restart).pack(side='right', padx=4)
        self._btn(top, '▶ Run', self.run_script, color='#89b4fa', fg='#11111b').pack(side='right', padx=4)

        form = tk.Frame(self, bg='#1e1e2e', padx=12, pady=10)
        form.pack(fill='x')

        tk.Label(form, text='① INPUT', bg='#1e1e2e', fg='#89b4fa', font=('Segoe UI', 9, 'bold')).pack(anchor='w', pady=(0, 4))
        r1 = self._row(form)
        self._label(r1, 'Folder with files:').pack(side='left')
        self._entry(r1, self.folder_var)
        self._btn(r1, 'Browse folder', self.choose_folder).pack(side='left', padx=(6,0))

        r2 = self._row(form)
        self._label(r2, 'Select file:').pack(side='left')
        self.file_combo = tk.OptionMenu(r2, self.score_var, '')
        self.file_combo.configure(bg='#313244', fg='#cdd6f4', relief='flat', activebackground='#45475a', highlightthickness=0)
        self.file_combo['menu'].configure(bg='#313244', fg='#cdd6f4')
        self.file_combo.pack(side='left', fill='x', expand=True)
        self._btn(r2, 'Refresh', self.refresh_file_list).pack(side='left', padx=(6,0))

        tk.Frame(form, bg='#313244', height=1).pack(fill='x', pady=8)

        tk.Label(form, text='② OPTIONS', bg='#1e1e2e', fg='#89b4fa', font=('Segoe UI', 9, 'bold')).pack(anchor='w', pady=(0, 4))
        r3 = self._row(form)
        self._label(r3, 'rsID mode:').pack(side='left')
        tk.Radiobutton(r3, text='Auto', value='auto', variable=self.rsid_mode_var, bg='#1e1e2e', fg='#cdd6f4', selectcolor='#313244', activebackground='#1e1e2e', activeforeground='#cdd6f4').pack(side='left')
        tk.Radiobutton(r3, text='Interactive', value='interactive', variable=self.rsid_mode_var, bg='#1e1e2e', fg='#cdd6f4', selectcolor='#313244', activebackground='#1e1e2e', activeforeground='#cdd6f4').pack(side='left', padx=(12,0))

        r4 = self._row(form)
        self._label(r4, 'Start index:').pack(side='left')
        self._entry(r4, self.start_idx_var, expand=False, width=10)
        tk.Label(r4, text='End index:', bg='#1e1e2e', fg='#cdd6f4').pack(side='left', padx=(14,6))
        self._entry(r4, self.end_idx_var, expand=False, width=10)
        tk.Label(r4, text='(empty = all rows)', bg='#1e1e2e', fg='#a6adc8').pack(side='left', padx=(12,0))

        r5 = self._row(form)
        tk.Checkbutton(r5, text='Create BIM file', variable=self.output_bim_var, bg='#1e1e2e', fg='#cdd6f4', selectcolor='#313244', activebackground='#1e1e2e', activeforeground='#cdd6f4').pack(side='left', padx=(0,18))
        tk.Checkbutton(r5, text='Keep checkpoint file', variable=self.keep_checkpoint_var, bg='#1e1e2e', fg='#cdd6f4', selectcolor='#313244', activebackground='#1e1e2e', activeforeground='#cdd6f4').pack(side='left', padx=(0,18))
        tk.Checkbutton(r5, text='Show PLINK command hint', variable=self.show_plink_hint_var, bg='#1e1e2e', fg='#cdd6f4', selectcolor='#313244', activebackground='#1e1e2e', activeforeground='#cdd6f4').pack(side='left')

        self.console = scrolledtext.ScrolledText(self, bg='#11111b', fg='#cdd6f4', insertbackground='white', font=('Consolas', 10), state='disabled', wrap='word')
        self.console.pack(fill='both', expand=True, padx=10, pady=(0,6))

        bottom = tk.Frame(self, bg='#1e1e2e', padx=10, pady=6)
        bottom.pack(fill='x')
        tk.Label(bottom, text='Interactive input:', bg='#1e1e2e', fg='#89b4fa').pack(side='left', padx=(0,8))
        self.entry = tk.Entry(bottom, bg='#313244', fg='#cdd6f4', insertbackground='white', font=('Consolas', 11), relief='flat')
        self.entry.pack(side='left', fill='x', expand=True, ipady=5)
        self.entry.bind('<Return>', self.send_input)
        self._btn(bottom, 'Send', self.send_input, color='#89b4fa', fg='#11111b').pack(side='left', padx=(6,0))

        tk.Label(self, textvariable=self.status, anchor='w', bg='#181825', fg='#a6adc8', padx=10, pady=4).pack(fill='x', side='bottom')

    def choose_folder(self):
        folder = filedialog.askdirectory(title='Select folder containing PGS files')
        if folder:
            self.folder_var.set(os.path.abspath(folder))
            self.refresh_file_list()
            self.status.set(f'Folder selected: {folder}')

    def refresh_file_list(self):
        folder = self.folder_var.get().strip()
        if not folder or not os.path.isdir(folder):
            messagebox.showwarning('No folder', 'Select a valid folder first.')
            return
        files = []
        skip_markers = ['_plink', '_checkpoint']
        skip_exts = {'.bim', '.bed', '.fam', '.log', '.py', '.nosex', '.profile'}
        for root, dirs, names in os.walk(folder):
            dirs[:] = [d for d in dirs if not d.startswith('.')]
            for name in names:
                low = name.lower()
                full = os.path.abspath(os.path.join(root, name))
                if any(m in low for m in skip_markers):
                    continue
                if low.endswith('.chain.gz'):
                    continue
                if low in {'license', 'requirements.txt'}:
                    continue
                ext = os.path.splitext(low)[1]
                if ext in skip_exts:
                    continue
                if low.endswith('.txt') or low.endswith('.tsv') or low.endswith('.csv') or low.endswith('.txt.gz') or low.endswith('.tsv.gz') or low.endswith('.csv.gz'):
                    files.append(full)
        files = sorted(set(files))
        menu = self.file_combo['menu']
        menu.delete(0, 'end')
        if not files:
            self.score_var.set('')
            messagebox.showinfo('No files found', f'No score files found in:\n{folder}')
            return
        for f in files:
            label = os.path.relpath(f, folder)
            menu.add_command(label=label, command=lambda v=f: self.score_var.set(v))
        self.score_var.set(files[0])
        self.status.set(f'{len(files)} score file(s) found.')

    def set_waiting(self, waiting):
        self.waiting_for_input = waiting
        self.status.set('⏳ Waiting for interactive input...' if waiting else 'Running...')

    def echo_input(self, value):
        self.console.configure(state='normal')
        self.console.insert(tk.END, f'>>> {value}\n')
        self.console.see(tk.END)
        self.console.configure(state='disabled')

    def send_input(self, event=None):
        value = self.entry.get().strip()
        if not value:
            return
        self.input_queue.put(value)
        self.entry.delete(0, tk.END)
        self.status.set('Interactive input sent.' if self.waiting_for_input else 'Input queued.')

    def restart(self):
        self.console.configure(state='normal')
        self.console.delete('1.0', tk.END)
        self.console.configure(state='disabled')
        while not self.input_queue.empty():
            try:
                self.input_queue.get_nowait()
            except Exception:
                break
        self.status.set('Ready.')

    def run_script(self):
        score_path = self.score_var.get().strip().strip('"')
        if not score_path:
            messagebox.showerror('Missing file', 'Select a PGS score file first.')
            return
        score_path = os.path.abspath(score_path)
        if not os.path.isfile(score_path):
            messagebox.showerror('File not found', f'File does not exist:\n{score_path}')
            return
        if self.start_idx_var.get().strip() and not self.start_idx_var.get().strip().isdigit():
            messagebox.showerror('Invalid input', 'Start index must be an integer.')
            return
        if self.end_idx_var.get().strip() and not self.end_idx_var.get().strip().isdigit():
            messagebox.showerror('Invalid input', 'End index must be an integer.')
            return

        cfg = RunConfig(
            score_file=score_path,
            output_bim=self.output_bim_var.get(),
            rsid_mode=self.rsid_mode_var.get(),
            start_idx=int(self.start_idx_var.get().strip() or '0'),
            end_idx=int(self.end_idx_var.get().strip()) if self.end_idx_var.get().strip() else None,
            keep_checkpoint=self.keep_checkpoint_var.get(),
            show_plink_hint=self.show_plink_hint_var.get(),
        )

        self.console.configure(state='normal')
        self.console.delete('1.0', tk.END)
        self.console.configure(state='disabled')
        self.status.set('Running...')

        def target():
            old_out, old_err, old_in = sys.stdout, sys.stderr, sys.stdin
            sys.stdout = GuiOutput(self.console)
            sys.stderr = GuiOutput(self.console)
            sys.stdin = GuiInput(self)
            try:
                global PLINK_AVAILABLE
                PLINK_AVAILABLE = check_plink_installed()
                if not setup_environment():
                    print('\n❌ Failed to set up the environment.')
                    return
                process_score_file(score_path, cfg)
                print('\n✅ Finished.')
            except Exception:
                traceback.print_exc()
            finally:
                sys.stdout, sys.stderr, sys.stdin = old_out, old_err, old_in
                self.after(0, lambda: self.status.set('Finished.'))

        self.worker = threading.Thread(target=target, daemon=True)
        self.worker.start()

if __name__ == '__main__':
    App().mainloop()
