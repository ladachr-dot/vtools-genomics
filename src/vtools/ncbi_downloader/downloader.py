from __future__ import annotations

import argparse
import json
import re
import sys
import time
import traceback
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def load_checkpoint(path: Path) -> dict:
    if not path.exists():
        return {"completed": [], "failed": [], "meta": {"processed": 0, "total": 0, "last_id": None}}
    return json.loads(path.read_text(encoding="utf-8"))


def save_checkpoint(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    low = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


def extract_rsid(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    m = re.search(r"(rs\d+)", str(value), flags=re.IGNORECASE)
    return m.group(1).lower() if m else None


class NCBIDownloader:
    def __init__(self, email: str | None = None, api_key: str | None = None, timeout_sec: int = 30):
        self.email = email
        self.api_key = api_key
        self.timeout_sec = timeout_sec
        self.session = requests.Session()
        retries = Retry(total=5, backoff_factor=1.2, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
        self.session.mount("https://", HTTPAdapter(max_retries=retries))

    def fetch_snp_docsum(self, rsid: str) -> dict:
        params = {"db": "snp", "id": rsid, "retmode": "json", "version": "2.0"}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        r = self.session.get(f"{NCBI_EUTILS}/esummary.fcgi", params=params, timeout=self.timeout_sec)
        r.raise_for_status()
        return r.json()


def parse_ncbi(payload: dict, rsid: str) -> dict:
    result = payload.get("result", {}) if isinstance(payload, dict) else {}
    block = None
    for key, value in result.items():
        if key == "uids":
            continue
        if isinstance(value, dict):
            block = value
            break
    if not block:
        return {"query_rsid": rsid, "ncbi_uid": None, "ncbi_title": None, "ncbi_chr": None, "ncbi_position": None, "ncbi_gene": None, "ncbi_taxid": None}
    genomic = (block.get("genomicinfo") or [{}])[0]
    genes = (block.get("genes") or [{}])[0]
    return {
        "query_rsid": rsid,
        "ncbi_uid": block.get("uid"),
        "ncbi_title": block.get("title") or block.get("name"),
        "ncbi_chr": genomic.get("chr") or genomic.get("chromosome"),
        "ncbi_position": genomic.get("chrpos") or genomic.get("position"),
        "ncbi_gene": genes.get("name") or genes.get("geneid"),
        "ncbi_taxid": block.get("taxid"),
    }


def read_input_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        return pd.read_csv(path, low_memory=False)


def enrich_file(input_path: str, output_path: str, checkpoint_path: str, id_column: str = "rsid", email: str | None = None, api_key: str | None = None, delay_sec: float = 0.34, batch_size: int = 50, progress=None):
    progress = progress or (lambda msg: None)
    src = Path(input_path).expanduser()
    out = Path(output_path).expanduser()
    chk = Path(checkpoint_path).expanduser()
    df = read_input_table(src)
    if id_column not in df.columns:
        guessed = pick_col(df, [id_column, "rsid", "riskAllele", "SNPS", "variantId"])
        if not guessed:
            raise ValueError(f"Didn't find ID column: {id_column}. Avaliable columns: {', '.join(map(str, df.columns))}")
        id_column = guessed
    df["rsid"] = df[id_column].apply(extract_rsid)
    work = df[df["rsid"].notna()].copy()
    if work.empty:
        raise ValueError("Didn't find rsID in the table.")
    unique_ids = work["rsid"].drop_duplicates().tolist()

    cp = load_checkpoint(chk)
    done = set(cp.get("completed", []))
    failed = set(cp.get("failed", []))
    pending = [x for x in unique_ids if x not in done]

    results = []
    if out.exists():
        try:
            old = pd.read_excel(out)
            keep_cols = [c for c in old.columns if c == "query_rsid" or c.startswith("ncbi_")]
            results = old[keep_cols].drop_duplicates().to_dict(orient="records") if keep_cols else []
        except Exception:
            results = []

    client = NCBIDownloader(email=email, api_key=api_key)
    total = len(unique_ids)
    processed = len(done)

    for idx, rsid in enumerate(pending, start=1):
        try:
            payload = client.fetch_snp_docsum(rsid)
            row = parse_ncbi(payload, rsid)
            results = [r for r in results if r.get("query_rsid") != rsid]
            results.append(row)
            done.add(rsid)
            failed.discard(rsid)
            processed += 1
            status = "ok"
        except Exception as e:
            failed.add(rsid)
            status = f"error: {e}"
        cp["completed"] = sorted(done)
        cp["failed"] = sorted(failed)
        cp["meta"] = {"processed": processed, "total": total, "last_id": rsid, "updated_at": time.strftime("%Y-%m-%d %H:%M:%S")}
        save_checkpoint(chk, cp)
        progress({"event": "item", "current": processed, "total": total, "rsid": rsid, "status": status})

        if idx % batch_size == 0 or idx == len(pending):
            enrich_df = pd.DataFrame(results)
            merged = work.merge(enrich_df, left_on="rsid", right_on="query_rsid", how="left")
            out.parent.mkdir(parents=True, exist_ok=True)
            merged.to_excel(out, index=False)
            progress({"event": "batch_saved", "rows": len(merged), "output": str(out)})
        time.sleep(delay_sec)

    if processed == total and not failed and chk.exists():
        chk.unlink()
        progress({"event": "checkpoint_deleted", "checkpoint": str(chk)})
    progress({"event": "done", "processed": processed, "total": total, "output": str(out)})


def build_parser():
    p = argparse.ArgumentParser(description="NCBI Downloader (CLI only) with checkpoint resume")
    p.add_argument("input_path", nargs="?", help="Path to TSV/CSV/XLSX. If omitted, you'll be asked for it interactively.")
    p.add_argument("--output", default=None)
    p.add_argument("--checkpoint", default=None)
    p.add_argument("--id-column", default=None)
    p.add_argument("--email", default=None)
    p.add_argument("--api-key", default=None)
    p.add_argument("--delay", type=float, default=0.34, help="Delay between requests, seconds")
    p.add_argument("--batch-size", type=int, default=50, help="Save output every N processed IDs")
    return p


def ask(prompt: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default else ""
    val = input(f"{prompt}{suffix}: ").strip()
    return val if val else (default or "")


def prompt_for_missing_args(args):
    """Interactively fill in anything the user didn't pass on the command line,
    so double-clicking the script (or running it with no args) still works."""
    if not args.input_path:
        print("=== NCBI Downloader (interactive mode) ===")
        print("Press Enter to accept the default shown in [brackets].\n")
        while not args.input_path:
            args.input_path = input("Path to input file (TSV/CSV/XLSX): ").strip()
            if not args.input_path:
                print("This field is required.")

    stem = Path(args.input_path).expanduser().stem
    if not args.output:
        args.output = ask("Output Excel path", str(Path.home() / f"{stem}_ncbi.xlsx"))
    if not args.checkpoint:
        args.checkpoint = ask("Checkpoint (resume) file path", str(Path.home() / f"{stem}.checkpoint.json"))
    if not args.id_column:
        args.id_column = ask("Column with rsID or variant text", "rsid")
    if args.email is None:
        args.email = ask("NCBI email (recommended, press Enter to skip)") or None
    if args.api_key is None:
        args.api_key = ask("NCBI API key (optional, press Enter to skip)") or None
    return args


def wait_before_exit():
    try:
        input("\nPress Enter to close this window...")
    except Exception:
        pass


def main():
    args = build_parser().parse_args()
    ran_interactively = args.input_path is None
    args = prompt_for_missing_args(args)

    enrich_file(
        input_path=args.input_path,
        output_path=args.output,
        checkpoint_path=args.checkpoint,
        id_column=args.id_column,
        email=args.email,
        api_key=args.api_key,
        delay_sec=args.delay,
        batch_size=args.batch_size,
        progress=lambda x: print(json.dumps(x, ensure_ascii=False)),
    )

    if ran_interactively:
        wait_before_exit()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted by user. Progress up to the last saved checkpoint/batch is kept.")
        wait_before_exit()
        raise SystemExit(1)
    except Exception as e:
        details = "".join(traceback.format_exception_only(type(e), e)).strip()
        print(f"\nNCBI Downloader failed.\n\n{details}", file=sys.stderr)
        wait_before_exit()
        raise SystemExit(1)
