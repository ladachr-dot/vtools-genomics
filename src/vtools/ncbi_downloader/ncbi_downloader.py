from __future__ import annotations

import argparse
import json
import re
import sys
import threading
import time
import traceback
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    TK_AVAILABLE = True
    TK_IMPORT_ERROR = None
except Exception as e:
    tk = None
    filedialog = messagebox = ttk = None
    TK_AVAILABLE = False
    TK_IMPORT_ERROR = e

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def show_startup_error(message: str) -> None:
    if TK_AVAILABLE:
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("NCBI Downloader error", message)
            root.destroy()
            return
        except Exception:
            pass
    print(message, file=sys.stderr)
    try:
        input("\nPress Enter to close...")
    except Exception:
        pass


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
            raise ValueError(f"Не нашла колонку с ID: {id_column}. Доступные колонки: {', '.join(map(str, df.columns))}")
        id_column = guessed
    df["rsid"] = df[id_column].apply(extract_rsid)
    work = df[df["rsid"].notna()].copy()
    if work.empty:
        raise ValueError("Не найдено ни одного rsID во входной таблице.")
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


if TK_AVAILABLE:
    class DownloaderApp(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("NCBI Downloader")
            self.geometry("900x680")
            self.queue = []
            self.worker = None
            self.build_ui()
            self.after(150, self.poll_queue)

        def build_ui(self):
            frame = ttk.Frame(self)
            frame.pack(fill="both", expand=True)
            pad = {"padx": 10, "pady": 6}
            header = "Загрузка данных из NCBI по rsID с докачкой.\nЗаполните поля ниже, затем нажмите Start / Resume download."
            ttk.Label(frame, text=header, justify="left").grid(row=0, column=0, columnspan=3, sticky="w", padx=10, pady=10)
            ttk.Label(frame, text="Input GWAS/PGS file").grid(row=1, column=0, sticky="w", **pad)
            self.input_var = tk.StringVar()
            ttk.Entry(frame, textvariable=self.input_var, width=74).grid(row=1, column=1, sticky="ew", **pad)
            ttk.Button(frame, text="Browse", command=self.pick_file).grid(row=1, column=2, **pad)
            ttk.Label(frame, text="Output Excel (.xlsx)").grid(row=2, column=0, sticky="w", **pad)
            self.output_var = tk.StringVar(value=str(Path.home() / "Desktop" / "study_ncbi.xlsx"))
            ttk.Entry(frame, textvariable=self.output_var, width=74).grid(row=2, column=1, sticky="ew", **pad)
            ttk.Label(frame, text="Resume checkpoint (.json)").grid(row=3, column=0, sticky="w", **pad)
            self.checkpoint_var = tk.StringVar(value=str(Path.home() / "Desktop" / "study_ncbi.checkpoint.json"))
            ttk.Entry(frame, textvariable=self.checkpoint_var, width=74).grid(row=3, column=1, sticky="ew", **pad)
            ttk.Label(frame, text="Column with rsID or variant text").grid(row=4, column=0, sticky="w", **pad)
            self.id_var = tk.StringVar(value="rsid")
            ttk.Entry(frame, textvariable=self.id_var, width=30).grid(row=4, column=1, sticky="w", **pad)
            ttk.Label(frame, text="NCBI email (recommended)").grid(row=5, column=0, sticky="w", **pad)
            self.email_var = tk.StringVar()
            ttk.Entry(frame, textvariable=self.email_var, width=40).grid(row=5, column=1, sticky="w", **pad)
            ttk.Label(frame, text="NCBI API key (optional)").grid(row=6, column=0, sticky="w", **pad)
            self.api_var = tk.StringVar()
            ttk.Entry(frame, textvariable=self.api_var, width=40, show="*").grid(row=6, column=1, sticky="w", **pad)
            self.progress = ttk.Progressbar(frame, mode="determinate")
            self.progress.grid(row=7, column=0, columnspan=3, sticky="ew", padx=10, pady=12)
            self.status_var = tk.StringVar(value="Ready")
            ttk.Label(frame, textvariable=self.status_var).grid(row=8, column=0, columnspan=3, sticky="w", padx=10)
            btns = ttk.Frame(frame)
            btns.grid(row=9, column=0, columnspan=3, sticky="w", padx=10, pady=4)
            ttk.Button(btns, text="Start / Resume download", command=self.start).pack(side="left", padx=4)
            ttk.Button(btns, text="What to write here", command=self.show_help).pack(side="left", padx=4)
            self.log = tk.Text(frame, height=24, wrap="word")
            self.log.grid(row=10, column=0, columnspan=3, sticky="nsew", padx=10, pady=10)
            frame.columnconfigure(1, weight=1)
            frame.rowconfigure(10, weight=1)

        def pick_file(self):
            path = filedialog.askopenfilename(filetypes=[("Genomics files", "*.tsv *.csv *.xlsx *.xls"), ("All files", "*.*")])
            if path:
                self.input_var.set(path)
                stem = Path(path).stem
                self.output_var.set(str(Path.home() / "Desktop" / f"{stem}_ncbi.xlsx"))
                self.checkpoint_var.set(str(Path.home() / "Desktop" / f"{stem}.checkpoint.json"))

        def show_help(self):
            message = (
                "Что писать в поля:\n\n"
                "1. Input GWAS/PGS file — путь к вашему входному файлу .tsv, .csv или .xlsx.\n"
                "2. Output Excel (.xlsx) — куда сохранить итоговую Excel-таблицу.\n"
                "3. Resume checkpoint (.json) — файл состояния для докачки.\n"
                "4. Column with rsID or variant text — обычно rsid, riskAllele, SNPS или variantId.\n"
                "5. NCBI email (recommended) — ваш email для Entrez.\n"
                "6. NCBI API key (optional) — ключ NCBI, если он у вас есть."
            )
            messagebox.showinfo("How to fill the form", message)

        def emit(self, payload: dict):
            self.queue.append(payload)

        def start(self):
            if self.worker and self.worker.is_alive():
                return
            args = dict(
                input_path=self.input_var.get(),
                output_path=self.output_var.get(),
                checkpoint_path=self.checkpoint_var.get(),
                id_column=self.id_var.get(),
                email=self.email_var.get() or None,
                api_key=self.api_var.get() or None,
                progress=self.emit,
            )
            self.worker = threading.Thread(target=self.run_job, args=(args,), daemon=True)
            self.worker.start()
            self.status_var.set("Running...")

        def run_job(self, kwargs):
            try:
                enrich_file(**kwargs)
            except Exception as e:
                self.emit({"event": "fatal", "message": str(e)})

        def poll_queue(self):
            while self.queue:
                item = self.queue.pop(0)
                if item.get("event") == "item":
                    total = max(int(item.get("total", 0)), 1)
                    current = int(item.get("current", 0))
                    self.progress.configure(maximum=total, value=current)
                    msg = f"{current}/{total} | {item.get('rsid')} | {item.get('status')}"
                    self.status_var.set(msg)
                    self.log.insert("end", msg + "\n")
                    self.log.see("end")
                elif item.get("event") == "batch_saved":
                    self.log.insert("end", f"Saved {item.get('rows')} rows -> {item.get('output')}\n")
                    self.log.see("end")
                elif item.get("event") == "checkpoint_deleted":
                    self.log.insert("end", f"Checkpoint deleted: {item.get('checkpoint')}\n")
                    self.log.see("end")
                elif item.get("event") == "done":
                    self.status_var.set(f"Done: {item.get('processed')}/{item.get('total')}")
                    self.log.insert("end", f"Completed. Output: {item.get('output')}\n")
                    self.log.see("end")
                elif item.get("event") == "fatal":
                    self.status_var.set("Failed")
                    messagebox.showerror("Error", item.get("message"))
            self.after(150, self.poll_queue)


def build_parser():
    p = argparse.ArgumentParser(description="Single-file NCBI Downloader with GUI and checkpoint resume")
    p.add_argument("input_path", nargs="?", help="Path to TSV/CSV/XLSX")
    p.add_argument("--output", default=str(Path.home() / "Desktop" / "study_ncbi.xlsx"))
    p.add_argument("--checkpoint", default=str(Path.home() / "Desktop" / "study_ncbi.checkpoint.json"))
    p.add_argument("--id-column", default="rsid")
    p.add_argument("--email", default=None)
    p.add_argument("--api-key", default=None)
    p.add_argument("--gui", action="store_true")
    return p


def main():
    args = build_parser().parse_args()
    if args.gui or not args.input_path:
        if not TK_AVAILABLE:
            raise RuntimeError(
                "Tkinter is not available in this Python environment. Run from terminal in CLI mode or install Python with Tk support."
            )
        app = DownloaderApp()
        app.mainloop()
        return
    enrich_file(
        input_path=args.input_path,
        output_path=args.output,
        checkpoint_path=args.checkpoint,
        id_column=args.id_column,
        email=args.email,
        api_key=args.api_key,
        progress=lambda x: print(json.dumps(x, ensure_ascii=False)),
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        details = "".join(traceback.format_exception_only(type(e), e)).strip()
        show_startup_error(f"NCBI Downloader could not start.\n\n{details}")
        raise SystemExit(1)
