"""Unit tests for pure helper functions in vtools.ncbi_downloader.downloader.

Covers rsID extraction and checkpoint load/save, which do not require any
network access to NCBI Entrez.
"""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from vtools.ncbi_downloader.downloader import extract_rsid, load_checkpoint, save_checkpoint


def test_extract_rsid_from_plain_value():
    assert extract_rsid("rs12345") == "rs12345"


def test_extract_rsid_case_insensitive_and_embedded():
    assert extract_rsid("Variant RS987 (risk)") == "rs987"


def test_extract_rsid_none_when_absent():
    assert extract_rsid("chr1:12345") is None


def test_extract_rsid_handles_none_and_nan():
    assert extract_rsid(None) is None


def test_checkpoint_roundtrip(tmp_path):
    chk = tmp_path / "resume.json"
    fresh = load_checkpoint(chk)
    assert fresh == {"completed": [], "failed": [], "meta": {"processed": 0, "total": 0, "last_id": None}}

    fresh["completed"] = ["rs1", "rs2"]
    save_checkpoint(chk, fresh)
    reloaded = load_checkpoint(chk)
    assert reloaded["completed"] == ["rs1", "rs2"]
