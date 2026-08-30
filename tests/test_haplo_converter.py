"""Unit tests for pure ISOGG<->YFull lookup logic in vtools.haplo_converter.converter.

These tests exercise the in-memory ISOGG_TO_YFULL mapping and the notation
detector directly, without touching the network or downloading the YFull
tree (only a local, throwaway SQLite file under tmp_path is created).
"""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from vtools.haplo_converter.converter import HaploDB, Converter


def _converter(tmp_path):
    db = HaploDB(db_path=str(tmp_path / "haplo_test.db"))
    return Converter(db)


def test_get_yfull_equivalent_exact_match(tmp_path):
    conv = _converter(tmp_path)
    assert conv.get_yfull_equivalent("R1a") == "R-Y482"


def test_get_yfull_equivalent_falls_back_to_parent(tmp_path):
    conv = _converter(tmp_path)
    assert conv.get_yfull_equivalent("R1a1a1b1a99") != ""


def test_get_yfull_equivalent_unknown_returns_empty(tmp_path):
    conv = _converter(tmp_path)
    assert conv.get_yfull_equivalent("ZZZNOTREAL") == ""


def test_detect_hg_notation_yfull(tmp_path):
    conv = _converter(tmp_path)
    assert conv.detect_hg_notation("R-M269") == "YFULL"


def test_detect_hg_notation_isogg(tmp_path):
    conv = _converter(tmp_path)
    assert conv.detect_hg_notation("R1b1a1b1a1b") == "ISOGG"


def test_detect_hg_notation_snp_like_is_not_isogg(tmp_path):
    conv = _converter(tmp_path)
    assert conv.detect_hg_notation("M269") is None
