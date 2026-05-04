"""
Tests for vtools-genomics core functions.
Covers: extract_rsid, parse_ncbi, pick_col, load_checkpoint (ncbi_downloader)

Run with:
    pytest tests/test_vtools.py -v
"""

import json
import pytest
import pandas as pd
from pathlib import Path

# --- import functions under test ---
from vtools.ncbi_downloader.ncbi_downloader import extract_rsid, parse_ncbi, pick_col, load_checkpoint


# ─────────────────────────────────────────────
# extract_rsid
# ─────────────────────────────────────────────

class TestExtractRsid:
    def test_plain_rsid(self):
        assert extract_rsid("rs12345") == "rs12345"

    def test_rsid_uppercase_normalized(self):
        # function should return lowercase
        assert extract_rsid("RS12345") == "rs12345"

    def test_rsid_embedded_in_text(self):
        assert extract_rsid("variant rs9988023 in region") == "rs9988023"

    def test_chr_pos_no_rsid(self):
        assert extract_rsid("1:204518n") is None

    def test_none_input(self):
        assert extract_rsid(None) is None

    def test_nan_input(self):
        assert extract_rsid(float("nan")) is None

    def test_empty_string(self):
        assert extract_rsid("") is None

    def test_numeric_input(self):
        # int with no rsID pattern
        assert extract_rsid(42) is None

    def test_rsid_with_leading_zeros(self):
        assert extract_rsid("rs0001") == "rs0001"


# ─────────────────────────────────────────────
# parse_ncbi
# ─────────────────────────────────────────────

class TestParseNcbi:
    def _make_payload(self, uid="123", chr="1", pos=12345, gene="BRCA1"):
        return {
            "result": {
                "uids": [uid],
                uid: {
                    "uid": uid,
                    "title": f"rs{uid}",
                    "taxid": 9606,
                    "genomicinfo": [{"chr": chr, "chrpos": pos}],
                    "genes": [{"name": gene}],
                }
            }
        }

    def test_normal_payload(self):
        payload = self._make_payload(uid="123", chr="1", pos=12345, gene="BRCA1")
        result = parse_ncbi(payload, "rs123")
        assert result["ncbi_chr"] == "1"
        assert result["ncbi_position"] == 12345
        assert result["ncbi_gene"] == "BRCA1"
        assert result["query_rsid"] == "rs123"

    def test_empty_payload(self):
        result = parse_ncbi({}, "rs999")
        assert result["query_rsid"] == "rs999"
        assert result["ncbi_chr"] is None
        assert result["ncbi_position"] is None

    def test_missing_result_key(self):
        result = parse_ncbi({"other": "data"}, "rs000")
        assert result["ncbi_uid"] is None

    def test_non_dict_payload(self):
        result = parse_ncbi(None, "rs111")
        assert result["query_rsid"] == "rs111"
        assert result["ncbi_gene"] is None

    def test_empty_genomicinfo(self):
        payload = {
            "result": {
                "uids": ["1"],
                "1": {
                    "uid": "1",
                    "title": "rs1",
                    "taxid": 9606,
                    "genomicinfo": [],
                    "genes": [],
                }
            }
        }
        result = parse_ncbi(payload, "rs1")
        assert result["ncbi_chr"] is None
        assert result["ncbi_gene"] is None

    def test_taxid_preserved(self):
        payload = self._make_payload(uid="7", chr="X", pos=999)
        result = parse_ncbi(payload, "rs7")
        assert result["ncbi_taxid"] == 9606


# ─────────────────────────────────────────────
# pick_col
# ─────────────────────────────────────────────

class TestPickCol:
    def setup_method(self):
        self.df = pd.DataFrame(columns=["RSID", "Effect_Allele", "beta", "CHR_POS"])

    def test_exact_case_insensitive_match(self):
        assert pick_col(self.df, ["rsid"]) == "RSID"

    def test_first_candidate_wins(self):
        # both rsid and chr_pos exist; should return rsid (first in list)
        result = pick_col(self.df, ["rsid", "chr_pos"])
        assert result == "RSID"

    def test_fallback_to_second_candidate(self):
        result = pick_col(self.df, ["snpid", "effect_allele"])
        assert result == "Effect_Allele"

    def test_no_match_returns_none(self):
        assert pick_col(self.df, ["nonexistent", "also_missing"]) is None

    def test_empty_candidates(self):
        assert pick_col(self.df, []) is None

    def test_empty_dataframe(self):
        assert pick_col(pd.DataFrame(), ["rsid"]) is None


# ─────────────────────────────────────────────
# load_checkpoint
# ─────────────────────────────────────────────

class TestLoadCheckpoint:
    def test_missing_file_returns_defaults(self, tmp_path):
        path = tmp_path / "nonexistent.json"
        result = load_checkpoint(path)
        assert result["completed"] == []
        assert result["failed"] == []
        assert "processed" in result["meta"]

    def test_loads_existing_checkpoint(self, tmp_path):
        data = {"completed": ["rs1", "rs2"], "failed": [], "meta": {"processed": 2, "total": 5, "last_id": "rs2"}}
        path = tmp_path / "cp.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        result = load_checkpoint(path)
        assert result["completed"] == ["rs1", "rs2"]
        assert result["meta"]["processed"] == 2

    def test_returns_dict(self, tmp_path):
        path = tmp_path / "empty.json"
        result = load_checkpoint(path)
        assert isinstance(result, dict)
