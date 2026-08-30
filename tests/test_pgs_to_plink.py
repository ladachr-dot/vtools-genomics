"""Unit tests for pure helper functions in vtools.pgs_to_plink.pgstoplink.

These tests only exercise functions that do not touch the network, PLINK,
or the filesystem outside of pytest's tmp_path fixture, so they can run in
CI without external dependencies (no NCBI calls, no PLINK binary, no
liftover chain files required).
"""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from vtools.pgs_to_plink.pgstoplink import (
    normalize_column_name,
    find_column,
    extract_genome_build,
)


def test_normalize_column_name_basic():
    assert normalize_column_name("Effect-Allele") == "effect_allele"
    assert normalize_column_name(" chr Position ") == "chr_position"
    assert normalize_column_name("BETA") == "beta"


def test_find_column_matches_synonym():
    columns = ["hm_rsid", "hm_effect_allele", "hm_beta"]
    synonyms = ["rsid", "rs_id", "snp", "variant_id", "hm_rsid"]
    assert find_column(columns, synonyms) == "hm_rsid"


def test_find_column_no_match_returns_none():
    columns = ["foo", "bar"]
    synonyms = ["rsid", "snp"]
    assert find_column(columns, synonyms) is None


def test_extract_genome_build_reads_header_comment(tmp_path):
    score_file = tmp_path / "score.txt"
    score_file.write_text(
        "#genome_build=GRCh38\n"
        "rsID\teffect_allele\tbeta\n"
        "rs1\tA\t0.1\n",
        encoding="utf-8",
    )
    assert extract_genome_build(str(score_file)) == "GRCh38"


def test_extract_genome_build_missing_header_is_unknown(tmp_path):
    score_file = tmp_path / "score.txt"
    score_file.write_text("rsID\teffect_allele\tbeta\nrs1\tA\t0.1\n", encoding="utf-8")
    assert extract_genome_build(str(score_file)) == "unknown"


def test_extract_genome_build_missing_file_is_unknown():
    assert extract_genome_build("/nonexistent/path/score.txt") == "unknown"
