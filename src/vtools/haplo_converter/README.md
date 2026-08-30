# PhyloResolve

**PhyloResolve** is a bidirectional ISOGG ↔ YFull haplogroup converter with a live SNP index, SQLite caching. It resolves Y-DNA haplogroup labels and SNP markers to their position in the YFull phylogenetic tree and vice versa.

---

## Features

- **Live SNP index** built from the YFull `current_tree.json` — no hard-coded SNP list
- **SQLite cache** with three indexes so the JSON is never re-parsed on subsequent runs
- **Auto-detection** of input notation (ISOGG: `R1b1a1`, YFull: `R-M269`)
- **Bidirectional conversion** ISOGG → YFull and YFull → ISOGG
- **Three-level SNP search**: exact match → SQL LIKE → fuzzy difflib
- **Haplogroup search with parent fallback**: climbs the tree until a match is found
- **Static fallback dictionary** (~300 entries) for top-level nodes (`R`, `I`, `J`, …) not covered by SNP markers
- **Detailed output statuses**: `matched_by_name_snp`, `matched_by_alternate_snp`, `ambiguous_multiple_yfull_matches`, `not_found_in_yfull_tree`, `missing_snp`
- **Excel output**: auto-sized columns + frozen header row (`freeze_panes`)
- **CLI mode** for scripted/batch use

---

## Requirements

```
Python 3.10+
requests
pandas        # for Excel conversion (Convert XLSX)
openpyxl      # for Excel conversion (Convert XLSX)
tkinter       # included in most Python distributions (for GUI)
```

Install dependencies:

```bash
pip install requests pandas openpyxl
```

---

## Quick Start

### CLI

Convert an ISOGG SNP Index workbook to Excel with YFull columns:

```bash
python phylo_resolve_en.py "SNP Index.xlsx"
```

Use a local tree file (no network required):

```bash
python phylo_resolve_en.py "SNP Index.xlsx" --tree current_tree.json --no-download
```

Specify output path and sheet name:

```bash
python phylo_resolve_en.py "SNP Index.xlsx" --sheet Human --output result.xlsx
```

---

## CLI Arguments

| Argument | Default | Description |
|---|---|---|
| `input` | — | Path to ISOGG SNP Index `.xlsx` or `.csv` |
| `--sheet` | `Human` | Excel sheet name |
| `--output`, `-o` | auto | Output `.xlsx` or `.csv` path |
| `--tree` | `current_tree.json` | Path to local YFull tree JSON |
| `--url` | GitHub raw URL | URL to download the YFull tree from |
| `--no-download` | off | Skip network download, use local tree only |
| `--force-download` | off | Re-download even if a local tree exists |
| `--db` | `haplo_pro.db` | Path to the SQLite database file |
| `--gui` | off | Launch the GUI regardless of other arguments |

---

## GUI Walkthrough

### Toolbar buttons

| Button | Action |
|---|---|
| 🌐 YFull (online) | Download `current_tree.json` from GitHub and import it |
| 📥 Load JSON | Import a local `current_tree.json` file |
| 📊 Load ISOGG CSV | Import an ISOGG SNP Index CSV into the database |
| 🔄 Convert CSV | Batch-convert a plain CSV/TXT file (any cell containing a haplogroup) |
| 📋 Convert XLSX | Convert an ISOGG SNP Index workbook to Excel with YFull columns |
| 🗑 Clear DB | Wipe all data from the SQLite database |
| ❓ Help | Show a quick in-app guide |

### Search box

Type any of the following and press **SEARCH** or hit **Enter**:

- ISOGG haplogroup: `R1b1a1`, `J2a1`, `N1c1a1`
- YFull haplogroup: `R-M269`, `I-DF29`, `N-M231`
- SNP name: `M269`, `Z280`, `CTS1211`, `FGC8628`

### Output fields

```
🧬 Query:          R1b1a1
   Haplogroup:     R1B1A1
   YFull equiv.:   R-L389
   Path:           /ROOT/R/R-M207/R-M173/R-L754/R-L761/R-L389
   Source:         STATIC_DICT
```

---

## File Conversion

### Convert CSV (`🔄 Convert CSV`)

Scans every cell of a plain CSV or TXT file. Any cell that looks like a haplogroup (ISOGG or YFull notation) is converted. The output is a new CSV with columns:

| Column | Description |
|---|---|
| `row` | Source row number |
| `column` | Source column number |
| `original` | Original value |
| `detected_format` | `ISOGG` or `YFULL` |
| `direction` | `ISOGG->YFULL` or `YFULL->ISOGG` |
| `converted` | Converted haplogroup label |

### Convert XLSX (`📋 Convert XLSX`)

Reads an **ISOGG SNP Index** workbook (expects columns `Subgroup Name`, `Name`, `Alternate Names`). Adds YFull columns at the left of the sheet and saves as a new `.xlsx`:

| Column | Description |
|---|---|
| `ISOGG_hg_clean` | Cleaned ISOGG haplogroup label |
| `Name_SNP_candidates` | SNP candidates extracted from the `Name` column |
| `Alternate_SNP_candidates` | SNP candidates extracted from `Alternate Names` |
| `matched_snp` | The SNP that was matched in the YFull tree |
| `YFull_hg` | Matched YFull haplogroup (e.g. `R-M269`) |
| `YFull_path` | Full path in the YFull tree |
| `status` | Match status (see below) |
| `all_yfull_matches` | All found matches, pipe-separated |
| `approx_YFull_if_not_found` | Best-effort guess when no match is found |

---

## Match Statuses

| Status | Meaning |
|---|---|
| `matched_by_name_snp` | Match found using the primary `Name` column SNP |
| `matched_by_alternate_snp` | Match found using an `Alternate Names` SNP |
| `ambiguous_multiple_yfull_matches` | SNP matched more than one YFull haplogroup |
| `not_found_in_yfull_tree` | SNP was extracted but is absent from the YFull tree |
| `missing_snp` | No SNP candidates could be extracted from the row |

---

## Architecture

```
phylo_resolve_en.py
├── HaploDB          — SQLite layer (3 indexes: name, hg, path)
├── Converter        — business logic
│   ├── normalize_for_db / detect_hg_notation
│   ├── import_yfull_json   (iterative stack + progress)
│   ├── import_csv
│   ├── search_by_haplogroup_fast
│   ├── convert_haplogroup_value
│   ├── convert_file  (CSV batch)
│   └── convert_excel (ISOGG SNP Index → xlsx)
├── App              — Tkinter GUI layer
└── main()           — CLI entry point
```

### Key design decisions

- **Hyphen is preserved** during normalization — `R-M269` is never mangled to `RM269`.
- **Static dictionary is a fallback**, not the primary source. The live YFull index is checked first; the dictionary only fills gaps for top-level nodes without a distinct SNP marker.
- **`search_by_haplogroup_fast` has a single definition** (with `looks_like_haplogroup` guard) — the duplicate from earlier versions is removed.
- **`import_yfull_json` uses abstractions** (`node_snps`, `node_children`, `split_yfull_snps`) that handle both `dict` and `list` root formats and different field names (`snps`, `SNPs`, `markers`).

---

## Notes

- The YFull tree (`current_tree.json`) can contain 200 000+ nodes. Initial import takes 30–120 seconds. Subsequent launches read from SQLite instantly.
- The SQLite database (`haplo_pro.db`) is created in the working directory. Delete it to force a full re-import.
- Excel conversion requires `pandas` and `openpyxl`. Search and CSV conversion work without them.
