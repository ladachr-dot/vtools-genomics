# Haplo Converter

`haplo_converter` is the Y-DNA haplogroup conversion module inside `vtools-genomics`. It helps convert between ISOGG-style labels (for example, `R1b1a1`) and YFull-style labels (for example, `R-M269`), search SNP and haplogroup entries, and build a local searchable cache from the public YFull tree JSON.
The module is implemented in `converter.py` and is exposed through the unified project CLI as:

```bash
vtools haplo-convert
```

This command launches an interactive console menu for loading data, searching the local database, and converting files.

## What it does

- Converts haplogroup labels between ISOGG and YFull notation when a mapping is available.
- Downloads the public YFull tree JSON and imports its SNP-to-haplogroup data into a local SQLite database.
- Searches by haplogroup, SNP, substring, or approximate SNP name using exact, `LIKE`, and fuzzy matching.
- Batch-converts plain `.csv` or `.txt` files by scanning cells for recognizable haplogroup labels and saving converted results to a new file.

## Current interface

The current version is a **console application**, not a desktop GUI. The `vtools` entry point calls `HaploCli().run()` from `converter.py`, which opens a numbered text menu in the terminal.

Main menu actions:

- Load YFull tree from the internet.
- Load a local JSON file.
- Load a CSV file into the database.
- Convert a text/CSV file.
- Search the database.
- Clear the database.

## Data sources

The module uses two main sources of information:

- The public YFull tree JSON, downloaded from `https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json`.
- Optional CSV data imported by the user, expected to contain at least `Subgroup Name` and `Name` columns.

At runtime, the module stores imported records in a local SQLite database called `haplo_pro.db`, and it also keeps a local copy of `current_tree.json` when downloaded.

## How matching works

The converter contains a built-in ISOGG-to-YFull mapping dictionary and also builds a searchable SNP index from the YFull JSON tree.

For haplogroup lookup, it tries to:

- Normalize the user input.
- Detect whether the label looks like ISOGG or YFull notation.
- Resolve exact matches first.
- Fall back to parent haplogroups when an exact mapping is unavailable.

For broader search, it can:

- Search exact SNP names in SQLite.
- Search partial matches in `name`, `hg`, `path`, and `source`.
- Use fuzzy matching with `difflib.get_close_matches()` for SNP-like queries.

## Installation

From the repository root:

```bash
pip install -e .
```

The project package requires Python 3.9+ in `pyproject.toml`, and the haplogroup converter itself relies on standard-library modules plus `requests` for YFull tree download.

## Usage

### Launch from the unified CLI

```bash
vtools haplo-convert
```

This opens the interactive console menu.

### Run the module directly

You can also run the module file itself:

```bash
python src/vtools/haplo_converter/converter.py
```

That starts the same interactive workflow implemented under `if __name__ == "__main__"`.

## Menu actions

### 1. Load YFull tree from the internet

This option downloads the latest `current_tree.json`, saves it locally next to the module, parses the tree, extracts SNP markers, and writes them into SQLite.

Use this when the database is empty or when you want to refresh the local YFull-derived index.

### 2. Load a local JSON file

This imports a local YFull-style JSON tree file into the SQLite database instead of downloading it.

Use this if you already have a saved `current_tree.json` or want to work offline after obtaining the file.

### 3. Load a CSV file

This imports a CSV file into the same SQLite cache. The parser looks for a header row containing `Subgroup Name` and `Name` within the first 100 rows of the file.

The imported records are normalized, stored in SQLite, and enriched with a resolved path when a matching YFull parent path can be inferred.

### 4. Convert a file

This mode scans every cell in a `.csv` or `.txt` file and tries to recognize haplogroup notation in each value.

When matches are found, the module writes a new file named like `<original>_converted.csv` or `<original>_converted.txt` and includes these columns:

- `row`
- `column`
- `original`
- `detected_format`
- `direction`
- `converted`

### 5. Search

Search accepts:

- Haplogroup labels such as `R1b1a1` or `J2`.
- YFull labels such as `R-M269`.
- SNP names such as `M269`, `Z280`, or `CTS1211`.
- Partial text for broader lookup.

For haplogroups, the tool uses a fast path with notation detection and parent fallback. For SNPs and general text, it searches the SQLite table directly and can fall back to fuzzy matching.

### 6. Clear database

This deletes all rows from the local SQLite table but does not remove the code itself.

Use it when you want to rebuild the local cache from scratch.

## Files created at runtime

The converter may create these runtime files in its working/module directory:

- `haplo_pro.db` — local SQLite cache.
- `current_tree.json` — downloaded YFull tree JSON.
- `*_converted.csv` or `*_converted.txt` — results of batch conversion.

These are generated artifacts and usually should not be committed to the repository.

## Example workflow

### Initialize the local cache

1. Run `vtools haplo-convert`.
2. Choose the option to load the YFull tree from the internet.
3. Wait until the import finishes.
4. Use the search option to test a known haplogroup or SNP.

### Convert a plain text or CSV file

1. Run `vtools haplo-convert`.
2. Choose `Convert file`.
3. Enter the path to a `.txt` or `.csv` file.
4. Review the generated `_converted` output file.

## Limitations

- The module is interactive and menu-driven; it does not currently expose a rich argument-based non-interactive CLI for all operations.
- YFull conversion quality depends on the built-in mapping dictionary and the available YFull JSON/SNP data.
- CSV import expects a specific structure and will fail if `Subgroup Name` and `Name` cannot be detected.
- The module stores state locally in SQLite, so results depend on what has already been imported into the local database.

## Relation to the project

`haplo_converter` is one of the three subtools in `vtools-genomics`, alongside the PGS-to-PLINK converter and the NCBI metadata downloader. The root CLI groups them under a single `vtools` command.
