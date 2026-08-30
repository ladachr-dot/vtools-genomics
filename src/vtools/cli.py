"""Unified command-line entry point for the vtools-genomics toolkit.

Exposes the three tools as subcommands of a single `vtools` command:

    vtools pgs-convert SCORE_FILE [OPTIONS]
    vtools ncbi-download INPUT_PATH [OPTIONS]
    vtools haplo-convert

Each subcommand wraps the corresponding module under `vtools.<tool>` and
keeps its original behaviour and options intact.
"""
from __future__ import annotations

import sys
from pathlib import Path

import click

SRC_ROOT = Path(__file__).resolve().parents[1]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from vtools.pgs_to_plink.pgstoplink import RunConfig, process_score_file, setup_environment, check_plink_installed
from vtools.ncbi_downloader.downloader import enrich_file
from vtools.haplo_converter.converter import Cli as HaploCli


@click.group()
def cli():
    """vtools - genomic data processing utilities (PGS/PLINK, NCBI, ISOGG/YFull)."""


@cli.command("pgs-convert")
@click.argument("score_file", type=click.Path(exists=True))
@click.option("--rsid-mode", type=click.Choice(["auto", "interactive"]), default="auto",
              help="How to resolve ambiguous rsID candidates.")
@click.option("--start-idx", type=int, default=0, help="First row index to process.")
@click.option("--end-idx", type=int, default=None, help="Last row index (exclusive) to process.")
@click.option("--bim", "output_bim", is_flag=True, default=False, help="Also create a .bim file.")
@click.option("--keep-checkpoint/--no-keep-checkpoint", default=True,
              help="Keep the checkpoint file after a successful run.")
@click.option("--plink-hint", "show_plink_hint", is_flag=True, default=False,
              help="Print an example PLINK command after finishing.")
def pgs_convert(score_file, rsid_mode, start_idx, end_idx, output_bim, keep_checkpoint, show_plink_hint):
    """Convert a PGS score file into PLINK-ready format, with optional rsID lookup and BIM creation."""
    score_path = str(Path(score_file).resolve())
    cfg = RunConfig(
        score_file=score_path,
        output_bim=output_bim,
        rsid_mode=rsid_mode,
        start_idx=start_idx,
        end_idx=end_idx,
        keep_checkpoint=keep_checkpoint,
        show_plink_hint=show_plink_hint,
    )
    if check_plink_installed():
        click.echo("plink is installed and available in PATH")
    else:
        click.echo("plink was not found in PATH (only needed to run the scoring step)")
    if not setup_environment():
        raise click.ClickException("Failed to set up the environment (dependencies/liftover chains).")
    process_score_file(score_path, cfg)


@cli.command("ncbi-download")
@click.argument("input_path", type=click.Path(exists=True))
@click.option("--output", default=None, help="Output Excel path (default: <input>_ncbi.xlsx next to the input).")
@click.option("--checkpoint", default=None, help="Checkpoint (resume) file path.")
@click.option("--id-column", default="rsid", help="Column with rsID or variant text.")
@click.option("--email", default=None, help="NCBI email (recommended).")
@click.option("--api-key", default=None, help="NCBI API key (optional).")
@click.option("--delay", type=float, default=0.34, help="Delay between batch saves, seconds.")
@click.option("--batch-size", type=int, default=50, help="Save output every N processed IDs.")
def ncbi_download(input_path, output, checkpoint, id_column, email, api_key, delay, batch_size):
    """Batch-enrich a TSV/CSV/XLSX table with NCBI rsID metadata (chromosome, position, gene)."""
    stem = Path(input_path).expanduser().stem
    output = output or str(Path(input_path).expanduser().parent / f"{stem}_ncbi.xlsx")
    checkpoint = checkpoint or str(Path(input_path).expanduser().parent / f"{stem}.checkpoint.json")
    enrich_file(
        input_path=input_path,
        output_path=output,
        checkpoint_path=checkpoint,
        id_column=id_column,
        email=email,
        api_key=api_key,
        delay_sec=delay,
        batch_size=batch_size,
        progress=lambda x: click.echo(str(x)),
    )


@cli.command("haplo-convert")
def haplo_convert():
    """Launch the interactive ISOGG <-> YFull haplogroup converter (console menu)."""
    HaploCli().run()


def main():
    cli()


if __name__ == "__main__":
    main()
