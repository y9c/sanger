#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Chromatogram File Utils.

wrap cfutils into cli app
- update in 20190405
"""

import click
import matplotlib as mpl
from pathlib import Path

mpl.use("Agg", force=True)

from cfutils.run import report_mutation
from cfutils.parser import parse_abi
from cfutils.tracks import join_tracks, split_track, export_tracks, slice_track
from cfutils.qc import read_metrics


@click.group()
@click.option("--debug/--no-debug", default=False)
def cli(debug):
    """Chromatogram File Utils."""
    if debug:
        click.echo("Debug mode is on")


# call mutation
@cli.command()
@click.option("--query", prompt="QUERY (abi file): ", help="Query file in abi format")
@click.option(
    "--subject",
    prompt="SUBJECT (fasta file): ",
    help="Subject file in fasta format as ref",
)
@click.option("--outdir", default=None, required=False, help="Output directory")
@click.option("--outbase", default=None, required=False, help="Output basename")
@click.option(
    "--aligned/--mutated",
    default=False,
    help="Report all aligned sites or mutation sites only",
)
@click.option(
    "--plot/--no-plot",
    default=False,
    help="Generate figure of mutation in chromatogram.",
)
def mut(query, subject, outdir, outbase, aligned, plot):
    """do mutation calling, then report in tsv and pdf."""
    report_mutation(
        query_ab1_file=query,
        subject_fasta_file=subject,
        output_dir=outdir,
        file_basename=outbase,
        report_all_sites=aligned,
        report_mut_plot=plot,
    )


@cli.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(exists=True))
def qc(files):
    """Report per-read QC metrics as a table."""
    header = ["file", "bases", "gc%", "N%", "meanQ", "minQ", "loQ%",
              "trim_start", "trim_end"]
    click.echo("\t".join(header))
    for f in files:
        rec = parse_abi(f)
        m = read_metrics(rec)
        row = [Path(f).name, int(m["n_bases"]), f"{m['gc_percent']:.1f}",
               f"{100*m['n_fraction']:.1f}", f"{m['mean_qual']:.1f}",
               int(m["min_qual"]), f"{100*m['low_qual_fraction']:.1f}",
               int(m["trim_start"]), int(m["trim_end"])]
        click.echo("\t".join(str(x) for x in row))


@cli.group()
def track():
    """Split, join, slice and slice chromatogram trace files."""


@track.command("join")
@click.argument("files", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("--outbase", "outbase", required=True, help="Output basename")
@click.option("--outdir", default=".", help="Output directory")
@click.option("--gap", default=24, type=int, help="Gap samples between reads")
def track_join(files, outbase, outdir, gap):
    """Join multiple ABI chromatograms end to end into one record."""
    records = [parse_abi(f) for f in files]
    joined = join_tracks(*records, name=outbase, gap=gap)
    paths = export_tracks([joined], outdir, fmt="npz")
    click.echo(f"Joined {len(files)} tracks into {len(joined.seq)} bases")
    for p in paths:
        click.echo(f"wrote {p}")


@track.command("split")
@click.argument("file", type=click.Path(exists=True))
@click.option("--cuts", required=True, help="Comma-separated 1-based cut positions")
@click.option("--outdir", default=".", help="Output directory")
@click.option("--fmt", default="npz", type=click.Choice(["npz", "tsv"]))
def track_split(file, cuts, outdir, fmt):
    """Split one ABI chromatogram into segments at the given cut positions."""
    record = parse_abi(file)
    cut_positions = [int(c) for c in cuts.split(",") if c.strip() != ""]
    pieces = split_track(record, cut_positions)
    paths = export_tracks(pieces, outdir, fmt=fmt)
    for p in paths:
        click.echo(f"wrote {p}")


@track.command("slice")
@click.argument("file", type=click.Path(exists=True))
@click.option("--start", required=True, type=int)
@click.option("--end", required=True, type=int)
@click.option("--outdir", default=".")
@click.option("--fmt", default="npz", type=click.Choice(["npz", "tsv"]))
def track_slice(file, start, end, outdir, fmt):
    """Extract a 1-based inclusive region [start, end] from a chromatogram."""
    record = parse_abi(file)
    seg = slice_track(record, start, end)
    paths = export_tracks([seg], outdir, fmt=fmt)
    for p in paths:
        click.echo(f"wrote {p}")
