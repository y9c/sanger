#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Chromatogram File Utils - command-line application.

The CLI uses rich-click for a styled, grouped help (themed sections) and is
organised around the package modules::

    cfutils mut            mutation calling & reporting
    cfutils qc             per-read quality-control metrics
    cfutils track          split / join / slice chromatogram trace files
    cfutils edit           trim, strip primers, reverse-complement
    cfutils basecall       re-call bases from raw four-channel traces
    cfutils assemble       reference-guided pileup & consensus
    cfutils analyze        sequence biology (translate, motifs, restriction)
    cfutils export         FASTA / VCF / JSON / batch summary
    cfutils plot           chromatogram rendering (+ features / DNA viewer)
"""

import importlib.metadata

import rich_click as click

try:
    __VERSION__ = importlib.metadata.version("cfutils")
except Exception:  # pragma: no cover - not installed as a package
    __VERSION__ = "0.0.0.dev62"

CTX = dict(help_option_names=["-h", "--help"])


class CfutilsGroup(click.RichGroup):
    """Root group that surfaces bad input/data errors cleanly.

    Library ``ValueError`` / ``RuntimeError`` (bad input or data) are turned
    into a concise ``ClickException`` message instead of a raw traceback.
    Normal help/abort control flow (``Exit`` / ``Abort``) still propagates.
    """

    def invoke(self, ctx):
        try:
            return super().invoke(ctx)
        except (click.exceptions.Exit, click.exceptions.Abort):
            raise
        except (ValueError, RuntimeError) as exc:
            raise click.ClickException(str(exc)) from exc


# --------------------------------------------------------------------------- #
#  Help / option theming (rich-click)
# --------------------------------------------------------------------------- #
click.rich_click.COMMAND_GROUPS = {
    "cfutils": [
        {"name": "Mutation & QC", "commands": ["mut", "qc"]},
        {"name": "Trace Editing", "commands": ["track", "edit"]},
        {"name": "Base Calling & Assembly", "commands": ["basecall", "assemble"]},
        {"name": "Sequence Analysis", "commands": ["analyze"]},
        {"name": "Export & Visualization", "commands": ["export", "plot"]},
    ],
}

click.rich_click.OPTION_GROUPS = {
    "cfutils mut": [
        {"name": "Inputs", "options": ["--query", "--subject"]},
        {"name": "Output", "options": ["--outdir", "--outbase"]},
        {"name": "Reporting", "options": ["--aligned/--mutated", "--plot/--no-plot"]},
    ],
    "cfutils edit trim": [
        {"name": "Trimming", "options": ["--cutoff"]},
        {"name": "Output", "options": ["--outdir"]},
    ],
    "cfutils track split": [
        {"name": "Split", "options": ["--cuts", "--outdir", "--fmt"]},
    ],
    "cfutils export vcf": [
        {"name": "Inputs", "options": ["--query", "--subject", "--min-qual"]},
    ],
}

click.rich_click.STYLE_OPTION = "bold green"
click.rich_click.STYLE_OPTION_DEFAULT = "dim"
click.rich_click.STYLE_COMMAND = "bold blue"


# --------------------------------------------------------------------------- #
#  Helpers
# --------------------------------------------------------------------------- #
def _emit(record, outdir):
    """Write a record to FASTA and echo the path."""
    from pathlib import Path

    from cfutils.export import to_fasta

    outdir = Path(outdir or ".")
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{record.name or 'out'}.fa"
    path.write_text(to_fasta(record))
    click.echo(f"✓ wrote {path}")


# --------------------------------------------------------------------------- #
#  Root group
# --------------------------------------------------------------------------- #
@click.group(
    cls=CfutilsGroup,
    invoke_without_command=False,
    help="cfutils — Sanger sequencing chromatogram analysis toolkit.",
    context_settings=CTX,
)
@click.version_option(__VERSION__, "--version", "-v")
def cli():
    """cfutils — Sanger sequencing chromatogram analysis toolkit."""


# --------------------------------------------------------------------------- #
#  Mut & QC
# --------------------------------------------------------------------------- #
@cli.command(
    help="Call mutations against a reference and report (tsv + optional pdf).",
    no_args_is_help=True,
    context_settings=CTX,
)
@click.option("--query", "-q", required=True, help="Query ABI file.")
@click.option("--subject", "-s", required=True, help="Reference FASTA file.")
@click.option("--outdir", "-o", default=None, help="Output directory")
@click.option("--outbase", "-b", default=None, help="Output basename")
@click.option(
    "--aligned/--mutated",
    default=False,
    help="Report all aligned sites (default: mutations only)",
)
@click.option(
    "--plot/--no-plot",
    default=False,
    help="Generate a chromatogram figure of the mutations",
)
def mut(query, subject, outdir, outbase, aligned, plot):
    """Call mutations against a reference and report (tsv + optional pdf)."""
    from cfutils.run import report_mutation

    report_mutation(
        query_ab1_file=query,
        subject_fasta_file=subject,
        output_dir=outdir,
        file_basename=outbase,
        report_all_sites=aligned,
        report_mut_plot=plot,
    )


@cli.command(
    help="Report per-read QC metrics as a table.",
    no_args_is_help=True,
    context_settings=CTX,
)
@click.argument("files", nargs=-1, required=True, type=click.Path(exists=True))
def qc(files):
    """Report per-read QC metrics as a table."""
    from pathlib import Path

    from cfutils.qc import read_metrics

    header = [
        "file",
        "bases",
        "gc%",
        "N%",
        "meanQ",
        "minQ",
        "loQ%",
        "trim_start",
        "trim_end",
    ]
    click.echo("\t".join(header))
    for f in files:
        m = read_metrics(__parse_abi(f))
        click.echo(
            "\t".join(
                str(x)
                for x in [
                    Path(f).name,
                    int(m["n_bases"]),
                    f"{m['gc_percent']:.1f}",
                    f"{100 * m['n_fraction']:.1f}",
                    f"{m['mean_qual']:.1f}",
                    int(m["min_qual"]),
                    f"{100 * m['low_qual_fraction']:.1f}",
                    int(m["trim_start"]),
                    int(m["trim_end"]),
                ]
            )
        )


# --------------------------------------------------------------------------- #
#  Trace Editing
# --------------------------------------------------------------------------- #
@cli.group(context_settings=CTX)
def track():
    """Split, join and slice chromatogram trace files."""


@track.command("join", no_args_is_help=True, context_settings=CTX)
@click.argument("files", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("--outbase", "-b", required=True, help="Output basename")
@click.option("--outdir", "-o", default=".", help="Output directory")
@click.option("--gap", "-g", default=24, type=int, help="Gap samples between reads")
def track_join(files, outbase, outdir, gap):
    """Join multiple ABI chromatograms end to end into one record."""
    from cfutils.tracks import export_tracks, join_tracks

    joined = join_tracks(*[__parse_abi(f) for f in files], name=outbase, gap=gap)
    click.echo(f"✓ Joined {len(files)} tracks into {len(joined.seq)} bases")
    for p in export_tracks([joined], outdir, fmt="npz"):
        click.echo(f"  wrote {p}")


@track.command("split", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option(
    "--cuts", "-c", required=True, help="Comma-separated 1-based cut positions"
)
@click.option("--outdir", "-o", default=".")
@click.option("--fmt", "-f", default="npz", type=click.Choice(["npz", "tsv"]))
def track_split(file, cuts, outdir, fmt):
    """Split one ABI chromatogram into segments at the given cut positions."""
    from cfutils.tracks import export_tracks, split_track

    cut_positions = [int(c) for c in cuts.split(",") if c.strip()]
    pieces = split_track(__parse_abi(file), cut_positions)
    for p in export_tracks(pieces, outdir, fmt=fmt):
        click.echo(f"✓ wrote {p}")


@track.command("slice", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--start", "-s", required=True, type=int)
@click.option("--end", "-e", required=True, type=int)
@click.option("--outdir", "-o", default=".")
@click.option("--fmt", "-f", default="npz", type=click.Choice(["npz", "tsv"]))
def track_slice(file, start, end, outdir, fmt):
    """Extract a 1-based inclusive region [start, end] from a chromatogram."""
    from cfutils.tracks import export_tracks, slice_track

    seg = slice_track(__parse_abi(file), start, end)
    for p in export_tracks([seg], outdir, fmt=fmt):
        click.echo(f"✓ wrote {p}")


@cli.group(context_settings=CTX)
def edit():
    """Trim reads and strip primers (keeps traces aligned)."""


@edit.command("trim", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--cutoff", "-c", default=0.05, type=float)
@click.option("--outdir", "-o", default=".")
def edit_trim(file, cutoff, outdir):
    """Mott quality-trim the read (keeps peak/trace axes aligned)."""
    from cfutils.transform import trim

    _emit(trim(__parse_abi(file), cutoff=cutoff), outdir)


@edit.command("trim-ends", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--min-qual", "-q", default=20, type=int)
@click.option("--outdir", "-o", default=".")
def edit_trim_ends(file, min_qual, outdir):
    """Hard-trim low-quality 5' and 3' ends."""
    from cfutils.transform import trim_ends

    _emit(trim_ends(__parse_abi(file), min_qual=min_qual), outdir)


@edit.command("strip-primers", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--forward", "-f", default="", help="5' primer sequence")
@click.option("--reverse", "-r", default=None, help="3' primer sequence")
@click.option("--outdir", "-o", default=".")
def edit_strip_primers(file, forward, reverse, outdir):
    """Remove primer sequences from the read ends."""
    from cfutils.transform import strip_primers

    _emit(strip_primers(__parse_abi(file), forward=forward, reverse=reverse), outdir)


@edit.command("reverse", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--outdir", "-o", default=".")
def edit_reverse(file, outdir):
    """Reverse-complement the whole chromatogram record."""
    from cfutils.transform import reverse_complement_record

    _emit(reverse_complement_record(__parse_abi(file)), outdir)


# --------------------------------------------------------------------------- #
#  Base Calling & Assembly
# --------------------------------------------------------------------------- #
@cli.group(context_settings=CTX)
def basecall():
    """Re-call bases from raw four-channel traces."""


@basecall.command("call", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option(
    "--hetero-threshold",
    "-r",
    default=0.45,
    type=float,
    help="second-peak ratio for ambiguity calling (0 disables)",
)
@click.option("--outdir", "-o", default=".")
def basecall_call(file, hetero_threshold, outdir):
    """Re-call bases from raw traces, print sequence, save FASTA."""
    from pathlib import Path

    from cfutils.basecaller import call_bases

    rec = __parse_abi(file, rescale=False)
    res = call_bases(rec, hetero_threshold=hetero_threshold)
    click.echo(f">{rec.name}  n={res.n_calls} ambiguous={res.n_ambiguous}")
    click.echo(res.sequence)
    outdir = Path(outdir or ".")
    outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / f"{rec.name}_recalled.fa"
    p.write_text(f">{rec.name}\n{res.sequence}\n")
    click.echo(f"✓ wrote {p}")


@basecall.command("hetero", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--min-ratio", "-r", default=0.45, type=float)
def basecall_hetero(file, min_ratio):
    """List mixed/heterozygous bases detected in the traces."""
    from cfutils.basecaller import call_bases

    rec = __parse_abi(file, rescale=False)
    res = call_bases(rec, hetero_threshold=min_ratio)
    click.echo("pos\tbase\t2nd\tratio\tQ")
    for c in res.calls:
        if c.is_ambiguous:
            click.echo(
                f"{c.position}\t{c.base}\t{c.second_base}"
                f"\t{c.second_ratio:.2f}\t{c.quality}"
            )


@cli.group(context_settings=CTX)
def assemble():
    """Reference-guided pileup and consensus calling."""


@assemble.command("consensus", no_args_is_help=True, context_settings=CTX)
@click.argument("reads", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("--reference", "-r", required=True, type=click.Path(exists=True))
@click.option("--min-cov", "-c", default=1, type=int)
@click.option("--min-freq", "-f", default=0.5, type=float)
@click.option("--quality", "-q", default=0, type=int, help="drop bases below Q")
def assemble_consensus(reads, reference, min_cov, min_freq, quality):
    """Pile up reads on a reference and print the consensus sequence."""
    from cfutils.assembly import consensus, pileup

    table = pileup(
        [__parse_abi(f) for f in reads],
        __parse_fasta(reference),
        quality_threshold=quality,
        min_cov=min_cov,
    )
    click.echo(consensus(table, min_freq=min_freq))


# --------------------------------------------------------------------------- #
#  Sequence Analysis
# --------------------------------------------------------------------------- #
@cli.group(context_settings=CTX)
def analyze():
    """Sequence-level biology on a called sequence."""


@analyze.command("translate", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--frame", "-f", default=1, type=int)
def analyze_translate(file, frame):
    """Translate the called sequence in a reading frame."""
    from cfutils.analysis import translate

    click.echo(translate(__parse_abi(file).seq, frame=frame))


@analyze.command("motif", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--motif", "-m", required=True)
@click.option("--both-strands/--forward-only", default=False)
def analyze_motif(file, motif, both_strands):
    """Find motif occurrence positions."""
    from cfutils.analysis import find_motifs

    click.echo(
        "\t".join(
            map(
                str,
                find_motifs(__parse_abi(file).seq, motif, both_strands=both_strands),
            )
        )
    )


@analyze.command("rest", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
def analyze_rest(file):
    """Scan for common restriction sites."""
    from cfutils.analysis import restriction_sites

    for name, pos in restriction_sites(__parse_abi(file).seq).items():
        if pos:
            click.echo(f"{name}\t{','.join(map(str, pos))}")


@analyze.command("gc", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--window", "-w", default=30, type=int)
@click.option("--step", "-s", default=1, type=int)
def analyze_gc(file, window, step):
    """Sliding GC% windows (position<TAB>gc%)."""
    from cfutils.analysis import gc_windows

    for pos, gc in gc_windows(__parse_abi(file).seq, window, step):
        click.echo(f"{pos}\t{gc}")


# --------------------------------------------------------------------------- #
#  Export & Visualization
# --------------------------------------------------------------------------- #
@cli.group(context_settings=CTX)
def export():
    """Export results to standard lab / pipeline formats."""


@export.command("fasta", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--outdir", "-o", default=".")
def export_fasta(file, outdir):
    """Write the called sequence as FASTA."""
    rec = __parse_abi(file)
    _emit(rec, outdir)


@export.command("vcf", no_args_is_help=True, context_settings=CTX)
@click.option("--query", "-q", required=True, type=click.Path(exists=True))
@click.option("--subject", "-s", required=True, type=click.Path(exists=True))
@click.option("--min-qual", "-m", default=0, type=int)
def export_vcf(query, subject, min_qual):
    """Emit variant calls against a reference as VCF."""
    from pathlib import Path

    from cfutils.align import call_mutations
    from cfutils.export import to_vcf

    sites = call_mutations(
        __parse_abi(query), __parse_fasta(subject), report_all_sites=True
    )
    click.echo(
        to_vcf(sites, reference_name=Path(subject).stem, min_qual=min_qual), nl=False
    )


@export.command("json", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
def export_json(file):
    """Emit a self-describing JSON QC/analysis report."""
    from cfutils.export import to_json

    click.echo(to_json(__parse_abi(file)))


@export.command("batch", no_args_is_help=True, context_settings=CTX)
@click.argument("files", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("--outdir", "-o", default=".")
@click.option("--fmt", "-f", default="csv", type=click.Choice(["csv", "tsv", "json"]))
def export_batch(files, outdir, fmt):
    """Write a QC summary table for many reads."""
    from cfutils.export import write_batch

    click.echo(
        "✓ wrote " + write_batch([__parse_abi(f) for f in files], outdir, fmt=fmt)
    )


@cli.group(context_settings=CTX)
def plot():
    """Render chromatograms (optionally with features / DNA viewer)."""


@plot.command("chrom", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--start", "-s", default=None, type=int)
@click.option("--end", "-e", default=None, type=int)
@click.option("--out", "-o", default="chrom.png")
def plot_chrom(file, start, end, out):
    """Plot a chromatogram region."""
    import matplotlib.pyplot as plt

    from cfutils.show import plot_chromatograph

    region = (start, end if end is not None else start + 30) if start else None
    fig, ax = plt.subplots(figsize=(16, 5))
    plot_chromatograph(__parse_abi(file), region=region, ax=ax)
    fig.savefig(out, bbox_inches="tight")
    click.echo(f"✓ wrote {out}")


@plot.command("features", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--start", "-s", default=None, type=int)
@click.option("--end", "-e", default=None, type=int)
@click.option("--out", "-o", default="features.png")
def plot_features(file, start, end, out):
    """Plot a chromatogram with attached feature annotations."""
    import matplotlib.pyplot as plt

    from cfutils.features import plot_features as _pf
    from cfutils.show import plot_chromatograph

    region = (start, end if end is not None else start + 30) if start else None
    fig, ax = plt.subplots(figsize=(16, 5))
    plot_chromatograph(__parse_abi(file), region=region, ax=ax)
    _pf(__parse_abi(file), ax)
    fig.savefig(out, bbox_inches="tight")
    click.echo(f"✓ wrote {out}")


@plot.command("dnaviewer", no_args_is_help=True, context_settings=CTX)
@click.argument("file", type=click.Path(exists=True))
@click.option("--start", "-s", default=None, type=int)
@click.option("--end", "-e", default=None, type=int)
@click.option("--out", "-o", default="dna.png")
def plot_dnaviewer(file, start, end, out):
    """Plot a chromatogram with a DNA Features Viewer map (optional dep)."""
    from cfutils.dnalink import plot_combined

    region = (start, end if end is not None else start + 30) if start else None
    fig, (ax_feat, ax_chrom) = plot_combined(__parse_abi(file), region=region)
    fig.savefig(out, bbox_inches="tight")
    click.echo(f"✓ wrote {out}")


# --------------------------------------------------------------------------- #
#  Internal parsing helpers (avoid circular import at module load)
# --------------------------------------------------------------------------- #
def __parse_abi(file, rescale=True):
    from cfutils.parser import parse_abi

    return parse_abi(file, rescale=rescale)


def __parse_fasta(file):
    from cfutils.parser import parse_fasta

    return parse_fasta(file)


if __name__ == "__main__":
    cli()
