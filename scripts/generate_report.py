#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Full end-to-end analysis + figure generation on the bundled real Sanger data.

Runs the whole cfutils pipeline on ``data/B5-M13R_B07.ab1`` vs
``data/ref.fa``, produces a set of publication-style figures under
``docs/figures/``, and writes:

* ``docs/analysis.md`` -- Jekyll (Just-the-Docs) markdown page embedding them.
* ``docs/report.html``  -- standalone, self-contained HTML report with the
  figures base64-embedded (portable, shareable without the repo).

Run:  python -m scripts.generate_report
"""

from __future__ import annotations

import base64
import os
from datetime import datetime
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cfutils.align import align_chromatograph, call_mutations
from cfutils.assembly import pileup, consensus, coverage
from cfutils.composite import side_by_side
from cfutils.features import ChromatogramFeature, plot_features
from cfutils.parser import parse_abi, parse_fasta
from cfutils.quality import QualityFilter
from cfutils.show import (
    annotate_mutation,
    center_region,
    highlight_base,
    plot_chromatograph,
)

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"
FIGDIR = HERE / "docs" / "figures"
QUERY_AB1 = DATA / "B5-M13R_B07.ab1"
SUBJECT_FA = DATA / "ref.fa"

plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight"})


# --------------------------------------------------------------------------- #
#  Analysis helpers
# --------------------------------------------------------------------------- #
def nuc_freq(seq: str) -> dict:
    return {b: seq.count(b) for b in "ACGT"}


def sliding_gc(seq: str, window: int = 30, step: int = 1):
    """Sliding GC% as (x_positions, gc_values)."""
    xs, gc = [], []
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i : i + window]
        xs.append(i + window / 2)
        gc.append(100.0 * (w.count("G") + w.count("C")) / window)
    return xs, gc


def gc_panel(ax, trace_x, peaks, seq, record, start_pos=None):
    """Minimal user panel function for the composite plot (GC%)."""
    n = len(seq)
    if n < 4 or peaks is None or len(peaks) == 0:
        ax.set_yticks([])
        return
    gc = [100.0 * (seq[max(0, i - 1):i + 2].count("G")
                   + seq[max(0, i - 1):i + 2].count("C")) / 3
          for i in range(n)]
    ax.fill_between(peaks[:n], gc, color="#2ca02c", alpha=0.35)
    ax.set_ylim(0, 100)
    ax.set_ylabel("GC %")


def qual_panel(ax, trace_x, peaks, seq, record, start_pos=None):
    """User panel drawing per-base quality aligned to the region trace."""
    qual = record.letter_annotations.get("phred_quality", [])
    if not qual:
        return
    start_pos = start_pos or 1
    lo = start_pos - 1
    q = np.asarray(qual[lo : lo + len(peaks)])
    ax.plot(peaks, q, color="#d62728", lw=1.2)
    ax.axhline(20, color="grey", ls="--", lw=0.8)
    ax.axhline(50, color="grey", ls=":", lw=0.8)
    ax.set_ylim(0, 65)
    ax.set_ylabel("Phred Q")


# --------------------------------------------------------------------------- #
#  Figures
# --------------------------------------------------------------------------- #
def fig_quality_profile(rec, out):
    """Per-base Phred quality across the whole read."""
    q = np.asarray(rec.letter_annotations["phred_quality"])
    x = np.arange(1, len(q) + 1)
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(x, q, color="#1f77b4", lw=1.1)
    ax.fill_between(x, q, 0, alpha=0.12, color="#1f77b4")
    ax.axhline(20, color="orange", ls="--", lw=1, label="low-quality threshold (Q20)")
    ax.axhline(50, color="grey", ls=":", lw=1, label="high-quality threshold (Q50)")
    ax.set_xlabel("Read position (1-based)")
    ax.set_ylabel("Phred quality")
    ax.set_ylim(0, np.max(q) + 5)
    ax.set_title(f"Per-base quality profile of {rec.name} (n={len(q)})")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    fig.savefig(out)
    plt.close(fig)


def fig_chromatogram_mutations(rec, mutations, region, out, title):
    """Chromatogram region with quality-filtered mutations highlighted."""
    fig, ax = plt.subplots(figsize=(15, 4.5))
    plot_chromatograph(rec, region=region, ax=ax)
    for mut in mutations:
        if region[0] - 6 <= mut.cf_pos <= region[1] + 6:
            highlight_base(mut.cf_pos, rec, ax, passed_filter=True)
            annotate_mutation(mut, rec, ax)
    ax.set_title(title, loc="left")
    fig.savefig(out)
    plt.close(fig)


def fig_composite(rec, out):
    """Side-by-side chromatogram + GC% + quality panels (composite API)."""
    region = (55, 85)
    fig, (ax_chrom, ax_gc, ax_qual) = plt.subplots(
        3, 1,
        figsize=(15, 6.5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1, 1]},
    )
    plot_chromatograph(rec, region=region, ax=ax_chrom)
    # build panels using the same resolved region x so everything aligns
    _panel_on_axes(ax_gc, rec, region, gc_panel)
    _panel_on_axes(ax_qual, rec, region, qual_panel)
    ax_chrom.set_title("Chromatogram + influence panels (Side-by-side / composite)", loc="left")
    fig.savefig(out)
    plt.close(fig)


def _panel_on_axes(ax, rec, region, func):
    from cfutils.composite import _trace_data
    seg_x, sel_peaks, sel_seq, start_pos, _e = _trace_data(rec, region)
    func(ax, seg_x, sel_peaks, sel_seq, rec, start_pos)


def fig_feature_overlay(rec, mutations, out):
    """Feature-annotation overlay on a chromatogram (features module)."""
    region = center_region(mutations[1].cf_pos, 35, len(rec))
    fig, ax = plt.subplots(figsize=(15, 4.8))
    plot_chromatograph(rec, region=region, ax=ax)
    feats = [
        ChromatogramFeature(start=region[0] + 1, end=region[1] - 1, strand=+1,
                            color="#7fbf7b", label="amplicon"),
    ]
    for mut in mutations[:4]:
        feats.append(
            ChromatogramFeature.from_sitepair(mut, color="#d62728"))
    plot_features(rec, ax, features=feats, show_legend=False)
    ax.set_title("Feature annotation overlay (ChromatogramFeature API)", loc="left")
    fig.savefig(out)
    plt.close(fig)


def fig_lollipop(rec, mutation_snps, all_sites, out):
    """Mutation lollipop plot: quality vs reference position, marking SNPs."""
    snp_pos = {s.ref_pos: s for s in mutation_snps}
    xs, ys = [], []
    for s in all_sites:
        if s.qual_site is None:
            continue
        xs.append(s.ref_pos)
        ys.append(s.qual_site)
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.scatter(xs, ys, s=10, color="#9ecae1", label="all aligned sites", zorder=2)
    for s in mutation_snps:
        ax.scatter([s.ref_pos], [s.qual_site], s=60, color="#d62728",
                   zorder=3, edgecolor="black")
        ax.vlines(s.ref_pos, 0, s.qual_site, color="#d62728", alpha=0.6, lw=1.5)
    ax.axhline(20, color="orange", ls="--")
    ax.set_xlabel("Reference position")
    ax.set_ylabel("Site Phred quality")
    ax.set_title("SNP lollipop: high-quality variants across the read",
                 loc="left")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    fig.savefig(out)
    plt.close(fig)


def fig_assembly(rec, ref, out):
    """Pileup depth + consensus along the reference (assembly module)."""
    table = pileup([rec], ref)
    cov = coverage(table)
    positions = np.array([p for p, _ in cov])
    depths = np.array([d for _, d in cov])
    cons = consensus(table)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(14, 5), sharex=True,
        gridspec_kw={"height_ratios": [1, 2]},
    )
    ax1.bar(positions, depths, width=1, color="#3182bd")
    ax1.set_ylabel("Read depth")
    ax1.set_ylim(0, max(1, depths.max() + 0.5))
    ax1.set_title("Reference-guided pileup & consensus (assembly module)",
                  loc="left")

    # consensus vs reference: colour bases that differ (including N-fill)
    refs = np.array([c.ref_base for p, c in sorted(table.columns.items())])
    cons_b = np.array(list(cons))
    colors = np.where(refs == cons_b, "#c6dbef", "#d95f0e")
    ax2.bar(np.arange(len(cons_b)), np.ones(len(cons_b)),
            color=colors, width=1)
    for i, b in enumerate(cons_b):
        ax2.text(i, 0.5, b, ha="center", va="center", fontsize=6,
                 family="monospace")
    ax2.set_yticks([])
    ax2.set_ylabel("consensus")
    ax2.set_xlabel("Reference position (sorted pileup columns)")
    ax2.set_xlim(-0.5, len(cons_b) - 0.5)
    ax2.set_title(f"consensus length {len(cons_b)}", fontsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_dna_viewer(rec, snps, out):
    """Combined DNA-features-viewer map + chromatogram (dnalink module)."""
    from cfutils.dnalink import plot_combined
    from cfutils.show import center_region
    feats = [ChromatogramFeature.from_sitepair(s, color="#d62728")
             for s in snps]
    # forward / reverse primer arrows and the amplicon
    feats += [
        ChromatogramFeature(start=1, end=len(rec), strand=+1, color="#7fbf7b",
                            label="amplicon", kind="region"),
        ChromatogramFeature(start=90, end=130, strand=+1, color="#1f77b4",
                            label="primer F", kind="primer"),
        ChromatogramFeature(start=190, end=230, strand=-1, color="#ff7f0e",
                            label="primer R", kind="primer"),
    ]
    region = center_region(215, 30, len(rec))  # around the SNP cluster
    fig, (ax_feat, ax_chrom) = plot_combined(
        rec, features=feats, region=region, plot_sequence=True)
    ax_feat.set_title(
        "DNA Features Viewer map + cfutils chromatogram (combined)", loc="left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_basecaller_compare(rec, out):
    """Compare vendor per-base quality vs the re-called basecaller quality."""
    from cfutils.parser import parse_abi as _p
    from cfutils.basecaller import call_bases
    raw = _p(str(QUERY_AB1), rescale=False)
    res = call_bases(raw)
    vend = np.asarray(rec.letter_annotations["phred_quality"])
    mine = np.asarray(res.qualities)
    n = min(len(vend), len(mine))
    x = np.arange(1, n + 1)
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(x, vend[:n], color="#1f77b4", lw=1.0, label="vendor quality")
    ax.plot(x, mine[:n], color="#d62728", lw=1.0, alpha=0.8, label="cfutils re-call")
    ax.axhline(20, color="orange", ls="--", lw=1)
    ax.set_xlabel("Read position")
    ax.set_ylabel("Phred quality")
    ax.set_title("Vendor quality vs cfutils base-caller quality", loc="left")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    fig.savefig(out)
    plt.close(fig)


def fig_base_composition(rec, out):
    """Nucleotide composition of the called read (analytical)."""
    freq = nuc_freq(rec.seq)
    labels = ["A", "C", "G", "T"]
    vals = [freq[l] for l in labels]
    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, vals, color=["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 2, f"{v}\n({100*v/len(rec.seq):.1f}%)",
                ha="center", fontsize=9)
    ax.set_ylabel("Count")
    ax.set_title("Called-base composition of the read", loc="left")
    ax.set_ylim(0, max(vals) * 1.2)
    fig.savefig(out)
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  Report writers
# --------------------------------------------------------------------------- #
def write_markdown(context, path):
    ctx = context
    md = f"""---
title: Real-data analysis walkthrough
layout: home
nav_order: 8
---

# Analysis of a real Sanger trace

This page walks through the full cfutils pipeline on the bundled real data
(`data/B5-M13R_B07.ab1` vs `data/ref.fa`), generated by
`python -m scripts.generate_report`.

## Sample summary

| Metric | Value |
|---|---|
| Read name | `{ctx['name']}` |
| Called bases | {ctx['n_bases']} |
| Mean Phred quality | {ctx['mean_qual']:.1f} |
| Min Phred quality | {ctx['min_qual']} |
| N fraction | {100*ctx['n_fraction']:.1f}% |
| Low-quality (<Q20) fraction | {100*ctx['low_qual_fraction']:.1f}% |
| Mott-trimmed interval | {ctx['trim_start']}–{ctx['trim_end']} ({ctx['trimmed_len']} bp) |
| Reads aligned to ref | {ctx['n_aligned']} |
| Base-caller re-call accuracy | {ctx['bc_acc']:.1f}% (n={ctx['bc_calls']}) |
| Mean re-called quality | {ctx['bc_mean_qual']:.1f} |
| Quality-filtered, real SNPs | {ctx['n_snps']} |

## Results

The read is high-quality overall (mean Q{ctx['mean_qual']:.0f}); the dye
primers degrade at the 3' end (see the quality profile).  A python re-call from
the raw four-channel traces reproduces the vendor base call with
**{ctx['bc_acc']:.1f}%** accuracy.  Comparing against the reference, **{ctx['n_snps']}**
high-confidence single-nucleotide variants are called:

| Ref pos | Change | Site Q | Local Q |
|---|---|---|---|
""" + "".join(
        f"| {s.ref_pos} | {s.ref_base}→{s.cf_base} | {s.qual_site} | {s.qual_local} |\n"
        for s in ctx['snps']
    ) + f"""
## Figures

### Per-base quality profile
![quality profile](figures/{ctx['figs']['quality']})
The 5' ~30 bp and 3' tail drop below Q20, while the interior is Q50+.

### Variant chromatograms
Mutations called by alignment are highlighted and annotated on the trace.

![double SNP at ref 60/61](figures/{ctx['figs']['chrom_at601']})
![SNP cluster ~169-177](figures/{ctx['figs']['chrom_cluster']})

### Side-by-side (composite) plotting
Chromatogram with GC% and per-base quality panels sharing the same x axis —
this is how external tools can plot their own signal alongside cfutils.

![composite](figures/{ctx['figs']['composite']})

### Feature overlay API
![feature overlay](figures/{ctx['figs']['features']})

### SNP lollipop
![lollipop](figures/{ctx['figs']['lollipop']})

### Pileup & consensus (assembly)
![assembly](figures/{ctx['figs']['assembly']})

### Base composition
![composition](figures/{ctx['figs']['composition']})

### DNA Features Viewer integration
![dna viewer](figures/{ctx['figs']['dna_viewer']})
The same chromatogram rendered together with a
[DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)
feature map (top) and the called sequence, sharing one x axis.

### Base-caller benchmark
![basecaller](figures/{ctx['figs']['basecall']})
A Python re-call of the raw four-channel traces reproduces the vendor base call
with **{ctx['bc_acc']:.1f}%** accuracy (mean Q{ctx['bc_mean_qual']:.0f}).

> Figures auto-generated; re-run `python -m scripts.generate_report` to refresh.
"""
    path.write_text(md)
    print(f"wrote {path}")


def write_html(context, path, figure_paths):
    rows = "".join(
        f"<tr><td>{s.ref_pos}</td><td>{s.ref_base}&rarr;{s.cf_base}</td>"
        f"<td>{s.qual_site}</td><td>{s.qual_local}</td></tr>"
        for s in context["snps"]
    )
    figures = ""
    labels = [
        ("quality", "Per-base quality profile"),
        ("chrom_at601", "Double SNP at ref 60/61"),
        ("chrom_cluster", "SNP cluster ~169-177"),
        ("composite", "Side-by-side composite"),
        ("features", "Feature annotation overlay"),
        ("lollipop", "SNP lollipop plot"),
        ("assembly", "Pileup & consensus"),
        ("composition", "Base composition"),
        ("dna_viewer", "DNA Features Viewer + chromatogram (combined)"),
        ("basecall", "Vendor quality vs cfutils base-caller quality"),
    ]
    for key, label in labels:
        b64 = base64.b64encode(figure_paths[key].read_bytes()).decode()
        figures += f'<figure><img src="data:image/png;base64,{b64}" alt="{label}"/>' \
                   f"<figcaption>{label}</figcaption></figure>\n"

    ctx = context
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>cfutils — real-data analysis report</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
      margin:2rem auto;max-width:1000px;color:#222;line-height:1.5}}
 h1,h2{{border-bottom:2px solid #eee;padding-bottom:.3rem}}
 table{{border-collapse:collapse;margin:1rem 0}}
 td,th{{border:1px solid #ddd;padding:.4rem .7rem}}
 figure{{margin:1.5rem 0;text-align:center}}
 img{{max-width:100%;border:1px solid #eee;border-radius:4px}}
 figcaption{{font-size:.85rem;color:#666;margin-top:.3rem}}
 code{{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px}}
</style></head><body>
<h1>cfutils · real-data analysis report</h1>
<p>Generated {ctx['date']} by <code>scripts/generate_report.py</code>.</p>
<h2>Sample</h2>
<table>
<tr><th>Read</th><td><code>{ctx['name']}</code></td></tr>
<tr><th>Called bases</th><td>{ctx['n_bases']}</td></tr>
<tr><th>Mean Phred quality</th><td>{ctx['mean_qual']:.1f}</td></tr>
<tr><th>Min Phred quality</th><td>{ctx['min_qual']}</td></tr>
<tr><th>N fraction</th><td>{100*ctx['n_fraction']:.1f}%</td></tr>
<tr><th>Low-quality (&lt;Q20)</th><td>{100*ctx['low_qual_fraction']:.1f}%</td></tr>
<tr><th>Mott-trimmed interval</th><td>{ctx['trim_start']}–{ctx['trim_end']} ({ctx['trimmed_len']} bp)</td></tr>
<tr><th>Aligned sites</th><td>{ctx['n_aligned']}</td></tr>
<tr><th>Raw ref&ne;cf sites</th><td>{ctx['raw_mut']}</td></tr>
<tr><th>Base-caller accuracy</th><td>{ctx['bc_acc']:.1f}% (n={ctx['bc_calls']})</td></tr>
<tr><th>High-confidence real SNPs</th><td>{ctx['n_snps']}</td></tr>
</table>
<h2>High-confidence variants</h2>
<table><tr><th>Ref pos</th><th>Change</th><th>Site Q</th><th>Local Q</th></tr>
{rows}</table>
<h2>Figures</h2>
{figures}
</body></html>"""
    path.write_text(html)
    print(f"wrote {path} ({len(html)//1024} KB)")


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)

    rec = parse_abi(str(QUERY_AB1))
    ref = parse_fasta(str(SUBJECT_FA))

    all_sites = align_chromatograph(rec, ref)
    all_calls = call_mutations(rec, ref, report_all_sites=True)
    raw_mut = [s for s in all_calls if s.ref_base != s.cf_base]
    qf = QualityFilter(min_base_qual=20, min_local_qual=20, flank_base_num=5)
    snps = [
        s for s in raw_mut
        if not {"N", "-"} & {s.ref_base, s.cf_base}
        and qf.passed(s)
    ]

    # ---- figures -----------------------------------------------------------
    figs = {}
    figs["quality"] = FIGDIR / "quality_profile.png"
    figs["chrom_at601"] = FIGDIR / "chromatogram_snp_60_61.png"
    figs["chrom_cluster"] = FIGDIR / "chromatogram_cluster_169_177.png"
    figs["composite"] = FIGDIR / "composite_side_by_side.png"
    figs["features"] = FIGDIR / "feature_overlay.png"
    figs["lollipop"] = FIGDIR / "snp_lollipop.png"
    figs["assembly"] = FIGDIR / "assembly_pileup.png"
    figs["composition"] = FIGDIR / "base_composition.png"
    figs["dna_viewer"] = FIGDIR / "dna_viewer_combined.png"
    figs["basecall"] = FIGDIR / "basecaller_quality.png"

    fig_quality_profile(rec, figs["quality"])
    fig_chromatogram_mutations(
        rec, raw_mut, center_region(snps[0].cf_pos, 20, len(rec)),
        figs["chrom_at601"], "SNPs ref60_61 highlighted",
    )
    fig_chromatogram_mutations(
        rec, raw_mut, center_region(215, 20, len(rec)),
        figs["chrom_cluster"], "SNP cluster ref169-177",
    )
    fig_composite(rec, figs["composite"])
    fig_feature_overlay(rec, snps, figs["features"])
    fig_lollipop(rec, snps, all_calls, figs["lollipop"])
    fig_assembly(rec, ref, figs["assembly"])
    fig_base_composition(rec, figs["composition"])
    fig_dna_viewer(rec, snps, figs["dna_viewer"])
    fig_basecaller_compare(rec, figs["basecall"])

    from collections import Counter
    from cfutils.qc import read_metrics
    from cfutils.basecaller import call_bases, basecaller_score
    from cfutils.parser import parse_abi as _parse_abi_raw

    qual = rec.letter_annotations["phred_quality"]
    mean_qual = float(np.mean(qual))
    qcm = read_metrics(rec)
    raw_rec = _parse_abi_raw(str(QUERY_AB1), rescale=False)
    bc = basecaller_score(call_bases(raw_rec), raw_rec.seq)

    context = {
        "name": rec.name,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_bases": len(rec),
        "mean_qual": mean_qual,
        "min_qual": int(qcm["min_qual"]),
        "n_fraction": qcm["n_fraction"],
        "trim_start": int(qcm["trim_start"]),
        "trim_end": int(qcm["trim_end"]),
        "trimmed_len": int(qcm["trimmed_len"]),
        "gc_percent": qcm["gc_percent"],
        "low_qual_fraction": qcm["low_qual_fraction"],
        "bc_calls": bc["n_calls"],
        "bc_acc": 100 * bc["accuracy"],
        "bc_mean_qual": bc["mean_quality"],
        "n_aligned": len(all_sites),
        "raw_mut": len(raw_mut),
        "n_snps": len(snps),
        "snps": snps,
        "figs": {k: v.name for k, v in figs.items()},
    }

    write_markdown(context, HERE / "docs" / "analysis.md")
    write_html(context, HERE / "docs" / "report.html", figs)

    print(f"\nSaved {len(figs)} figures to {FIGDIR}")


if __name__ == "__main__":
    main()
