#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate the README example figures from the bundled real Sanger data.

The examples are written against the **high-level :class:`Chromatogram`
object** (see ``cfutils.chromatogram``), so the images match the README usage.

Outputs under ``examples/``:

* ``mutation_call.png``    - chromatogram with called variants highlighted
* ``quality_profile.png``  - per-base Phred quality + trimming
* ``side_by_side.png``     - chromatogram + GC% + quality panels
* ``feature_overlay.png``  - feature annotation (primers/amplicon/SNPs)
* ``basecall_hetero.png``  - re-called bases with mixed/ambiguous sites
* ``dna_features.png``     - DNA Features Viewer map + chromatogram
* ``assembly.png``         - pileup depth + consensus vs a reference

Run:  python -m scripts.make_readme_examples
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cfutils import Chromatogram, style
from cfutils.align import call_mutations
from cfutils.composite import _trace_data
from cfutils.features import ChromatogramFeature, plot_features
from cfutils.parser import parse_fasta
from cfutils.show import annotate_mutation, center_region, highlight_base

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"
OUT = HERE / "examples"

style.apply()


def _panel(ax, cg, region, func):
    seg_x, sel_peaks, sel_seq, start, _e = _trace_data(cg.to_record, region)
    func(ax, sel_peaks, sel_seq, cg, start)


def gc_panel(ax, peaks, seq, cg, start):
    n = len(seq)
    arr = [
        100.0
        * (
            seq[max(0, i - 1) : i + 2].count("G")
            + seq[max(0, i - 1) : i + 2].count("C")
        )
        / 3
        for i in range(n)
    ]
    ax.fill_between(peaks[:n], arr, color="#2ca02c", alpha=0.35)
    ax.set_ylim(0, 100)
    ax.set_ylabel("GC%")


def qual_panel(ax, peaks, seq, cg, start):
    q = cg.quality
    lo = (start or 1) - 1
    ax.plot(peaks, q[lo : lo + len(peaks)], color="#d62728", lw=1.2)
    ax.axhline(20, color="grey", ls="--", lw=0.8)
    ax.set_ylim(0, 65)
    ax.set_ylabel("Phred Q")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cg = Chromatogram.from_abi(str(DATA / "B5-M13R_B07.ab1"))
    ref = parse_fasta(str(DATA / "ref.fa"))

    all_calls = call_mutations(cg.to_record, ref, report_all_sites=True)
    mutations = [s for s in all_calls if s.ref_base != s.cf_base]
    snps = [
        s
        for s in mutations
        if "N" not in (s.ref_base + s.cf_base) and "-" not in (s.ref_base + s.cf_base)
    ]

    # 1. mutation call (Chromatogram.plot + annotate)
    region = center_region(snps[0].cf_pos, 20, cg.length)
    fig, ax = plt.subplots(figsize=(15, 4.2))
    cg.plot(region=region, ax=ax)
    for mut in mutations:
        if region[0] - 6 <= mut.cf_pos <= region[1] + 6:
            highlight_base(mut.cf_pos, cg.to_record, ax, passed_filter=True)
            annotate_mutation(mut, cg.to_record, ax)
    snp = snps[0]
    ax.set_title(
        f"Mutation calling ({snp.ref_base}{snp.ref_pos}{snp.cf_base} highlighted)",
        loc="left",
    )
    fig.savefig(OUT / "mutation_call.png")
    plt.close(fig)

    # 2. quality profile + trimming (Chromatogram quality/qc)
    q = cg.quality
    m = cg.qc()
    fig, ax = plt.subplots(figsize=(15, 3.6))
    ax.plot(np.arange(1, len(q) + 1), q, color="#1f77b4", lw=1.1)
    ax.fill_between(np.arange(1, len(q) + 1), q, 0, alpha=0.12, color="#1f77b4")
    ax.axhline(20, color="orange", ls="--", lw=1)
    ax.axvspan(
        m["trim_start"],
        m["trim_end"],
        color="green",
        alpha=0.10,
        label="Mott-trimmed region",
    )
    ax.set_xlabel("Read position")
    ax.set_ylabel("Phred Q")
    ax.set_title(
        "Quality profile (mean Q=%.1f, CRL=%d)" % (cg.mean_quality, m["crl"]),
        loc="left",
    )
    ax.legend(loc="lower right")
    fig.savefig(OUT / "quality_profile.png")
    plt.close(fig)

    # 3. side-by-side composite (Chromatogram traces via cg.plot)
    region = (55, 85)
    fig, (axc, axg, axq) = plt.subplots(
        3, 1, figsize=(15, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1, 1]}
    )
    cg.plot(region=region, ax=axc)
    _panel(axg, cg, region, gc_panel)
    _panel(axq, cg, region, qual_panel)
    axc.set_title("Chromatogram + GC% + quality (side-by-side)", loc="left")
    fig.savefig(OUT / "side_by_side.png")
    plt.close(fig)

    # 4. feature overlay (ChromatogramFeature)
    r = center_region(215, 32, cg.length)
    fig, ax = plt.subplots(figsize=(15, 4.5))
    cg.plot(region=r, ax=ax)
    feats = [ChromatogramFeature.from_sitepair(s, color="#d62728") for s in snps]
    feats.append(
        ChromatogramFeature(
            start=r[0] + 1, end=r[1] - 1, strand=+1, color="#7fbf7b", label="amplicon"
        )
    )
    plot_features(cg.to_record, ax, features=feats, show_legend=False)
    ax.set_title("Feature overlay (primers/amplicon/SNPs)", loc="left")
    fig.savefig(OUT / "feature_overlay.png")
    plt.close(fig)

    # 5. basecall + heterozygotes (Chromatogram.basecall)
    res = cg.basecall()
    region = (1, 40)
    fig, ax = plt.subplots(figsize=(15, 4.2))
    plot_rec = cg.to_record
    from cfutils.show import plot_chromatograph

    plot_chromatograph(plot_rec, region=region, ax=ax)
    peaks = plot_rec.annotations["peak positions"]
    for c in res.calls:
        if c.is_ambiguous and region[0] <= c.position <= region[1]:
            # map the 1-based called position to its rescaled trace x
            ax.text(
                peaks[c.position - 1],
                1.02,
                f"{c.base}",
                color="purple",
                fontsize="x-large",
                fontweight="bold",
                ha="center",
            )
    ax.set_title("Re-called bases (ambiguous/mixed sites in purple)", loc="left")
    fig.savefig(OUT / "basecall_hetero.png")
    plt.close(fig)

    # 6. DNA Features Viewer combined
    from cfutils.dnalink import plot_combined

    feats = [ChromatogramFeature.from_sitepair(s, color="#d62728") for s in snps]
    feats.append(
        ChromatogramFeature(
            start=1, end=cg.length, strand=+1, color="#7fbf7b", label="amplicon"
        )
    )
    fig, (axf, axc2) = plot_combined(
        cg.to_record,
        features=feats,
        region=center_region(215, 25, cg.length),
        plot_sequence=True,
    )
    axf.set_title("DNA Features Viewer map + chromatogram", loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "dna_features.png")
    plt.close(fig)

    # 7. assembly (pileup depth + consensus)
    from cfutils.assembly import consensus, coverage, pileup

    table = pileup([cg.to_record], ref)
    cov = coverage(table)
    positions = np.array([p for p, _ in cov])
    depths = np.array([d for _, d in cov])
    cons = consensus(table)
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(15, 4.5), sharex=True, gridspec_kw={"height_ratios": [1, 2]}
    )
    ax1.bar(positions, depths, width=1, color="#3182bd")
    ax1.set_ylabel("depth")
    ax1.set_title("Reference-guided pileup & consensus", loc="left")
    refs = [c.ref_base for p, c in sorted(table.columns.items())]
    ax2.bar(
        np.arange(len(cons)),
        np.ones(len(cons)),
        color=["#d95f0e" if a != b else "#c6dbef" for a, b in zip(refs, cons)],
        width=1,
    )
    ax2.set_yticks([])
    ax2.set_xlabel("reference position")
    fig.tight_layout()
    fig.savefig(OUT / "assembly.png")
    plt.close(fig)

    print(f"wrote {len(list(OUT.glob('*.png')))} example figures to {OUT}")


if __name__ == "__main__":
    main()
