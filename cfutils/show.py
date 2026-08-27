#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""
Chromatogram File Utils.

show alignment with matplotlib

author:     Fabio Zanini
date:       09/12/13
content:    Plot functions for Sanger chromatographs.
modified:   By Ye Chang in 2018-05-14
"""

from collections import defaultdict
from typing import Optional, Tuple

from ._mpl import mpl, plt, Axes, HAVE_MPL, require_matplotlib
from .align import SitePair, align_chromatograph
from .parser import SeqRecord  # Import the custom SeqRecord class
from .utils import get_logger, reverse_complement

LOGGER = get_logger(__name__)


def plot_chromatograph(
    seq: SeqRecord,
    region: Optional[Tuple[int, int]] = None,
    ax: Optional[Axes] = None,
    color_map: Optional[dict] = None,
    show_bases: bool = True,
    show_positions: bool = True,
    show_rc: bool = False,
) -> Axes:
    """
    Plot Sanger chromatograph.

    region: include both start and end (1-based)
    """
    require_matplotlib()
    if ax is None:
        ax = plt.gca()
        # _, ax = plt.subplots(1, 1, figsize=(16, 6))

    if seq is None:
        return ax

    if region is None:
        # turn into 0 based for better indexing
        region_start, region_end = 0, len(seq)
    else:
        region_start = max(region[0], 0)
        region_end = min(region[1], len(seq) - 1)

    _colors = defaultdict(lambda: "purple", {"A": "g", "C": "b", "G": "k", "T": "r"})
    if color_map is not None:
        _colors.update(color_map)

    # Get signals
    peaks = seq.annotations["peak positions"]
    trace_x = seq.annotations["trace_x"]
    traces_y = [seq.annotations["channel " + str(i)] for i in range(1, 5)]
    bases = seq.annotations["channels"]

    xlim_left, xlim_right = peaks[region_start] - 1, peaks[region_end] + 0.5

    # Ensure seq is treated as a string
    sequence_str = seq.seq

    # subset peak and sequence
    # TODO: this might fix the bug
    peak_start = peaks[0]
    peak_zip = [
        (p, s)
        for i, (p, s) in enumerate(zip(peaks, sequence_str))
        if region_start <= i <= region_end
    ]
    peaks, sequence_str = list(zip(*peak_zip))

    # subset trace_x and traces_y together
    trace_zip = [
        (x + peak_start, *ys)
        for x, *ys in zip(trace_x, *traces_y)
        if xlim_left <= x <= xlim_right
    ]
    if not trace_zip:
        return ax
    trace_x, *traces_y = list(zip(*trace_zip))

    # Plot traces
    trmax = max(map(max, traces_y))
    for base in bases:
        chanel_index = bases.index(base)
        trace_y = [1.0 * ci / trmax for ci in traces_y[chanel_index]]
        if show_rc:
            base = reverse_complement(base)
        ax.plot(trace_x, trace_y, color=_colors[base], lw=2, label=base)
        ax.fill_between(trace_x, 0, trace_y, facecolor=_colors[base], alpha=0.125)

    # Plot bases at peak positions
    if show_bases:
        for i, peak in enumerate(peaks):
            b = reverse_complement(sequence_str[i]) if show_rc else sequence_str[i]
            ax.text(
                peak,
                -0.11,
                b,
                color=_colors[b],
                va="center",
                ha="center",
                alpha=0.66,
                fontsize="x-large",
                fontweight="bold",
            )
        ax.set_ylim(bottom=-0.15, top=1.05)
    else:
        ax.set_ylim(bottom=-0.05, top=1.05)

    #  peaks[0] - max(2, 0.02 * (peaks[-1] - peaks[0])),
    #  right=peaks[-1] + max(2, 0.02 * (peaks[-1] - peaks[0])),
    ax.set_xlim(xlim_left + 0.5, xlim_right)

    if show_positions:
        ax.set_xticks(peaks)
        ax.set_xticklabels(list(range(region_start + 1, region_end + 2)))
    else:
        ax.set_xticks([])

    if show_rc:
        ax.invert_xaxis()

    # hide y axis
    ax.set_yticklabels([])
    ax.get_yaxis().set_visible(False)
    # hide border
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # hide grid
    ax.grid(False)
    # set legend
    ax.legend(loc="upper left", bbox_to_anchor=(0.95, 0.99))
    return ax


def show_reference(
    query_record: SeqRecord,
    subject_record: SeqRecord,
    ax: Axes,
    ref_central: Optional[int] = None,
) -> Axes:
    """
    show the reference of the chromatograph.

    design: if location is not proviode, do the alignment first
    @param seq: input SeqRecord of ref
    """

    sitepairs = align_chromatograph(query_record, subject_record)
    sitepairs_indexing = {s.cf_pos: s for s in sitepairs}
    cf_sites = [int(i.get_text()) for i in ax.get_xticklabels()]
    matched_sitepairs = [sitepairs_indexing[pos] for pos in cf_sites]
    for i, peak in enumerate(ax.get_xticks()):
        ax.text(
            peak,
            1.05,
            matched_sitepairs[i].ref_base,
            color="dimgrey",
            va="bottom",
            ha="center",
            alpha=0.85,
            fontsize="xx-large",
            fontweight="bold",
            clip_on=False,
        )
        if ref_central is not None:
            ref_pos = matched_sitepairs[i].ref_pos - ref_central
        else:
            ref_pos = matched_sitepairs[i].ref_pos
        ax.text(
            peak,
            1.12,
            ref_pos,
            color="dimgrey",
            va="bottom",
            ha="center",
            alpha=0.85,
            fontsize="medium",
            fontweight="normal",
            clip_on=False,
        )
    return ax


def center_region(
    pos: int, half_window: int, seq_length: int
) -> Tuple[int, int]:
    """Return a 1-based ``(start, end)`` region of width ``2*half_window``
    centred on ``pos``, clamped to the sequence bounds.

    This addresses the "selected base is not in the middle" issue by letting
    callers build a plot window that always keeps base ``pos`` centred when
    possible.
    """
    start = pos - half_window
    end = pos + half_window
    # shift window back into bounds rather than cropping it unevenly
    if start < 1:
        end += 1 - start
        start = 1
    if end > seq_length:
        start -= end - seq_length
        end = seq_length
    start = max(1, start)
    return int(start), int(end)


def highlight_base(
    pos_highlight: int, seq: SeqRecord, ax: Axes, passed_filter=True
) -> Axes:
    """
    Highlight the area around a peak with a rectangle.
    """
    require_matplotlib()
    peaks = seq.annotations["peak positions"]
    peak = peaks[pos_highlight - 1]

    xmin, xmax = ax.get_xlim()
    if not xmin <= peak < xmax:
        raise ValueError("peak not within plot bounds")

    if pos_highlight == 1:
        xmin = xmin
    else:
        xmin = 0.5 * (peaks[pos_highlight - 1] + peaks[pos_highlight - 2])

    if pos_highlight == len(peaks):
        xmax = xmax
    else:
        xmax = 0.5 * (peaks[pos_highlight - 1] + peaks[pos_highlight])
    ymin, ymax = ax.get_ylim()

    if passed_filter:
        fcolor = "yellow"
    else:
        fcolor = "grey"
    rec = mpl.patches.Rectangle(
        (xmin, ymin),
        (xmax - xmin),
        (ymax - ymin),
        edgecolor="none",
        facecolor=fcolor,
        alpha=0.3,
    )
    ax.add_patch(rec)
    return ax


def annotate_mutation(mut: SitePair, seq: SeqRecord, ax) -> Axes:
    """
    Annotate mutation pattern chromatograph position.
    """
    peaks = seq.annotations["peak positions"]
    peak = peaks[mut.cf_pos - 1]
    ax.text(
        peak,
        0.99,
        f"{mut.ref_base}{mut.ref_pos}{mut.cf_base}",
        color="c",
        fontsize="large",
        fontweight="bold",
        rotation=45,
        ha="center",
        va="center",
    )
    return ax
