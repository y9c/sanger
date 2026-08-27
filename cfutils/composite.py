#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Composite / side-by-side plotting.

The core cfutils workflow plots a single chromatogram.  A recurring request is
to show *another* tool's output (GC content, coverage, a second read, a custom
signal computed by the user's own analysis) right next to the chromatogram,
sharing the same x-axis so features line up.

This module exposes a small, dependency-light API:

* :func:`side_by_side` -- chromatogram (with optional feature overlay) plus an
  arbitrary user-provided panel function that draws on a shared-x axis.

The panel function receives ``(ax, trace_x, peaks, seq, query_record)`` and is
free to plot anything (e.g. a ``fill_between`` of a per-sample signal, a bar
chart at peak positions, or another tool's trace) -- because it draws on the
same axes objects, positions are automatically aligned with the chromatogram.
"""

from __future__ import annotations

from typing import Callable, List, Optional, Tuple

from ._mpl import plt, Axes, Figure, require_matplotlib
from .features import plot_features
from .parser import SeqRecord
from .show import plot_chromatograph

__all__ = ["side_by_side", "add_panel"]


def _trace_data(record: SeqRecord, region: Optional[Tuple[int, int]] = None):
    """Resolve trace arrays for a (optional) region."""
    peaks = record.annotations["peak positions"]
    trace_x = record.annotations.get("trace_x",
                                     list(range(len(peaks))))
    n = len(peaks)

    if region is None:
        start1, end1 = 1, n
    else:
        start1 = max(1, min(region[0], n))
        end1 = max(start1, min(region[1], n))

    lo, hi = start1 - 1, end1  # py index slice
    sel_peaks = peaks[lo:hi]
    if not sel_peaks:
        raise ValueError("empty region for panel")

    x0, x1 = float(sel_peaks[0]) - 2, float(sel_peaks[-1]) + 2
    mask = [x0 <= x <= x1 for x in trace_x]
    seg_x = [x for x, m in zip(trace_x, mask) if m]
    return seg_x, sel_peaks, record.seq[lo:hi], start1, end1


def side_by_side(
    query: SeqRecord,
    panel_func: Callable[[Axes, List[float], List[float], str, SeqRecord, int], None],
    region: Optional[Tuple[int, int]] = None,
    features: Optional[list] = None,
    show_features: bool = True,
    panel_share_x: bool = True,
    height_ratios: Tuple[int, int] = (4, 1),
    figure_width: float = 16,
) -> Tuple[Figure, Tuple[Axes, Axes]]:
    """Plot a chromatogram and a user panel side by side on a shared axis.

    Args:
        query: parsed chromatogram :class:`SeqRecord`.
        panel_func: ``f(ax, trace_x, peaks, seq, record, start_pos)`` drawing
            the second panel.  ``trace_x``/``peaks`` are the (region-resolved)
            trace x coordinates, ``seq`` the region's bases, ``record`` the
            unchanged query record, and ``start_pos`` the 1-based position of
            the first base in the region (so per-base data can be placed at the
            correct trace coordinate).
        region: optional 1-based ``(start, end)`` region to restrict plotting.
        features: optional explicit feature list passed to :func:`plot_features`.
        show_features: whether to overlay attached/listed features on the
            chromatogram panel.
        panel_share_x: share the x-axis between panels so positions align.
        height_ratios: relative heights of the chromatogram vs the panel.
        figure_width: width of the figure in inches.

    Returns:
        ``(fig, (ax_chrom, ax_panel))``.
    """
    require_matplotlib()
    fig, (ax_chrom, ax_panel) = plt.subplots(
        2, 1,
        figsize=(figure_width, 6),
        sharex=panel_share_x,
        gridspec_kw={"height_ratios": list(height_ratios)},
    )

    plot_chromatograph(query, region=region, ax=ax_chrom)
    if show_features:
        plot_features(query, ax_chrom, features=features)

    seg_x, sel_peaks, sel_seq, start_pos, _end = _trace_data(query, region)
    panel_func(ax_panel, seg_x, sel_peaks, sel_seq, query, start_pos)
    if panel_share_x:
        _link_xlims(ax_chrom, ax_panel)
        ax_chrom.tick_params(labelbottom=False)
    return fig, (ax_chrom, ax_panel)


def add_panel(fig: Figure, ax_chrom: Axes,
              panel_func: Callable[[Axes, List[float], List[float], str, SeqRecord, int], None],
              query: SeqRecord, region: Optional[Tuple[int, int]] = None) -> Axes:
    """Add a second panel beneath an existing chromatogram axes (late binding).

    Convenience when you already have a chromatogram plotted and just want an
    aligned extra panel.
    """
    require_matplotlib()
    from matplotlib.gridspec import GridSpec

    gs = ax_chrom.get_subplotspec().get_gridspec()
    nrows = gs.get_geometry()[0] + 1
    # create a new grid below
    gs2 = GridSpec(nrows, 1, height_ratios=[3] * (nrows - 1) + [1])
    ax_panel = fig.add_subplot(gs2[nrows - 1, 0], sharex=ax_chrom)
    seg_x, sel_peaks, sel_seq, start_pos, _e = _trace_data(query, region)
    panel_func(ax_panel, seg_x, sel_peaks, sel_seq, query, start_pos)
    ax_chrom.tick_params(labelbottom=False)
    return ax_panel


def _link_xlims(ax_top: Axes, ax_bottom: Axes) -> None:
    """Best-effort x-limit alignment so peaks line up between panels."""
    xl = ax_top.get_xlim()
    ax_bottom.set_xlim(xl)
    if ax_bottom.get_xticks() is None or len(ax_bottom.get_xticklabels()) == 0:
        pass
