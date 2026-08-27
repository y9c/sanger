#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot cfutils chromatograms together with DNA Features Viewer.

`DNA Features Viewer <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer>`_
renders rich feature maps (genes, primers, mutation sites) with automatic label
handling.  This module bridges the two libraries so you can show a cfutils
Sanger chromatogram and a DNA-features-viewer feature track in the same figure,
sharing a common x axis.

``dna_features_viewer`` is an **optional** dependency; importing this module does
not require it, but calling any function here does.

Typical usage::

    from cfutils.dnalink import plot_combined, to_graphic_record
    fig, (ax_feat, ax_chrom) = plot_combined(query_record, region=(55, 90))
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional, Tuple, Union

from ._mpl import plt, Axes, Figure, require_matplotlib
from .features import ChromatogramFeature, iter_features
from .parser import SeqRecord
from .show import plot_chromatograph

if TYPE_CHECKING:
    from .features import ChromatogramFeature

__all__ = [
    "to_graphic_features",
    "to_graphic_record",
    "plot_combined",
]


def _dfv():
    """Lazily import dna_features_viewer (optional dependency)."""
    try:
        import dna_features_viewer as dfv
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "dna_features_viewer is required for cfutils.dnalink; "
            "install it with `pip install dna_features_viewer`."
        ) from exc
    return dfv


def to_graphic_features(
    features: List[ChromatogramFeature], length: int,
) -> List:
    """Translate cfutils :class:`ChromatogramFeature` to dfv ``GraphicFeature``.

    Coordinates are converted from 1-based to dfv's 0-based convention.
    """
    _ = _dfv()
    dfv = _dfv()
    from dna_features_viewer import GraphicFeature

    out = []
    for f in features:
        start = max(0, f.start - 1)
        end = min(length, f.end)
        if end <= start:
            end = start + 1  # dfv needs a non-empty span
        out.append(GraphicFeature(
            start=start, end=end, strand=f.strand, color=f.color,
            label=f.label,
        ))
    return out


def to_graphic_record(
    record: SeqRecord, features: Optional[List[ChromatogramFeature]] = None,
):
    """Build a dfv ``GraphicRecord`` from a cfutils record + its features.

    When ``features`` is omitted the record's attached features are used.
    """
    _ = _dfv()
    from dna_features_viewer import GraphicRecord

    feats = list(features) if features is not None else list(iter_features(record))
    gfs = to_graphic_features(feats, len(record))
    # carry the called sequence through so dfv can draw it if requested
    return GraphicRecord(
        sequence_length=len(record),
        features=gfs,
        sequence=record.seq,
    )


def plot_combined(
    record: SeqRecord,
    features: Optional[List[ChromatogramFeature]] = None,
    region: Optional[Tuple[int, int]] = None,
    ax: Optional[Axes] = None,
    figure_width: float = 16,
    plot_sequence: bool = False,
    with_ruler: bool = True,
) -> Tuple[Figure, Tuple[Axes, Axes]]:
    """Plot a DNA-features-viewer feature map above a cfutils chromatogram.

    The top panel is the dfv feature track (base coordinates), the bottom panel
    is the chromatogram.  When plotting a region, both panels are cropped to the
    same 1-based ``[start, end]`` so the feature and traces line up.

    Args:
        record: parsed chromatogram :class:`SeqRecord`.
        features: dfv-style features; defaults to the record's attached ones.
        region: optional 1-based ``(start, end)`` region.
        ax: optional existing axis to draw the feature track into (for advanced
            compositing); when given a single ``ax`` is returned.
        figure_width: figure width in inches.
        plot_sequence: also draw the nucleotide sequence under the feature map.
        with_ruler: draw the coordinate ruler.

    Returns:
        ``(fig, (ax_feat, ax_chrom))`` or, when ``ax`` is supplied, ``ax``.
    """
    dfv = _dfv()
    from dna_features_viewer import GraphicRecord

    require_matplotlib()
    seq_len = len(record)
    if region is not None:
        start1, end1 = max(1, region[0]), min(seq_len, region[1])
    else:
        start1, end1 = 1, seq_len

    feats = list(features) if features is not None else list(iter_features(record))
    gfs = to_graphic_features(feats, seq_len)
    gfr = GraphicRecord(sequence_length=seq_len, features=gfs, sequence=record.seq)
    # crop the record to the selected region (0-based dfv coords)
    view = gfr.crop((start1 - 1, end1))

    if ax is not None:
        view.plot(ax=ax, with_ruler=with_ruler)
        if plot_sequence:
            view.plot_sequence(ax)
        return ax

    fig, (ax_feat, ax_chrom) = plt.subplots(
        2, 1, figsize=(figure_width, 6),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 3]},
    )
    view.plot(ax=ax_feat, with_ruler=with_ruler)
    if plot_sequence:
        view.plot_sequence(ax_feat)

    # dfv maps 0-based base positions; put the chromatogram on the same axis
    plot_chromatograph(record, region=(start1, end1), ax=ax_chrom)
    ax_feat.set_xlim(*ax_chrom.get_xlim())
    ax_feat.tick_params(labelbottom=False)
    return fig, (ax_feat, ax_chrom)
