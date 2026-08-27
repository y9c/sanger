#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Chromatogram feature annotation & overlay API.

Inspired by the DNA Features Viewer object model (data decoupled from
rendering), this module lets any tool attach annotation features (mutations,
primers, genes, quality dropout) to a parsed Sanger :class:`SeqRecord` and
plot them side by side with the chromatogram.

Typical usage::

    from cfutils.features import ChromatogramFeature, plot_features

    feat = ChromatogramFeature(start=42, end=60, strand=+1,
                               color="#ff8888", label="primer F")
    plot_features(query_record, ax, features=[feat])
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Iterable, List, Optional

from ._mpl import Axes, mpl, require_matplotlib

if TYPE_CHECKING:
    from .parser import SeqRecord

__all__ = [
    "ChromatogramFeature",
    "add_feature",
    "iter_features",
    "plot_features",
    "peak_to_x",
]


@dataclass
class ChromatogramFeature:
    """A single annotation feature over a chromatogram sequence.

    ``start``/``end`` are 1-based coordinates along the (called) sequence.
    ``strand`` controls the arrow direction (``+1`` right, ``-1`` left).
    """

    start: int
    end: int
    strand: int = +1
    color: str = "#cccccc"
    label: Optional[str] = None
    kind: str = "feature"
    #: optional extra payload carried through (e.g. a mutation SitePair)
    data: object = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if self.end < self.start:
            raise ValueError("feature end must be >= start")
        if self.strand not in (-1, 0, 1):
            raise ValueError("strand must be -1, 0 or +1")

    @classmethod
    def from_sitepair(
        cls, site, color: str = "#e0573b", label_fmt: str = "{ref}{ref_pos}{cf}"
    ):
        """Build a single-base feature point around an aligned mutation site."""
        return cls(
            start=site.cf_pos,
            end=site.cf_pos,
            strand=0,
            color=color,
            label=label_fmt.format(
                ref=site.ref_base, ref_pos=site.ref_pos, cf=site.cf_base
            ),
            kind="mutation",
            data=site,
        )


#: annotations key used to store features attached to a SeqRecord
FEATURES_ANNOTATION_KEY = "features"


def add_feature(record: "SeqRecord", feature: ChromatogramFeature) -> None:
    """Attach a feature to a parsed record (mutates its annotations)."""
    feats = record.annotations.setdefault(FEATURES_ANNOTATION_KEY, [])
    feats.append(feature)


def iter_features(record: "SeqRecord") -> Iterable[ChromatogramFeature]:
    """Yield all features attached to a record (empty if none)."""
    return iter(record.annotations.get(FEATURES_ANNOTATION_KEY, []))


def peak_to_x(record: "SeqRecord", positions: Iterable[int]) -> List[float]:
    """Map 1-based sequence positions to x-axis (trace) coordinates.

    Uses the recorded peak positions; peaks outside the trace range are
    clipped to the nearest valid peak x.
    """
    peaks = record.annotations["peak positions"]
    if not peaks:
        return []
    out = []
    for pos in positions:
        idx = max(0, min(pos - 1, len(peaks) - 1))
        out.append(peaks[idx])
    return out


def _plot_arrow(ax, x0, x1, y, color, lw, alpha):
    """Small axis-aligned arrow used for strand-aware features."""
    from matplotlib.patches import FancyArrow

    if x1 == x0:
        return  # degenerate single-point feature; caller draws a marker
    head = 0.6 * min(abs(x1 - x0), 3.0)
    return FancyArrow(
        x0,
        y,
        x1 - x0,
        0.0,
        width=lw,
        head_width=0.30,
        head_length=head,
        length_includes_head=True,
        facecolor=color,
        edgecolor="none",
        alpha=alpha,
    )


def plot_features(
    record: "SeqRecord",
    ax: Axes,
    features: Optional[Iterable[ChromatogramFeature]] = None,
    y_bottom: float = 1.15,
    band_height: float = 0.30,
    show_legend: bool = True,
    alpha: float = 0.85,
) -> Axes:
    """Render feature arrows/boxes above the chromatogram on ``ax``.

    Coordinates are mapped through the recorded peak positions, so features
    stay aligned with the trace even after trimming or rescaling.

    Args:
        record: parsed chromatogram :class:`SeqRecord`.
        ax: target matplotlib axes (usually the chromatogram axes).
        features: optional explicit list; defaults to attached features.
        y_bottom: vertical offset (in data units) of the feature band.
        band_height: height of the feature band.
        show_legend: whether to draw a legend for labelled features.
        alpha: transparency of the coloured features.
    """
    require_matplotlib()
    if features is None:
        features = list(iter_features(record))
    features = list(features)
    if not features:
        return ax

    band_y = y_bottom + band_height / 2.0
    lw = band_height * 0.55
    drawn_labels = []

    for feat in features:
        x0, x1 = peak_to_x(record, (feat.start, feat.end))
        y_correction = band_height * 0.16

        if feat.strand != 0 and x1 != x0:
            arrow = _plot_arrow(ax, x0, x1, band_y, feat.color, lw * 0.5, alpha)
            if arrow is not None:
                ax.add_patch(arrow)
            direction = -1 if feat.strand < 0 else 0
            x0, x1 = x0 + direction * 0.4, x1 + direction * 0.4
            band_y = y_bottom + band_height / 2.0 + y_correction
        else:
            # neutral / point feature -> a vertical marker
            height = band_height * 0.8
            ax.plot(
                [x0, x0],
                [y_bottom + 0.05, y_bottom + height],
                color=feat.color,
                lw=2.2,
                alpha=alpha,
                solid_capstyle="round",
            )
            ax.plot(x0, band_y, marker="|", ms=14, color=feat.color, alpha=alpha)

        if feat.label:
            ax.text(
                x0 + 0.15 * (x1 - x0) if feat.strand else x0,
                band_y + 0.02,
                feat.label,
                color=feat.color,
                fontsize="small",
                fontweight="bold",
                va="bottom",
                ha="left" if (feat.strand >= 0 or x1 == x0) else "right",
                clip_on=False,
            )
            drawn_labels.append((feat.label, feat.color))

    # keep the reference band clear
    ax.set_ylim(top=max(ax.get_ylim()[1], y_bottom + band_height + 0.05))

    if show_legend and drawn_labels:
        handles = [
            mpl.lines.Line2D([], [], color=c, lw=3, label=label, alpha=alpha)
            for label, c in drawn_labels
        ]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.95, 0.99))
    return ax


def annotate_features(record: "SeqRecord", ax: Axes, **kwargs) -> Axes:
    """Alias for :func:`plot_features` using attached features only."""
    return plot_features(record, ax, **kwargs)
