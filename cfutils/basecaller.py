#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Base-calling from raw four-channel Sanger traces.

ABI files store both the machine-called sequence (``PBAS2``) and the raw
analysed four colour-channel traces (``DATA9``-``DATA12``).  This module
implements a dependency-light (numpy-only) peak-driven base caller that can
re-derive a called sequence and quality directly from the traces.

This is *not* a substitute for the vendor's base-caller (which uses the
analysis records and machine models), but it is a fast, transparent reference
implementation that is useful for cross-checking, re-processing legacy traces,
and teaching how peak heights map to base calls.

Algorithm
---------
1. Baseline-correct each channel (subtract a running median).
2. Detect per-channel peaks (local maxima above a threshold, with a minimum
   spacing).
3. Cluster peaks across channels within a small x-window -> one base position.
4. Call the base of the strongest channel; derive a Phred-like quality from
   the ratio of the strongest peak to the second-strongest channel.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from .parser import SeqRecord

__all__ = [
    "BaseCall",
    "BaseCallResult",
    "detect_peaks",
    "call_bases",
    "basecaller_score",
]

_CHANNEL_BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


def _channel_traces(record: SeqRecord) -> Tuple[np.ndarray, List[str]]:
    """Return (n_channels x n_samples, channel_base_order) arrays."""
    traces = np.array(
        [record.annotations["channel " + str(i)] for i in range(1, 5)],
        dtype=float,
    )
    channel_order = list(record.annotations.get("channels", "ACGT"))
    return traces, channel_order


def _baseline_correct(trace: np.ndarray, window: int = 51) -> np.ndarray:
    """Subtract a running median baseline to remove dye bleed / drift."""
    if trace.size == 0:
        return trace
    import numpy as np

    n = trace.size
    half = window // 2
    baseline = np.empty_like(trace)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        baseline[i] = np.median(trace[lo:hi])
    return trace - baseline


def _local_maxima(signal: np.ndarray, min_distance: int = 4,
                  threshold: float = 0.0) -> np.ndarray:
    """Indices of local maxima separated by >= min_distance samples."""
    if signal.size == 0:
        return np.array([], dtype=int)
    # candidate local maxima: greater than both immediate neighbours
    mask = (signal[1:-1] > signal[:-2]) & (signal[1:-1] > signal[2:])
    candidates = np.flatnonzero(mask) + 1
    if candidates.size == 0:
        return np.array([], dtype=int)

    # greedy suppression of candidates too close together, keep the taller
    order = candidates[np.argsort(signal[candidates])[::-1]]
    kept = []
    for idx in order:
        if all(abs(int(idx) - k) >= min_distance for k in kept):
            kept.append(int(idx))
    kept = np.array(sorted(kept), dtype=int)
    return kept[signal[kept] > threshold]


@dataclass
class BaseCall:
    """One called base with supporting evidence."""

    position: int          # 1-based position in the called sequence
    base: str
    trace_x: float         # x-sample coordinate of the winning peak
    confidence: float      # 0..1, margin of the winning channel
    quality: int           # approximate Phred-like quality
    second_base: Optional[str] = None
    second_ratio: Optional[float] = None


@dataclass
class BaseCallResult:
    """Result of calling bases on one chromatogram."""

    calls: List[BaseCall]

    @property
    def sequence(self) -> str:
        return "".join(c.base for c in self.calls)

    @property
    def qualities(self) -> List[int]:
        return [c.quality for c in self.calls]

    @property
    def n_calls(self) -> int:
        return len(self.calls)

    def accuracy(self, reference_seq: str) -> float:
        """Fraction of calls matching a reference (same-length) sequence."""
        ref = str(reference_seq).upper()
        mine = self.sequence.upper()
        if not mine:
            return 0.0
        m = min(len(ref), len(mine))
        if m == 0:
            return 0.0
        return sum(a == b for a, b in zip(ref[:m], mine[:m])) / m


def _quality_from_margin(margin: float) -> int:
    """Map a confidence margin to a Phred-like quality score."""
    if margin <= 0:
        return 0
    # margin == 1 -> perfectly dominant; 0.5 -> coin flip
    if margin >= 0.85:
        return 60
    if margin >= 0.7:
        return 40
    if margin >= 0.55:
        return 20
    return 5 + int(30 * margin)


def detect_peaks(record: SeqRecord, min_distance: int = 5,
                 threshold_quantile: float = 0.2) -> List[Tuple[float, int]]:
    """Detect trace peaks per channel.

    Returns a flat list of ``(x_sample, channel_index)``, grouped loosely in
    x order. ``channel_index`` refers to the position in
    ``record.annotations["channels"]``.
    """
    traces, order = _channel_traces(record)
    peaks: List[Tuple[float, int]] = []
    for ch, trace in enumerate(traces):
        corrected = _baseline_correct(trace)
        thresh = np.quantile(corrected, threshold_quantile)
        for x in _local_maxima(corrected, min_distance=min_distance,
                               threshold=thresh):
            peaks.append((float(x), ch))
    peaks.sort()
    return peaks


def call_bases(record: SeqRecord, min_distance: int = 5,
               cluster_dx: int = 6, threshold_quantile: float = 0.2,
               max_calls: Optional[int] = None,
               use_peak_positions: bool = True) -> BaseCallResult:
    """Call bases from the raw traces of a chromatogram.

    Two modes:

    * ``use_peak_positions=True`` (default, recommended): use the vendor peak
      positions recorded in the ABI (``PLOC2``) as the reference peak sites and
      re-call the base at each from the four channels.  This is fast and
      accurate, and yields a per-base dye-peak evidence quality.  Requires the
      record parsed with ``parse_abi(..., rescale=False)``.
    * ``use_peak_positions=False``: detect peaks de novo from the traces.

    Args:
        record: parsed chromatogram whose four channel traces are used.
        min_distance: de novo peak spacing (ignored in peak-position mode).
        cluster_dx: de novo cross-channel clustering radius.
        threshold_quantile: baseline threshold quantile for detection.
        max_calls: cap the number of base positions (for speed on long reads).
        use_peak_positions: use the vendor peak positions when present.

    Returns:
        A :class:`BaseCallResult`.
    """
    traces, order = _channel_traces(record)
    n_channels, n_samples = traces.shape
    corrected = np.array([_baseline_correct(t) for t in traces])

    if use_peak_positions and len(record.annotations.get("peak positions", [])):
        peak_x = np.asarray(record.annotations["peak positions"], dtype=float)
        # plus the origin offset if the trace was truncated (rescaled mode)
        x0 = float(np.asarray(record.annotations["trace_x"], dtype=float)[0]) \
            if "trace_x" in record.annotations and len(record.annotations["trace_x"]) else 0.0
        ref_positions = peak_x - x0
        positions = ref_positions[::max(1, min_distance // 5)] if False else ref_positions
        if max_calls:
            positions = positions[:max_calls]
        x_vals = positions
    else:
        # de novo: cluster detected peaks across channels
        all_peaks = []  # (x, ch, height)
        for ch, t in enumerate(corrected):
            thresh = np.quantile(t, threshold_quantile)
            for x in _local_maxima(t, min_distance=min_distance, threshold=thresh):
                all_peaks.append((float(x), ch, float(t[x])))
        all_peaks.sort(key=lambda p: p[0])
        clusters: List[List[Tuple[int, float, float]]] = []
        for _x, ch, h in all_peaks:
            if clusters and _x - clusters[-1][-1][1] <= cluster_dx:
                clusters[-1].append((ch, _x, h))
            else:
                clusters.append([(ch, _x, h)])
        if max_calls:
            clusters = clusters[:max_calls]
        x_vals = [float(np.mean([p[1] for p in c])) for c in clusters]

    calls: List[BaseCall] = []
    for i, x in enumerate(x_vals, start=1):
        idx = int(round(x))
        if idx < 0:
            idx = 0
        if idx >= n_samples:
            idx = n_samples - 1
        # sample each channel in a small window around the peak
        window = max(1, min_distance // 2)
        heights = []
        for ch in range(n_channels):
            lo = max(0, idx - window)
            hi = min(n_samples, idx + window + 1)
            seg = corrected[ch, lo:hi]
            heights.append(float(seg.max()) if seg.size else 0.0)
        heights = np.asarray(heights)
        peak_i = int(np.argmax(heights))
        base = order[peak_i] if peak_i < len(order) else "N"
        # margin = how much the winner beats the runner-up
        sorted_heights = np.sort(heights)[::-1]
        top, second = sorted_heights[0], sorted_heights[1] if sorted_heights.size > 1 else 0.0
        second_ratio = (second / top) if top > 0 else 0.0
        margin = 1.0 - second_ratio
        calls.append(BaseCall(
            position=i,
            base=base if top > 0 else "N",
            trace_x=float(x),
            confidence=margin,
            quality=_quality_from_margin(margin),
            second_base=None,
            second_ratio=second_ratio,
        ))

    return BaseCallResult(calls=calls)


def basecaller_score(result: BaseCallResult, reference_seq: str) -> dict:
    """Return a small summary dict comparing calls against a reference."""
    acc = result.accuracy(reference_seq)
    return {
        "n_calls": result.n_calls,
        "accuracy": acc,
        "mean_quality": float(np.mean(result.qualities)) if result.qualities else 0.0,
        "n_gt_q20": sum(q >= 20 for q in result.qualities),
    }
