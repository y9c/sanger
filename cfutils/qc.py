#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Per-read quality-control (QC) metrics for Sanger chromatograms.

Provides the numbers a bench scientist / pipeline needs to judge read quality:
read length, base composition, N fraction, mean / worst Phred quality,
low-quality fraction, Mott-trimmed length, and peak resolution.  Also a batch
:func:`summarize` that produces a table (list of dicts) for many reads.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, Iterable, List, Tuple

import numpy as np

if TYPE_CHECKING:
    from .parser import SeqRecord

__all__ = [
    "read_metrics",
    "trimmed_bounds",
    "summarize",
    "LOW_QUAL",
]

#: conventional low-quality Phred threshold
LOW_QUAL = 20


def _mott_bounds(record: "SeqRecord", cutoff: float = 0.05,
                 segment: int = 20) -> Tuple[int, int]:
    """1-based start/end of the high-quality segment (Mott trimming).

    Returns ``(1, len)`` when the read is too short to trim.
    """
    if len(record) <= segment:
        return 1, len(record)
    qual = record.letter_annotations["phred_quality"]
    score = [cutoff - (10 ** (q / -10.0)) for q in qual]
    cummul = [0]
    start = 0
    trimmed_start = 0
    for i in range(1, len(score)):
        v = cummul[-1] + score[i]
        if v < 0:
            cummul.append(0)
        else:
            cummul.append(v)
            if not start:
                trimmed_start = i
                start = True
    if not cummul:
        return 1, len(record)
    trimmed_end = cummul.index(max(cummul))
    return max(1, trimmed_start + 1), max(trimmed_end, 1)


def trimmed_bounds(record: "SeqRecord", cutoff: float = 0.05,
                   segment: int = 20) -> Tuple[int, int]:
    """Public helper returning 1-based high-quality segment bounds."""
    return _mott_bounds(record, cutoff=cutoff, segment=segment)


def read_metrics(record: "SeqRecord") -> Dict[str, float]:
    """Compute QC metrics for a single chromatogram record.

    Returns a dict with keys::

        n_bases, gc_percent, n_n, n_fraction, mean_qual, min_qual,
        low_qual_fraction, trim_start, trim_end, trimmed_len,
        n_peaks, peak_resolution
    """
    seq = record.seq
    n = len(seq)
    qual = np.asarray(record.letter_annotations.get("phred_quality", []),
                      dtype=float)
    peaks = np.asarray(record.annotations.get("peak positions", []), dtype=float)

    n_n = sum(a == "N" for a in seq.upper())
    trim_start, trim_end = _mott_bounds(record)

    metrics: Dict[str, float] = {
        "n_bases": float(n),
        "gc_percent": 100.0 * (seq.count("G") + seq.count("C")) / n if n else 0.0,
        "n_n": float(n_n),
        "n_fraction": n_n / n if n else 0.0,
        "mean_qual": float(qual.mean()) if qual.size else 0.0,
        "min_qual": float(qual.min()) if qual.size else 0.0,
        "low_qual_fraction": float((qual < LOW_QUAL).mean()) if qual.size else 0.0,
        "trim_start": float(trim_start),
        "trim_end": float(trim_end),
        "trimmed_len": float(trim_end - trim_start + 1),
        "n_peaks": float(len(peaks)),
    }
    if len(peaks) > 2:
        diffs = np.diff(peaks)
        metrics["peak_resolution"] = float(np.mean(diffs))
    else:
        metrics["peak_resolution"] = 0.0
    return metrics


def summarize(records: Iterable["SeqRecord"]) -> List[Dict[str, float]]:
    """Return a list of metric dicts (one per record) for batch QC."""
    return [read_metrics(rec) for rec in records]
