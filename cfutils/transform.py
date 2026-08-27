#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Transformations on whole chromatogram records.

These operate on the record as a unit -- not just the string -- so the trace
arrays, peak axis and quality stay consistent with the sequence.  This
includes:

* :func:`trim` -- Mott quality trimming that keeps traces/peaks aligned (fixing
  the long-standing "chromatogram switches position after trim" issue).
* :func:`reverse_complement_record` -- reverse/complement a whole record,
  including its traces and peaks (Snapgene-style reversed view).

Both reuse the axis-preserving slicing in :mod:`cfutils.tracks` so the
sequence, quality, peaks and channels never drift apart.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Tuple

import numpy as np

from .qc import trimmed_bounds
from .tracks import slice_track

if TYPE_CHECKING:
    from .parser import SeqRecord

__all__ = ["trim", "reverse_complement_record", "NCHANNELS"]


def trim(record: "SeqRecord", cutoff: float = 0.05, segment: int = 20,
         name: Optional[str] = None) -> "SeqRecord":
    """Return a quality-trimmed copy of the record.

    Richard Mott's trimming algorithm finds the contiguous segment of the read
    with the highest cumulative base quality.  Unlike the legacy
    ``_abi_trim`` (which only sliced the sequence and left the trace/peak axes
    untouched), this reuses axis-preserving slicing so peak positions, channel
    traces and ``trace_x`` all stay aligned with the trimmed sequence.

    Args:
        record: the chromatogram to trim.
        cutoff: base-score cutoff for the Mott algorithm (0.05 default).
        segment: minimum sequence length before trimming is attempted.
        name: optional name for the trimmed record.
    """
    start, end = trimmed_bounds(record, cutoff=cutoff, segment=segment)
    if start <= 1 and end >= len(record):
        return record
    trimmed = slice_track(record, start, end,
                          name=name or f"{record.name or 'trimmed'}_trimmed")
    return trimmed


def reverse_complement_record(record: "SeqRecord",
                              name: Optional[str] = None) -> "SeqRecord":
    """Reverse-complement a whole chromatogram record.

    The sequence and quality are reversed/complemented, and the trace + peak
    axis are mirrored so a reverse-strand view stays correctly aligned.

    Args:
        record: chromatogram to reverse-complement.
        name: optional name for the output record.
    """
    comp = str.maketrans("ACGTN", "TGCAN")
    seq = record.seq.translate(comp)[::-1]
    # quality reverses with the sequence (same per-base ordering)
    qual = list(record.letter_annotations.get("phred_quality", []))[::-1]

    channels = np.array(
        [record.annotations["channel " + str(i)] for i in range(1, 5)],
        dtype=float,
    )
    peaks = np.asarray(record.annotations.get("peak positions", []), dtype=float)
    trace_x = np.asarray(record.annotations.get("trace_x", []), dtype=float)

    new = record.__class__(
        seq,
        id=record.id,
        name=name or f"{record.name or 'rc'}_rc",
        description=record.description,
        annotations=dict(record.annotations),
        letter_annotations={"phred_quality": qual},
    )
    # mirror the traces and peaks about the trace axis so they point left
    if trace_x.size and channels.size:
        xmax = float(trace_x.max())
        new.annotations["trace_x"] = (xmax - trace_x).tolist()
        # reverse the sample order within each channel
        for i in range(1, 5):
            new.annotations["channel " + str(i)] = channels[i - 1][::-1].tolist()
        new.annotations["peak positions"] = (xmax - peaks).tolist()
    return new


NCHANNELS = 4
