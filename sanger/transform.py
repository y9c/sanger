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

Both reuse the axis-preserving slicing in :mod:`sanger.tracks` so the
sequence, quality, peaks and channels never drift apart.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

from .qc import trimmed_bounds
from .tracks import slice_track

if TYPE_CHECKING:
    from .parser import SeqRecord

__all__ = [
    "trim",
    "trim_ends",
    "trim_leading_ns",
    "strip_primers",
    "reverse_complement_record",
    "NCHANNELS",
]


def trim(
    record: "SeqRecord",
    cutoff: float = 0.05,
    segment: int = 20,
    name: Optional[str] = None,
) -> "SeqRecord":
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
    trimmed = slice_track(
        record, start, end, name=name or f"{record.name or 'trimmed'}_trimmed"
    )
    return trimmed


def trim_ends(
    record: "SeqRecord", min_qual: int = 20, name: Optional[str] = None
) -> "SeqRecord":
    """Hard-trim the low-quality 5' and 3' ends of a read.

    The trimmed region is the longest interior span whose bases are all at or
    above ``min_qual`` (a common, simple lab criterion).  Peak/trace axes stay
    aligned via the same slicing used by :func:`trim`.

    Args:
        record: the chromatogram to trim.
        min_qual: minimum per-base quality to keep.
        name: optional output name.
    """
    qual = record.letter_annotations.get("phred_quality", [])
    if not qual or len(qual) != len(record):
        return record
    # walk in from both ends while quality < min_qual
    start = 0
    while start < len(qual) and qual[start] < min_qual:
        start += 1
    end = len(qual) - 1
    while end >= start and qual[end] < min_qual:
        end -= 1
    if end < start:
        return slice_track(record, 1, 1, name=name or f"{record.name or 'tr'}_tr")
    return slice_track(
        record, start + 1, end + 1, name=name or f"{record.name or 'tr'}_trimmed"
    )


def trim_leading_ns(record: "SeqRecord", name: Optional[str] = None) -> "SeqRecord":
    """Remove the unreliable leading run of ``N``/low-quality bases.

    The first 20-40 bases of a Sanger read are typically poorly resolved
    (short products migrate erratically) and are called as ``N``.  This trims
    the leading ``N`` run (and trailing one), keeping peak/trace axes aligned.
    """
    seq = record.seq.upper()
    n = len(seq)
    start = 0
    while start < n and seq[start] == "N":
        start += 1
    end = n
    while end > start and seq[end - 1] == "N":
        end -= 1
    if start == 0 and end == n:
        return record
    return slice_track(
        record, start + 1, end, name=name or f"{record.name or 'ns'}_noNs"
    )


def strip_primers(
    record: "SeqRecord",
    forward: str = "",
    reverse: Optional[str] = None,
    name: Optional[str] = None,
) -> "SeqRecord":
    """Remove primer sequences from the 5' (and optionally 3') read ends.

    The primer is matched at the read start (``forward``) and, if given, the
    ``reverse`` primer at the 3' end.  Exact matching only; if a primer is not
    found the corresponding end is left untouched.  Traces/peaks stay aligned.

    Args:
        record: the chromatogram to strip.
        forward: 5' primer sequence (matched and removed if present).
        reverse: 3' primer sequence (matched and removed if present).
        name: optional output name.
    """
    seq = record.seq.upper()
    start = 1
    if forward:
        f = forward.upper()
        if seq.startswith(f):
            start = len(f) + 1
        else:
            # tolerate the primer being reverse-complemented on the read
            rc = reverse_complement_record(record).seq.upper()
            if rc.startswith(f):
                pass  # user should reverse-complement the record first
    end = len(record)
    if reverse:
        r = reverse.upper()
        if seq.endswith(r):
            end = len(record) - len(r)
    if start == 1 and end == len(record):
        return record
    return slice_track(
        record, start, end, name=name or f"{record.name or 's'}_primers_removed"
    )


def reverse_complement_record(
    record: "SeqRecord", name: Optional[str] = None
) -> "SeqRecord":
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
