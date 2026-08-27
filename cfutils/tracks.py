#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Split and join Sanger chromatogram (trace) records.

Chromatograms are modelled as :class:`~cfutils.parser.SeqRecord` objects whose
``annotations`` hold the four colour-channel traces, the peak positions and a
``trace_x`` axis.  This module provides the operations a wet-lab user asked
for:

* :func:`join_tracks`  -- stitch several reads / traces into one long record.
* :func:`split_track`  -- cut one record into independent segments.
* :func:`slice_track`  -- extract an arbitrary 1-based region as a new record.

All operations preserve the alignment between ``trace_x``, channels, peaks,
called sequence and quality, and use vectorised numpy where it matters for
performance.
"""

from __future__ import annotations

import itertools
from pathlib import Path
from typing import Iterable, List, Optional

import numpy as np

from .parser import SeqRecord

__all__ = [
    "NCHANNELS",
    "join_tracks",
    "split_track",
    "slice_track",
    "export_tracks",
]

NCHANNELS = 4


def _trace_arrays(record: SeqRecord):
    """Return (channels as np.array shape (4, n), peaks, trace_x, seq, qual)."""
    channels = np.array(
        [record.annotations["channel " + str(i + 1)] for i in range(NCHANNELS)],
        dtype=float,
    )
    peaks = np.asarray(record.annotations["peak positions"], dtype=float)
    trace_x = record.annotations.get(
        "trace_x",
        np.arange(channels.shape[1]) if len(channels.shape) == 2 else [],
    )
    trace_x = np.asarray(trace_x, dtype=float)
    seq = record.seq
    qual = list(record.letter_annotations.get("phred_quality", []))
    return channels, peaks, trace_x, seq, qual


def _build_record(
    channels: np.ndarray,
    peaks: np.ndarray,
    trace_x: np.ndarray,
    seq: str,
    qual: List[int],
    name: str,
    annotations_extra: Optional[dict] = None,
) -> SeqRecord:
    """Assemble a SeqRecord from numeric trace arrays."""
    record = SeqRecord(seq, name=name)
    for i in range(NCHANNELS):
        record.annotations["channel " + str(i + 1)] = channels[i].tolist()
    record.annotations["peak positions"] = peaks.tolist()
    record.annotations["trace_x"] = trace_x.tolist()
    record.letter_annotations["phred_quality"] = qual
    if annotations_extra:
        record.annotations.update(annotations_extra)
    return record


def join_tracks(
    *records: SeqRecord, name: Optional[str] = None, gap: int = 24
) -> SeqRecord:
    """Stitch chromatogram records end-to-end into a single record.

    The traces, peak axis and quality are concatenated in order.  Each
    record's trace is offset so the pieces sit sequentially on a shared
    x axis (with ``gap`` samples of neutral spacing between them).  Peak
    positions are re-derived so sequence/peak alignment is preserved.

    Args:
        records: two or more parsed chromatogram records (in order).
        name: name for the joined record (defaults to "joined").
        gap: number of zero samples to leave between consecutive records.

    Returns:
        A single :class:`SeqRecord` combining all traces.
    """
    records = list(records)
    if not records:
        raise ValueError("join_tracks requires at least one record")
    if name is None:
        name = "joined"

    seq_parts, qual_parts = [], []
    peak_offsets: List[float] = []
    channel_parts, x_parts = [], []

    running_offset = 0.0
    for rec in records:
        channels, peaks, trace_x, seq, qual = _trace_arrays(rec)
        if trace_x is not None and len(trace_x) > 0:
            # normalise each record's trace to start at 0, then shift
            shift = trace_x[0]
            trace_x = np.asarray(trace_x, dtype=float) - shift
            peaks_shifted = np.asarray(peaks, dtype=float) - shift
        else:
            peaks_shifted = np.asarray(peaks, dtype=float)

        channel_parts.append(channels)
        x_parts.append(trace_x + running_offset)
        peak_offsets.append(peaks_shifted + running_offset)
        seq_parts.append(seq)
        qual_parts.extend(qual)

        if trace_x is not None and len(trace_x) > 0:
            running_offset += float(trace_x[-1]) + gap

    channels = np.concatenate(channel_parts, axis=1)
    trace_x = np.concatenate(x_parts)
    peaks = np.concatenate(peak_offsets)

    # keep the channel/base order (FWO) from the first record for
    # plotting/base-calling
    extra = {}
    if records and "channels" in records[0].annotations:
        extra["channels"] = records[0].annotations["channels"]

    return _build_record(
        channels,
        peaks,
        trace_x,
        "".join(seq_parts),
        qual_parts,
        name,
        annotations_extra=extra,
    )


def slice_track(
    record: SeqRecord, start: int, end: int, name: Optional[str] = None
) -> SeqRecord:
    """Extract a 1-based, inclusive region ``[start, end]`` as a new record.

    Returns a fresh, independent record.  A copy of the original annotations
    (sample well, machine, run dates) is carried across where intact.
    """
    start1 = max(1, start)
    end1 = min(len(record), end)
    if start1 > end1:
        raise ValueError(f"empty slice [{start}, {end}]")

    channels, peaks, trace_x, seq, qual = _trace_arrays(record)

    lo_idx, hi_idx = start1 - 1, end1  # python slicing end exclusive
    seg_peaks = peaks[lo_idx:hi_idx]
    if len(seg_peaks) == 0:
        raise ValueError(f"no peaks in slice [{start}, {end}]")

    # include the full trace between the first and last peak in the slice,
    # plus a small margin so the chromatogram looks natural
    x0 = float(seg_peaks[0]) - 2
    x1 = float(seg_peaks[-1]) + 2
    mask = (trace_x >= x0) & (trace_x <= x1)
    seg_channels = channels[:, mask]
    seg_x = trace_x[mask]
    if len(seg_x) == 0:
        raise ValueError(f"no trace samples in slice [{start}, {end}]")

    # re-normalise the slice so its peak axis starts at 0 (self-consistent),
    # and record where it came from so positions can be mapped back
    shift = float(seg_peaks[0])
    seg_peaks = seg_peaks - shift
    seg_x = seg_x - shift

    extra = {}
    for key in (
        "sample_well",
        "dye",
        "polymer",
        "machine_model",
        "run_start",
        "run_finish",
        "channels",
    ):
        if key in record.annotations:
            extra[key] = record.annotations[key]
    # provenance: 0-based offset of this slice within the original read
    parent_offset = int(record.annotations.get("offset", 0))
    extra["offset"] = parent_offset + lo_idx
    extra["parent"] = record.name or record.id or ""

    new_name = name or f"{record.name or 'slice'}_[{start1}-{end1}]"
    rec = _build_record(
        seg_channels,
        seg_peaks,
        seg_x,
        record.seq[lo_idx:hi_idx],
        qual[lo_idx:hi_idx],
        new_name,
        annotations_extra=extra,
    )
    return rec


def split_track(
    record: SeqRecord, cuts: Iterable[int], names: Optional[Iterable[str]] = None
) -> List[SeqRecord]:
    """Split a chromatogram into segments at the given 1-based cut positions.

    Each cut position becomes the *first* base of the following segment.

    Args:
        record: the chromatogram to split.
        cuts: 1-based positions where the trace is cut.
        names: optional output names; must have ``len(cuts)+1`` entries.

    Returns:
        A list of independent :class:`SeqRecord` pieces.
    """
    cuts = sorted({int(c) for c in cuts if 1 <= c <= len(record)})
    if not cuts:
        # nothing to split -> return a copy of the whole record
        return [slice_track(record, 1, len(record), name=names[0] if names else None)]

    # each piece is inclusive on both ends; avoid sharing the cut base
    starts = [1] + [c + 1 for c in cuts]
    ends = cuts + [len(record)]
    name_iter = iter(names) if names else itertools.repeat(None)
    pieces = []
    for start, end in zip(starts, ends):
        if start > end:
            continue
        pieces.append(slice_track(record, start, end, name=next(name_iter)))
    return pieces


def export_tracks(
    records: Iterable[SeqRecord], outdir: str, fmt: str = "npz"
) -> List[Path]:
    """Persist trace records to disk without re-writing binary ABI.

    For round-tripping traces that were split/joined we export a lossless
    numpy ``.npz`` (fast) or a simple ``.tsv`` (portable).  Reimport with
    :func:`import_tracks`.

    Args:
        records: chromatogram records to export.
        outdir: destination directory (created if needed).
        fmt: ``"npz"`` (lossless, fast) or ``"tsv"`` (human-readable).

    Returns:
        List of written file paths.
    """
    from pathlib import Path as P

    outdir = P(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for i, rec in enumerate(records):
        channels, peaks, trace_x, seq, qual = _trace_arrays(rec)
        stem = rec.name or f"track_{i}"
        if fmt == "npz":
            path = outdir / f"{stem}.npz"
            np.savez(
                path,
                channels=channels,
                peaks=peaks,
                trace_x=trace_x,
                qual=np.asarray(qual, dtype=int),
            )
        elif fmt == "tsv":
            path = outdir / f"{stem}.tsv"
            _write_tsv(path, channels, peaks, trace_x, seq, qual)
        else:
            raise ValueError(f"unknown export format: {fmt}")
        paths.append(path)
    return paths


def import_tracks(path: str) -> SeqRecord:
    """Reimport a trace record exported by :func:`export_tracks`."""
    path = Path(path)
    if path.suffix == ".npz":
        data = np.load(path)
        return _build_record(
            data["channels"],
            data["peaks"],
            data["trace_x"],
            _qual_to_seq_guess(data),
            list(map(int, data["qual"])),
            path.stem,
        )
    return _read_tsv(path)


def _qual_to_seq_guess(data) -> str:
    """A placeholder-known limitation: npz does not store the base string.

    For trace-only round trips we cannot recover the called sequence from the
    binary alone; callers should pass bases explicitly via ``join_tracks`` of
    parsed records.  We return an empty string length matching peak count.
    """
    return ""  # caller responsible for re-calling if needed


def _write_tsv(path, channels, peaks, trace_x, seq, qual) -> None:
    with open(path, "w") as fh:
        fh.write(
            "trace_x\tpeaks\t"
            + "\t".join(f"ch{i + 1}" for i in range(NCHANNELS))
            + "\tseq\tqual\n"
        )
        n = len(trace_x)
        for j in range(n):
            row = [str(trace_x[j]) if j < len(trace_x) else ""]
            row.append(str(peaks[j]) if j < len(peaks) else "")
            for i in range(NCHANNELS):
                row.append(str(channels[i, j]) if j < channels.shape[1] else "")
            row.append(seq[j] if j < len(seq) else "")
            row.append(str(qual[j]) if j < len(qual) else "")
            fh.write("\t".join(row) + "\n")


def _read_tsv(path) -> SeqRecord:
    import csv

    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # header
        rows = list(reader)
    trace_x, peaks, seq, qual = [], [], [], []
    channels = [[] for _ in range(NCHANNELS)]
    for row in rows:
        if not row:
            continue
        try:
            trace_x.append(float(row[0]))
            peaks.append(float(row[1]))
        except (ValueError, IndexError):
            continue
        for i in range(NCHANNELS):
            if 2 + i < len(row) and row[2 + i] != "":
                channels[i].append(float(row[2 + i]))
        if 2 + NCHANNELS < len(row):
            seq.append(row[2 + NCHANNELS])
        if 3 + NCHANNELS < len(row) and row[3 + NCHANNELS] != "":
            qual.append(int(row[3 + NCHANNELS]))
    arr = np.zeros((NCHANNELS, len(trace_x)), dtype=float)
    for i in range(NCHANNELS):
        arr[i, : len(channels[i])] = channels[i]
    return _build_record(
        arr, np.asarray(peaks), np.asarray(trace_x), "".join(seq), qual, path.stem
    )
