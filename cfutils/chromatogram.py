#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
High-level :class:`Chromatogram` domain object.

The low-level API returns lightweight ``SeqRecord`` objects whose data lives in
ad-hoc ``annotations`` / ``letter_annotations`` dicts (traces, peak positions,
quality, channel order).  This wrapper centralises those accessors into typed
properties and bundles the common operations -- QC, base-calling, mutation
calling, trimming, plotting and exporting -- behind one object.

It is additive and non-breaking: ``parse_abi`` still returns a ``SeqRecord``;
use :func:`Chromatogram.from_abi` (or :meth:`Chromatogram.from_record`) to get
the richer object.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from .parser import SeqRecord, parse_abi

__all__ = ["Chromatogram"]

_NCHANNELS = 4


class Chromatogram:
    """A parsed Sanger chromatogram with typed accessors and operations.

    Wraps a :class:`~cfutils.parser.SeqRecord` and exposes the trace data as
    convenient properties, plus methods that delegate to the analysis modules.
    """

    def __init__(self, record: SeqRecord, path: Optional[str] = None):
        if not isinstance(record, SeqRecord):
            raise TypeError("Chromatogram wraps a SeqRecord")
        self.record = record
        self.path = str(path) if path else None

    # ------------------------------------------------------------------ #
    #  Constructors
    # ------------------------------------------------------------------ #
    @classmethod
    def from_abi(cls, path: str, rescale: bool = True) -> "Chromatogram":
        """Parse an ABI file into a :class:`Chromatogram`."""
        return cls(parse_abi(path, rescale=rescale), path=path)

    @classmethod
    def from_record(cls, record: SeqRecord) -> "Chromatogram":
        """Wrap an already-parsed :class:`SeqRecord`."""
        return cls(record)

    # ------------------------------------------------------------------ #
    #  Accessors
    # ------------------------------------------------------------------ #
    @property
    def name(self) -> str:
        return self.record.name

    @property
    def id(self) -> str:
        return self.record.id

    @property
    def sequence(self) -> str:
        """The called base sequence."""
        return self.record.seq

    @property
    def length(self) -> int:
        return len(self.record)

    @property
    def quality(self) -> np.ndarray:
        """Per-base Phred quality as a float array."""
        return np.asarray(
            self.record.letter_annotations.get("phred_quality", []), dtype=float
        )

    @property
    def peaks(self) -> np.ndarray:
        """Peak x-positions."""
        return np.asarray(
            self.record.annotations.get("peak positions", []), dtype=float
        )

    @property
    def trace_x(self) -> np.ndarray:
        return np.asarray(self.record.annotations.get("trace_x", []), dtype=float)

    @property
    def channels(self) -> str:
        """Channel/base order (e.g. ``GATC``)."""
        return str(self.record.annotations.get("channels", "ACGT"))

    @property
    def traces(self) -> np.ndarray:
        """(4, n_samples) array of the four dye channels."""
        return np.array(
            [self.record.annotations["channel " + str(i)] for i in range(1, 5)],
            dtype=float,
        )

    @property
    def to_record(self) -> SeqRecord:
        """The wrapped low-level :class:`SeqRecord`."""
        return self.record

    # ------------------------------------------------------------------ #
    #  Derived metrics
    # ------------------------------------------------------------------ #
    @property
    def gc_percent(self) -> float:
        s = self.sequence
        return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0

    @property
    def mean_quality(self) -> float:
        q = self.quality
        return float(q.mean()) if q.size else 0.0

    def qc(self) -> dict:
        """Full QC metrics (see :mod:`cfutils.qc`)."""
        from .qc import (
            continuous_read_length,
            noise_metric,
            read_metrics,
            signal_intensity,
        )

        m = read_metrics(self.record)
        m["crl"] = continuous_read_length(self.record)
        m["signal_intensity"] = signal_intensity(self.record)
        m["snr"] = noise_metric(self.record)
        return m

    # ------------------------------------------------------------------ #
    #  Operations (delegate to the analysis modules)
    # ------------------------------------------------------------------ #
    def slice(self, start: int, end: int, name: Optional[str] = None) -> "Chromatogram":
        """Return a sub-region ``[start, end]`` as a new Chromatogram."""
        from .tracks import slice_track

        return Chromatogram(
            slice_track(self.record, start, end, name=name), path=self.path
        )

    def trim(
        self, cutoff: float = 0.05, segment: int = 20, name: Optional[str] = None
    ) -> "Chromatogram":
        """Mott quality-trim (keeps traces/peaks aligned)."""
        from .transform import trim

        return Chromatogram(
            trim(self.record, cutoff=cutoff, segment=segment, name=name), path=self.path
        )

    def trim_leading_ns(self) -> "Chromatogram":
        from .transform import trim_leading_ns

        return Chromatogram(trim_leading_ns(self.record), path=self.path)

    def reverse_complement(self) -> "Chromatogram":
        from .transform import reverse_complement_record

        return Chromatogram(reverse_complement_record(self.record), path=self.path)

    def basecall(self, hetero_threshold: float = 0.45):
        """Re-call bases from the raw (unrescaled) four-channel traces."""
        from .basecaller import call_bases

        raw = self.record
        # base/peak analysis needs raw trace coordinates; re-parse if possible
        if self.path:
            raw = parse_abi(self.path, rescale=False)
        return call_bases(raw, hetero_threshold=hetero_threshold)

    def read_quality(self):
        """Alias for :meth:`basecall` quality summary."""
        return self.basecall()

    def call_mutations(self, reference: SeqRecord, report_all: bool = False):
        """Call variants against a reference.  Returns a list of SitePairs."""
        from .align import call_mutations

        return call_mutations(self.record, reference, report_all_sites=report_all)

    def detect_orientation(self, reference: SeqRecord) -> int:
        """Return +1 (forward) or -1 (reverse-complement) vs ``reference``."""
        from .align import detect_orientation

        rec = self.record
        if self.path:
            rec = self.record
        return detect_orientation(rec, reference)

    def analyze(
        self,
        kind: str = "translate",
        motif: Optional[str] = None,
        frame: int = 1,
        both_strands: bool = False,
    ):
        """Sequence-level analysis (translate / motif / restriction / gc)."""
        from .analysis import find_motifs, gc_windows, restriction_sites, translate

        seq = self.sequence
        if kind == "translate":
            return translate(seq, frame=frame)
        if kind == "motif":
            if not motif:
                raise ValueError("kind='motif' requires a motif sequence")
            return find_motifs(seq, motif, both_strands=both_strands)
        if kind == "restriction":
            return {k: v for k, v in restriction_sites(seq).items() if v}
        if kind == "gc":
            return gc_windows(seq)
        raise ValueError(f"unknown analysis: {kind!r}")

    def to_fasta(self, width: int = 60) -> str:
        from .export import to_fasta

        return to_fasta(self.record, name=self.name, width=width)

    def to_vcf(self, reference: SeqRecord, min_qual: int = 0) -> str:
        from .align import call_mutations
        from .export import to_vcf

        sites = call_mutations(self.record, reference, report_all_sites=True)
        return to_vcf(sites, reference_name=reference.name or "ref", min_qual=min_qual)

    def export(self, outdir: str, fmt: str = "fasta") -> str:

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        dest = outdir / f"{self.name or 'seq'}.fa"
        dest.write_text(self.to_fasta())
        return str(dest)

    def plot(self, region: Optional[Tuple[int, int]] = None, ax=None):
        """Render the chromatogram (and return ``(fig, ax)`` or ``ax``)."""
        from ._mpl import require_matplotlib

        require_matplotlib()
        from .show import plot_chromatograph

        if ax is None:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(16, 5))
            plot_chromatograph(self.record, region=region, ax=ax)
            return fig, ax
        plot_chromatograph(self.record, region=region, ax=ax)
        return ax

    # ------------------------------------------------------------------ #
    #  dunders
    # ------------------------------------------------------------------ #
    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        return self.sequence

    def __repr__(self) -> str:
        return f"Chromatogram('{self.name}', len={self.length})"
