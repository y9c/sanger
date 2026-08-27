#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Assembly & consensus calling for multiple Sanger reads.

Sanger sequencing usually yields a handful of reads for a region.  This module
implements a lightweight reference-guided "pileup": every read is aligned to a
reference, per-reference-position base counts are aggregated, and a consensus
sequence (with quality) is derived -- the foundation for assembling overlapping
reads into one answer.

The alignment core is the C-accelerated ``ssw`` library (already used by
:mod:`cfutils.align`), so pileup scales well to many/long reads.

Typical usage::

    from cfutils.assembly import pileup, consensus
    reads = [parse_abi(f) for f in ["a.ab1", "b.ab1"]]
    ref   = parse_fasta("ref.fa")
    table = pileup(reads, ref)
    cons  = consensus(table)
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Tuple

from .align import align_chromatograph
from .parser import SeqRecord

if TYPE_CHECKING:
    pass

__all__ = [
    "PileupTable",
    "PileupColumn",
    "pileup",
    "consensus",
    "coverage",
]

#: bases we consider unambiguous for consensus; everything else is N
_CANONICAL = ("A", "C", "G", "T")


@dataclass
class PileupColumn:
    """Aggregated observations at one reference position.

    ``counts`` maps base (or ``"-"`` for a deletion) to the number of reads
    supporting it.  ``n_reads`` is the number of reads covering the position.
    """

    ref_pos: int
    ref_base: str
    counts: Counter = field(default_factory=Counter)
    n_reads: int = 0

    @property
    def depth(self) -> int:
        return sum(self.counts.values())

    def consensus_base(self, min_freq: float = 0.5,
                       min_count: int = 1) -> str:
        """Most-supported canonical base, or ``ref_base`` if no majority.

        A deletion (``"-"``) wins only when it is the clear plurality and the
        favourite canonical base falls below ``min_freq``.
        """
        if not self.counts:
            return self.ref_base
        total = self.depth or 1
        best, best_n = self.counts.most_common(1)[0]
        # only trust unambiguous bases; treat ambiguous as absent
        if best == "-":
            if best_n >= min_count and best_n / total >= min_freq:
                return "-"
            return self.ref_base
        if best in _CANONICAL and best_n >= min_count:
            return best
        return self.ref_base


@dataclass
class PileupTable:
    """Reference-position indexed pileup of aligned reads."""

    columns: Dict[int, PileupColumn] = field(default_factory=dict)

    def __getitem__(self, ref_pos: int) -> PileupColumn:
        return self.columns[ref_pos]

    def __len__(self) -> int:
        return len(self.columns)

    def iter_columns(self) -> List[PileupColumn]:
        return [self.columns[p] for p in sorted(self.columns)]

    def consensus_sequence(self, min_freq: float = 0.5,
                           min_count: int = 1) -> str:
        """Return the consensus string, in reference order (deletions kept as gaps)."""
        return "".join(
            self.columns[p].consensus_base(min_freq, min_count)
            for p in sorted(self.columns)
        )


def _accumulate(columns: Dict[int, PileupColumn], site, read_base: str,
                seen: set, key: tuple) -> None:
    """Record one read's base at a reference position (once per read/position)."""
    if key in seen:
        return
    seen.add(key)
    col = columns.setdefault(site.ref_pos, None)
    if col is None:
        col = PileupColumn(ref_pos=site.ref_pos, ref_base=site.ref_base)
        columns[site.ref_pos] = col
    col.counts[read_base] += 1
    col.n_reads += 1


def pileup(
    reads: Iterable[SeqRecord],
    reference: SeqRecord,
    quality_threshold: int = 0,
    min_cov: int = 1,
) -> PileupTable:
    """Aggregate aligned reads into a reference-indexed pileup.

    Args:
        reads: parsed chromatogram records to pile up onto the reference.
        reference: the reference :class:`SeqRecord` (usually a ``parse_fasta``).
        quality_threshold: drop an aligned base when its site quality is below
            this (0 disables quality filtering).
        min_cov: only keep reference columns covered by at least this many reads.

    Returns:
        A :class:`PileupTable`.
    """
    reads = list(reads)
    n_total = len(reads)
    if n_total == 0:
        raise ValueError("pileup requires at least one read")

    columns: Dict[int, PileupColumn] = {}

    for read_idx, read in enumerate(reads):
        seen: set = set()
        for site in align_chromatograph(read, reference):
            ref_base = site.ref_base
            read_base = site.cf_base
            if read_base == "-":
                base = "-"
            else:
                if quality_threshold and site.qual_site is not None \
                        and site.qual_site < quality_threshold:
                    continue  # low-quality base excluded
                base = read_base if read_base in _CANONICAL else read_base
            _accumulate(columns, site, base, seen, (read_idx, site.ref_pos))

    if min_cov > 1:
        columns = {p: c for p, c in columns.items() if c.n_reads >= min_cov}

    return PileupTable(columns=columns)


def consensus(
    table: PileupTable, min_freq: float = 0.5, min_count: int = 1
) -> str:
    """Consensus string from a :class:`PileupTable`."""
    return table.consensus_sequence(min_freq=min_freq, min_count=min_count)


def coverage(table: PileupTable) -> List[Tuple[int, int]]:
    """Return ``(ref_pos, depth)`` pairs for every pileup column, in order."""
    return [(c.ref_pos, c.depth) for c in table.iter_columns()]
