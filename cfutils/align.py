#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-05-27 20:19


"""align query sequence with ref.

Use 1-based for all the position
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from .parser import SeqRecord
from .utils import get_logger

LOGGER = get_logger(__name__)

try:  # ssw is a compiled extension; fall back to a pure-Python aligner
    import ssw

    _HAVE_SSW = True
except Exception:  # pragma: no cover - ssw not installed
    ssw = None
    _HAVE_SSW = False

try:  # our own Cython Smith-Waterman accelerator (self-contained)
    from ._swalign import sw_align as _cy_swalign

    _HAVE_CY = True
except Exception:  # pragma: no cover - cython extension not built
    _cy_swalign = None
    _HAVE_CY = False


class _PyAlignment:
    """Minimal stand-in for the ``ssw`` alignment result object."""

    def __init__(self, reference_begin, query_begin, query_aligned, ref_aligned):
        self.reference_begin = reference_begin
        self.query_begin = query_begin
        # ``alignment`` is a (query, scores, reference) triple so that
        # ``zip(*alignment.alignment)`` yields per-column base triples.
        self.alignment = (query_aligned, [0] * len(query_aligned), ref_aligned)


def _sw_align(reference: str, query: str, match=2, mismatch=-1, gap=-1):
    """Pure-Python/NumPy Smith-Waterman local alignment of ``query`` onto ``reference``.

    Returns a :class:`_PyAlignment`.  Used only when the compiled ``ssw``
    extension is unavailable, so cfutils stays dependency-light.  The DP fill
    is vectorised per reference row with NumPy (linear-gap recurrence), giving
    near-native speed on reads of a few kilobases.
    """
    n, m = len(reference), len(query)
    if n == 0 or m == 0:
        return _PyAlignment(0, 0, query, reference)
    # encode bases to small ints for fast equality comparisons
    code = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
    ref_i = np.fromiter(
        (code.get(c.upper(), 4) for c in reference), dtype=np.int8, count=n
    )
    qry_i = np.fromiter((code.get(c.upper(), 4) for c in query), dtype=np.int8, count=m)
    s_row = (ref_i[:, None] == qry_i[None, :]).astype(np.int8)

    jn = np.arange(1, m + 1)
    gap_j = gap * jn  # gap penalty accumulated along a row
    # tr: 0 stop, 1 diag, 2 up(ref gap), 3 left(query gap)
    tr = np.zeros((n + 1, m + 1), dtype=np.uint8)
    prev = np.zeros(m + 1, dtype=np.int32)
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        s = s_row[i - 1]
        up = prev[1:] + gap
        diag = prev[:-1] + np.where(s == 1, match, mismatch)
        base = np.maximum(up, diag)
        base = np.maximum(base, 0)
        # left propagation: f[j] = max(_j<=j, base[_j] + gap*(j-_j))
        f = np.maximum.accumulate(base - gap_j) + gap_j
        cur = np.zeros(m + 1, dtype=np.int32)
        cur[1:] = f.astype(np.int32)
        # choose move: 0 if f<=0, else 1 diag, 2 up, 3 left
        mv = np.zeros(m, dtype=np.uint8)
        diag_won = f == diag
        up_won = (~diag_won) & (f == up)
        mv[(f > 0) & diag_won] = 1
        mv[(f > 0) & up_won] = 2
        mv[(f > 0) & (~diag_won) & (~up_won)] = 3
        tr[i, 1:] = mv
        prev = cur
        mi = int(np.argmax(f))
        if f[mi] > best:
            best, bi, bj = int(f[mi]), i, mi + 1

    # traceback
    i, j = bi, bj
    q_al, r_al = [], []
    while i > 0 and j > 0:
        move = tr[i, j]
        if move == 0:
            break
        if move == 1:
            q_al.append(query[j - 1])
            r_al.append(reference[i - 1])
            i -= 1
            j -= 1
        elif move == 2:
            q_al.append("-")
            r_al.append(reference[i - 1])
            i -= 1
        else:
            q_al.append(query[j - 1])
            r_al.append("-")
            j -= 1
    q_al = "".join(reversed(q_al))
    r_al = "".join(reversed(r_al))
    return _PyAlignment(i, j, q_al, r_al)


def _align(reference: str, query: str):
    """Dispatch to the fastest available aligner.

    Priority: our Cython SW (self-contained, C speed) -> ssw (compiled, if
    installed) -> pure-Python/NumPy Smith-Waterman.
    """
    if _HAVE_CY:
        bi, bj, q_al, r_al = _cy_swalign(reference, query)
        return _PyAlignment(bi, bj, q_al, r_al)
    if _HAVE_SSW:
        return ssw.Aligner().align(reference=reference, query=query)
    return _sw_align(reference, query)


@dataclass
class SitePair:
    """Object for storing align pair at mutation site."""

    ref_pos: int
    ref_base: str
    cf_pos: int
    cf_base: str
    qual_site: Optional[int] = None
    qual_local: Optional[int] = None

    def __repr__(self):
        return f"{self.ref_base}({self.ref_pos})->{self.cf_base}({self.cf_pos})"


def run_align(reference: str, query: str) -> List[SitePair]:
    """Align query sequence with reference sequence.

    Args:
        reference (str): The reference sequence.
        query (str): The query sequence.

    Returns:
        List[SitePair]: A list of SitePair objects representing alignment.
    """
    # normalise case so that lowercase reference/query bases do not create
    # spurious "mutations" from case-only differences
    reference = str(reference).upper()
    query = str(query).upper()
    alignment = _align(reference, query)
    results = []
    # begin positions are 0-based; cfutils uses 1-based positions throughout
    query_pos = alignment.query_begin + 1
    ref_pos = alignment.reference_begin + 1
    for query_base, _, ref_base in zip(*alignment.alignment):
        results.append(
            SitePair(
                ref_pos=ref_pos,
                ref_base=ref_base,
                cf_pos=query_pos,
                cf_base=query_base,
            )
        )
        if query_base != "-":
            query_pos += 1
        if ref_base != "-":
            ref_pos += 1
    return results


def get_quality(pos: int, query_record: SeqRecord, flank_base_num=0) -> Tuple[int, int]:
    """get quality of site and local region.

    change flank_base_num to number gt 0 to get mean qual within region
    """
    qual = query_record.letter_annotations["phred_quality"]
    qual_site = qual[pos - 1]
    qual_flank = qual[
        max(0, pos - 1 - flank_base_num) : min(len(qual), pos + flank_base_num)
    ]
    qual_local = int(sum(qual_flank) / len(qual_flank))
    return qual_site, qual_local


def align_chromatograph(
    query_record: SeqRecord, subject_record: SeqRecord
) -> List[SitePair]:
    """run align.

    @return: list of SitePair about all sites
    """
    sitepairs = run_align(
        reference=str(subject_record.seq), query=str(query_record.seq)
    )
    LOGGER.info(f"{query_record.name}: Total aligned number: {len(sitepairs)}")
    for site in sitepairs:
        site.qual_site, site.qual_local = get_quality(
            site.cf_pos, query_record, flank_base_num=5
        )
        LOGGER.debug(f"{site}\tlocal:{site.qual_local}\tsite:{site.qual_site}")
    return sitepairs


def call_mutations(
    query_record: SeqRecord,
    subject_record: SeqRecord,
    report_all_sites: bool = False,
) -> List[SitePair]:
    """run align and call mutations.

    @return: list of SitePair about mutation sites
    """
    sitepairs = align_chromatograph(query_record, subject_record)
    mutations = []
    for site in sitepairs:
        if report_all_sites:
            mutations.append(site)
            LOGGER.debug(f"Site ({site}) is reported!")
        else:
            if site.ref_base != site.cf_base:
                mutations.append(site)
                LOGGER.debug(f"Site ({site}) is with mutation!")
    if not report_all_sites:
        LOGGER.info(f"{query_record.name}: Total mutation number: {len(mutations)}")
    return mutations


def detect_orientation(
    query_record: SeqRecord, subject_record: SeqRecord, window: int = 200
) -> int:
    """Decide whether a read is forward or reverse relative to the reference.

    Compares the score of aligning the read's 5' segment to the reference in
    the forward orientation against the reverse-complemented orientation and
    returns ``+1`` (forward) or ``-1`` (reverse-complement).  Usefully,
    ``rc == -1`` means you should :func:`cfutils.transform.reverse_complement_record`
    before analysis.  Returns ``+1`` on tie/low-signal.
    """
    from .utils import reverse_complement

    ref = str(subject_record.seq).upper()
    q_head = str(query_record.seq)[:window].upper()
    rc = reverse_complement(q_head)
    fwd = run_align(ref, q_head)
    rev = run_align(ref, rc)
    fwd_score = sum(1 for s in fwd if s.cf_base != "-" and s.ref_base != "-")
    rev_score = sum(1 for s in rev if s.cf_base != "-" and s.ref_base != "-")
    return +1 if fwd_score >= rev_score else -1
