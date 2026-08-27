#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Quality-score handling and mutation filtering.

Sanger base quality (Phred) values are carried per base in the record's
``letter_annotations["phred_quality"]``. This module centralises the quality
logic that used to live inside the plotting/reporting code so it can be
reused, tested and tuned independently.

Default thresholds reflect typical guidance for Sanger mutation calling:
site quality (the called base) ``>= 20`` and local (window-averaged) quality
``>= 20``.
"""

from __future__ import annotations

from dataclasses import dataclass
from statistics import mean
from typing import TYPE_CHECKING, Iterable, List, Optional, Tuple

if TYPE_CHECKING:
    from .align import SitePair
    from .parser import SeqRecord

__all__ = [
    "QualityFilter",
    "DEFAULT_BASE_QUAL",
    "DEFAULT_LOCAL_QUAL",
    "site_qualities",
    "passed_filter",
]

#: module-level defaults (kept for backward compatibility with old hard-coded 50/20)
DEFAULT_BASE_QUAL = 20
DEFAULT_LOCAL_QUAL = 20


@dataclass
class QualityFilter:
    """Tunable quality thresholds for Sanger mutation calling.

    ``flank_base_num`` controls the window used to compute local quality:
    0 means only the single called base, >0 averages over ``2*n+1`` bases
    centred on the site.
    """

    min_base_qual: int = DEFAULT_BASE_QUAL
    min_local_qual: int = DEFAULT_LOCAL_QUAL
    flank_base_num: int = 0

    def passed(self, site: "SitePair") -> bool:
        """Return True if the site passes both base and local quality gates."""
        if site.qual_site is None or site.qual_local is None:
            return site.qual_site is None and site.qual_local is None or (
                site.qual_site is not None
                and site.qual_site >= self.min_base_qual
                and (site.qual_local is None
                     or site.qual_local >= self.min_local_qual)
            )
        return (
            site.qual_site >= self.min_base_qual
            and site.qual_local >= self.min_local_qual
        )

    def filter(
        self, sites: Iterable["SitePair"]
    ) -> List["SitePair"]:
        """Return only the sites passing the quality filter."""
        return [s for s in sites if self.passed(s)]


def site_qualities(
    query_record: "SeqRecord", pos: int, flank_base_num: int = 0
) -> Tuple[int, int]:
    """Return ``(site_quality, local_quality)`` for a 1-based position.

    ``site_quality`` is the Phred quality of the called base.
    ``local_quality`` is the average quality over ``2*flank_base_num+1``
    centred bases (or just the site when ``flank_base_num == 0``).
    """
    qual = query_record.letter_annotations["phred_quality"]
    idx = pos - 1
    site_q = qual[idx]
    lo = max(0, idx - flank_base_num)
    hi = min(len(qual), idx + flank_base_num + 1)
    local_q = int(mean(qual[lo:hi]))
    return site_q, local_q


def passed_filter(site: "SitePair", min_base_qual: int = DEFAULT_BASE_QUAL,
                  min_local_qual: int = DEFAULT_LOCAL_QUAL) -> bool:
    """Functional helper, mirroring ``QualityFilter.passed`` defaults."""
    return QualityFilter(min_base_qual, min_local_qual).passed(site)
