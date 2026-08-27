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

import ssw
from .parser import SeqRecord

from .utils import get_logger

LOGGER = get_logger(__name__)


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
    aligner = ssw.Aligner()
    alignment = aligner.align(reference=reference, query=query)
    results = []
    query_pos = alignment.query_begin
    ref_pos = alignment.reference_begin
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
