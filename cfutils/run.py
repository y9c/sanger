#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Chromatogram File Utils.

do some wrap functions
"""

import os
from datetime import datetime
from pathlib import Path
from typing import List

from cfutils.align import call_mutations
from cfutils.parser import SeqRecord, parse_abi, parse_fasta
from cfutils.quality import QualityFilter
from cfutils.show import annotate_mutation, highlight_base, plot_chromatograph

from ._mpl import plt, require_matplotlib
from .utils import get_logger

LOGGER = get_logger(__name__)


def do_mutation_showing(
    query_record: SeqRecord,
    mutations: List,
    output_fig_file: str,
    quality_filter: QualityFilter = None,
) -> None:
    """report mutations in pdf format."""
    require_matplotlib()
    if quality_filter is None:
        quality_filter = QualityFilter(min_base_qual=20, min_local_qual=20)

    mutations = sorted(mutations, key=lambda x: x.cf_pos)
    flanking_size = 6
    windows_size = 30
    mutation_windows = []
    start_pos = max(1, mutations[0].cf_pos - flanking_size)
    mutation_region = []
    for idx, mut in enumerate(mutations):
        if mut.cf_pos + flanking_size <= start_pos + windows_size:
            mutation_region.append(mut)
        else:
            mutation_windows.append(mutation_region)
            start_pos = max(1, mutations[idx].cf_pos - flanking_size)
            mutation_region = [mut]
    mutation_windows.append(mutation_region)

    fig, axes = plt.subplots(
        len(mutation_windows), figsize=(20, 5 * len(mutation_windows))
    )
    for idx, mutation_region in enumerate(mutation_windows):
        if len(mutation_windows) == 1:
            ax = axes
        else:
            ax = axes[idx]
        region_start = max(1, mutation_region[0].cf_pos - flanking_size)
        plot_chromatograph(
            query_record,
            region=(region_start, region_start + windows_size),
            ax=ax,
        )
        for mut in mutation_region:
            base_passed = quality_filter.passed(mut)
            highlight_base(mut.cf_pos, query_record, ax, passed_filter=base_passed)
            annotate_mutation(mut, query_record, ax)
    fig.savefig(output_fig_file, bbox_inches="tight")


def report_mutation(
    query_ab1_file,
    subject_fasta_file,
    output_dir=None,
    file_basename=None,
    report_all_sites=False,
    report_mut_plot=False,
):
    """reprot mutation within region."""
    if output_dir is None:
        output_dir = os.path.join(
            os.getcwd(),
            "CFresult_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
        )
    os.makedirs(output_dir, exist_ok=True)

    if file_basename is None:
        file_basename = (
            Path(query_ab1_file).stem + "_vs_" + Path(subject_fasta_file).stem
        )

    query_record = parse_abi(query_ab1_file)
    subject_record = parse_fasta(subject_fasta_file)

    sites = call_mutations(
        query_record,
        subject_record,
        report_all_sites=report_all_sites,
    )
    # save mutation / alignment to tsv file
    with open(os.path.join(output_dir, file_basename + ".tsv"), "w") as f_mut:
        header = [
            "RefLocation",
            "RefBase",
            "CfLocation",
            "CfBase",
            "SiteQual",
            "LocalQual",
        ]
        f_mut.write("\t".join(header) + "\n")
        for site in sites:
            f_mut.write(
                f"{site.ref_pos}\t{site.ref_base}\t{site.cf_pos}\t{site.cf_base}\t{site.qual_site}\t{site.qual_local}\n"
            )

    # do forget to filter mutation for plot
    if report_all_sites:
        mutations = [s for s in sites if s.ref_base != s.cf_base]
        LOGGER.info(f"{query_record.name}: Mutation number for plot: {len(mutations)}")
    else:
        mutations = sites

    # show mutation in pdf file
    if mutations and report_mut_plot:
        output_fig_file = os.path.join(output_dir, file_basename + ".pdf")
        do_mutation_showing(query_record, mutations, output_fig_file)
