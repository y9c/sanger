"""Chromatogram File Utils — Sanger sequencing analysis toolkit.

Modules
-------
parser      : ABI/FASTA parsing into lightweight SeqRecord objects.
align       : reference-guided alignment & mutation calling (ssw-backed).
show        : matplotlib rendering of chromatograms and mutations.
features    : feature-annotation overlay API (dna-features-viewer inspired).
quality     : quality-score computation and mutation filtering.
tracks      : split / join / slice chromatogram trace files.
assembly    : multi-read pileup and consensus calling.
composite   : side-by-side plotting with external tools.
run / cli   : high-level reporting and the command-line interface.
"""

from .parser import SeqRecord, parse_abi, parse_fasta
from .align import SitePair, align_chromatograph, call_mutations, run_align
from .show import (
    plot_chromatograph,
    show_reference,
    highlight_base,
    annotate_mutation,
    center_region,
)
from .features import ChromatogramFeature, add_feature, plot_features
from .quality import QualityFilter, site_qualities
from .tracks import join_tracks, split_track, slice_track
from .assembly import pileup, consensus, coverage
from .composite import side_by_side, add_panel
from .basecaller import call_bases, detect_peaks
from .qc import read_metrics, summarize, trimmed_bounds
from .transform import trim, reverse_complement_record
from .dnalink import to_graphic_record, plot_combined

__version__ = "0.0.0.dev62"

__all__ = [
    "SeqRecord", "parse_abi", "parse_fasta",
    "SitePair", "align_chromatograph", "call_mutations", "run_align",
    "plot_chromatograph", "show_reference", "highlight_base",
    "annotate_mutation", "center_region",
    "ChromatogramFeature", "add_feature", "plot_features",
    "QualityFilter", "site_qualities",
    "join_tracks", "split_track", "slice_track",
    "pileup", "consensus", "coverage",
    "side_by_side", "add_panel",
    "call_bases", "detect_peaks",
    "read_metrics", "summarize", "trimmed_bounds",
    "trim", "reverse_complement_record",
    "to_graphic_record", "plot_combined",
    "__version__",
]
