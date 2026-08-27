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

from .align import (
    SitePair,
    align_chromatograph,
    call_mutations,
    detect_orientation,
    run_align,
)
from .analysis import find_motifs, gc_windows, restriction_sites, translate
from .assembly import consensus, coverage, pileup
from .basecaller import call_bases, detect_peaks
from .chromatogram import Chromatogram
from .composite import add_panel, side_by_side
from .dnalink import plot_combined, to_graphic_record
from .export import batch_summary, to_fasta, to_json, to_vcf, write_batch
from .features import ChromatogramFeature, add_feature, plot_features
from .parser import SeqRecord, parse_abi, parse_fasta
from .qc import (
    continuous_read_length,
    noise_metric,
    read_metrics,
    signal_intensity,
    summarize,
    trimmed_bounds,
)
from .quality import QualityFilter, site_qualities
from .show import (
    annotate_mutation,
    center_region,
    highlight_base,
    plot_chromatograph,
    show_reference,
)
from .tracks import join_tracks, slice_track, split_track
from .transform import (
    reverse_complement_record,
    strip_primers,
    trim,
    trim_ends,
    trim_leading_ns,
)

__version__ = "0.0.0.dev62"

__all__ = [
    "SeqRecord",
    "parse_abi",
    "parse_fasta",
    "SitePair",
    "align_chromatograph",
    "call_mutations",
    "run_align",
    "plot_chromatograph",
    "show_reference",
    "highlight_base",
    "annotate_mutation",
    "center_region",
    "ChromatogramFeature",
    "add_feature",
    "plot_features",
    "QualityFilter",
    "site_qualities",
    "join_tracks",
    "split_track",
    "slice_track",
    "pileup",
    "consensus",
    "coverage",
    "side_by_side",
    "add_panel",
    "call_bases",
    "detect_peaks",
    "read_metrics",
    "summarize",
    "trimmed_bounds",
    "continuous_read_length",
    "signal_intensity",
    "noise_metric",
    "trim",
    "trim_ends",
    "trim_leading_ns",
    "strip_primers",
    "reverse_complement_record",
    "to_graphic_record",
    "plot_combined",
    "translate",
    "find_motifs",
    "restriction_sites",
    "gc_windows",
    "to_fasta",
    "to_vcf",
    "to_json",
    "batch_summary",
    "write_batch",
    "detect_orientation",
    "Chromatogram",
    "__version__",
]
