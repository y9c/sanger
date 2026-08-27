#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Model Context Protocol (MCP) server for sanger.

Exposes sanger as agent-friendly tools so LLM agents / MCP clients can analyse
Sanger sequencing data programmatically.  The server runs locally and reads
files from the host, so it is intended for local (or sandboxed) use.

Run with::

    sanger-mcp                       # stdio transport
    python -m sanger.mcp_server       # identical

Tools (each returns JSON-serialisable data):

* ``read_chromatogram``   - summary of an ABI read
* ``qc_metrics``          - full per-read QC metrics (incl. CRL, SNR)
* ``call_mutations``      - variants vs a reference (SNPs/indels)
* ``re_call_bases``       - re-call bases from raw traces + heterozygotes
* ``analyze_sequence``    - translate / motifs / restriction sites / GC
* ``trim_read``           - quality-trim a read
* ``export_sequence``     - write FASTA (or VCF with a reference)
* ``plot_chromatogram``   - render a PNG of a chromatogram region

Optional extras: ``pip install "sanger[agent]"`` (i.e. ``mcp>=2``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

try:
    from mcp.server.mcpserver import MCPServer

    _HAVE_MCP = True
except Exception:  # pragma: no cover - mcp not installed
    MCPServer = None
    _HAVE_MCP = False

from .parser import parse_abi, parse_fasta

mcp = MCPServer("sanger")


def _require_plot():
    from ._mpl import require_matplotlib

    require_matplotlib()


# --------------------------------------------------------------------------- #
#  I/O / analysis tools
# --------------------------------------------------------------------------- #
def read_chromatogram(path: str) -> dict:
    """Parse an ABI chromatogram and return a concise summary."""
    rec = parse_abi(path)
    from .qc import continuous_read_length, noise_metric, read_metrics, signal_intensity

    m = read_metrics(rec)
    return {
        "name": rec.name,
        "id": rec.id,
        "length": len(rec),
        "sequence": rec.seq,
        "gc_percent": round(m["gc_percent"], 2),
        "mean_qual": round(m["mean_qual"], 2),
        "min_qual": int(m["min_qual"]),
        "n_fraction": round(m["n_fraction"], 4),
        "low_qual_fraction": round(m["low_qual_fraction"], 4),
        "trim_interval": [int(m["trim_start"]), int(m["trim_end"])],
        "crl": continuous_read_length(rec),
        "signal_intensity": round(signal_intensity(rec), 1),
        "snr": round(noise_metric(rec), 1),
    }


def qc_metrics(path: str) -> dict:
    """Return the full QC metric dict for a chromatogram."""
    from .qc import continuous_read_length, noise_metric, read_metrics, signal_intensity

    rec = parse_abi(path)
    m = read_metrics(rec)
    m["crl"] = continuous_read_length(rec)
    m["signal_intensity"] = signal_intensity(rec)
    m["snr"] = noise_metric(rec)
    return {k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()}


def call_mutations(
    query_ab1: str,
    subject_fasta: str,
    report_all: bool = False,
    min_base_qual: int = 20,
    min_local_qual: int = 20,
) -> dict:
    """Call variants of ``query_ab1`` against ``subject_fasta``.

    Returns a dict with ``n_aligned``, ``variants`` (list of
    {ref_pos,ref_base,cf_base,quality}) and ``confident`` variants passing the
    quality filter.  Indels are flagged with ``-`` in ref_base/cf_base.
    """
    from .align import align_chromatograph
    from .align import call_mutations as cm
    from .quality import QualityFilter

    q = parse_abi(query_ab1)
    r = parse_fasta(subject_fasta)
    all_sites = align_chromatograph(q, r)
    sites = cm(q, r, report_all_sites=True)
    qf = QualityFilter(min_base_qual=min_base_qual, min_local_qual=min_local_qual)
    var = []
    for s in sites:
        if s.ref_base == s.cf_base:
            continue
        var.append(
            {
                "ref_pos": s.ref_pos,
                "ref_base": s.ref_base,
                "cf_pos": s.cf_pos,
                "cf_base": s.cf_base,
                "quality": s.qual_site,
                "local_qual": s.qual_local,
                "passed": qf.passed(s),
            }
        )
    return {
        "query": q.name,
        "reference": r.name,
        "n_aligned": len(all_sites),
        "n_variants": len(var),
        "n_confident": sum(1 for v in var if v["passed"]),
        "variants": var,
    }


def re_call_bases(path: str, hetero_threshold: float = 0.45) -> dict:
    """Re-call bases from the raw four-channel traces of ``path``.

    Returns the called sequence, mean quality, and any heterozygous/mixed-base
    sites detected.
    """
    from .basecaller import call_bases

    rec = parse_abi(path, rescale=False)
    res = call_bases(rec, hetero_threshold=hetero_threshold)
    return {
        "name": rec.name,
        "n_calls": res.n_calls,
        "n_ambiguous": res.n_ambiguous,
        "mean_quality": round(float(sum(res.qualities) / len(res.qualities)), 1)
        if res.qualities
        else 0.0,
        "sequence": res.sequence,
        "heterozygotes": [
            {"pos": p, "major": a, "minor": b, "minor_fraction": f}
            for p, a, b, f in res.heterozygotes()
        ],
    }


def analyze_sequence(
    path: str,
    analysis: str = "translate",
    motif: Optional[str] = None,
    frame: int = 1,
    both_strands: bool = False,
) -> dict:
    """Run a sequence-level analysis on the called bases.

    ``analysis`` is one of ``translate``, ``motif``, ``restriction``, ``gc``.
    """
    from .analysis import find_motifs, gc_windows, restriction_sites, translate

    seq = parse_abi(path).seq
    if analysis == "translate":
        return {"result": translate(seq, frame=frame)}
    if analysis == "motif":
        if not motif:
            raise ValueError("analysis='motif' requires a --motif sequence")
        return {"positions": find_motifs(seq, motif, both_strands=both_strands)}
    if analysis == "restriction":
        return {k: v for k, v in restriction_sites(seq).items() if v}
    if analysis == "gc":
        return {"windows": [{"pos": p, "gc": g} for p, g in gc_windows(seq)]}
    raise ValueError(f"unknown analysis: {analysis!r}")


def trim_read(
    path: str,
    mode: str = "mott",
    cutoff: float = 0.05,
    min_qual: int = 20,
    strip_ns: bool = True,
) -> dict:
    """Quality-trim a read and return the trimmed sequence.

    ``mode`` is ``mott`` (cumulative-score segment), ``ends`` (hard-trim low-
    quality ends) or ``none``.  Optionally reports length before/after.
    """
    from .transform import trim, trim_ends, trim_leading_ns

    rec = parse_abi(path)
    before = len(rec)
    if mode == "mott":
        out = trim(rec, cutoff=cutoff)
    elif mode == "ends":
        out = trim_ends(rec, min_qual=min_qual)
    elif mode == "none":
        out = rec
    else:
        raise ValueError(f"unknown trim mode: {mode!r}")
    if strip_ns and len(out) > 2:
        out = trim_leading_ns(out)
    return {
        "name": out.name,
        "length_before": before,
        "length_after": len(out),
        "sequence": out.seq,
    }


def export_sequence(
    path: str, outdir: str, fmt: str = "fasta", reference: Optional[str] = None
) -> dict:
    """Write the called sequence (or VCF vs ``reference``) to ``outdir``."""
    from .align import call_mutations as cm
    from .export import to_fasta, to_vcf

    rec = parse_abi(path)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    if fmt == "fasta":
        dest = out / f"{rec.name or 'seq'}.fa"
        dest.write_text(to_fasta(rec))
    elif fmt == "vcf":
        if not reference:
            raise ValueError("export_sequence fmt='vcf' requires a --reference")
        ref = parse_fasta(reference)
        sites = cm(rec, ref, report_all_sites=True)
        dest = out / f"{rec.name or 'seq'}.vcf"
        dest.write_text(to_vcf(sites, reference_name=ref.name))
    else:
        raise ValueError(f"unknown export fmt: {fmt!r}")
    return {"written": str(dest)}


def plot_chromatogram(
    path: str,
    out: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    plot_features: bool = False,
) -> dict:
    """Render a PNG of a chromatogram region to ``out`` (returns its path)."""
    _require_plot()
    import matplotlib.pyplot as plt

    from .show import plot_chromatograph

    rec = parse_abi(path)
    region = (start, end if end is not None else start + 30) if start else None
    fig, ax = plt.subplots(figsize=(16, 5))
    plot_chromatograph(rec, region=region, ax=ax)
    if plot_features:
        from .features import plot_features as pf

        pf(rec, ax)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return {"written": str(Path(out).resolve())}


# --------------------------------------------------------------------------- #
#  Register tools with the MCP server
# --------------------------------------------------------------------------- #
if _HAVE_MCP:
    mcp.tool()(read_chromatogram)
    mcp.tool()(qc_metrics)
    mcp.tool()(call_mutations)
    mcp.tool()(re_call_bases)
    mcp.tool()(analyze_sequence)
    mcp.tool()(trim_read)
    mcp.tool()(export_sequence)
    mcp.tool()(plot_chromatogram)


def main():
    """Run the MCP server over stdio (the default for local agents)."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
