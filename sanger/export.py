#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Export Sanger results to standard lab / pipeline formats.

Covers the formats a wet-lab or downstream pipeline actually needs:

* **FASTA**  -- the called / trimmed / consensus sequence.
* **VCF**    -- variant calls against a reference (the interchange standard).
* **JSON**   -- a self-describing per-read report (QC + basecall + variants).
* **CSV/TSV**-- a batch summary table for many reads.

The module is symmetric with :mod:`sanger.parser` so records round-trip.
"""

from __future__ import annotations

import json
from typing import Dict, Iterable, List, Optional, Sequence

from .parser import SeqRecord

__all__ = [
    "to_fasta",
    "to_vcf",
    "to_json",
    "batch_summary",
    "write_batch",
    "MISSING",
]

MISSING = "."


def to_fasta(record: SeqRecord, name: Optional[str] = None, width: int = 60) -> str:
    """Return a FASTA string of the record's called sequence (wrapped)."""
    name = name or record.name or record.id or "seq"
    seq = record.seq
    lines = [f">{name}"]
    for i in range(0, len(seq), width):
        lines.append(seq[i : i + width])
    return "\n".join(lines) + "\n"


def to_vcf(
    variants: Sequence,
    sample_id: str = "SAMPLE",
    reference_name: str = "ref",
    min_qual: int = 0,
) -> str:
    """Render variant sites as a VCF 4.2 string.

    ``variants`` are :class:`~sanger.align.SitePair` (or any object exposing
    ``ref_pos``, ``ref_base``, ``cf_base``, ``qual_site``).  Only sites where
    both REF and ALT are canonical and different are emitted; indels (``-``)
    are skipped.
    """
    header = [
        "##fileformat=VCFv4.2",
        "##source=sanger",
        f"##reference={reference_name}",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    rows = []
    for v in variants:
        ref, alt = (v.ref_base or "").upper(), (v.cf_base or "").upper()
        if ref in ("-", "N", "") or alt in ("-", "N", ""):
            continue
        if ref == alt:
            continue
        qual = v.qual_site if v.qual_site is not None else MISSING
        if qual != MISSING and min_qual and qual < min_qual:
            continue
        rows.append(f"{reference_name}\t{v.ref_pos}\t.\t{ref}\t{alt}\t{qual}\t.\t.")
    return "\n".join(header + rows) + "\n"


def to_json(record: SeqRecord, extra: Optional[dict] = None, indent: int = 2) -> str:
    """Return a self-describing JSON report for a single record.

    Includes identifiers, sequence composition, quality summary and any
    attached features/mutations supplied via ``extra``.
    """
    from .qc import read_metrics

    m = read_metrics(record)
    feats = [
        {
            "start": f.start,
            "end": f.end,
            "strand": f.strand,
            "color": f.color,
            "label": f.label,
            "kind": f.kind,
        }
        for f in record.annotations.get("features", [])
    ]
    payload = {
        "name": record.name,
        "id": record.id,
        "sequence": record.seq,
        "length": len(record),
        "qc": m,
        "features": feats,
    }
    if extra:
        payload.update(extra)
    return json.dumps(payload, indent=indent)


def batch_summary(
    records: Iterable[SeqRecord], include_bases: bool = False
) -> List[Dict]:
    """Return a list of per-read summary dicts (the QC table rows)."""
    from .qc import read_metrics

    rows = []
    for rec in records:
        m = read_metrics(rec)
        row = {
            "name": rec.name,
            "id": rec.id,
            "bases": int(m["n_bases"]),
            "gc_percent": round(m["gc_percent"], 1),
            "n_fraction": round(m["n_fraction"], 3),
            "mean_qual": round(m["mean_qual"], 1),
            "min_qual": int(m["min_qual"]),
            "low_qual_fraction": round(m["low_qual_fraction"], 3),
            "trim_start": int(m["trim_start"]),
            "trim_end": int(m["trim_end"]),
            "trimmed_len": int(m["trimmed_len"]),
        }
        if include_bases:
            row["sequence"] = rec.seq
        rows.append(row)
    return rows


def write_batch(
    records: Iterable[SeqRecord], outdir: str, fmt: str = "csv", name: str = "summary"
) -> str:
    """Write a batch QC summary of many reads.

    Args:
        records: chromatogram records (or a list of files to be parsed).
        outdir: destination directory (created if needed).
        fmt: ``csv``, ``tsv`` or ``json``.
        name: output basename (without extension).

    Returns:
        Path of the written file.
    """
    import csv
    from pathlib import Path

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = batch_summary(records)
    keys = list(rows[0].keys()) if rows else []

    if fmt == "json":
        path = outdir / f"{name}.json"
        path.write_text(json.dumps(rows, indent=2))
    elif fmt in ("csv", "tsv"):
        delim = "," if fmt == "csv" else "\t"
        path = outdir / f"{name}.{fmt}"
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=keys, delimiter=delim)
            writer.writeheader()
            writer.writerows(rows)
    else:
        raise ValueError(f"unknown batch format: {fmt}")
    return str(path)
