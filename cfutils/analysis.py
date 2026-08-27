#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sequence-level biology helpers for Sanger reads.

Beyond the chromatogram mechanics, labs routinely need simple per-sequence
computations on the called DNA: translate a coding region, find restriction
sites / primer motifs, compute GC windows, and detect a poly-A / repeat tail.
This module is deliberately small and dependency-free.

All coordinates returned are 1-based, consistent with the rest of cfutils.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

__all__ = [
    "CODON_TABLE",
    "translate",
    "find_motifs",
    "restriction_sites",
    "RE_BASIC",
    "gc_windows",
    "reverse_complement",
]

#: standard genetic code (U/C/T are all accepted for T)
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

#: a small built-in set of common type-II restriction enzymes:
#: name -> (recognition site, cut position counting from the site's start).
RE_BASIC = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "XhoI": ("CTCGAG", 1),
    "SalI": ("GTCGAC", 1),
    "NotI": ("GCGGCCGC", 2),
    "KpnI": ("GGTACC", 1),
    "PstI": ("CTGCAG", 1),
    "NcoI": ("CCATGG", 1),
    "SacI": ("GAGCTC", 1),
}


def translate(seq: str, frame: int = 1, stop_as="X") -> str:
    """Translate a DNA sequence in the given 1-based reading frame.

    ``*`` stop codons are emitted as ``stop_as`` (default ``X``).  Incomplete
    final codons and codons containing ``N`` are emitted as ``stop_as``.
    """
    frame = ((frame - 1) % 3) + 1
    s = seq.upper().replace("U", "T")
    start = frame - 1
    protein = []
    for i in range(start, len(s) - 2, 3):
        codon = s[i : i + 3]
        if "N" in codon:
            protein.append(stop_as)
        else:
            aa = CODON_TABLE.get(codon, stop_as)
            protein.append(stop_as if aa == "*" else aa)
    return "".join(protein)


def find_motifs(seq: str, motif: str, both_strands: bool = False) -> List[int]:
    """Return 1-based start positions where ``motif`` occurs in ``seq``.

    When ``both_strands`` is True, occurrences on the reverse strand are also
    reported — i.e. positions where ``reverse_complement(motif)`` appears in
    the forward sequence (the reverse strand then *reads* as ``motif``).
    """
    s = seq.upper()
    motif = motif.upper()
    hits = _find_all(s, motif)
    if both_strands:
        hits += _find_all(s, reverse_complement(motif))
    return sorted(set(hits))


def restriction_sites(
    seq: str,
    enzymes: Optional[Dict[str, Tuple[str, int]]] = None,
    both_strands: bool = True,
) -> Dict[str, List[int]]:
    """Scan a sequence for restriction sites.

    ``enzymes`` maps enzyme names to ``(recognition_site, cut_offset)``; it
    defaults to :data:`RE_BASIC`.  Returned values are 1-based positions of the
    recognition site start on the forward strand (a site on the reverse strand
    is reported at its forward-strand occurrence of the site sequence).

    Returns a dict keyed by enzyme name, values are lists of positions.
    """
    enzymes = enzymes or RE_BASIC
    s = seq.upper()
    out: Dict[str, List[int]] = {}
    for name, (site, cut) in enzymes.items():
        site = site.upper()
        fwd = _find_all(s, site)
        if both_strands:
            fwd += _find_all(s, reverse_complement(site))
        out[name] = sorted(set(fwd))
    return out


def _find_all(text: str, sub: str) -> List[int]:
    return [
        i + 1 for i in range(len(text) - len(sub) + 1) if text[i : i + len(sub)] == sub
    ]


def gc_windows(seq: str, window: int = 30, step: int = 1) -> List[Tuple[int, float]]:
    """Sliding GC% windows.

    Returns ``(center_position, gc_percent)`` pairs, positions 1-based.
    """
    s = seq.upper()
    out = []
    for i in range(0, len(s) - window + 1, step):
        w = s[i : i + window]
        gc = 100.0 * (w.count("G") + w.count("C")) / window
        out.append((i + window // 2 + 1, round(gc, 2)))
    return out


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string (handles N and IUPAC ambiguity codes)."""
    comp = str.maketrans(
        {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "R": "Y",
            "Y": "R",
            "S": "S",
            "W": "W",
            "K": "M",
            "M": "K",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
            "a": "t",
            "t": "a",
            "c": "g",
            "g": "c",
            "r": "y",
            "y": "r",
            "s": "s",
            "w": "w",
            "k": "m",
            "m": "k",
            "b": "v",
            "v": "b",
            "d": "h",
            "h": "d",
            "n": "n",
        }
    )
    return seq.translate(comp)[::-1]
