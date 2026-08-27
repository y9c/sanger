[![Pypi Releases](https://img.shields.io/pypi/v/sanger.svg)](https://pypi.python.org/pypi/sanger)
[![Downloads](https://static.pepy.tech/badge/sanger)](https://pepy.tech/project/sanger)

# sanger

A lightweight toolkit for **Sanger sequencing** data: chromatogram
visualization, alignment & mutation calling, quality control, base-calling,
trimming, assembly and export — usable as a CLI, a Python library, or an
MCP server for LLM agents.

## Quick start

```bash
pip install sanger
sanger mut -q read.ab1 -s ref.fa -o out --plot
```

```python
from sanger import Chromatogram, parse_fasta

cg = Chromatogram.from_abi("./data/B5-M13R_B07.ab1")
print(cg.length, cg.mean_quality, cg.gc_percent)   # 1141 50.2 52.0
ref = parse_fasta("./data/ref.fa")
print(cg.qc())                                    # QC metrics (CRL, SNR)
print(cg.to_vcf(ref))                             # variant calling -> VCF
cg.plot(region=(55, 90))                          # render a region
```

## Examples

The gallery below is produced from the bundled real ABI sample
(`data/B5-M13R_B07.ab1` vs `data/ref.fa`) by `python -m scripts.make_readme_examples`.

**Mutation calling** — the SNP `T61A` is highlighted and annotated.

![Mutation calling](https://raw.githubusercontent.com/y9c/sanger/master/examples/mutation_call.png)

**Quality control** — per-base Phred quality, CRL and Mott trimming.

![Quality profile](https://raw.githubusercontent.com/y9c/sanger/master/examples/quality_profile.png)

**Side-by-side panels** — chromatogram + GC% + quality on a shared axis.

![Side-by-side](https://raw.githubusercontent.com/y9c/sanger/master/examples/side_by_side.png)

**Feature overlay** — primers, amplicon and SNPs on the trace.

![Feature overlay](https://raw.githubusercontent.com/y9c/sanger/master/examples/feature_overlay.png)

**Re-called bases** — mixed/heterozygous sites are marked (M/W/K).

![Re-called bases](https://raw.githubusercontent.com/y9c/sanger/master/examples/basecall_hetero.png)

**Assembly** — pileup depth and consensus against a reference.

![Assembly](https://raw.githubusercontent.com/y9c/sanger/master/examples/assembly.png)

## Installation

The default install pulls only `numpy`, `click` and `rich-click` — no C
compiler and no plotting library required. Parsing, QC, alignment (a bundled
Cython Smith-Waterman with a NumPy fallback), analysis and export all work out
of the box.

```bash
pip install sanger
```

Optional extras:

```bash
pip install "sanger[plot]"     # matplotlib -> chromatogram figures
pip install "sanger[viewer]"   # DNA Features Viewer integration
pip install "sanger[agent]"    # MCP server for LLM agents
pip install "sanger[all]"      # everything
```

The bundled Cython Smith-Waterman accelerator (`sanger._swalign`,
self-contained, no `ssw` dependency) is compiled automatically when a C
compiler is present at build time; otherwise the NumPy fallback is used.

From source:

```bash
git clone git@github.com:y9c/sanger.git
cd sanger
make init       # install dependencies
make test       # run the test-suite
```

## Command-line interface

Built on `rich-click`, with themed command groups (`sanger --help`):

```text
sanger mut            mutation calling & reporting
sanger qc             per-read quality-control metrics
sanger track          split / join / slice chromatogram trace files
sanger edit           trim, strip primers, reverse-complement
sanger basecall       re-call bases from raw four-channel traces
sanger assemble       reference-guided pileup & consensus
sanger analyze        sequence biology (translate, motifs, restriction)
sanger export         FASTA / VCF / JSON / batch summary
sanger plot           chromatogram rendering (+ features / DNA viewer)
```

Examples:

```bash
sanger mut -q read.ab1 -s ref.fa -o out --plot        # mutation report + figure
sanger qc r1.ab1 r2.ab1                               # QC table
sanger track split read.ab1 -c 20,40 -f tsv           # split traces
sanger edit trim read.ab1 -c 0.05 -o out              # Mott quality trim
sanger edit strip-primers read.ab1 -f AAAA -r CCCA    # primer removal
sanger basecall call read.ab1 -r 0.45                 # re-call bases
sanger basecall hetero read.ab1                       # mixed/heterozygous sites
sanger assemble consensus a.ab1 b.ab1 -r ref.fa       # reference-guided consensus
sanger analyze rest read.ab1                          # restriction sites
sanger analyze translate read.ab1 -f 1                # protein translation
sanger export vcf -q read.ab1 -s ref.fa               # VCF of variants
sanger export batch *.ab1 -o out -f csv               # batch QC table
sanger plot dnaviewer read.ab1 --start 50 --end 100   # with DNA Features Viewer
```

## Python API

The high-level [`Chromatogram`](#chromatogram-object) object is the easiest way
to work with the toolkit; the low-level modules remain available for custom work.

```python
from sanger import Chromatogram, parse_fasta

cg = Chromatogram.from_abi("./data/B5-M13R_B07.ab1")
ref = parse_fasta("./data/ref.fa")
```

<details>
<summary>Common operations</summary>

**Mutation calling & quality filtering**

```python
from sanger.quality import QualityFilter

snps = cg.call_mutations(ref)
confident = QualityFilter(min_base_qual=20, min_local_qual=20).filter(snps)
print([f"{s.ref_base}{s.ref_pos}{s.cf_base}" for s in confident])
```

**QC metrics (incl. continuous read length)**

```python
m = cg.qc()
print(m["mean_qual"], m["trim_start"], m["trim_end"], m["crl"], m["snr"])
```

**Re-call bases from raw traces + mixed/heterozygous sites**

```python
res = cg.basecall()
print(res.sequence)
for pos, major, minor, frac in res.heterozygotes(min_ratio=0.2):
    print(pos, major, minor, frac)
```

**Trim / reverse-complement / orientation**

```python
trimmed = cg.trim().trim_leading_ns()     # Mott trim then drop leading Ns
rc      = cg.reverse_complement()
ori     = cg.detect_orientation(ref)      # +1 forward, -1 reverse-complement
```

**Sequence-level analysis**

```python
print(cg.analyze("translate", frame=1))            # protein translation
print(cg.analyze("restriction"))                    # {'EcoRI': [27], ...}
print(cg.analyze("motif", motif="AATT"))            # motif positions
```

**Export**

```python
print(cg.to_fasta())
print(cg.to_vcf(ref))
cg.export("out")                                    # write FASTA to disk
```

**Feature overlay (for external tools)**

```python
from sanger import ChromatogramFeature
from sanger.features import plot_features

feat = ChromatogramFeature(start=90, end=130, strand=+1, label="primer F")
fig, ax = cg.plot(region=(80, 140))
plot_features(cg.to_record, ax, features=[feat])
```

**Side-by-side with another tool's output (shared x-axis)**

```python
from sanger.composite import side_by_side

def my_panel(ax, trace_x, peaks, seq, record, start=None):
    ax.bar(range(len(seq)), [1.0] * len(seq), color="0.66")
    ax.set_ylabel("my tool signal")

fig, (ax_chrom, ax_panel) = side_by_side(cg.to_record, my_panel, region=(10, 40))
```

**Consensus / assembly from many reads**

```python
from sanger import parse_abi
from sanger.assembly import pileup, consensus

reads = [parse_abi(f) for f in ["a.ab1", "b.ab1", "c.ab1"]]
table = pileup(reads, ref, quality_threshold=20)
print(consensus(table))
```

**Plot together with DNA Features Viewer** (needs `sanger[viewer]`)

```python
from sanger import ChromatogramFeature
from sanger.dnalink import plot_combined

feats = [ChromatogramFeature(start=90, end=130, strand=+1, label="primer F")]
fig, (ax_feat, ax_chrom) = plot_combined(cg.to_record, features=feats, region=(55, 90))
```

</details>

## Chromatogram object

A terse, idiomatic workflow in one object:

```python
from sanger import Chromatogram

cg = Chromatogram.from_abi("./data/B5-M13R_B07.ab1")
cg.length, cg.mean_quality, cg.gc_percent, cg.channels   # 1141 50.2 52.0 GATC
cg.qc()                            # QC metrics (CRL, SNR)
cg.basecall()                      # re-call bases from raw traces
cg.call_mutations(ref)             # variant calling
cg.trim().trim_leading_ns()        # quality trimming
cg.reverse_complement()            # reverse-strand view
cg.analyze("restriction")          # restriction-site scan
cg.plot(region=(55, 90))           # render a region
cg.to_fasta(), cg.to_vcf(ref)      # export
cg.export("out")                   # write to disk
```

## Agent / MCP

sanger ships a [Model Context Protocol](https://modelcontextprotocol.io) server
so LLM agents and MCP clients can call the toolkit as tools:

```bash
pip install "sanger[agent]"     # adds mcp>=2
sanger-mcp                        # run the MCP server over stdio
python -m sanger.mcp_server       # identical
```

| Tool | Purpose |
|---|---|
| `read_chromatogram(path)` | parse an ABI → summary (length, GC%, quality, CRL) |
| `qc_metrics(path)` | full per-read QC metrics (incl. CRL, signal, SNR) |
| `call_mutations(query ab1, subject fa)` | variants vs a reference (SNPs/indels) |
| `re_call_bases(path)` | re-call bases from raw traces + heterozygotes |
| `analyze_sequence(path, kind)` | translate / motifs / restriction / GC |
| `trim_read(path, mode)` | quality-trim a read |
| `export_sequence(path, outdir)` | write FASTA or VCF |
| `plot_chromatogram(path, out, start, end)` | render a chromatogram PNG |

Register it in an MCP client's config, e.g.:

```json
{ "mcpServers": { "sanger": { "command": "sanger-mcp" } } }
```

## ChangeLog

- Replace the external `ssw` aligner with a bundled Cython Smith-Waterman
  (+ NumPy fallback) — no external alignment dependency.
- Reverse-complement the chromatogram file (inspired by Snapgene).
- Add per-base quality filtering, continuous read length (CRL), signal/SNR QC.
- Add base-calling from raw four-channel traces + heterozygote (mixed-base) calling.
- Add feature overlay API and DNA Features Viewer integration.
- Add split / join / slice trace operations with provenance (`offset`).
- Add reference-guided pileup & consensus (assembly module).
- Add FASTA / VCF / JSON / batch export.
- Fix false-positive mutations from lowercase reference bases (case-insensitive).
- Build the CLI with rich-click (themed groups, short options, `--version`).

## TODO

- [x] call mutation by alignment and plot Chromatogram graphic
- [x] add a doc
- [x] change x-axis by peak location
- [x] fix bug that chromatogram switches position after trim
- [x] wrap as a CLI app
- [x] return quality score in output (quality module + report)
- [x] fix selected base not being centred (center_region)
- [x] fix `plot_chromatograph` rendering bug (full-region bounds, channels)
- [x] add projection/assembly (assembly module)
- [x] preserve trimmed-origin positions when slicing/joining (offset provenance)
