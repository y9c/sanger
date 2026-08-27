[![Readthedocs](https://readthedocs.org/projects/cfutils/badge/?version=latest)](https://cfutils.readthedocs.io/en/latest/?badge=latest)
[![Pypi Releases](https://img.shields.io/pypi/v/cfutils.svg)](https://pypi.python.org/pypi/cfutils)
[![Downloads](https://static.pepy.tech/badge/cfutils)](https://pepy.tech/project/cfutils)

**Chromatogram File Utils**

For Sanger sequencing data visualizing, alignment, mutation calling, and trimming etc.

## Demo

![plot chromatogram with mutation](https://raw.githubusercontent.com/y9c/cfutils/master/data/plot.png)

> command to generate the demo above

```bash
cfutils mut --query ./data/B5-M13R_B07.ab1 --subject ./data/ref.fa --outdir ./data/ --plot
```

## How to use?

- You can have mutation detection and visualization in one step using the command line.

```bash
cfutils mut --help
```

- You can also integrate the result matplotlib figures and use it as a python module.

An example:

```python
import matplotlib.pyplot as plt
import numpy as np

from cfutils.parser import parse_abi
from cfutils.show import plot_chromatograph

seq = parse_abi("./data/B5-M13R_B07.ab1")
peaks = seq.annotations["peak positions"][100:131]

fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
plot_chromatograph(
    seq,
    region=(100, 130),
    ax=axes[0],
    show_bases=True,
    show_positions=True,
    color_map=dict(zip("ATGC", ["C0", "C2", "C1", "C4"])),
)
axes[1].bar(peaks, np.random.randn(len(peaks)), color="0.66")
plt.show()
```

![plot chromatogram in_matplotlib](https://raw.githubusercontent.com/y9c/cfutils/master/data/matplotlib_example.png)

## How to install?

### lightweight core

The default install pulls only `numpy`, `click` and `rich-click` — no C
compiler and no plotting library required:

```bash
pip install cfutils
```

Parsing, QC, alignment (a NumPy Smith-Waterman fallback is built in), analysis
and export all work out of the box.

### optional extras

```bash
pip install "cfutils[plot]"     # matplotlib -> chromatogram figures
pip install "cfutils[align]"    # ssw (C) -> faster alignment
pip install "cfutils[viewer]"   # DNA Features Viewer integration
pip install "cfutils[all]"      # everything
```

The bundled Cython Smith-Waterman accelerator (`cfutils._swalign`) is compiled
automatically when a C compiler is present at build time and used
transparently; otherwise the NumPy fallback is used.

### from source

```bash
git clone git@github.com:y9c/cfutils.git
cd cfutils
make init       # install dependencies
make test       # run the test-suite
```

## ChangeLog

- Reverse completement the chromatogram file. (Inspired by Snapgene)
- build as python package for pypi
- fix bug that highlighting wrong base
- replace blastn with buildin python aligner
- add `features` module: dna-features-viewer style feature overlay API
- add `quality` module: tunable quality-score filtering (was hard-coded in plot)
- add `tracks` module: split / join / slice chromatogram trace files
- add `assembly` module: multi-read pileup + consensus calling
- add `composite` module: side-by-side plotting with external tools
- add `center_region` + fix `highlight_base` end-boundary bug in `show`
- add `basecaller` module: re-call bases from raw four-channel traces
- add `qc` module: per-read QC metrics + batch summaries
- add `transform` module: public Mott `trim` (keeps traces/peaks aligned) + `reverse_complement_record`
- add `utils.normalize_ambiguity` (IUPAC handling)
- fix false-positive mutations from lowercase reference bases (case-insensitive compare)
- add `dnalink` module: plot cfutils chromatograms together with DNA Features Viewer
- expose `parse_abi(rescale=...)`; base/peak analysis requires `rescale=False`
- add `analysis` module: translate, motif/restriction-site scan, sliding GC
- add `export` module: FASTA / VCF / JSON / batch summary
- add heterozygote (mixed-base) calling to `basecaller`
- add `trim_ends` / `strip_primers` to `transform`
- rich-click CLI UX: themed command groups, option groups, short options, `--version`

## Command-line interface

The CLI is built on `rich-click` and organised into themed groups (see
`cfutils --help`):

```text
cfutils mut            mutation calling & reporting
cfutils qc             per-read quality-control metrics
cfutils track          split / join / slice chromatogram trace files
cfutils edit           trim, strip primers, reverse-complement
cfutils basecall       re-call bases from raw four-channel traces
cfutils assemble       reference-guided pileup & consensus
cfutils analyze        sequence biology (translate, motifs, restriction)
cfutils export         FASTA / VCF / JSON / batch summary
cfutils plot           chromatogram rendering (+ features / DNA viewer)
```

Examples:

```bash
cfutils mut -q read.ab1 -s ref.fa -o out --plot        # mutation report + figure
cfutils qc r1.ab1 r2.ab1                                # QC table
cfutils track split read.ab1 -c 20,40 -f tsv            # split traces
cfutils edit trim read.ab1 -c 0.05 -o out               # Mott quality trim
cfutils edit strip-primers read.ab1 -f AAAA -r CCCA     # primer removal
cfutils basecall call read.ab1 -r 0.45                  # re-call bases
cfutils basecall hetero read.ab1                        # mixed/heterozygous sites
cfutils assemble consensus a.ab1 b.ab1 -r ref.fa        # reference-guided consensus
cfutils analyze rest read.ab1                           # restriction sites
cfutils analyze translate read.ab1 -f 1                 # protein translation
cfutils export vcf -q read.ab1 -s ref.fa                # VCF of variants
cfutils export batch *.ab1 -o out -f csv                # batch QC table
cfutils plot dnaviewer read.ab1 --start 50 --end 100    # with DNA Features Viewer
```

## New analysis modules

### split / join trace files
```bash
cfutils track join a.ab1 b.ab1 --outbase joined --outdir out/
cfutils track split a.ab1 --cuts 20,40 --outdir out/
cfutils track slice a.ab1 --start 1 --end 200 --outdir out/
```

### quality filtering (was hard-coded in the plotter)
```python
from cfutils.align import call_mutations
from cfutils.quality import QualityFilter
sites = call_mutations(query, reference)
passed = QualityFilter(min_base_qual=20, min_local_qual=20).filter(sites)
```

### feature overlay (compatible with external tools)
```python
from cfutils.parser import parse_abi
from cfutils.features import ChromatogramFeature, plot_features
from cfutils.show import plot_chromatograph
import matplotlib.pyplot as plt

rec = parse_abi("./data/B5-M13R_B07.ab1")
fig, ax = plt.subplots(figsize=(16, 5))
feat = ChromatogramFeature(start=100, end=131, strand=+1,
                           color="#ff8888", label="M13R primer")
plot_chromatograph(rec, region=(90, 140), ax=ax)
plot_features(rec, ax, features=[feat])
plt.show()
```

### side-by-side with another tool's output
```python
from cfutils.composite import side_by_side
from cfutils.parser import parse_abi

def my_panel(ax, trace_x, peaks, seq, record):
    ax.bar(range(len(seq)), [1.0]*len(seq), color="0.66")
    ax.set_ylabel("my tool signal")

rec = parse_abi("./data/B5-M13R_B07.ab1")
fig, (ax_chrom, ax_panel) = side_by_side(rec, my_panel, region=(10, 40))
```

### consensus / assembly from many reads
```python
from cfutils.parser import parse_abi, parse_fasta
from cfutils.assembly import pileup, consensus
reads = [parse_abi(f) for f in ["a.ab1", "b.ab1", "c.ab1"]]
ref = parse_fasta("ref.fa")
table = pileup(reads, ref, quality_threshold=20)
print(consensus(table))
```

### plot together with DNA Features Viewer
```python
from cfutils.parser import parse_abi
from cfutils.features import ChromatogramFeature
from cfutils.dnalink import plot_combined, to_graphic_record

rec = parse_abi("./data/B5-M13R_B07.ab1")
feats = [ChromatogramFeature(start=90, end=130, strand=+1, label="primer F")]
fig, (ax_feat, ax_chrom) = plot_combined(rec, features=feats, region=(55, 90))
```
Requires the optional extra: `pip install cfutils[viewer]` (i.e. `dna-features-viewer`).

### re-call bases from raw traces
```python
from cfutils.parser import parse_abi
from cfutils.basecaller import call_bases, basecaller_score
rec = parse_abi("./data/B5-M13R_B07.ab1", rescale=False)  # raw traces
res = call_bases(rec)
print(res.sequence)
print(basecaller_score(res, rec.seq))
```

### trim / reverse-complement a whole record
```python
from cfutils.transform import trim, reverse_complement_record
short = trim(rec)                 # Mott quality trim, keeps peaks/trace aligned
rc    = reverse_complement_record(rec)
```

### QC metrics (incl. continuous read length)
```python
from cfutils.qc import read_metrics, continuous_read_length, noise_metric, signal_intensity
m = read_metrics(rec)
print(m["mean_qual"], m["trim_start"], m["trim_end"])
print("CRL:", continuous_read_length(rec))   # longest run with 20-base avg Q >= 20
print("signal:", signal_intensity(rec), "SNR:", noise_metric(rec))
```

### automatic orientation detection
```python
from cfutils.align import detect_orientation
ori = detect_orientation(rec, ref)           # +1 forward, -1 reverse-complement
```

### mixed / heterozygous bases
```python
from cfutils.parser import parse_abi
from cfutils.basecaller import call_bases
res = call_bases(parse_abi("./data/B5-M13R_B07.ab1", rescale=False))
print(res.sequence)
for pos, major, minor, frac in res.heterozygotes(min_ratio=0.2):
    print(pos, major, minor, frac)           # e.g. 4 M C 0.74
```

### export to standard formats
```python
from cfutils.export import to_fasta, to_vcf, to_json, write_batch
print(to_fasta(rec))                                   # FASTA
print(to_vcf(mutations, reference_name="ref"))          # VCF of variants
print(to_json(rec))                                     # self-describing JSON
write_batch([rec1, rec2], "out", fmt="csv")             # batch QC table
```

## Performance & acceleration

cfutils keeps the common path fast without heavy dependencies:

* **Base-calling / trace analysis** — the running-median baseline correction is
  vectorised with NumPy (sliding-window median), giving ~12x over a naive loop.
* **Alignment** — dispatch order is: bundled **Cython Smith-Waterman**
  (`cfutils._swalign`, self-contained) -> **ssw** (if installed) -> **NumPy
  Smith-Waterman**.  All three share the same interface.
* **Parsing** — the ABI reader already unpacks channel/quality arrays with a
  single C `struct.unpack` per record (~7 ms per `.ab1`).

Benchmarks on the bundled 1141 bp sample: parse ~7 ms; full read↔reference
alignment ~15 ms; re-call bases ~35 ms; mutation report + figure ~0.5 s.

## TODO

- [x] call mutation by alignment and plot Chromatogram graphic
- [x] add a doc
- [x] change xaxis by peak location
- [ ] fix bug that chromatogram switch pos after trim (partially addressed via tracks/composite)
- [x] wrap as a cli app
- [x] return quality score in output (quality module + report)
- [x] fix issue that selected base is not in the middle (center_region)
- [ ] fix plot_chromatograph rendering bug (further validation needed)
- [x] add projection feature to make align and assemble possible (assembly module)
- [ ] preserve trimmed-origin positions when slicing/joining so ref coords stay stable
