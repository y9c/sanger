# API Reference

## sanger.run.report_mutation

```python
def report_mutation(
    query_ab1_file: str,
    subject_fasta_file: str,
    output_dir: Optional[str] = None,
    file_basename: Optional[str] = None,
    report_all_sites: bool = False,
    report_mut_plot: bool = False,
) -> None:
    """Report mutation within region, optionally generate plot and TSV output."""
```

- **query_ab1_file**: Path to the ABI file (Sanger sequencing data)
- **subject_fasta_file**: Path to the reference FASTA file
- **output_dir**: Output directory for results (optional)
- **file_basename**: Output file basename (optional)
- **report_all_sites**: If True, report all aligned sites; otherwise, only mutation sites
- **report_mut_plot**: If True, generate a PDF plot of the mutation region

---

## sanger.parser.parse_abi

```python
def parse_abi(path: str) -> SeqRecord:
    """Parse an ABI file and return a SeqRecord object."""
```

---

## sanger.show.plot_chromatograph

```python
def plot_chromatograph(
    seq: SeqRecord,
    region: Optional[Tuple[int, int]] = None,
    ax: Optional[matplotlib.axes.Axes] = None,
    show_bases: bool = True,
    show_positions: bool = True,
    color_map: Optional[dict] = None,
) -> matplotlib.axes.Axes:
    """Plot a chromatogram for a given sequence region."""
```

---

For more details, see the source code and docstrings in each module.
