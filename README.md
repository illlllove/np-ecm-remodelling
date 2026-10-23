# Nucleus pulposus ECM remodelling: figure reproduction

This repository contains the Python code and processed inputs used to redraw seven main figures and 24 supplementary figures for a study of fibrocartilage-like extracellular matrix remodelling in human intervertebral discs.

The study distinguishes the relative abundance of remodelling-state cells, expression within those cells, and the signal measured in bulk tissue. These are different biological quantities: the presence of a reproducible ECM programme does not imply a uniform increase in its cell fraction or expression across degeneration cohorts.

## Scope

This is a figure-reproduction repository. It contains frozen donor-level summaries, marker-expression tables, UMAP coordinates, palettes and figure legends. It does not contain the original raw expression matrices or the complete pipelines for quality control, clustering, statistical inference and upstream multi-omics analyses. Running the scripts redraws existing results; it does not refit the original analyses.

The supplied inputs reproduce the canonical V6 main and supplementary figures. No result-dependent gene selection, label transfer, cell filtering or new statistical testing is performed by these drawing scripts.

## Contents

- `Figure 1.py` through `Figure 7.py`: main figures.
- `Supplementary Figure S1.py` through `Supplementary Figure S24.py`: supplementary figures.
- `Inputs/`: 61 frozen plotting inputs.
- `Inputs.json`: input checksums and source-relative provenance.
- `requirements.txt` and `Environment.json`: exact Python dependency versions.
- `Run figures.py`: run one, several or all figures.
- `SHA256.txt`: checksums for this repository snapshot.

## Installation and use

The reference environment is Python 3.12.14 on Windows x86-64. Install the locked dependencies in an isolated Python environment:

```text
python -m pip install -r requirements.txt
```

Run commands from the repository directory:

```text
python "Run figures.py" M1
python "Run figures.py" M2 M4 M6
python "Run figures.py" S23
python "Run figures.py" all
```

An individual file can also be run directly:

```text
python "Figure 1.py"
```

Each figure file contains its own drawing functions. It reads only its packaged plotting inputs and writes PDF, 300 dpi PNG and temporary files to `Outputs/`. Output and cache files are not tracked by Git.

PDFium is used to render PNG files. Raster antialiasing can differ from other PDF readers; the underlying PDF drawing and numerical results are unchanged.

## Portable Windows archive

The complete distribution additionally contains a Windows x86-64 Python runtime and its installed dependencies. That large archive is distributed separately rather than committed to this repository. Once extracted, `Run figures.cmd` in the complete distribution runs without requiring the original project directory or a separate Python installation.

The repository is maintained at https://github.com/illlllove/np-ecm-remodelling.
The archival DOI will be added after the Zenodo record is registered.

## Data provenance and reuse

Source-relative paths and SHA-256 checksums are recorded in `Inputs.json`. These files are processed plotting inputs derived from public research datasets and completed analyses, not newly collected individual-level clinical data. Raw-data accessions and cohort roles should be consulted in the associated manuscript.

Original project code and documentation are licensed under MIT; see `LICENSE`.
Data-source conditions and third-party software licenses are separate and are not replaced by the project code license; see `Reuse.txt`.

## Citation

Authors, in order: Liang Wenyue; Yao Hongping; Tu Wang; Zhai Yu.
Version: 1.0.0. See `CITATION.cff` for machine-readable citation metadata.
The Zenodo DOI will be added once the archived record has been published.
