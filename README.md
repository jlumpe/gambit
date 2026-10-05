# GAMBIT

[![Build Status](https://github.com/jlumpe/gambit/actions/workflows/ci.yml/badge.svg)](https://github.com/jlumpe/gambit/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/gambit-genomics/badge/?version=latest)](https://gambit-genomics.readthedocs.io/en/latest/?badge=latest)
<!-- [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/gambit/README.html) -->
[![Bioconda](https://img.shields.io/conda/vn/bioconda/gambit)](https://bioconda.github.io/recipes/gambit/README.html)
[![PyPI](https://img.shields.io/pypi/v/gambit)](https://pypi.org/project/gambit/)
[![Python versions](https://img.shields.io/pypi/pyversions/gambit)](https://pypi.org/project/gambit/)

GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking) is a tool for rapid taxonomic identification of microbial pathogens.
It uses an efficient genomic distance metric along with a curated database of approximately 50,000 reference genomes (derived from NCBI
[RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)) to identify genome assemblies from across the Bacterial kingdom in seconds.

See below for basic installation and usage instructions, or check out the
[documentation](https://gambit-genomics.readthedocs.io/en/latest) for more detailed information and
a basic tutorial.

## About

Copyright © 2016-2026 Jared Lumpe

GAMBIT has been a personal project of mine for many years. Although there have been numerous
contributors to the publication, it is not a product of any lab or institution.

GAMBIT is provided as free software under the terms of the [AGPLv3 license](https://github.com/jlumpe/gambit/blob/main/LICENSE).
It is not covered by any type of software patent.

### Publication

Lumpe J, Gumbleton L, Gorzalski A, Libuit K, Varghese V, et al. (2023) GAMBIT (Genomic Approximation
Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome
sequencing of bacterial isolates for clinical identification. PLOS ONE 18(2): e0277575.
[https://doi.org/10.1371/journal.pone.0277575](https://doi.org/10.1371/journal.pone.0277575)

See [jlumpe/gambit-publication](https://github.com/jlumpe/gambit-publication) for a reproducible
workflow to generate all analyses and figures in the paper.

### Contact

Please contact Jared Lumpe at [jared@jaredlumpe.com](mailto:jared@jaredlumpe.com) with any questions or feedback.

## Installation



### Conda

Install from Bioconda:

```bash
conda install -c conda-forge -c bioconda gambit
```



### Pixi

Install the `gambit` command globally using
[Pixi](https://pixi.prefix.dev/latest/global_tools/introduction/):

```bash
pixi global install -c conda-forge -c bioconda gambit
```

To add GAMBIT to an existing Pixi workspace instead, use `pixi add gambit` (after ensuring the
`conda-forge` and `bioconda` channels are added to the workspace).

### Pip

```bash
pip install gambit
```

Pre-built wheels are only provided for Linux (x86_64) and CPython 3.9-3.14. On other platforms pip
will attempt to build from the source distribution (see next section).

### From source

```bash
git clone https://github.com/jlumpe/gambit.git
cd gambit
pip install .  # or "pip install -e ." for an editable install
```

Requires a C compiler with OpenMP support. Not supported on macOS (Apple clang lacks `-fopenmp`) or
Windows. macOS users should install using Conda or Pixi instead, and Windows users should use
[WSL](https://learn.microsoft.com/en-us/windows/wsl/).

## Database files

GAMBIT requires a reference database for taxonomic identification. Download the database files and
place them in a directory of your choice:

- [gambit-refseq-curated-1.0.gdb](https://storage.googleapis.com/jlumpe-gambit/public/databases/refseq-curated/1.0/gambit-refseq-curated-1.0.gdb)
- [gambit-refseq-curated-1.0.gs](https://storage.googleapis.com/jlumpe-gambit/public/databases/refseq-curated/1.0/gambit-refseq-curated-1.0.gs)



## Basic usage

```
gambit [-d /path/to/database/] query [-o results.csv] genome1.fasta genome2.fasta ...
```

Positional arguments are one or more FASTA files containing query genome assemblies. You must
provide the path to the directory containing the database files using either the `-d` option
(*before* the `query` subcommand) or by setting the `GAMBIT_DB_PATH` environment variable.

See the documentation for additional details on the command line interface and description of the output.
