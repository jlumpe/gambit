# GAMBIT
[![Build Status](https://github.com/hesslab-gambit/gambit/actions/workflows/ci.yml/badge.svg)](https://github.com/hesslab-gambit/gambit/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/hesslab-gambit/badge/?version=latest)](https://hesslab-gambit.readthedocs.io/en/latest/?badge=latest)

GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking) is a tool for rapid taxonomic identification of microbial pathogens.
It uses an extremely efficient genomic distance metric along with a curated database of over 50,000 reference genomes (derived from NCBI [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/))
to identify query genomes within seconds.

Developed by Jared Lumpe in collaboration with the David Hess lab at Santa Clara University. Publication coming soon.

See the [documentation](https://hesslab-gambit.readthedocs.io/en/latest/?badge=stable) for more
detailed information on the tool and how to use it.


## Installation

### Building from source

Install build dependencies:

    pip install setuptools numpy cython

Clone the repository:

    git clone https://github.com/hesslab-gambit/gambit

Build and install:

    cd gambit
    pip install .


### Download reference database

Download the following files and place them in a directory of your choice:

* [gambit-genomes-1.0b1-210719.db](https://storage.googleapis.com/hesslab-gambit-public/databases/refseq-curated/1.0-beta/gambit-genomes-1.0b1-210719.db)
* [gambit-signatures-1.0b1-210719.h5](https://storage.googleapis.com/hesslab-gambit-public/databases/refseq-curated/1.0-beta/gambit-signatures-1.0b1-210719.h5)

(Note - these are marked as "beta" but little is likely to change in the upcoming 1.0 release).


## Usage

    gambit [-d DB_DIR] query [-f {csv|json}] [-o OUTFILE] GENOMES*

* `GENOMES` are one or more files containing query genomes. Currently only FASTA format is supported.
  They must be assembled but may consist of multiple contigs. Support for unassembled raw reads in
  FASTQ format is in development.
* `DB_DIR` is the directory containing the database files. Alternatively you may set the
  `GAMBIT_DB_PATH` environment variable in order to avoid typing this each time.
* `OUTFILE` is the file to write results to.
* `-f` sets the output format (defaults to `csv`).

See the documentation for additional details on the command line interface and description of the output.


## Contact

For questions regarding usage of the software itself, please contact Jared Lumpe at [mjlumpe@gmail.com](mailto:mjlumpe@gmail.com).
All other questions should be directed to David Hess at [dchess@scu.edu](mailto:dchess@scu.edu).
