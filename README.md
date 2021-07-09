# GAMBIT
[![Build Status](https://github.com/hesslab-gambit/gambit/actions/workflows/ci.yml/badge.svg)](https://github.com/hesslab-gambit/gambit/actions/workflows/ci.yml)

GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking) is a tool for rapid taxonomic identification of microbial pathogens.
It uses an extremely efficient genomic distance metric along with a curated database of over 50,000 reference genomes (derived from NCBI [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/))
to identify query genomes within seconds.

Developed by Jared Lumpe in collaboration with the David Hess lab at Santa Clara University.


## Installation

### Building from source

Install build dependencies:

    pip install setuptools numpy cython

Clone the repository:

    git clone https://github.com/hesslab-gambit/gambit

Build and install:

    cd gambit
    pip install .


### Database files

TODO


## Usage

    gambit [OPTIONS] query [--csv | --json] [-o OUTPUT] GENOMES*

Query genomes must be assembled but may consist of multiple contigs. Currently only FASTA format is supported.
Support for unassembled raw reads (FASTQ format) is in development.


## Contact

For questions regarding usage of the software itself, please contact Jared Lumpe at [mjlumpe@gmail.com](mailto:mjlumpe@gmail.com).
All other questions should be directed to David Hess at [dchess@scu.edu](mailto:dchess@scu.edu).
