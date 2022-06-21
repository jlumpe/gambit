.. _install-page:

Installation and Setup
**********************


Install from Bioconda
=====================

The recommended way to install the tool is through the conda package manager (available
`here <https://docs.conda.io/en/latest/miniconda.html>`_)::

    conda install -c bioconda gambit


Install from source
===================

Installing from source requires the ``cython`` package as well as a C compiler be installed on your
system. Clone the repository and navigate to the directory, and then run::

    pip install .

Or do an editable development install with::

    pip install -e .


Database files
==============

Download the following files and place them in a directory of your choice:

* `gambit-genomes-1.0b1-210719.db <https://storage.googleapis.com/hesslab-gambit-public/databases/refseq-curated/1.0-beta/gambit-genomes-1.0b1-210719.db>`_
* `gambit-signatures-1.0b1-210719.h5 <https://storage.googleapis.com/hesslab-gambit-public/databases/refseq-curated/1.0-beta/gambit-signatures-1.0b1-210719.h5>`_
