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

Download files for the latest database release from the :ref:`Database Releases` page and place them
in a directory of your choice. The directory should not contain any other files with the same
extension.
