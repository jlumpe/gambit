.. _install-page:

Installation and Setup
**********************

Python package and command-line tool
====================================

Install from Bioconda
---------------------

The recommended way to install the tool is through the `Conda`_ package manager (I recommend the
`Miniconda`_ distribution) from the `Bioconda`_ channel::

    conda install -c bioconda gambit


.. _Conda: https://www.anaconda.com/products/distribution
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _Bioconda: https://bioconda.github.io/


Install from source
-------------------

Installing from source requires the ``cython`` package as well as a C compiler be installed on your
system. Clone the repository and navigate to the directory, and then run::

    pip install .

Or do an editable development install with::

    pip install -e .


.. _install-db:

Database files
==============

You will need a GAMBIT reference database (consisting of one ``.gdb`` and one ``.gs`` file) to
perform taxonomic classification using the ``gambit query`` command. Download files for the latest
database release from the :ref:`Database Releases` page and place them in a directory of your
choice. The directory should not contain any other files with the same extensions.
