.. _install-page:

Installation and Setup
**********************

Python package and command-line tool
====================================

Conda
-----

Install with the `Conda`_ package manager from the `Bioconda`_ channel::

    conda install -c conda-forge -c bioconda gambit

.. _Conda: https://docs.conda.io/
.. _Bioconda: https://bioconda.github.io/


Pixi
----

Install the ``gambit`` command globally using `Pixi`_::

    pixi global install -c conda-forge -c bioconda gambit

To add GAMBIT to an existing Pixi workspace instead, use ``pixi add gambit`` (after ensuring the
``conda-forge`` and ``bioconda`` channels are added to the workspace).

.. _Pixi: https://pixi.prefix.dev/latest/global_tools/introduction/


Pip
---

Install from `PyPI`_::

    pip install gambit

Pre-built wheels are only provided for Linux (x86_64) and CPython 3.12-3.14. On other platforms pip
will attempt to build from the source distribution (see :ref:`install-source`).

.. _PyPI: https://pypi.org/project/gambit/


.. _install-source:

From source
-----------

Clone the repository and install with pip::

    git clone https://github.com/jlumpe/gambit.git
    cd gambit
    pip install .  # or "pip install -e ." for an editable install

Requires a C compiler with OpenMP support. Not supported on macOS (Apple clang lacks ``-fopenmp``)
or Windows. macOS users should install using Conda or Pixi instead, and Windows users should use
`WSL`_.

.. _WSL: https://learn.microsoft.com/en-us/windows/wsl/


.. _install-db:

Database files
==============

You will need a GAMBIT reference database (consisting of one ``.gdb`` and one ``.gs`` file) to
perform taxonomic classification using the ``gambit query`` command. Download files for the latest
database release from the :ref:`Database Releases` page and place them in a directory of your
choice. The directory should not contain any other files with the same extensions.
