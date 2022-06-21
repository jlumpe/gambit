Tutorial
********


Before starting, make sure you have followed the instructions in :ref:`install-page`.

Also see the :ref:`cli-page` page for complete documentation of all GAMBIT subcommands
and options.


Telling GAMBIT where the database files are
===========================================

Classification of unknown genomes requires the GAMBIT reference database files. You can let GAMBIT
know which directory contains these files in one of two ways:

Via the command line
--------------------

The first is to explicitly pass it via the command line using the ``--db`` option, like so::

    gambit --db /path/to/database/ COMMAND ...

Note that this option must appear immediately after ``gambit`` and before the command name.

Using an environment variable
-----------------------------

The second is to use the ``GAMBIT_DB_PATH`` environment variable. This can be done by
running the following command at the beginning of your shell session::

    export GAMBIT_DB_PATH="/path/to/database/"

Alternatively, you can add this line to your ``.bashrc`` to have it apply to all future sessions
(make sure to restart your current session after doing so).


Genome input
============

Genome assemblies used as input must be in FASTA format, optionally compressed with gzip.

Most commands accept a list of genome files as positional arguments, e.g.::

    gambit COMMAND [OPTIONS] genome1.fasta genome2.fasta ...

Alternatively you can use the ``-l`` option to provide a text file containing the genome file
names/paths, one per line. The paths in this file are considered relative to the directory given by
the ``--ldir`` option, or the current working directory if not given. Note that the ``gambit dist``
command has different names for these options because there are two lists of genomes to specify.
See the :ref:`cli-page` page for more complete information.


Predicting taxonomy of unknown genomes
======================================

TODO


Calculating GAMBIT distances
============================

TODO


Pre-computing k-mer signatures
==============================

TODO


Creating relatedness trees
==========================

TODO

