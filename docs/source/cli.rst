.. _cli-page:

Command Line Interface
**********************


Root command group
==================

.. program:: gambit

::

   gambit [OPTIONS] COMMAND [ARGS]...

Some top-level options are set at the root command group, and should be specified `before` the name
of the subcommand to run.

Options
-------

.. option:: -d, --db DIR

    Path to the directory containing reference database files. Required by most subcommands.
    As an alternative you can specify the database location with the :envvar:`GAMBIT_DB_PATH`
    environment variable.

Environment variables
---------------------

.. envvar:: GAMBIT_DB_PATH

    Alternative to :option:`-d` for specifying path to database.


Querying the database
=====================

query
-----

.. program:: gambit query

::

    gambit query [OPTIONS] (-s SIGFILE | -l LIST | GENOMES...)


Predict taxonomy of microbial samples from genome sequences.

``GENOMES`` are one or more FASTA files containing assembled query genomes. Alternatively

a file containing pre-calculated signatures may be used with the ``--sigfile`` option. The
reference database must be specified from the root command group.


Options
.......


.. option:: -s, --sigfile FILE

    Path to file containing query signatures.

.. option:: -o, --output FILE

   File to write output to. If omitted will write to stdout.

.. option:: -f, --outfmt {csv|json|archive}

   Results format (see next section).


Result Formats
--------------

CSV
...

A .csv file with one row per query. Contains the following columns:

* ``query.name`` - Name of query.
* ``query.path`` - Path to query file, if any.
* ``predicted.name`` - Name of predicted taxon.
* ``predicted.rank`` - Rank of predicted taxon.
* ``predicted.ncbi_id`` - ID of predicted taxon in NCBI taxonomy database.
* ``predicted.threshold`` - Distance threshold of predicted taxon.
* ``closest.distance`` - Distance to closest genome.
* ``closest.description`` - Description of closest genome.


JSON
....

A machine-readable format meant to be used in pipelines.

.. todo::
   Document schema


Archive
.......

A more verbose JSON-based format used for testing and development.


Generating and inspecting k-mer signatures
==========================================

signatures info
---------------

.. program:: gambit signatures info

::

   gambit signatures info [OPTIONS] FILE


Print information about a GAMBIT signatures file. Defaults to a basic human-readable format.


Options
.......

.. option:: -j, --json

   Print information in JSON format. Includes more information than standard output.

.. option:: -p, --pretty

   Prettify JSON output to make it more human-readable.

.. option:: -i, --ids

   Print IDs of all signatures in file.


signatures create
-----------------

.. program:: gambit signatures create

::

   gambit signatures create [OPTIONS] GENOMES

Calculate GAMBIT signatures of ``GENOMES`` and write to file.

The ``-k`` and ``--prefix`` options may be omitted if a reference database is specified through the
root command group, in which case the parameters of the database will be used.


Options
.......

.. option:: -o, --output FILE

   Path to write file to (required).

.. option:: -k INTEGER

   Length of k-mers to find (does not include length of prefix).

.. option:: -p, --prefix STRING

   K-mer prefix to match, a non-empty string of DNA nucleotide codes.

.. option:: -i, --ids FILE

   File containing IDs to assign to signatures in file metadata. Should contain one ID per line.

.. option:: -m, --meta-json FILE

   JSON file containing metadata to attach to file.

   .. todo::
      Document schema
