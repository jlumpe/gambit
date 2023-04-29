.. _cli-page:

Command Line Interface
**********************

Genome assembly files accepted by the CLI must be in FASTA format, optionally compressed with gzip.


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

    Path to the directory containing reference database files. Required by the :ref:`query <query-cmd>`
    subcommand.
    As an alternative you can specify the database location with the :envvar:`GAMBIT_DB_PATH`
    environment variable.

Environment variables
---------------------

.. envvar:: GAMBIT_DB_PATH

    Alternative to :option:`-d` for specifying path to database.


Querying the database
=====================

.. _query-cmd:

"query" command
---------------

.. program:: gambit query

::

    gambit query [OPTIONS] (-s SIGFILE | -l LISTFILE | GENOMES...)

Predict taxonomy of microbial samples from genome sequences.

The reference database must be specified from the root command group.

Query genomes
.............

Query genomes can be specified using one of the following methods:

* Give paths of one or more genome files as positional arguments.
* Use the ``-l`` option to specify a text file containing paths of the genome files.
* Use the ``-s`` option to use a signatures file created with the
  `signatures create <signatures-create-cmd_>`_ command.

.. option:: -l LISTFILE

   File containing paths to genomes, one per line.

.. option:: --ldir DIRECTORY

   Parent directory of paths in file given by ``-l`` option.

.. option:: -s, --sigfile FILE

   A genome signatures file.


Additional Options
..................

.. option:: -o, --output FILE

   File to write output to. If omitted will write to stdout.

.. option:: -f, --outfmt {csv|json|archive}

   Results format (see next section).

.. option:: --progress / --no-progress

   Show/don't show progress meter.

.. option:: -c, --cores INT

   Number of CPU cores to use.


.. _query-result-formats:

Result Formats
--------------

CSV
...

A .csv file with one row per query. Contains the following columns:

- ``query`` - Query genome file name (minus extension).
- ``predicted`` - Predicted taxon.

  - ``predicted.name`` - Name of taxon.
  - ``predicted.rank`` - Taxonomic rank (genus, species, etc.).
  - ``predicted.ncbi_id`` - Numeric ID in NCBI taxonomy database, if any.
  - ``predicted.threshold`` - Classification threshold.

- ``closest`` - Reference genome closest to query.

  - ``closest.distance`` - Distance to closest genome.
  - ``closest.decription`` - Text description.

- ``next`` - Next most specific taxon for which the classification threshold was not met.

  - ``next.name``
  - ``next.rank``
  - ``next.ncbi_id``
  - ``next.threshold``


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

.. _signatures-info-cmd:

"signatures info" command
-------------------------

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


.. _signatures-create-cmd:

"signatures create" command
---------------------------

.. program:: gambit signatures create

::

   gambit signatures create [OPTIONS] -o OUTFILE (-l LISTFILE | GENOMES...)

Calculate GAMBIT signatures of a set of genomes and write to a binary file.


Input/output
............

.. option:: -l LISTFILE

   File containing paths to genomes, one per line.

.. option:: --ldir DIRECTORY

   Parent directory of paths in file given by ``-l`` option.

.. option:: -o, --output FILE

   Path to write file to (required).

K-mer parameters
................

.. option:: -k INTEGER

   Length of k-mers to find (does not include length of prefix). Default is 11.

.. option:: -p, --prefix STRING

   K-mer prefix to match, a non-empty string of DNA nucleotide codes. Default is ATGAC.

Metadata
........

.. option:: -i, --ids FILE

   File containing IDs to assign to signatures in file metadata. Should contain one ID per line.
   If omitted will use file names stripped of extensions.

.. option:: -m, --meta-json FILE

   JSON file containing metadata to attach to file.

   .. todo::
      Document metadata schema

Additional Options
..................

.. option:: --progress / --no-progress

   Show/don't show progress meter.

.. option:: -c, --cores INT

   Number of CPU cores to use.


Calculating genomic distances
=============================

"dist" command
--------------

.. program:: gambit dist

::

   gambit dist [OPTIONS] -o OUTFILE
       (-q GENOME... | --ql LISTFILE | --qs SIGFILE)
       (-r GENOME... | --rl LISTFILE | --rs SIGFILE | --square | --use-db)

Calculate pairwise distances between a set of query genomes and a set of reference genomes.
Output is a .csv file. If using ``--qs`` along with ``--rs`` or ``-use-db``, the k-mer parameters
of the query signature file must match the reference parameters.

Query genomes
.............

.. option:: -q GENOME

    Path to a single genome file. May be used multiple times.

.. option:: --ql LISTFILE

   File containing paths of genome files, one per line.

.. option:: --qdir DIRECTORY

   Parent directory of paths in file given by ``--ql`` option.

.. option:: --qs SIGFILE

  A genome signatures file.

Reference genomes
.................

.. option:: -r GENOME

    Path to a single genome file. May be used multiple times.

.. option:: --rl LISTFILE

   File containing paths of genome files, one per line.

.. option:: --rdir DIRECTORY

   Parent directory of paths in file given by ``--rl`` option.

.. option:: --rs SIGFILE

  A genome signatures file.

.. option:: -s, --square

   Use same genomes as the query.

.. option:: -d, --use-db

   Use all genomes in reference database.

Output
......

.. option:: -o FILE

   File to write output to. Required.

K-mer parameters
................

Only allowed if query and reference genomes do not come from precomputed signature files.

.. option:: -k INTEGER

   Length of k-mers to find (does not include length of prefix). Default is 11.

.. option:: -p, --prefix STRING

   K-mer prefix to match, a non-empty string of DNA nucleotide codes. Default is ATGAC.

Additional options
..................

.. option:: --progress / --no-progress

   Show/don't show progress meter.

.. option:: -c, --cores INT

   Number of CPU cores to use.


Creating relatedness trees
==========================

"gambit tree" command
---------------------

.. program:: gambit tree

::

   gambit tree [OPTIONS] (-l LISTFILE | -s SIGFILE | GENOMES...)

Estimate a relatedness tree for a set of genomes and output in Newick format.

Input/output
............

.. option:: -l LISTFILE

   File containing paths of genome files, one per line.

.. option:: --ldir DIRECTORY

   Parent directory of paths in file given by ``-l`` option.

.. option:: -s, --sigfile SIGFILE

  A genome signatures file.

.. option:: -o FILE

   File to write output to. If omitted will write to stdout.

.. todo::

   Allow using a distance matrix calculated using ``gambit dist``.

K-mer parameters
................

Not allowed if the ``-s/--sigfile`` option was used.

.. option:: -k INTEGER

   Length of k-mers to find (does not include length of prefix). Default is 11.

.. option:: -p, --prefix STRING

   K-mer prefix to match, a non-empty string of DNA nucleotide codes. Default is ATGAC.

Additional options
..................

.. option:: --progress / --no-progress

   Show/don't show progress meter.

.. option:: -c, --cores INT

   Number of CPU cores to use.
