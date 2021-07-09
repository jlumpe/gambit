Command Line Interface
**********************

Root command group
==================

.. program:: gambit

::

   gambit [OPTIONS] COMMAND [ARGS]...


Options
-------

.. option:: -d DB_DIR

    Path to directory containing GAMBIT database files. Must contain exactly one ``.db`` and one
    ``.h5`` file. Required by most subcommands. As an alternative you can specify the database
    location with the :envvar:`GAMBIT_DB_PATH` environment variable.


Environment
-----------

.. envvar:: GAMBIT_DB_PATH

    Alternative to :option:`-d` for specifying path to database.


Commands
========

query
-----

.. program:: gambit query

::

   gambit query [OPTIONS] FILES...


Predict taxonomy of microbial samples from genome sequences.

Files must contain assembled genome sequences, but may have multiple contigs.


.. option:: -o, --output OUTFILE

   File to write output to. If omitted will write to stdout.

.. option:: -s, --seqfmt {fasta}

   Format of genome sequence files. Currently only FASTA is supported.

.. option:: -f, --outfmt {json|csv}

   Output format.


Output Formats
==============

JSON
----

TODO

CSV
---

(Not yet implemented)
