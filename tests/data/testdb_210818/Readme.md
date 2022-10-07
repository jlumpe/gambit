# testdb_210818

This directory contains data for the `testdb_210818` test database. Unlike the
previously used `testdb_210126`, all files are small enough to be added in
version control so nothing needs to be downloaded from an external data
repository. To use this database from the CLI, just pass this directory with the
`--db` parameter at the root level.


## Files

* `ref-genomes.gdb` - reference genomes metadata.
* `ref-signatures.gs` - reference genome signatures.
* `queries/`
  * `queries.csv` - table listing all query files and expected results.
  * `genomes/` - contains query genome files in FASTA format.
  * `query-signatures.gs` - precalculated signatures for query genomes.
* `results/` - pre-calculated results using query files in `queries`.
* `generate-results.py` - script which generates result files in `results/`.
  Verifies against expected result attributes in `queries.csv`.
