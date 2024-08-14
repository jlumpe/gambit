# testdb_210818

This directory contains data for the `testdb_210818` test database. Unlike the
previously used `testdb_210126`, all files are small enough to be added in
version control so nothing needs to be downloaded from an external data
repository. To use this database from the CLI, just pass this directory with the
`--db` parameter at the root level.


## Files

* `ref-genomes.gdb`: reference genomes metadata.
* `ref-signatures.gs`: reference genome signatures.
* `ref-genomes.csv`: CSV file of basic reference genome properties (sort of redundant with `ref-genomes.gdb`).
* `ref-genomes/`: contains reference genome files in FASTA format.
* `queries/`
  * `queries.csv`: table listing all query files and expected results.
  * `genomes/`: contains query genome files in FASTA format.
  * `query-signatures.gs`: precalculated signatures for query genomes.
* `results/`: pre-calculated results using query files in `queries`, exported in the "archive" JSON
  format. Two sets of results, one with strict mode enabled and one without. These are used to
  reconsitute the `gambit.query.QueryResults` instances using `gambit.results.ResultsArchiveReader`.
* `generate-results.py`: script which generates result files in `results/`. This will need to be
  re-run if the query results object changes structure or if the "archive" JSON format changes.
  Results are verified against contents of `queries.csv` before exporting.


### Query genome properties

`queries.csv` contains information on expected results for each query genome. This should stay
constant even if the exported files change format in future releases.

Contains the following columns:

- `name`: File name.
- `predicted`: Name of predicted taxon in strict mode, or empty if no prediction.
- `primary`: Description of primary genome match in strict mode, or empty if no prediction.
- `closest`: Description of closest genome match.
- `warnings`: Whether warnings should be generated in strict mode.

In non-strict mode, the primary match will the set to the closest match.
