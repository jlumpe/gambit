# Changelog


## 1.0.1

* Misc
  * Better error reporting when database files cannot be found (in CLI and API).
  * Minor documentation updates.


## 1.0.0

* New features
	* `tree` command for generating hierarchical clustering trees from distance matrices.
* General
	* Preferred extensions for genome database files and signatures files have been changed from
	  `.db` and `.h5` to `.gdb` and `.gs`.
* Performance improvements
	* Use process-based parallelism by default for parsing multiple sequence files (much faster).
	* Speed up `gambit dist` with `-s` option applied.
* CLI
	* Strip directory and extension from input file IDs. This applies to CSV output for querying and
	  distance calculation and IDs in generated signature files.
	* `-k` and `--prefix` parameters now default to values used RefSeq database.
	* Add option to specify number of cores to use.
	* Add option to disable progress bar printing.


## 0.5.1

Minor edits to project README and metadata.


## 0.5.0

* New features
	* `gambit dist` command for calculating distance matrices.
* CLI
	* Sequence file input
		* Explicitly restrict input to FASTA format only.
		* Files may be gzipped.
		* Read input file lists from text files.
	* Minor changes to options of subcommands in `signatures` group.
* API
	* `gambit.db` subpackage:
		* Database-loading funcs moved to class methods of `ReferenceDatabase`.
		* Additional taxonomy tree methods.
		* Some additional internal reorganization/refactoring.


## 0.4.0

* New features
	* Result reporting
		* Results include list of closest reference genomes. This is only reported in JSON-based
		  output formats.
		* New "next_taxon" attribute, indicating the next most specific taxon for which the
		  threshold was not met.
* CLI
	* `signatures info` subcommand uses current reference DB by default.
* Documentation
	* Some improvements to API docs.
* API and internals
	* `calc_signature()` function can take multiple sequences as input.
	* Remove `calc_signature_parse()` function.
	* Refactoring
		* Rename `GAMBITDatabase` -> `ReferenceDatabase`, `gambit.db.gambitdb` -> `.refdb`
		* Rename `gambit.signatures` -> `gambit.sigs`.
		* Merge `gambit.sigs.array`, `gambit.sigs.meta` -> `gambit.sigs.base`
		* Rename `gambit.io.export` -> `gambit.results`
		* Move generic sequence code from `gambit.kmers` to `gambit.seq`.
		* Merge `gambit.io.seq` -> `gambit.seq`.
		* Rename `load_database*` funcs -> `load_db*`.
		* Move `gambit.io.json` -> `gambit.util.json`, `gambit.io.util` -> `gambit.util.io`,
		  remove `gambit.io`.
		* Moved some other stuff between modules.
	* Improvements to `gambit.sigs.hdf5.HDF5Signatures`
		* Improvements to `.create()` method.
		* Support compression.
	* Format-independent functions for reading/writing signature data.
	* `jaccarddist_pairwise()` function.
	* Add more tree-based methods to `Taxon`.
	* `gambit.metric` changes
		* `jaccarddist_array` and `jaccarddist_matrix` functions now accept any sequence type (e.g.
		  `list`) for the `refs` argument, but with diminished performance.


## 0.3.0

* CLI updates
	* `gambit query` now accepts query signatures from a signature file.
	* New command group `gambit signatures` with `info` and `create` subcommands.
	* New `debug` command group (hidden).
* Performance enhancements
	* Signature calculation for multiple sequence files can be run in parallel.
	* Signature calculation with large `k` much faster.
	* Benchmarks for signature calculation.
* Documentation
	* Installation instructions
	* More complete CLI docs
* API and internals
	* Major refactor to `gambit.kmers` and `gambit.signatures`
	* `find_kmers()` renamed to `calc_signature()` and moved to `gambit.signatures.calc`, related
		functions also renamed and moved.
	* Refactored k-mer search into new `find_kmers()` function, which finds locations of prefix
		matches in sequence.
	* Several other classes and functions moved from `gambit.kmers` to `gambit.signatures` submodules.
	* Rearrangement of stuff within `gambit.signatures`.
	* Added required `kmerspec` attribute to `AbstractKmerArray`.
	* Renamed some `KmerSpec` attributes
	* Rename `gambit.kmers.reverse_complement()` -> `revcomp()`
	* Refactor of Jaccard functions
	* Removed `_sparse` from function names
	* Array and matrix functions now calculate distance only, renamed from `jaccard_*` to `jaccarddist_*`
* New features
	* Most functions which take DNA sequences now accept `str`, `bytes`, or `Bio.Seq.Seq`.
	* Convert signatures between compatible `KmerSpec`s.
	* `HDF5Signatures` `close()` method and context manager.
* Other
	* Updated Cython `kmers` code.
	* Many updates/improvements to tests.


## 0.2.2

* Replace `testdb_210126` with `testdb_210818`. Small enough to include all files, including reference signatures and query sequences, in version control.
* Store pre-calculated query results for tests.
* Some other minor test improvements and bug fixes.


## 0.2.1

* Added license
* Add `setuptools` to runtime dependencies
* Minor docstring edits


## 0.2.0

* User-facing
	* Rework JSON results format to be simpler and hide internal details
	* Add CSV results format (default)
	* Display progress while querying
	* Increase query performance and decrease memory usage
* Internal
	* Major redesign of `GAMBITDatabase` and query funcs
	* `GAMBITDatabase` stores indices of reference signatures instead of loading them all up front
	* Read reference signatures in chunks when calculating distance matrix
	* Maintain reference to SQLAlchemy Session object on `GAMBITDatabase`.
	* `strict` classification parameter
	* Enables new behavior of finding and reconciling all matching taxa
	* Defaults to off, which results old behavior of using only closest match
	* Fix bug in `consensus_taxon()` and add tests
	* Flexible, generic progress monitoring API
	* Add to long-running functions like querying, distance matrix calculation, and k-mer finding
	* "archive" export format for saving full result data.
	* Lots of test improvements
	* More type annotations
	* Update API docs


## 0.1.0

Initial version.
