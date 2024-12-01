#!/usr/bin/env python3

"""
Runs query with testdb_210818 database and its included query sequences and saves results in archive
format.

This script should be re-run whenever the expected results change.
"""

import sys
from pathlib import Path

from gambit.query import QueryParams, QueryResults, query_parse
from gambit.results import ResultsArchiveWriter
from gambit.util.misc import zip_strict
from gambit.util.io import FilePath


THISDIR = Path(__file__).parent
ROOTDIR = THISDIR.parent.parent.parent

sys.path.insert(0, str(ROOTDIR))
from tests.testdb import TestDB, TestQueryGenome
from tests.results import check_results as check_results_base


PARAMS = {
	'non_strict': QueryParams(classify_strict=False, report_closest=10),
	'strict': QueryParams(classify_strict=True, report_closest=10),
}


def check_results(queries: list[TestQueryGenome], query_files: list[FilePath], results: QueryResults):
	"""Check query results object against queries.csv table before exporting."""

	strict = results.params.classify_strict

	for query, query_file, item in zip_strict(queries, query_files, results.items):

		clsresult = item.classifier_result
		predicted = clsresult.predicted_taxon

		assert item.file == Path(query_file)

		# Check if warnings expected (only if in strict mode)
		assert bool(clsresult.warnings) == (strict and query['warnings'])

		# Closest
		assert clsresult.closest_match.genome.description == query['closest']

		# Predicted taxon
		if query['predicted']:
			assert predicted is not None
			assert clsresult.primary_match is not None

			if strict:
				assert predicted.name == query['predicted']
				assert clsresult.primary_match.genome.description == query['primary']

		else:
			assert predicted is None
			assert clsresult.primary_match is None
			assert item.report_taxon is None


def main():
	testdb = TestDB(THISDIR)
	db = testdb.refdb
	query_files = testdb.get_query_files(relative=True)

	writer = ResultsArchiveWriter(pretty=True)

	for label, params in PARAMS.items():
		print('Running query:', label)
		results = query_parse(db, query_files, params)
		check_results_base(results)
		check_results(testdb.query_genomes, query_files, results)

		with open(f'results/{label}.json', 'wt') as f:
			writer.export(f, results)

		print('done!\n\n')


if __name__ == '__main__':
	main()
