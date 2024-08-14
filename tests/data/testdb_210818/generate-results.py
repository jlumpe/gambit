#!/usr/bin/env python3

"""
Runs query with testdb_210818 database and its included query sequences and saves results in archive
format.

This script should be re-run whenever the expected results change.
"""

import sys
from pathlib import Path

from gambit.seq import SequenceFile
from gambit.db import reportable_taxon
from gambit.query import QueryParams, QueryResults, query_parse
from gambit.results import ResultsArchiveWriter
from gambit.util.misc import zip_strict


THISDIR = Path(__file__).parent
ROOTDIR = THISDIR.parent.parent.parent

sys.path.insert(0, str(ROOTDIR))
from tests.testdb import TestDB, TestQueryGenome


PARAMS = {
	'non_strict': QueryParams(classify_strict=False, report_closest=10),
	'strict': QueryParams(classify_strict=True, report_closest=10),
}


def check_results(queries: list[TestQueryGenome], query_files: list[SequenceFile], results: QueryResults):
	"""Check query results object against queries.csv table before exporting."""

	strict = results.params.classify_strict

	for query, query_file, item in zip_strict(queries, query_files, results.items):
		warnings = []

		clsresult = item.classifier_result
		predicted = clsresult.predicted_taxon

		assert item.input.file == query_file

		# No errors
		assert clsresult.success
		assert clsresult.error is None

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
				assert clsresult.primary_match == clsresult.closest_match
				assert predicted is clsresult.primary_match.matched_taxon

			assert item.report_taxon is reportable_taxon(predicted)

		else:
			assert predicted is None
			assert clsresult.primary_match is None
			assert item.report_taxon is None

		# Closest matches
		assert len(item.closest_genomes) == results.params.report_closest
		assert item.closest_genomes[0] == clsresult.closest_match

		for i in range(1, results.params.report_closest):
			assert item.closest_genomes[i].distance >= item.closest_genomes[i-1].distance

		# Next taxon
		nt = clsresult.next_taxon
		if nt is None:
			# Predicted should be most specific possible
			assert clsresult.closest_match.matched_taxon == clsresult.closest_match.genome.taxon

		else:
			assert nt.distance_threshold is not None
			assert nt.distance_threshold < clsresult.closest_match.distance

			# This should hold true as long as the primary match is the closest match, just warn if
			# it fails.
			if predicted is not None:
				if predicted not in nt.ancestors():
					warnings.append(f'Next taxon {nt.name} not a descendant of predicted taxon {predicted.name}')

		# Display warnings
		for w in warnings:
			print(f'[Query "{query["name"]}"]:', w, file=sys.stderr)


def main():
	testdb = TestDB(THISDIR)
	db = testdb.refdb
	query_files = testdb.get_query_files(relative=True)

	writer = ResultsArchiveWriter(pretty=True)

	for label, params in PARAMS.items():
		print('Running query:', label)
		results = query_parse(db, query_files, params)
		check_results(testdb.query_genomes, query_files, results)

		with open(f'results/{label}.json', 'wt') as f:
			writer.export(f, results)

		print('done!\n\n')


if __name__ == '__main__':
	main()
