#!/usr/bin/env python3

"""
Runs query with testdb_210818 database and its included query sequences and saves results in archive
format.

This script should be re-run whenever the expected results change.
"""

from pathlib import Path
from csv import DictReader

from gambit.io.seq import SequenceFile
from gambit.db import load_db_from_dir
from gambit.query import QueryParams, query_parse
from gambit.io.results.archive import ResultsArchiveWriter
from gambit.util.misc import zip_strict


PARAMS = {
	'non_strict': QueryParams(classify_strict=False),
	'strict': QueryParams(classify_strict=True),
}


def load_query_data():
	with open('queries/queries.csv', newline='') as f:
		rows = list(DictReader(f))

	genomes_dir = Path('queries/genomes')

	for row in rows:
		row['warnings'] = row['warnings'].lower() == 'true'
		row['file'] = SequenceFile(
			path=genomes_dir / (row['name'] + '.fasta'),
			format='fasta',
		)

	return rows


def check_results(queries, results):

	for query, item in zip_strict(queries, results.items):
		clsresult = item.classifier_result
		predicted = clsresult.predicted_taxon

		assert item.input.file == query['file']
		assert clsresult.success
		assert clsresult.error is None

		if results.params.classify_strict:
			if query['predicted']:
				assert predicted is not None
				assert predicted.name == query['predicted']
				assert clsresult.primary_match is not None
				assert clsresult.primary_match.genome.description == query['primary']
				assert item.report_taxon is (predicted if predicted.report else predicted.parent)

			else:
				assert predicted is None
				assert clsresult.primary_match is None
				assert item.report_taxon is None

			assert clsresult.closest_match.genome.description == query['closest']
			assert bool(clsresult.warnings) == query['warnings']

		else:
			if query['predicted']:
				assert clsresult.primary_match == clsresult.closest_match
				assert predicted is clsresult.primary_match.matched_taxon
				assert item.report_taxon is (predicted if predicted.report else predicted.parent)

			else:
				assert predicted is None
				assert clsresult.primary_match is None
				assert item.report_taxon is None

			assert clsresult.closest_match.genome.description == query['closest']


def main():
	queries = load_query_data()
	query_files = [query['file'] for query in queries]
	db = load_db_from_dir('')

	writer = ResultsArchiveWriter(pretty=True)

	for label, params in PARAMS.items():
		results = query_parse(db, query_files, params)
		check_results(queries, results)

		with open(f'results/{label}.json', 'wt') as f:
			writer.export(f, results)


if __name__ == '__main__':
	main()
