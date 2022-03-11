#!/usr/bin/env python3

"""
Runs query with testdb_210818 database and its included query sequences and saves results in archive
format.

This script should be re-run whenever the expected results change.
"""

import sys
from pathlib import Path
from csv import DictReader

from gambit.seq import SequenceFile
from gambit.db import ReferenceDatabase
from gambit.db.models import reportable_taxon
from gambit.query import QueryParams, query_parse
from gambit.results.archive import ResultsArchiveWriter
from gambit.util.misc import zip_strict


PARAMS = {
	'non_strict': QueryParams(classify_strict=False, report_closest=10),
	'strict': QueryParams(classify_strict=True, report_closest=10),
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
	strict = results.params.classify_strict

	for query, item in zip_strict(queries, results.items):
		warnings = []

		clsresult = item.classifier_result
		predicted = clsresult.predicted_taxon

		assert item.input.file == query['file']

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
		assert item.closest_genomes[0].genome.description == query['closest']

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
	queries = load_query_data()
	query_files = [query['file'] for query in queries]
	db = ReferenceDatabase.load_from_dir('.')

	writer = ResultsArchiveWriter(pretty=True)

	for label, params in PARAMS.items():
		results = query_parse(db, query_files, params)
		check_results(queries, results)

		with open(f'results/{label}.json', 'wt') as f:
			writer.export(f, results)


if __name__ == '__main__':
	main()
