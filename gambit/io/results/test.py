"""Funcs for testing exported data."""

import csv
import json
from typing import TextIO
from pathlib import Path

import numpy as np

from gambit.io.json import to_json
from gambit.query import QueryResults
from gambit.util.misc import zip_strict


def cmp_json_attrs(data, obj, attrnames):
	for attr in attrnames:
		assert data[attr] == getattr(obj, attr)


def check_json_results(file: TextIO,
                       results: QueryResults,
                       strict: bool = False,
                       ):
	"""Check exported JSON data matches the given results object.

	Parameters
	----------
	file
		Opened results file.
	results
		Query results to check against.
	strict
		If True, expect that ``data`` was exported from the exact same ``results`` object. Otherwise
		expect results from a separate query run with the same inputs.

	Raises
	------
	AssertionError
		If any of the checks fail.
	"""

	data = json.load(file)

	assert len(data['items']) == len(results.items)
	# assert data['params'] == to_json(results.params)
	cmp_json_attrs(data['genomeset'], results.genomeset, ['id', 'key', 'version', 'name', 'description'])
	assert data['signaturesmeta'] == to_json(results.signaturesmeta)
	assert data['gambit_version'] == results.gambit_version
	assert data['extra'] == results.extra

	if strict:
		assert data['timestamp'] == to_json(results.timestamp)

	for item, item_data in zip(results.items, data['items']):
		query = item_data['query']
		assert query['name'] == item.input.label

		if item.input.file is None:
			assert query['path'] is None
			assert query['format'] is None

		else:
			assert query['format'] == item.input.file.format

			if strict:
				assert query['path'] == str(item.input.file.path)
			else:
				assert Path(query['path']).name == item.input.file.path.name

		predicted_data = item_data['predicted_taxon']
		if item.report_taxon is None:
			assert predicted_data is None
		else:
			cmp_json_attrs(predicted_data, item.report_taxon, ['id', 'key', 'name', 'ncbi_id', 'rank'])
			assert np.isclose(predicted_data['distance_threshold'], item.report_taxon.distance_threshold)

		closest = item.classifier_result.closest_match
		assert np.isclose(item_data['closest_genome_distance'], item.classifier_result.closest_match.distance)

		cmp_json_attrs(
			item_data['closest_genome'],
			closest.genome,
			['key', 'description', 'organism', 'ncbi_db', 'ncbi_id', 'genbank_acc', 'refseq_acc'],
		)
		assert item_data['closest_genome']['id'] == closest.genome.genome_id

		for taxon, taxon_data in zip_strict(closest.genome.taxon.ancestors(True), item_data['closest_genome']['taxonomy']):
			cmp_json_attrs(taxon_data, taxon, ['id', 'key', 'name', 'ncbi_id', 'rank', 'distance_threshold'])


def check_csv_results(file: TextIO,
                      results: QueryResults,
                      strict: bool = False,
                      ):
	"""Check exported CSV data matches the given results object.

	Parameters
	----------
	file
		Opened results file.
	results
		Query results to check against.
	strict
		If True, expect that ``data`` was exported from the exact same ``results`` object. Otherwise
		expect results from a separate query run with the same inputs.

	Raises
	------
	AssertionError
		If any of the checks fail.
	"""

	rows = list(csv.DictReader(file))
	assert len(rows) == len(results.items)

	for item, row in zip(results.items, rows):
		assert row['query.name'] == item.input.label

		if item.input.file is None:
			assert row['query.path'] == ''
		elif strict:
			assert row['query.path'] == str(item.input.file.path)
		else:
			assert Path(row['query.path']).name == item.input.file.path.name

		if item.report_taxon is None:
			assert row['predicted.name'] == ''
			assert row['predicted.rank'] == ''
			assert row['predicted.ncbi_id'] == ''
			assert row['predicted.threshold'] == ''
		else:
			assert row['predicted.name'] == item.report_taxon.name
			assert row['predicted.rank'] == item.report_taxon.rank
			assert row['predicted.ncbi_id'] == str(item.report_taxon.ncbi_id or '')
			assert np.isclose(float(row['predicted.threshold']), item.report_taxon.distance_threshold)

		assert np.isclose(float(row['closest.distance']), item.classifier_result.closest_match.distance)
		assert row['closest.description'] == item.classifier_result.closest_match.genome.description
