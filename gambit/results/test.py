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

def cmp_taxon_json(taxon_data, taxon):
	if taxon is None:
		assert taxon_data is None
	else:
		assert taxon_data is not None
		cmp_json_attrs(taxon_data, taxon, ['id', 'key', 'name', 'ncbi_id', 'rank', 'distance_threshold'])

def cmp_annnotatedgenome_json(genome_data, genome):
	assert genome_data['id'] == genome.genome_id
	cmp_json_attrs(
		genome_data,
		genome,
		['key', 'description', 'organism', 'ncbi_db', 'ncbi_id', 'genbank_acc', 'refseq_acc'],
	)
	for taxon_data, taxon in zip_strict(genome_data['taxonomy'], genome.taxon.ancestors(True)):
		cmp_taxon_json(taxon_data, taxon)

def cmp_genomematch_json(match_data, match):
	assert np.isclose(match_data['distance'], match.distance)
	cmp_annnotatedgenome_json(match_data['genome'], match.genome)

	assert (match_data['matched_taxon'] is None) == (match.matched_taxon is None)
	if match.matched_taxon is not None:
		cmp_taxon_json(match_data['matched_taxon'], match.matched_taxon)

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

		# Predicted taxon
		predicted_data = item_data['predicted_taxon']
		cmp_taxon_json(predicted_data, item.report_taxon)
		if item.report_taxon is not None:
			assert np.isclose(predicted_data['distance_threshold'], item.report_taxon.distance_threshold)

		# Next taxon
		cmp_taxon_json(item_data['next_taxon'], item.classifier_result.next_taxon)

		# Closest genomes
		for match, match_data in zip_strict(item.closest_genomes, item_data['closest_genomes']):
			cmp_genomematch_json(match_data, match)


def cmp_csv_taxon(row, taxon, prefix):

	if taxon is None:
		assert row[prefix + '.name'] == ''
		assert row[prefix + '.rank'] == ''
		assert row[prefix + '.ncbi_id'] == ''
		assert row[prefix + '.threshold'] == ''
	else:
		assert row[prefix + '.name'] == taxon.name
		assert row[prefix + '.rank'] == taxon.rank
		assert row[prefix + '.ncbi_id'] == str(taxon.ncbi_id or '')
		assert np.isclose(float(row[prefix + '.threshold']), taxon.distance_threshold)


def check_csv_results(file: TextIO, results: QueryResults, strict: bool = False):
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

		cmp_csv_taxon(row, item.report_taxon, 'predicted')
		cmp_csv_taxon(row, item.classifier_result.next_taxon, 'next')

		closest = item.closest_genomes[0]
		assert np.isclose(float(row['closest.distance']), closest.distance)
		assert row['closest.description'] == closest.genome.description
