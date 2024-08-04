"""Funcs for testing exported data."""

import csv
import json
from typing import TextIO, Any, Iterable, Optional
from pathlib import Path

import numpy as np

from gambit.util.json import to_json
from gambit.query import QueryResults, QueryResultItem
from gambit.classify import GenomeMatch, ClassifierResult
from gambit.util.misc import zip_strict
from gambit.db.models import AnnotatedGenome, Taxon


def compare_genome_matches(match1: Optional[GenomeMatch], match2: Optional[GenomeMatch]):
	"""Assert two ``GenomeMatch`` instances are equal.

	The values for the ``distance`` attribute are only checked for approximate equality, to support
	instances where one was loaded from a results archive (saving and loading a float in JSON is
	lossy).

	Also allows one or both values to be None.
	"""
	if match1 is None or match2 is None:
		assert match1 is None and match2 is None
		return

	assert match1.genome == match2.genome
	assert match1.matched_taxon == match2.matched_taxon
	assert np.isclose(match1.distance, match2.distance)


def compare_classifier_results(result1: ClassifierResult, result2: ClassifierResult) -> bool:
	"""Assert two ``ClassifierResult`` instances are equal."""
	assert result1.success == result2.success
	assert result1.predicted_taxon == result2.predicted_taxon
	compare_genome_matches(result1.primary_match, result2.primary_match)
	compare_genome_matches(result1.closest_match, result2.closest_match)
	assert result1.next_taxon == result2.next_taxon
	assert set(result1.warnings) == set(result2.warnings)
	assert result1.error == result2.error


def compare_result_items(item1: QueryResultItem, item2: QueryResultItem) -> bool:
	"""Assert two ``QueryResultItem`` instances are equal.

	Does not compare the value of the ``input`` attributes.
	"""
	assert item1.report_taxon == item2.report_taxon
	compare_classifier_results(item1.classifier_result, item2.classifier_result)

	assert len(item1.closest_genomes) == len(item2.closest_genomes)
	for m1, m2 in zip(item1.closest_genomes, item2.closest_genomes):
		compare_genome_matches(m1, m2)


def cmp_json_attrs(data: dict[str, Any], obj, attrnames: Iterable[str]):
	"""Assert JSON data values equals object attribute values for the given keys/names."""

	for attr in attrnames:
		assert data[attr] == getattr(obj, attr)


def cmp_taxon_json(data: dict[str, Any], taxon: Optional[Taxon]):
	"""Assert Taxon instance matches data in JSON export."""

	if taxon is None:
		assert data is None

	else:
		assert data is not None
		cmp_json_attrs(data, taxon, ['id', 'key', 'name', 'ncbi_id', 'rank'])
		if taxon.distance_threshold is None:
			assert data['distance_threshold'] is None
		else:
			assert data['distance_threshold'] is not None
			assert np.isclose(data['distance_threshold'], taxon.distance_threshold)


def cmp_annnotatedgenome_json(data: dict[str, Any], genome: AnnotatedGenome):
	"""Assert AnnotatedGenome instance matches data in JSON export."""

	assert data['id'] == genome.genome_id
	cmp_json_attrs(
		data,
		genome,
		['key', 'description', 'organism', 'ncbi_db', 'ncbi_id', 'genbank_acc', 'refseq_acc'],
	)
	for taxon_data, taxon in zip_strict(data['taxonomy'], genome.taxon.ancestors(True)):
		cmp_taxon_json(taxon_data, taxon)


def cmp_genomematch_json(data, match: GenomeMatch):
	"""Assert GenomeMatch instance matches data in JSON export."""

	assert np.isclose(data['distance'], match.distance)
	cmp_annnotatedgenome_json(data['genome'], match.genome)

	cmp_taxon_json(data['matched_taxon'], match.matched_taxon)


def check_json_results(file: TextIO,
                       results: QueryResults,
                       strict: bool = False,
                       ):
	"""Assert exported JSON data matches the given results object.

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
	cmp_json_attrs(data['genomeset'], results.genomeset, ['id', 'key', 'version', 'name', 'description'])
	assert data['signaturesmeta'] == to_json(results.signaturesmeta)

	if strict:
		assert data['timestamp'] == to_json(results.timestamp)
		assert data['gambit_version'] == results.gambit_version
		assert data['extra'] == results.extra

	for item, item_data in zip(results.items, data['items']):

		# Compare data['query'] <-> item.input
		query = item_data['query']
		assert query['name'] == item.input.label

		if item.input.file is None:
			assert query['path'] is None
			assert query['format'] is None

		else:
			assert query['format'] == item.input.file.format

			# Check path matches exactly if strict mode, otherwise just file name
			if strict:
				assert query['path'] == str(item.input.file.path)
			else:
				assert Path(query['path']).name == item.input.file.path.name

		# Predicted/next taxon
		cmp_taxon_json(item_data['predicted_taxon'], item.report_taxon)
		cmp_taxon_json(item_data['next_taxon'], item.classifier_result.next_taxon)

		# Closest genomes
		assert len(item_data['closest_genomes']) == len(item.closest_genomes)
		for match, match_data in zip_strict(item.closest_genomes, item_data['closest_genomes']):
			cmp_genomematch_json(match_data, match)


def cmp_csv_taxon(row: dict[str, str], taxon: Optional[Taxon], prefix: str):

	if taxon is None:
		assert row[prefix + '.name'] == ''
		assert row[prefix + '.rank'] == ''
		assert row[prefix + '.ncbi_id'] == ''
		assert row[prefix + '.threshold'] == ''

	else:
		assert row[prefix + '.name'] == taxon.name
		assert row[prefix + '.rank'] == taxon.rank
		assert row[prefix + '.ncbi_id'] == str(taxon.ncbi_id or '')

		dt = row[prefix + '.threshold']
		if taxon.distance_threshold is None:
			assert dt == ''
		else:
			assert np.isclose(float(dt), taxon.distance_threshold)


def check_csv_results(file: TextIO, results: QueryResults, strict: bool = False):
	"""Assert exported CSV data matches the given results object.

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
		assert row['query'] == item.input.label

		cmp_csv_taxon(row, item.report_taxon, 'predicted')
		cmp_csv_taxon(row, item.classifier_result.next_taxon, 'next')

		closest = item.closest_genomes[0]
		assert np.isclose(float(row['closest.distance']), closest.distance)
		assert row['closest.description'] == closest.genome.description
