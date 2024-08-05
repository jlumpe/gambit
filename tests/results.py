"""Helper code for tests related to the QueryResults class or exported result data."""

import csv
import json
from typing import TextIO, Any, Iterable, Optional
from pathlib import Path
from warnings import warn

import numpy as np

from gambit.util.json import to_json
from gambit.query import QueryResults, QueryResultItem, QueryParams
from gambit.classify import GenomeMatch, ClassifierResult
from gambit.util.misc import zip_strict
from gambit.db.models import AnnotatedGenome, Taxon, reportable_taxon


def check_results(results: QueryResults, warnings: bool = True):
	"""Check invariants on query results object."""

	assert results.params is not None

	for item in results.items:
		check_result_item(item, results.params, warnings=warnings)


def check_result_item(item: QueryResultItem, params: QueryParams, warnings: bool = True):
	"""Check invariants on successful query result item."""

	clsresult = item.classifier_result
	predicted = clsresult.predicted_taxon

	# No errors
	assert clsresult.success
	assert clsresult.error is None

	# Predicted taxon
	if predicted is not None:
		assert clsresult.primary_match is not None

		if not params.classify_strict:
			assert clsresult.primary_match == clsresult.closest_match
			assert predicted is clsresult.primary_match.matched_taxon

		assert item.report_taxon is reportable_taxon(predicted)

	else:
		assert clsresult.primary_match is None
		assert item.report_taxon is None

	# Closest matches
	assert len(item.closest_genomes) == params.report_closest
	assert item.closest_genomes[0] == clsresult.closest_match

	# Check closest_genomes is sorted by distance
	for i in range(1, params.report_closest):
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
				if warnings:
					warn(
						f'[Query {item.label}]: '
						f'next taxon {nt.name} not a descendant of predicted taxon {predicted.name}'
					)


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


def compare_classifier_results(result1: ClassifierResult, result2: ClassifierResult):
	"""Assert two ``ClassifierResult`` instances are equal."""
	assert result1.success == result2.success
	assert result1.predicted_taxon == result2.predicted_taxon
	compare_genome_matches(result1.primary_match, result2.primary_match)
	compare_genome_matches(result1.closest_match, result2.closest_match)
	assert result1.next_taxon == result2.next_taxon
	assert set(result1.warnings) == set(result2.warnings)
	assert result1.error == result2.error


def compare_result_items(item1: QueryResultItem, item2: QueryResultItem):
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

		# Compare data['query'] <-> item.label / item.file
		query = item_data['query']
		assert query['name'] == item.label

		if item.file is None:
			assert query['path'] is None

		else:
			# Check path matches exactly if strict mode, otherwise just file name
			if strict:
				assert query['path'] == str(item.file)
			else:
				assert Path(query['path']).name == item.file.name

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
		assert row['query'] == item.label

		cmp_csv_taxon(row, item.report_taxon, 'predicted')
		cmp_csv_taxon(row, item.classifier_result.next_taxon, 'next')

		closest = item.closest_genomes[0]
		assert np.isclose(float(row['closest.distance']), closest.distance)
		assert row['closest.description'] == closest.genome.description
