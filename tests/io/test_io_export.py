from io import StringIO
import json
from csv import DictReader

import pytest
import numpy as np

from gambit.query import QueryResults, QueryResultItem, QueryInput, QueryParams
from gambit.classify import ClassifierResult, GenomeMatch
from gambit.db import ReferenceGenomeSet, Genome
from gambit.signatures import SignaturesMeta
from gambit.io.seq import SequenceFile
from gambit.io.json import to_json
from gambit.util.misc import zip_strict
from gambit.io.export.json import JSONResultsExporter
from gambit.io.export.csv import CSVResultsExporter
from gambit.io.export.archive import ResultsArchiveReader, ResultsArchiveWriter


@pytest.fixture()
def session(testdb_session):
	return testdb_session()


@pytest.fixture()
def results(session):
	"""Create a fake QueryResults object."""

	gset = session.query(ReferenceGenomeSet).one()

	# Taxa to use as matches
	taxa = gset.taxa.filter_by(rank='subspecies').order_by('id').limit(20).all()

	# Make classifier results
	classifier_results = []
	for i, taxon in enumerate(taxa):
		# Pretend some of these are in the NCBI database
		if i % 2:
			taxon.ncbi_id = 1000 + i

		genome = taxon.genomes.first()
		assert genome is not None

		match = GenomeMatch(
			genome=genome,
			distance=(i + 1) / 100,
			matched_taxon=taxon,
		)

		classifier_results.append(ClassifierResult(
			success=True,
			predicted_taxon=taxon,
			primary_match=match,
			closest_match=match,
		))

	# Add one bad result
	failed_genome = gset.genomes.join(Genome).order_by(Genome.key).first()
	assert failed_genome is not None
	classifier_results.append(ClassifierResult(
		success=False,
		predicted_taxon=None,
		primary_match=None,
		closest_match=GenomeMatch(
			genome=failed_genome,
			distance=.99,
			matched_taxon=None,
		),
		warnings=['One warning', 'Two warning'],
		error='Error message',
	))

	# Make result items
	items = []
	for i, cr in enumerate(classifier_results):
		predicted = cr.predicted_taxon
		items.append(QueryResultItem(
			input=QueryInput(f'query-{i}', SequenceFile(f'query-{i}.fasta', 'fasta')),
			classifier_result=cr,
			report_taxon=None if predicted is None else predicted.parent if i % 4 == 0 else predicted,
		))

	return QueryResults(
		items=items,
		params=QueryParams(chunksize=1234, classify_strict=True),
		genomeset=gset,
		signaturesmeta=SignaturesMeta(
			id='test',
			name='Test signatures',
			version='1.0',
			id_attr='key',
		),
		extra=dict(foo=1),
	)


def cmp_json_attrs(data, obj, attrnames):
	for attr in attrnames:
		assert data[attr] == getattr(obj, attr)

def test_json(results):
	"""Test JSONResultsExporter."""

	buf = StringIO()
	exporter = JSONResultsExporter()
	exporter.export(buf, results)

	buf.seek(0)
	data = json.load(buf)

	assert len(data['items']) == len(results.items)
	# assert data['params'] == to_json(results.params)
	cmp_json_attrs(data['genomeset'], results.genomeset, ['id', 'key', 'version', 'name', 'description'])
	assert data['signaturesmeta'] == to_json(results.signaturesmeta)
	assert data['gambit_version'] == results.gambit_version
	assert data['timestamp'] == to_json(results.timestamp)
	assert data['extra'] == results.extra

	for item, item_data in zip(results.items, data['items']):
		query = item_data['query']
		assert query['name'] == item.input.label
		assert query['path'] == str(item.input.file.path)
		assert query['format'] == item.input.file.format

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


def test_csv(results):
	"""Test CSVResultsExporter."""

	buf = StringIO()
	exporter = CSVResultsExporter()
	exporter.export(buf, results)

	buf.seek(0)
	rows = list(DictReader(buf, **exporter.format_opts))
	assert len(rows) == len(results.items)

	for item, row in zip(results.items, rows):
		assert row['query.name'] == item.input.label
		assert row['query.path'] == str(item.input.file.path)

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


def test_results_archive(session, results):
	"""Test ResultArchiveWriter/Reader."""

	buf = StringIO()
	writer = ResultsArchiveWriter()
	writer.export(buf, results)

	reader = ResultsArchiveReader(session)
	buf.seek(0)
	results2 = reader.read(buf)

	assert results2 == results
