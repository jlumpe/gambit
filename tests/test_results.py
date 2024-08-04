from io import StringIO

import pytest

from gambit.query import QueryResults, QueryResultItem, QueryInput, QueryParams
from gambit.classify import ClassifierResult, GenomeMatch
from gambit.db import ReferenceGenomeSet, Genome
from gambit.sigs import SignaturesMeta
from gambit.seq import SequenceFile
from gambit.results.json import JSONResultsExporter
from gambit.results.csv import CSVResultsExporter
from gambit.results.archive import ResultsArchiveReader, ResultsArchiveWriter
from .results import check_json_results, check_csv_results


def export_to_buffer(results: QueryResults, exporter) -> StringIO:
	"""Export query results to a `StringIO` buffer."""
	buf = StringIO()
	exporter.export(buf, results)
	buf.seek(0)
	return buf


@pytest.fixture()
def session(testdb):
	return testdb.Session()


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
			closest_genomes=[cr.closest_match],
		))

	# Set one input file to None
	items[-1].input.file = None

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


def test_json(results: QueryResults):
	"""Test JSONResultsExporter."""
	exporter = JSONResultsExporter()
	buf = export_to_buffer(results, exporter)
	check_json_results(buf, results, strict=True)


def test_csv(results: QueryResults):
	"""Test CSVResultsExporter."""
	exporter = CSVResultsExporter()
	buf = export_to_buffer(results, exporter)
	check_csv_results(buf, results, strict=True)


def test_results_archive(session, results: QueryResults):
	"""Test ResultArchiveWriter/Reader."""
	writer = ResultsArchiveWriter()
	buf = export_to_buffer(results, writer)
	reader = ResultsArchiveReader(session)
	results2 = reader.read(buf)

	assert results2 == results
