from io import StringIO

import pytest

from gambit.export.archive import ResultsArchiveReader, ResultsArchiveWriter
from gambit.query import QueryResults, QueryResultItem, QueryInput, QueryParams
from gambit.classify import ClassifierResult, GenomeMatch
from gambit.db import ReferenceGenomeSet, Taxon, Genome
from gambit.signatures import SignaturesMeta
from gambit.io.seq import SequenceFile


@pytest.fixture()
def session(testdb_session):
	return testdb_session()


@pytest.fixture()
def results(session):
	"""Create a fake QueryResults object."""

	gset = session.query(ReferenceGenomeSet).one()

	taxa = gset.taxa.filter_by(rank='subspecies').order_by('id').limit(20).all()
	classifier_results = []

	for i, taxon in enumerate(taxa):
		genome = taxon.genomes.first()
		assert genome is not None

		match = GenomeMatch(
			genome=taxon.genomes.first(),
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

	items = []

	for i, cr in enumerate(classifier_results):
		items.append(QueryResultItem(
			input=QueryInput(f'query-{i}', SequenceFile(f'query-{i}.fasta', 'fasta')),
			classifier_result=cr,
			report_taxon=cr.predicted_taxon,
		))

	return QueryResults(
		items=items,
		params=QueryParams(),
		genomeset=gset,
		signaturesmeta=SignaturesMeta(
			id='test',
			name='Test signatures',
			version='1.0',
			id_attr='key',
		),
		extra=dict(foo=1),
	)


def test_result_archive(session, results):

	buf = StringIO()

	writer = ResultsArchiveWriter(install_info=True)
	writer.export(buf, results)

	reader = ResultsArchiveReader(session)
	buf.seek(0)
	results2 = reader.read(buf)

	assert 'install_info' in results2.extra
	results2.extra.pop('install_info')

	assert results2 == results
