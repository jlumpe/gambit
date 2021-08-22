"""Test the gambit.query module."""

import pytest

from gambit.query import QueryInput, QueryParams, query_parse
from gambit.io.seq import SequenceFile
from gambit.util.misc import zip_strict
from gambit import __version__ as GAMBIT_VERSION


class TestQueryInput:
	"""Test QueryInput class."""

	def test_convert(self):
		file = SequenceFile('path/to/file.fa', 'fasta')
		qi = QueryInput('foo', file)

		assert QueryInput.convert(qi) is qi
		assert QueryInput.convert('foo') == QueryInput('foo', None)
		assert QueryInput.convert(file) == QueryInput('file.fa', file)

		with pytest.raises(TypeError):
			QueryInput.convert(3.4)


@pytest.mark.parametrize('classify_strict', [False, True])
def test_query_python(testdb, testdb_queries, classify_strict):
	"""Run a full query using the Python API."""

	query_files = [item['file'] for item in testdb_queries]
	params = QueryParams(classify_strict=classify_strict)

	results = query_parse(testdb, query_files, params)

	assert results.params == params
	assert results.genomeset == testdb.genomeset
	assert results.signaturesmeta == testdb.signatures.meta
	assert results.gambit_version == GAMBIT_VERSION

	for query, item in zip_strict(testdb_queries, results.items):
		clsresult = item.classifier_result
		predicted = clsresult.predicted_taxon

		assert item.input.file == query['file']
		assert clsresult.success
		assert clsresult.error is None

		if classify_strict:
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
