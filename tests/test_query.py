"""Test the gambit.query module."""

import pytest

from gambit.query import QueryInput, query_parse, compare_result_items
from gambit.seq import SequenceFile
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


def test_query_python(testdb, testdb_queries, testdb_results):
	"""Run a full query using the Python API."""
	ref_results = testdb_results
	params = ref_results.params
	query_files = [item['file'] for item in testdb_queries]

	results = query_parse(testdb, query_files, params)

	assert results.params == params
	assert results.genomeset == testdb.genomeset
	assert results.signaturesmeta == testdb.signatures.meta
	assert results.gambit_version == GAMBIT_VERSION

	for file, item, ref_item in zip_strict(query_files, results.items, ref_results.items):
		assert item.input.file == file
		compare_result_items(item, ref_item)
