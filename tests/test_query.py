"""Test the gambit.query module."""

import pytest

from gambit.query import QueryInput, QueryResults, query, query_parse
from gambit.seq import SequenceFile
from gambit.util.misc import zip_strict
from gambit import __version__ as GAMBIT_VERSION

from .testdb import TestDB
from .results import compare_result_items


class TestQueryInput:
	"""Test QueryInput class."""

	def test_convert(self):
		file = SequenceFile('path/to/file.fa', 'fasta')
		qi = QueryInput('foo', file)

		assert QueryInput.convert(qi) is qi
		assert QueryInput.convert('foo') == QueryInput('foo', None)
		assert QueryInput.convert(file) == QueryInput(str(file.path), file)

		with pytest.raises(TypeError):
			QueryInput.convert(3.4)


@pytest.mark.parametrize('strict', [False, True])
class TestQuery:
	"""Run a full query using the Python API."""

	def check_results(self, results: QueryResults, ref_results: QueryResults):

		assert results.params == ref_results.params
		assert results.genomeset == ref_results.genomeset
		assert results.signaturesmeta == ref_results.signaturesmeta
		assert results.gambit_version == GAMBIT_VERSION

		for item, ref_item in zip_strict(results.items, ref_results.items):
			compare_result_items(item, ref_item)

	def test_query(self, testdb: TestDB, strict: bool):
		"""Test the query() function."""

		ref_results = testdb.get_query_results(strict)
		params = ref_results.params
		query_sigs = testdb.query_signatures

		results = query(testdb.refdb, query_sigs, params)
		self.check_results(results, ref_results)

		for sigid, item in zip_strict(query_sigs.ids, results.items):
			assert item.input.file is None
			# assert item.input.label == sigid

	def test_query_parse(self, testdb: TestDB, strict: bool):
		"""Test the query_parse() function."""

		ref_results = testdb.get_query_results(strict)
		params = ref_results.params
		query_files = testdb.get_query_files()

		results = query_parse(testdb.refdb, query_files, params)
		self.check_results(results, ref_results)

		for file, item in zip_strict(query_files, results.items):
			assert item.input.file == file
			assert item.input.label == str(file.path)
