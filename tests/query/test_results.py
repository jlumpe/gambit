"""Test gambit.query.results module."""

import pytest

from gambit.query.results import QueryInput
from gambit.io.seq import SequenceFile


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
