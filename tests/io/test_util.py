"""Test gambit.io.util"""

import io

import pytest
import numpy as np

from gambit.io import util


class TestOpenCompressed:
	"""Test gambit.io.util.open_compressed."""

	@pytest.fixture(scope='class')
	def text_data(self):
		"""Random printable characters encoded as ASCII."""
		random = np.random.RandomState()
		return random.randint(32, 128, size=1000, dtype='b').tobytes()

	@pytest.fixture(scope='class', params=list(util.COMPRESSED_OPENERS))
	def compression(self, request):
		"""Compression method string."""
		return request.param

	@pytest.fixture(params=['w', 'wt', 'wb'])
	def text_file(self, request, text_data, compression, tmpdir):
		"""Path to file with text_data written to it using open_compressed."""

		file = tmpdir.join('chars.txt')
		mode = request.param

		if mode[-1] != 'b':
			to_write = text_data.decode('ascii')

		else:
			to_write = text_data

		with util.open_compressed(compression, file.strpath, mode) as fobj:
			fobj.write(to_write)

		return file

	@pytest.mark.parametrize('mode,binary', [
		('r', False),
		('rt', False),
		('rb', True),
	])
	def test_read(self, mode, binary, text_data, compression, text_file):
		"""Check that the file is readable and its contents match what was written."""

		with util.open_compressed(compression, text_file.strpath, mode) as fobj:
			contents = fobj.read()

		if binary:
			assert isinstance(contents, bytes)
			assert contents == text_data

		else:
			assert isinstance(contents, str)
			assert contents == text_data.decode('ascii')


class TestClosingIterator:
	"""Test the ClosingIterator class."""

	# Number of lines in text buffer
	NLINES = 100

	@pytest.fixture
	def fobj(self):
		"""Text buffer with a number in each line."""

		from io import StringIO

		# Write some lines to a buffer
		buf = StringIO()

		for i in range(self.NLINES):
			buf.write('{}\n'.format(i))

		buf.seek(0)

		return buf

	@pytest.fixture
	def iterator(self, fobj):
		"""Instance of ClosingIterator using fobj."""

		# Iterable which reads from fobj on every iteration."""
		iterable = (int(line.strip()) for line in fobj)

		# Create ClosingIterator object
		return util.ClosingIterator(iterable, fobj)

	def test_close_on_finish(self, iterator, fobj):
		"""Check that the stream gets closed when the iterator runs out."""

		assert not fobj.closed and not iterator.closed

		for val in iterator:
			assert not fobj.closed and not iterator.closed

		assert fobj.closed and iterator.closed

	def test_context(self, iterator, fobj):

		with iterator as rval:

			# Check __exit__() returns itself
			assert rval is iterator

			assert not fobj.closed and not iterator.closed

		# Should close after context exited
		assert fobj.closed and iterator.closed

	def test_close_method(self, iterator, fobj):
		"""Test the close() method."""

		assert not fobj.closed and not iterator.closed

		iterator.close()

		assert fobj.closed and iterator.closed
