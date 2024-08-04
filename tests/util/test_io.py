"""Test gambit.util.io."""

from pathlib import Path

import pytest
import numpy as np

import gambit.util.io as ioutil


class TestOpenCompressed:
	"""Test open_compressed()"""

	@pytest.fixture(scope='class')
	def text_data(self):
		"""Random printable characters encoded as ASCII."""
		random = np.random.RandomState()
		return random.randint(32, 128, size=1000, dtype='b').tobytes()

	@pytest.fixture(scope='class', params=['none', 'gzip'])
	def compression(self, request):
		"""Compression method string."""
		return request.param

	@pytest.fixture()
	def text_file(self, text_data: bytes, compression: str, tmp_path: Path):
		"""Path to file with text_data written to it using open_compressed."""

		file = tmp_path / 'chars.txt'

		with ioutil.open_compressed(file, 'wb', compression) as fobj:
			fobj.write(text_data)

		return file

	@pytest.mark.parametrize('binary', [True, False])
	@pytest.mark.parametrize('auto', [True, False])
	def test_read(self, binary: bool, auto: bool, text_data: bytes, text_file: Path, compression: str):
		"""Test we can read the file in both binary and text mode."""

		mode = 'rb' if binary else 'rt'

		with ioutil.open_compressed(text_file, mode, 'auto' if auto else compression) as fobj:
			contents = fobj.read()

		if binary:
			assert isinstance(contents, bytes)
			assert contents == text_data

		else:
			assert isinstance(contents, str)
			assert contents == text_data.decode('ascii')

	@pytest.mark.parametrize('write_mode', ['w', 'a', 'x'])
	@pytest.mark.parametrize('binary', [True, False])
	def test_write(self, write_mode: str, binary: bool, text_data: bytes, compression: str, tmp_path: Path):
		"""
		Test writing data using the w, a, and x modes.

		TODO - these are all identical when the file doesn't exist, test behavior when it does.
		"""

		file = tmp_path / 'chars.txt'
		mode = write_mode + ('b' if binary else 't')
		to_write = text_data if binary else text_data.decode('ascii')

		with ioutil.open_compressed(file, mode, compression) as fobj:
			fobj.write(to_write)

		with ioutil.open_compressed(file, 'rb', compression) as f:
			contents = f.read()

		assert contents == text_data


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
		return ioutil.ClosingIterator(iterable, fobj)

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


class TestMaybeOpen:
	"""Test the maybe_open() function."""

	@pytest.mark.parametrize('pathtype', [str, Path])
	def test_with_path(self, tmpdir, pathtype):
		path = pathtype(tmpdir / 'test.txt')

		with ioutil.maybe_open(path, 'w') as f:
			assert not f.closed
			assert f.mode == 'w'
			f.write('test')

		assert f.closed

		with ioutil.maybe_open(path, 'r') as f:
			assert not f.closed
			assert f.mode == 'r'
			assert f.read() == 'test'

		assert f.closed

	def test_with_fileobj(self, tmpdir):
		with open(tmpdir / 'test.txt', 'w') as f:
			with ioutil.maybe_open(f) as f2:
				assert f2 is f

			assert not f.closed


@pytest.mark.parametrize('strip', [False, True])
@pytest.mark.parametrize('skip_empty', [False, True])
def test_read_lines(strip, skip_empty, tmp_path):
	"""Test the read_lines() function."""
	file = tmp_path / 'test.txt'
	lines = list(map(str, range(20)))

	with open(file, 'w') as f:
		for i, line in enumerate(lines):
			# Extra random whitespace
			f.write(f'\t{line}  ' if strip else line)
			f.write('\n')
			if skip_empty and i % 3 == 0:
				f.write('  \n' if strip else '\n')

	lines2 = list(ioutil.read_lines(file, strip, skip_empty))
	assert lines2 == lines
