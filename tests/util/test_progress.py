from contextlib import contextmanager
from unittest.mock import patch
from string import ascii_letters

import pytest

from gambit.util.progress import default_progress_cls, progress_config, get_progress, iter_progress, \
	capture_progress, check_progress
from gambit.util.progress import NullProgressMeter, TestProgressMeter, TqdmProgressMeter, ClickProgressMeter


@contextmanager
def no_tqdm():
	"""Context manager which makes tqdm not importable even if installed."""

	# Setting value of the 'tqdm' key to None in sys.modules will raise a ModuleNotFound error
	# on import, even if the package is installed.
	with patch.dict('sys.modules', tqdm=None):
		yield


class TestDefaultProgressCls:
	"""Test the default_progress_cls() function and get_progress(True)."""

	def test_no_tqdm(self):
		"""Test with tqdm package not available."""
		with no_tqdm():
			with pytest.warns(UserWarning):
				cls = default_progress_cls()
			assert cls is NullProgressMeter

	def test_tqdm(self):
		pytest.importorskip('tqdm')
		assert default_progress_cls() is TqdmProgressMeter


class TestProgressConfigFunc:
	"""Test the progress_config() function."""

	def test_null(self):
		"""Test passing None and False as argument."""
		for arg in [None, False]:
			config = progress_config(arg, foo=1)
			assert config.callable == NullProgressMeter.create
			assert config.kw == dict(foo=1)

	@pytest.mark.parametrize('with_tqdm', [False, True])
	def test_true(self, with_tqdm):
		"""Test passing True as argument."""
		if with_tqdm:
			pytest.importorskip('tqdm')  # Skip if tqdm not available.
			config = progress_config(True, foo=1)
			assert config.callable == TqdmProgressMeter.create
			assert config.kw == dict(foo=1)

		else:
			with no_tqdm():
				with pytest.warns(UserWarning):
					config = progress_config(True, foo=1)
				assert config.callable == NullProgressMeter.create
				assert config.kw == dict(foo=1)

	def test_cls(self):
		"""Test passing AbstractProgressMeter subclass as argument."""
		for cls in [NullProgressMeter, TestProgressMeter]:
			config = progress_config(cls, foo=1)
			assert config.callable == cls.create
			assert config.kw == dict(foo=1)

	def test_str(self):
		for key, cls in [('tqdm', TqdmProgressMeter), ('click', ClickProgressMeter)]:
			config = progress_config(key, foo=1)
			assert config.callable == cls.create
			assert config.kw == dict(foo=1)

	def test_factory(self):
		"""Test passing a factory function as argument."""

		def factory(total, *, initial=None, **kw):
			return TestProgressMeter.create(total, initial=initial, foo=1, **kw)

		config = progress_config(factory, foo=1)
		assert config.callable == factory
		assert config.kw == dict(foo=1)

	def test_progressconfig(self):
		"""Test passing a factory function as argument."""

		config = TestProgressMeter.config(foo=1, bar=2)
		config2 = progress_config(config, bar=20, baz=3)

		assert config2.callable == TestProgressMeter.create
		assert config2.kw == dict(foo=1, bar=20, baz=3)

	def test_invalid(selfj):
		with pytest.raises(TypeError):
			get_progress(0, 100)


class TestGetProgress:
	"""Test the get_progress() function."""

	@pytest.fixture()
	def total(self):
		return 100

	@pytest.fixture(params=[0, 10])
	def initial(self, request):
		return request.param

	def test_null(self, total, initial):
		"""Test passing None and False as argument."""
		for arg in [None, False]:
			assert isinstance(get_progress(arg, total, initial=initial), NullProgressMeter)

	@pytest.mark.parametrize('with_tqdm', [False, True])
	def test_true(self, total, initial, with_tqdm):
		"""Test passing True as argument."""
		if with_tqdm:
			pytest.importorskip('tqdm')  # Skip if tqdm not available.
			meter = get_progress(True, total, initial=initial)
			assert isinstance(meter, TqdmProgressMeter)
			assert meter.total == total
			assert meter.n == initial

		else:
			with no_tqdm():
				with pytest.warns(UserWarning):
					meter = get_progress(True, total, initial=initial)
				assert isinstance(meter, NullProgressMeter)

	def test_cls(self, total, initial):
		"""Test passing AbstractProgressMeter subclass as argument."""
		for cls in [NullProgressMeter, TestProgressMeter]:
			meter = get_progress(cls, total, initial=initial)
			assert isinstance(meter, cls)

			if cls is not NullProgressMeter:
				assert meter.total == total
				assert meter.n == initial

	def test_str(self, total, initial):
		meter = get_progress('click', total, initial=initial)
		assert isinstance(meter, ClickProgressMeter)
		assert meter.total == total
		assert meter.n == initial

	def test_factory(self, total, initial):
		"""Test passing a factory function as argument."""

		def factory(total, *, initial=None, **kw):
			return TestProgressMeter.create(total, initial=initial, foo=1, **kw)

		meter = get_progress(factory, total, initial=initial, bar=2)

		assert isinstance(meter, TestProgressMeter)
		assert meter.total == total
		assert meter.n == initial
		assert meter.kw == dict(foo=1, bar=2)

	def test_progressconfig(self, total, initial):
		"""Test passing a factory function as argument."""

		config = TestProgressMeter.config(foo=1)
		meter = get_progress(config, total, initial=initial, bar=2)

		assert isinstance(meter, TestProgressMeter)
		assert meter.total == total
		assert meter.n == initial
		assert meter.kw == dict(foo=1, bar=2)

	def test_invalid(selfj):
		with pytest.raises(TypeError):
			get_progress(0, 100)


def test_capture_progress():
	"""Test the capture_progress() function."""
	config, instances = capture_progress(TestProgressMeter)
	assert instances == []

	instance1 = config.create(100)
	assert instances == [instance1]

	instance2 = config.create(100)
	assert instances == [instance1, instance2]

	instance3 = config.create(100)
	assert instances == [instance1, instance2, instance3]


@pytest.mark.parametrize('pass_total', [False, True])
@pytest.mark.parametrize('abort_early', [False, True])
def test_iter_progress(pass_total, abort_early):
	"""Test the iter_progress() function."""
	items = ascii_letters
	abort_at = 10

	if pass_total:
		iterable = iter(items)
		total = len(items)
	else:
		iterable = items
		total = None

	with iter_progress(iterable, TestProgressMeter, total=total, foo=1) as itr:
		assert isinstance(itr.meter, TestProgressMeter)
		assert itr.meter.total == len(items)
		assert itr.meter.kw == dict(foo=1)
		assert itr.meter.n == 0
		assert not itr.meter.closed

		for i, val in enumerate(itr):
			assert val == items[i]
			assert itr.meter.n == i
			assert not itr.meter.closed

			if abort_early and i == abort_at:
				break

		if abort_early:
			assert i == abort_at
			assert itr.meter.n == abort_at
			assert not itr.meter.closed
		else:
			assert i == len(items) - 1
			assert itr.meter.n == len(items)
			assert itr.meter.closed

	assert itr.meter.closed  # Always closed after exiting context


def test_check_progress():
	"""Test the check_progress function."""

	with check_progress() as pconf:
		with get_progress(pconf, 100) as meter:
			meter.moveto(100)

	with check_progress(total=100) as pconf:
		with get_progress(pconf, 100) as meter:
			meter.moveto(100)

	with check_progress(check_closed=False) as pconf:
		meter = get_progress(pconf, 100)
		meter.moveto(100)

	# Not completed
	with pytest.raises(AssertionError):
		with check_progress() as pconf:
			with get_progress(pconf, 100) as meter:
				meter.moveto(99)

	# Wrong total
	with pytest.raises(AssertionError):
		with check_progress(total=100) as pconf:
			with get_progress(pconf, 10) as meter:
				meter.moveto(10)

	# Not closed
	with pytest.raises(AssertionError):
		with check_progress() as pconf:
			meter = get_progress(pconf, 100)
			meter.moveto(100)

	# Not instantiated
	with pytest.raises(AssertionError):
		with check_progress():
			pass

	# Instantiated multiple times
	with pytest.raises(AssertionError):
		with check_progress() as pconf:
			with get_progress(pconf, 100) as meter1:
				meter1.moveto(100)
			with get_progress(pconf, 100) as meter2:
				meter2.moveto(100)


def test_NullProgressMeter():
	"""Test the NullProgressMeter class."""

	# All methods are no-ops so just test we can call interface funcs with no errors.
	meter = NullProgressMeter()
	meter.increment()
	meter.increment(10)
	meter.moveto(100)
	meter.close()

	# Accepts standard arguments but ignores them
	assert isinstance(meter.create(100), NullProgressMeter)
	assert isinstance(meter.create(100, initial=10), NullProgressMeter)
	assert isinstance(meter.create(100, foo=10), NullProgressMeter)


class TestTestProgressMeter:
	"""Test the TestProgressMeter class."""

	def test_basic(self):
		kw = dict(foo=1, bar=True)
		pbar = TestProgressMeter(100, **kw)

		assert pbar.total == 100
		assert pbar.n == 0
		assert not pbar.closed
		assert pbar.kw == kw

		pbar.increment()
		assert pbar.n == 1

		pbar.increment(10)
		assert pbar.n == 11

		pbar.increment(-1)
		assert pbar.n == 10

		pbar.moveto(50)
		assert pbar.n == 50

		pbar.moveto(40)
		assert pbar.n == 40

		with pytest.raises(ValueError):
			pbar.increment(100)

		with pytest.raises(ValueError):
			pbar.increment(-100)

		with pytest.raises(ValueError):
			pbar.moveto(101)

		with pytest.raises(ValueError):
			pbar.moveto(-1)

		pbar.close()
		assert pbar.closed

		with pytest.raises(RuntimeError):
			pbar.increment()

		with pytest.raises(RuntimeError):
			pbar.moveto(100)

	def test_no_allow_decrement(self):
		pbar = TestProgressMeter(100, allow_decrement=False)

		# Moving forward
		pbar.increment()
		pbar.increment(0)
		pbar.increment(10)
		pbar.moveto(50)

		# Moving backward
		with pytest.raises(ValueError):
			pbar.increment(-1)

		with pytest.raises(ValueError):
			pbar.moveto(40)


class TestClickProgressMeter:
	"""Test the ClickProgressMeter class."""
	# TODO


class TestTqdmProgressMeter:
	"""Test the TqdmProgressMeter class."""
	# TODO
