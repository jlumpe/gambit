"""Test gambit.io.json."""

from pathlib import Path
from datetime import date, datetime

import pytest
import numpy as np
from attr import attrs, attrib

import gambit.io.json as mjson


def roundtrip(obj, cls=None, checktype=True):
	"""Convert to JSON string and back again."""
	if cls is None:
		cls = type(obj)

	encoded = mjson.dumps(obj)
	obj2 = mjson.loads(encoded, cls)

	if checktype:
		assert isinstance(obj2, cls)

	return obj2


class TestBuiltins:
	"""Test conversion of builtin types."""

	def test_path(self):
		p = Path('foo/bar/baz')
		assert roundtrip(p) == p

	def test_date(self):
		d = date.today()
		assert roundtrip(d) == d

	def test_datetime(self):
		dt = datetime.now()
		assert roundtrip(dt) == dt


class TestNumpy:
	"""Test serialization of numpy types."""

	def test_scalar(self):
		assert mjson.dumps(np.int64(1)) == '1'
		assert mjson.dumps(np.float64(1)) == '1.0'

	def test_array(self):
		# TODO
		pass


@pytest.mark.parametrize('custom_to', [False, True])
@pytest.mark.parametrize('custom_from', [False, True])
def test_jsonable(custom_to, custom_from):
	"""Test the Jsonable mixin class."""

	@attrs()
	class TestCls(mjson.Jsonable):
		a: int = attrib()
		b: str = attrib()

	if custom_to:
		TestCls.__to_json__ = lambda self: dict(a=self.a + 1, b=self.b.upper())

	if custom_from:
		TestCls.__from_json__ = classmethod(lambda cls, data: cls(data['a'] - 1, data['b'].lower()))

	# Nested as attribute of another class
	@attrs()
	class TestCls2:
		x: TestCls = attrib()
		y: int = attrib()

	obj1 = TestCls(3, 'foo')
	obj2 = TestCls2(obj1, 10)

	data_standard = dict(a=3, b='foo')
	data_custom = dict(a=4, b='FOO')

	data_to = data_custom if custom_to else data_standard
	assert mjson.to_json(obj1) == data_to
	assert mjson.to_json(obj2) == dict(x=data_to, y=10)

	data_from = data_custom if custom_from else data_standard
	assert mjson.from_json(data_from, TestCls) == obj1
	assert mjson.from_json(dict(x=data_from, y=10), TestCls2) == obj2
