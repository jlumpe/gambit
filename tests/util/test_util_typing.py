"""Test the gambit.util.typing submodule."""

import typing
from typing import Union, Optional

from gambit.util.typing import is_union, union_types, is_optional, unwrap_optional


def test_is_union():
	"""Test the is_union() function."""
	assert is_union(Union[int, str])
	assert is_union(Union[int, str, bool])
	assert is_union(Optional[int])
	assert not is_union(Union)
	assert not is_union(Optional)
	assert not is_union(int)
	assert not is_union(None)
	assert not is_union(typing.List)
	assert not is_union(typing.Any)


def test_union_types():
	"""Test the union_types() function."""
	assert union_types(Union[int, str]) == (int, str)
	assert union_types(Union[int, str, bool]) == (int, str, bool)
	assert union_types(Optional[int]) == (int, type(None))


def test_is_optional():
	"""Test the is_optional() function."""
	assert is_optional(Optional[int])
	assert is_optional(Union[int, None])
	assert is_optional(Union[None, int])
	assert not is_optional(Union[int, str])
	assert not is_optional(Union)
	assert not is_optional(Optional)
	assert not is_optional(int)
	assert not is_optional(None)
	assert not is_optional(type(None))
	assert not is_optional(typing.Any)


def test_unwrap_optional():
	"""Test the unwrap_optional() function."""
	assert unwrap_optional(Optional[int]) is int
	assert unwrap_optional(Union[int, None]) is int
	assert unwrap_optional(Union[None, int]) is int
