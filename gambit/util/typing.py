"""Utilities based on the built-in ``typing`` module."""

import typing
from typing import Union, Any


def is_union(T) -> bool:
	"""Check if a type annotation is a *parameterized* :class:`typing.Union`.

	Parameters
	----------
	T
		Result of ``Union[A, B, ...]``.
	"""
	return isinstance(T, typing._GenericAlias) and T.__origin__ is typing.Union


def union_types(T) -> tuple:
	"""Get the types from a parameterized :class:`typing.Union`.

	Parameters
	----------
	T
		Result of ``Union[A, B, ...]``.
	"""
	return T.__args__


def is_optional(T) -> bool:
	"""Check if a parametrized union type is equivalent to one returned by :data:`typing.Optional`."""
	if not is_union(T):
		return False
	types = union_types(T)
	return len(types) == 2 and type(None) in types


def unwrap_optional(u):
	"""Get ``T`` from ``typing.Optional[T]``."""
	for T in union_types(u):
		if T is not type(None):
			return T

	raise ValueError(f'Not an Optional type: {u!r}')
