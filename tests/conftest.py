from pathlib import Path

import numpy as np
import pytest
from sqlalchemy import create_engine

from .testdb import TestDB


@pytest.fixture(scope='session')
def test_data():
	"""The directory containing test data."""
	return Path(__file__).parent / 'data'


@pytest.fixture(autouse=True)
def raise_numpy_errors():
	"""Raise exceptions for all Numpy errors in all tests.

	NOTE: this doesn't affect operations with Numpy scalars and so has limited usefulness.
	"""

	old_settings = np.seterr(all='raise')

	yield

	np.seterr(**old_settings)  # Not really necessary


@pytest.fixture(scope='session')
def make_empty_db():
	"""Function which creates an empty in-memory-database with initialized schema."""
	from gambit.db.models import Base

	def empty_db_factory():
		engine = create_engine('sqlite:///:memory:')
		Base.metadata.create_all(engine)
		return engine

	return empty_db_factory


@pytest.fixture(scope='session')
def testdb(test_data):
	"""Object which facilitates access to testdb_210818 data.

	This cleans things up a bit from the way it was before, which was a bunch of separate fixtures
	with session scope named "testdb_*".
	"""
	return TestDB(test_data / 'testdb_210818')
