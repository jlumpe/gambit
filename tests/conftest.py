from pathlib import Path
from urllib.request import urlretrieve
from warnings import warn

import numpy as np
import pytest
from sqlalchemy import create_engine

from .testdb import TestDB


REFSEQ_SIGNATURES_URL = (
	'https://storage.googleapis.com/jlumpe-gambit/public/databases/refseq-curated/1.0/'
	'gambit-refseq-curated-1.0.gs'
)


@pytest.fixture(scope='session')
def test_data():
	"""The directory containing test data."""
	return Path(__file__).parent / 'data'


@pytest.fixture(scope='session')
def refseq_db_signatures_file(test_data) -> Path:
	"""Path to the full curated RefSeq database signatures file.

	Downloads the file from the Database Releases page if it is not already cached locally.
	"""
	filename = REFSEQ_SIGNATURES_URL.split('/')[-1]
	path = (test_data / 'cache' / filename).resolve()

	if not path.is_file():
		warn(f'{path} does not exist, downloading...')
		path.parent.mkdir(parents=True, exist_ok=True)
		tmp_path = path.with_name(path.name + '.tmp')
		urlretrieve(REFSEQ_SIGNATURES_URL, tmp_path)
		tmp_path.rename(path)

	return path


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
