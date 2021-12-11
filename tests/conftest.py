from pathlib import Path
from csv import DictReader
import sqlite3

import numpy as np
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


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
def testdb_files(test_data):
	"""Paths to testdb_210818 files."""
	root = test_data / 'testdb_210818'
	return dict(
		root=root,
		ref_genomes=root / 'testdb_210818-genomes.db',
		ref_signatures=root / 'testdb_210818-signatures.h5',
		queries_table=root / 'queries/queries.csv',
		query_genomes=root / 'queries/genomes/',
		query_signatures=root / 'queries/query-signatures.h5',
		results=root / 'results/',
	)

@pytest.fixture(scope='session')
def testdb_engine(testdb_files):
	"""SQLAlchemy engine connected to test database."""
	return create_engine('sqlite:///' + str(testdb_files['ref_genomes']))

@pytest.fixture(scope='session')
def testdb_session(testdb_engine):
	"""Factory function which creates a new session for the test database."""
	from gambit.db.sqla import ReadOnlySession
	return sessionmaker(testdb_engine, class_=ReadOnlySession)

@pytest.fixture(scope='session')
def testdb_copy(testdb_files):
	"""Factory function which creates an in-memory copy of the test database."""

	def make_testdb_copy():
		src = sqlite3.connect(str(testdb_files['ref_genomes']))
		memory = sqlite3.connect(':memory:')
		src.backup(memory)
		engine = create_engine('sqlite://', creator=lambda: memory)
		return sessionmaker(engine)()

	return make_testdb_copy

@pytest.fixture(scope='session')
def testdb_signatures(testdb_files):
	"""K-mer signatures for testdb reference genomes."""
	from gambit.sigs.hdf5 import HDF5Signatures

	return HDF5Signatures.open(testdb_files['ref_signatures'])

@pytest.fixture(scope='session')
def testdb_query_signatures(testdb_files):
	"""K-mer signatures for testdb query genomes."""
	from gambit.sigs.hdf5 import HDF5Signatures

	return HDF5Signatures.open(testdb_files['query_signatures'])

@pytest.fixture(scope='session')
def testdb(testdb_session, testdb_signatures):
	"""Full ReferenceDatabase object for test db."""
	from gambit.db import ReferenceDatabase, ReferenceGenomeSet

	session = testdb_session()
	gset = session.query(ReferenceGenomeSet).one()
	return ReferenceDatabase(gset, testdb_signatures)

@pytest.fixture(scope='session')
def testdb_queries(testdb, testdb_files):
	"""Query files and their expected results."""
	from gambit.seq import SequenceFile

	genomes_dir = testdb_files['query_genomes']

	with open(testdb_files['queries_table'], newline='') as f:
		rows = list(DictReader(f))

	for row in rows:
		row['warnings'] = row['warnings'].lower() == 'true'
		row['file'] = SequenceFile(
			path=genomes_dir / (row['name'] + '.fasta'),
			format='fasta',
		)

	return rows

@pytest.fixture(scope='session', params=['non_strict', 'strict'])
def testdb_results(request, testdb_files, testdb_session):
	"""Pre-calculated query results.

	Use a yield statement here instead of a return, we want to keep a reference
	to the session object until teardown or else it may be garbage collected,
	which would render any ORM instances in the results object invalid.
	"""
	from gambit.io.results.archive import ResultsArchiveReader

	session = testdb_session()
	reader = ResultsArchiveReader(session)

	path = testdb_files['results'] / (request.param + '.json')
	yield reader.read(path)
