from pathlib import Path

import numpy as np
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from gambit.db.models import Base as models_base
from gambit.db.sqla import ReadOnlySession


@pytest.fixture(scope='session')
def test_data():
	"""The directory containing test data."""
	return Path(__file__).parent / 'data'


@pytest.fixture(autouse=True)
def raise_numpy_errors():
	"""Raise exceptions for all Numpy errors in all tests."""

	old_settings = np.seterr(all='raise')

	yield

	np.seterr(**old_settings)  # Not really necessary


@pytest.fixture(scope='session')
def make_empty_db():
	"""Function which creates an empty in-memory-database with initialized schema."""
	def empty_db_factory():
		engine = create_engine('sqlite:///:memory:')
		models_base.metadata.create_all(engine)
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
	)

@pytest.fixture(scope='session')
def testdb_engine(testdb_files):
	"""SQLAlchemy engine connected to test database."""
	return create_engine('sqlite:///' + str(testdb_files['ref_genomes']))

@pytest.fixture(scope='session')
def testdb_session(testdb_engine):
	"""Factory function which creates a new session for the test database."""
	return sessionmaker(testdb_engine, class_=ReadOnlySession)

@pytest.fixture(scope='session')
def testdb_copy(testdb_files):
	"""Factory function which creates an in-memory copy of the test database."""
	import sqlite3

	def make_testdb_copy():
		src = sqlite3.connect(str(testdb_files['ref_genomes']))
		memory = sqlite3.connect(':memory:')
		src.backup(memory)
		engine = create_engine('sqlite://', creator=lambda: memory)
		return sessionmaker(engine)()

	return make_testdb_copy

@pytest.fixture(scope='module')
def testdb_signatures(testdb_files):
	"""K-mer signatures for test genomes."""
	from gambit.signatures.hdf5 import HDF5Signatures

	return HDF5Signatures.open(testdb_files['ref_signatures'])

@pytest.fixture(scope='module')
def testdb(testdb_session, testdb_signatures):
	"""Full GAMBITDatabase object for test db."""
	from gambit.db import GAMBITDatabase, ReferenceGenomeSet

	session = testdb_session()
	gset = session.query(ReferenceGenomeSet).one()
	return GAMBITDatabase(gset, testdb_signatures)

@pytest.fixture(scope='module')
def testdb_queries(testdb, testdb_files):
	"""Query files and their expected results."""
	from gambit.io.seq import SequenceFile
	from csv import DictReader

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
