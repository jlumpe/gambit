"""Access test database data.
"""

from typing import Callable
from pathlib import Path
from types import SimpleNamespace
from csv import DictReader
import sqlite3
import gzip

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from gambit.seq import SequenceFile
from gambit.sigs import load_signatures
from gambit.db import ReferenceDatabase, ReadOnlySession, only_genomeset
from gambit.results.archive import ResultsArchiveReader


class LazyAttribute:
	"""Descriptor which initializes a property value the first time it is used."""

	def __init__(self, initializer: Callable, value_attr: str):
		self.initializer = initializer
		self.value_attr = value_attr
		self.__doc__ = initializer.__doc__

	def __get__(self, instance, owner=None):
		if instance is None:
			return self

		try:
			return instance.__dict__[self.value_attr]
		except KeyError:
			pass

		value = self.initializer(instance)
		setattr(instance, self.value_attr, value)
		return value


def lazy(f: Callable) -> LazyAttribute:
	attr = '_' + f.__name__
	return LazyAttribute(f, attr)


class TestDB:
	"""Object which provides access to test database resources.

	Many attributes are "lazy", meaning they are not initialized until first used. This is similar
	to how it would work if the attributes were separate Pytest fixtures.
	"""

	def __init__(self, root):
		root = Path(root)
		self.paths = SimpleNamespace(
			root=root,
			ref_genomes=root / 'ref-genomes.gdb',
			ref_signatures=root / 'ref-signatures.gs',
			refs_table=root / 'ref-genomes.csv',
			ref_genomes_dir=root / 'ref-genomes/',
			queries_table=root / 'queries/queries.csv',
			query_genomes_dir=root / 'queries/genomes/',
			query_signatures=root / 'queries/query-signatures.gs',
			results=root / 'results/',
		)

	@lazy
	def engine(self):
		"""SQLAlchemy engine connected to genome database."""
		return create_engine(f'sqlite:///{self.paths.ref_genomes}')

	@lazy
	def Session(self):
		"""Sessionmaker for the reference genome database."""
		return sessionmaker(self.engine, class_=ReadOnlySession)

	def copy_session(self):
		"""Create an in-memory copy of the test database."""
		src = sqlite3.connect(str(self.paths.ref_genomes))
		memory = sqlite3.connect(':memory:')
		src.backup(memory)
		engine = create_engine('sqlite://', creator=lambda: memory)
		return sessionmaker(engine)()

	@lazy
	def ref_signatures(self):
		"""K-mer signatures for reference genomes."""
		return load_signatures(self.paths.ref_signatures)

	@lazy
	def query_signatures(self):
		"""K-mer signatures for query genomes."""
		return load_signatures(self.paths.query_signatures)

	@lazy
	def kmerspec(self):
		return self.ref_signatures.kmerspec

	@lazy
	def refdb(self):
		"""Full ReferenceDatabase object."""
		session = self.Session()
		gset = only_genomeset(session)
		return ReferenceDatabase(gset, self.ref_signatures)

	@classmethod
	def _add_file_cols(cls, genomes_dir, row):
		row['file'] = SequenceFile(
			path=genomes_dir / (row['name'] + '.fasta'),
			format='fasta',
		)
		row['file_gz'] = SequenceFile(
			path=genomes_dir / (row['name'] + '.fasta.gz'),
			format='fasta',
			compression='gzip',
		)

	@lazy
	def query_genomes(self):
		"""Query genomes and their expected results."""

		with open(self.paths.queries_table, newline='') as f:
			rows = list(DictReader(f))

		for row in rows:
			row['warnings'] = row['warnings'].lower() == 'true'
			self._add_file_cols(self.paths.query_genomes_dir, row)

		return rows

	@lazy
	def ref_genomes(self):
		"""Reference genomes and their attributes."""

		with open(self.paths.refs_table, newline='') as f:
			rows = list(DictReader(f))

		for row in rows:
			self._add_file_cols(self.paths.ref_genomes_dir, row)

		return rows

	@classmethod
	def _ensure_gz(cls, items):
		"""Ensure gzipped versions of the query/ref files are available.

		These aren't added to version control, so they are created the first time they are needed.
		"""
		for item in items:
			dst = item['file_gz'].path
			if dst.is_file():
				continue

			with open(item['file'].path) as f:
				content = f.read()

			with gzip.open(dst, 'wt') as f:
				f.write(content)

	@classmethod
	def _get_genome_files(cls, items, gzipped):
		if gzipped:
			col = 'file_gz'
			cls._ensure_gz(items)
		else:
			col = 'file'
		return [q[col] for q in items]

	def get_query_files(self, gzipped: bool=False):
		return self._get_genome_files(self.query_genomes, gzipped)

	def get_ref_files(self, gzipped: bool=False):
		return self._get_genome_files(self.ref_genomes, gzipped)

	def get_query_results(self, strict: bool, session=None):
		"""Pre-calculated query results."""
		if session is None:
			session = self.refdb.session
		reader = ResultsArchiveReader(session)
		path = self.paths.results / ('strict.json' if strict else 'non_strict.json')
		return reader.read(path)
