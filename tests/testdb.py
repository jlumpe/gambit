"""Access test database data.
"""

from typing import Callable, TypeVar, Generic, Any, overload, TypedDict
from pathlib import Path
from dataclasses import dataclass
from csv import DictReader
import sqlite3
import gzip

from sqlalchemy import create_engine

from gambit.kmers import KmerSpec
from gambit.sigs import load_signatures, AnnotatedSignatures
from gambit.db import ReferenceDatabase, only_genomeset, file_sessionmaker, default_sessionmaker
from gambit.results import ResultsArchiveReader
from gambit.query import QueryResults
from gambit.util.io import FilePath


T = TypeVar('T')


class LazyAttribute(Generic[T]):
	"""Descriptor which initializes a property value the first time it is used."""

	def __init__(self, initializer: Callable[[Any], T], value_attr: str):
		self.initializer = initializer
		self.value_attr = value_attr
		self.__doc__ = initializer.__doc__

	@overload
	def __get__(self, instance: None, owner=None) -> 'LazyAttribute[T]':
		pass

	@overload
	def __get__(self, instance, owner=None) -> T:
		pass

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


def lazy(f: Callable[[Any], T]) -> LazyAttribute[T]:
	attr = '_' + f.__name__
	return LazyAttribute(f, attr)


@dataclass
class TestDBPaths:
	root: Path
	# Reference genomes .gdb file
	ref_genomes: Path
	# Reference genomes .gs file
	ref_signatures: Path
	# Reference genomes .csv
	refs_table: Path
	# Directory containing reference genome FASTA files
	ref_genomes_dir: Path
	# queries.csv
	queries_table: Path
	# Directory containing query genome FASTA files
	query_genomes_dir: Path
	# Query genomes .gs file
	query_signatures: Path
	# Directory containing QueryResults exports in archive format.
	results: Path


class TestQueryGenome(TypedDict):
	name: str
	predicted: str
	primary: str
	closest: str
	warnings: bool


class TestRefGenome(TypedDict):
	name: str
	key: str
	taxon: str


class TestDB:
	"""Object which provides access to test database resources.

	Many attributes are "lazy", meaning they are not initialized until first used. This is similar
	to how it would work if the attributes were separate Pytest fixtures.
	"""

	paths: TestDBPaths

	# Prevent pytest interpreting as containing test methods
	__test__ = False

	def __init__(self, root: FilePath):
		root = Path(root)
		self.paths = TestDBPaths(
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
	def Session(self):
		"""Sessionmaker for the reference genome database."""
		return file_sessionmaker(self.paths.ref_genomes)

	def copy_session(self):
		"""Create an in-memory copy of the test database."""
		src = sqlite3.connect(str(self.paths.ref_genomes))
		memory = sqlite3.connect(':memory:')
		src.backup(memory)
		engine = create_engine('sqlite://', creator=lambda: memory)
		return default_sessionmaker(engine)()

	@lazy
	def ref_signatures(self) -> AnnotatedSignatures:
		"""K-mer signatures for reference genomes."""
		return load_signatures(self.paths.ref_signatures)  # type: ignore

	@lazy
	def query_signatures(self) -> AnnotatedSignatures:
		"""K-mer signatures for query genomes."""
		return load_signatures(self.paths.query_signatures)  # type: ignore

	@lazy
	def kmerspec(self) -> KmerSpec:
		return self.ref_signatures.kmerspec  # type: ignore

	@lazy
	def refdb(self) -> ReferenceDatabase:
		"""Full ReferenceDatabase object."""
		session = self.Session()
		gset = only_genomeset(session)
		return ReferenceDatabase(gset, self.ref_signatures)

	@lazy
	def query_genomes(self) -> list[TestQueryGenome]:
		"""Query genomes and their expected results."""

		with open(self.paths.queries_table, newline='') as f:
			rows = list(DictReader(f))

		for row in rows:
			# Convert "warnings" column to bool
			row['warnings'] = row['warnings'].lower() == 'true'

		return rows  # type: ignore

	@lazy
	def ref_genomes(self) -> list[TestRefGenome]:
		"""Reference genomes and their attributes."""

		with open(self.paths.refs_table, newline='') as f:
			rows = list(DictReader(f))

		return rows  # type: ignore

	@classmethod
	def _ensure_gz(cls, file: Path, file_gz: Path):
		"""Ensure gzipped version of the query/ref file is available.

		These aren't added to version control, so they are created the first time they are needed.
		"""
		if file_gz.is_file():
			return

		with open(file) as f:
			content = f.read()

		with gzip.open(file_gz, 'wt') as f:
			f.write(content)

	def _get_genome_files(self, base: Path, names: list[str], gzipped: bool, relative: bool) -> list[Path]:
		base2 = base.relative_to(self.paths.root) if relative else base

		files = []

		for name in names:
			fname = name + '.fasta'

			if gzipped:
				fname_gz = fname + '.gz'
				self._ensure_gz(base / fname, base / fname_gz)
				path = base2 / fname_gz
			else:
				path = base2 / fname

			files.append(path)

		return files

	def get_query_files(self, gzipped: bool = False, relative: bool = False) -> list[Path]:
		return self._get_genome_files(
			self.paths.query_genomes_dir,
			[genome['name'] for genome in self.query_genomes],
			gzipped=gzipped,
			relative=relative,
		)

	def get_ref_files(self, gzipped: bool = False, relative: bool = False) -> list[Path]:
		return self._get_genome_files(
			self.paths.ref_genomes_dir,
			[genome['name'] for genome in self.ref_genomes],
			gzipped=gzipped,
			relative=relative,
		)

	def get_query_results(self, strict: bool, session=None) -> QueryResults:
		"""Pre-calculated query results."""
		if session is None:
			session = self.refdb.session
		reader = ResultsArchiveReader(session)
		path = self.paths.results / ('strict.json' if strict else 'non_strict.json')
		return reader.read(path)
