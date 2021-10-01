from typing import Optional, List, Dict, Any, Sequence
from pathlib import Path

import click
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

from gambit.db import locate_db_files, GAMBITDatabase
from gambit.db.models import only_genomeset
from gambit.db.sqla import ReadOnlySession
from gambit.sigs.hdf5 import HDF5Signatures
from gambit.sigs.meta import ReferenceSignatures
from gambit.io.seq import SequenceFile


class CLIContext:
	"""Click context object for GAMBIT CLI.

	Loads reference database data lazily the first time it is requested.

	Currently a single option (or environment variable) is used to specify the location of the
	database files, in the future options may be added to specify the reference genomes SQLite
	file and genome signatures file separately. Class methods treat them as being independent.

	Attributes
	----------
	root_context
		Click context object from root command group.
	db_path
		Path to directory containing database files, specified in root command group.
	has_genomes
		Whether reference genome metadata is available.
	has_signatures
		Whether reference signatures are available.
	has_database
		Whether reference genome metadata and reference signatures are both available.
	engine
		SQLAlchemy engine connecting to genomes database.
	Session
		SQLAlchemy session maker for genomes database.
	signatures
		Reference genome signatures.
	"""
	root_context: click.Context
	db_path: Optional[Path]
	has_genomes: bool
	has_signatures: bool
	has_database: bool
	engine: Optional[Engine]
	Session: Optional[sessionmaker]
	signatures: Optional[ReferenceSignatures]

	def __init__(self, root_context: click.Context):
		"""
		Parameters
		----------
		root_context
			Click context object from root command group.
		"""
		self.root_context = root_context

		db_path = root_context.params['db_path']
		self.db_path = None if db_path is None else Path(db_path)

		self._db_found = False
		self._has_genomes = None
		self._has_signatures = None
		self._signatures_path = None

		self._engine = None
		self._Session = None
		self._signatures = None

	def _find_db(self):
		"""Find database files."""
		if self._db_found:
			return

		if self.db_path is None:
			self._has_genomes = self._has_signatures = False

		else:
			self._has_genomes = self._has_signatures = True
			self._genomes_path, self._signatures_path = locate_db_files(self.db_path)

		self._db_found = True

	@property
	def has_genomes(self):
		if not self._db_found:
			self._find_db()
		return self._has_genomes

	@property
	def has_signatures(self):
		if not self._db_found:
			self._find_db()
		return self._has_signatures

	@property
	def has_database(self):
		return self.has_genomes and self.has_signatures

	def require_database(self):
		"""Raise an exception if genome metadata and signatures are not available."""
		if not self.has_database:
			raise click.ClickException('Must supply path to database directory.')

	def require_genomes(self):
		"""Raise an exception if genome metadata is not available."""
		self.require_database()

	def require_signatures(self):
		"""Raise an exception if signatures are not available."""
		self.require_database()

	def _init_genomes(self):
		if self._engine is not None or not self.has_genomes:
			return

		self._engine = create_engine(f'sqlite:///{self._genomes_path}')
		self._Session = sessionmaker(self.engine, class_=ReadOnlySession)

	@property
	def engine(self):
		self._init_genomes()
		return self._engine

	@property
	def Session(self):
		self._init_genomes()
		return self._Session

	@property
	def signatures(self):
		if self._signatures is None and self.has_signatures:
			self._signatures = HDF5Signatures.open(self._signatures_path)

		return self._signatures

	def get_database(self) -> GAMBITDatabase:
		"""Get reference database object."""
		self.require_database()
		session = self.Session()
		gset = only_genomeset(session)
		return GAMBITDatabase(gset, self.signatures)


def seqfmt_option():
	return click.option(
		'-s', '--seqfmt',
		type=click.Choice(['fasta']),
		default='fasta',
		help='Format of sequence files. Currently only FASTA is supported.',
	)

def genome_files_arg():
	return click.argument(
		'files',
		nargs=-1,
		type=click.Path(exists=True, dir_okay=False),
		metavar='GENOMES...',
	)

def seq_file_params():
	"""Decorator which adds sequence file parameters to command."""

	def decorator(f):
		seqfmt_dec = seqfmt_option()
		genomes_dec = genome_files_arg()
		return seqfmt_dec(genomes_dec(f))

	return decorator

def get_seq_files(params: Dict[str, Any]) -> List[SequenceFile]:
	"""Get list of sequence files from command parameters."""
	return SequenceFile.from_paths(params['files'], params['seqfmt'])


def print_table(rows: Sequence[Sequence], colsep: str = ' ', left: str = '', right: str = ''):
	"""Print a basic table."""

	echo = lambda s: click.echo(s, nl=False)

	rows = [list(map(str, row)) for row in rows]
	ncol = max(map(len, rows))

	widths = [0] * ncol
	for row in rows:
		for i, val in enumerate(row):
			widths[i] = max(widths[i], len(val))

	for row in rows:
		echo(left)

		for i, val in enumerate(row):
			echo(val.ljust(widths[i]))

			if i < ncol - 1:
				echo(colsep)

		echo(right)
		echo('\n')
