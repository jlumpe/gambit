from typing import Optional, List, Dict, Any
from pathlib import Path

import click
from attr import attrs, attrib
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

from gambit.db import locate_db_files, ReferenceGenomeSet
from gambit.db.sqla import ReadOnlySession
from gambit.signatures.hdf5 import HDF5Signatures
from gambit.io.seq import SequenceFile


@attrs
class CLIContext:
	"""Click context object for GAMBIT CLI.

	Attributes
	----------
	db_path
		Path to directory containing database files, specified in root command group.
	"""
	db_path: Optional[Path] = attrib(default=None)

	def __attrs_post_init__(self):

		self._db_found = False
		self._genomes_path = None
		self._signatures_path = None
		self._engine = None
		self._session = None
		self._Session = None
		self._gset = None
		self._signatures = None

	def _require_db(self):
		if self._db_found:
			return

		if self.db_path is None:
			raise click.ClickException('Must supply path to database directory.')

		self._genomes_path, self._signatures_path = locate_db_files(self.db_path)
		self._db_found = True

	def _require_genomes(self):
		if self._engine is not None:
			return

		self._require_db()

		self._engine = create_engine(f'sqlite:///{self._genomes_path}')
		self._Session = sessionmaker(self.engine(), class_=ReadOnlySession)
		self._session = self._Session()

	def engine(self) -> Engine:
		"""SQLAlchemy engine connecting to database."""
		self._require_genomes()
		return self._engine

	def session(self) -> ReadOnlySession:
		"""Create a new SQLAlchemy session for the database."""
		self._require_genomes()
		return self._session

	def genomeset(self) -> ReferenceGenomeSet:
		if self._gset is None:
			self._require_genomes()
			self._gset = self._session.query(ReferenceGenomeSet).one()

		return self._gset

	def signatures(self) -> HDF5Signatures:
		if self._signatures is None:
			self._require_db()
			self._signatures = HDF5Signatures.open(self._signatures_path)

		return self._signatures


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
