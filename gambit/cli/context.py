from typing import Optional
from pathlib import Path

import click
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

from gambit.db.fromfile import locate_db_files
from gambit.db.models import ReferenceGenomeSet
from gambit.db.sqla import ReadOnlySession
from gambit.signatures.hdf5 import HDF5Signatures


class CLIContext:
	"""Click context object for GAMBIT CLI.

	Attributes
	----------
	db_path
		Path to directory containing database files, specified in root command group.
	"""
	db_path: Optional[Path]

	def __init__(self, db_path):
		self.db_path = db_path

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
