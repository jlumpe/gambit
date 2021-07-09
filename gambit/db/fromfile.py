"""Load database data given file paths."""

import os
from pathlib import Path
from typing import Tuple

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import MultipleResultsFound, NoResultFound

from .models import ReferenceGenomeSet
from .gambitdb import GAMBITDatabase
from .sqla import ReadOnlySession
from gambit.io.util import FilePath


def open_genomeset(path: FilePath, session_cls=ReadOnlySession) -> ReferenceGenomeSet:
	"""Open an SQLite file containing GAMBIT reference genomes and get its :class:`ReferenceGenomeSet`.

	Parameters
	----------
	path
		Path to database file.
	session_cls
		:class:`sqlalchemy.orm.Session` subclass to use.

	Raises
	------
	RuntimeError
		If the database file does not contain a single genome set.
	"""
	engine = create_engine(f'sqlite:///{os.fspath(path)}')
	session = sessionmaker(engine, class_=session_cls)()

	try:
		return session.query(ReferenceGenomeSet).one()
	except MultipleResultsFound as e:
		raise RuntimeError('Database file contains multiple genome sets.') from e
	except NoResultFound as e:
		raise RuntimeError('Database file contains no genome sets.') from e


def locate_db_files(path: FilePath) -> Tuple[Path, Path]:
	"""Locate an SQLite genome database file and HDF5 signatures file in a directory.

	Files are located by extension, ``.db`` for SQLite file and ``.h5`` for signatures file.
	Does not look in subdirectories.

	Parameters
	----------
	path
		Path to directory to look within.

	Returns
	-------
		Paths to genomes database file and signatures file.

	Raises
	------
	RuntimeError
		If files could not be located or if multiple files with the same extension exist in the
		directory.
	"""
	path = Path(path)

	genomes_matches = list(path.glob('*.db'))
	if len(genomes_matches) == 0:
		raise RuntimeError(f'No genome database (.db) files found in directory {path}')
	if len(genomes_matches) > 1:
		raise RuntimeError(f'Multiple genome database (.db) files found in directory {path}')

	signatures_matches = list(path.glob('*.h5'))
	if len(signatures_matches) == 0:
		raise RuntimeError(f'No signature (.h5) files found in directory {path}')
	if len(signatures_matches) > 1:
		raise RuntimeError(f'Multiple signature (.h5) files found in directory {path}')

	return genomes_matches[0], signatures_matches[0]


def load_database(genomes_file: FilePath, signatures_file: FilePath) -> GAMBITDatabase:
	"""Load complete database given paths to SQLite genomes database file and HDF5 signatures file."""
	from gambit.signatures.hdf5 import HDF5Signatures
	gset = open_genomeset(genomes_file)
	sigs = HDF5Signatures.open(signatures_file)
	return GAMBITDatabase(gset, sigs)


def load_database_from_dir(path: FilePath) -> GAMBITDatabase:
	"""
	Load complete database given directory containing SQLite genomes database file and HDF5
	signatures file.

	See :func:`.locate_db_files` for how these files are located within the directory.

	Raises
	------
	RuntimeError
		If files cannot be located in directory.
	"""
	genomes_file, signatures_file = locate_db_files(path)
	return load_database(genomes_file, signatures_file)
