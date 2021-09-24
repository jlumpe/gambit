"""Load database data given file paths."""

import os
from pathlib import Path
from typing import Tuple

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .models import only_genomeset
from .gambitdb import GAMBITDatabase
from .sqla import ReadOnlySession
from gambit.io.util import FilePath


def file_sessionmaker(path: FilePath, class_=ReadOnlySession, **kw) -> sessionmaker:
	"""Get an SQLAlchemy ``sessionmaker`` for an sqlite database file.

	Parameters
	----------
	path
		Path to database file.
	class_
		SQLAlchemy ``Session`` subclass to use. Defaults to :class:`gambit.db.sqla.ReadOnlySession`.
	\\**kw
		Additional keyword arguments to :class:`sqlalchemy.orm.sessionmaker`.
	"""
	engine = create_engine(f'sqlite:///{os.fspath(path)}')
	return sessionmaker(engine, class_=class_, **kw)


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

	session = file_sessionmaker(genomes_file)()
	gset = only_genomeset(session)
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
