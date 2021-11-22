from pathlib import Path
from typing import Tuple, Sequence

from sqlalchemy.orm import object_session

from .models import ReferenceGenomeSet, AnnotatedGenome, genomes_by_id_subset, only_genomeset
from .sqla import file_sessionmaker
from gambit.sigs.meta import ReferenceSignatures
from gambit.io.util import FilePath


class ReferenceDatabase:
	"""Object containing reference genomes, their k-mer signatures, and associated data.

	This is all that is needed at runtime to run queries.

	Attributes
	----------
	genomeset
		Genome set containing reference genomes.
	genomes
		List of reference genomes.
	signatures
		K-mer signatures for each genome. A subtype of ``ReferenceSignatures``, so contains metadata
		on signatures as well as the signatures themselves. Type may represent signatures stored on
		disk (e.g. :class:`HDF5Signatures`) instead of in memory. OK to contain additional
		signatures not corresponding to any genome in ``genomes``.
	sig_indices
		Index of signature in ``signatures`` corresponding to each genome in ``genomes``.
		In sorted order to improve performance when iterating over them (improve locality if in
		memory and avoid seeking if in file).
	session
		The SQLAlchemy session ``genomeset`` and the elements of ``genomes`` belong to.
		It is important to keep a reference to this, just having references to the ORM objects
		themselves is not enough to keep the session from being garbage collected.

	Parameters
	----------
	genomeset
	signatures
	"""
	genomeset: ReferenceGenomeSet
	genomes: Sequence[AnnotatedGenome]
	signatures: ReferenceSignatures
	sig_indices: Sequence[int]

	def __init__(self, genomeset: ReferenceGenomeSet, signatures: ReferenceSignatures):
		self.genomeset = genomeset
		self.signatures = signatures
		self.session = object_session(genomeset)

		id_attr = signatures.meta.id_attr
		if id_attr is None:
			raise TypeError('id_attr field of signatures metadata cannot be None')

		self.genomes, self.sig_indices = genomes_by_id_subset(genomeset, id_attr, signatures.ids)

		n = genomeset.genomes.count()
		if len(self.genomes) != n:
			missing = n - len(self.genomes)
			raise ValueError(f'{missing} of {n} genomes not matched to signature IDs. Is the id_attr attribute of the signatures metadata correct?')


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


def load_database(genomes_file: FilePath, signatures_file: FilePath) -> ReferenceDatabase:
	"""Load complete database given paths to SQLite genomes database file and HDF5 signatures file."""
	from gambit.sigs.hdf5 import HDF5Signatures

	session = file_sessionmaker(genomes_file)()
	gset = only_genomeset(session)
	sigs = HDF5Signatures.open(signatures_file)
	return ReferenceDatabase(gset, sigs)


def load_database_from_dir(path: FilePath) -> ReferenceDatabase:
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
