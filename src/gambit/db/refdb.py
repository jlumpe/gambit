from pathlib import Path
from typing import Tuple, Sequence, Union, List, Dict, Optional, Any

from sqlalchemy.orm import object_session, Session
from sqlalchemy.orm.attributes import InstrumentedAttribute

from .models import ReferenceGenomeSet, AnnotatedGenome, Genome, only_genomeset
from .sqla import file_sessionmaker
from gambit.sigs.base import ReferenceSignatures, load_signatures
from gambit.util.io import FilePath


# Type alias for argument specifying genome id attribute
GenomeAttr = Union[str, InstrumentedAttribute]


class DatabaseLoadError(Exception):
	"""Raised when there is a problem loading a database.

	Attributes
	----------
	msg
	directory
		Directory we're attempting to load from.
	genomes_file
	signatures_file
	"""

	msg: str
	directory: Optional[Path]
	genomes_file: Optional[Path]
	signatures_file: Optional[Path]

	def __init__(self, msg, directory=None, genomes_file=None, signatures_file=None):
		self.msg = msg
		self.directory = directory
		self.genomes_file = genomes_file
		self.signatures_file = signatures_file


def load_genomeset(db_file: FilePath) -> Tuple[Session, ReferenceGenomeSet]:
	"""Get the only :class:`gambit.db.models.ReferenceGenomeSet` from a genomes database file."""
	session = file_sessionmaker(db_file)()
	gset = only_genomeset(session)
	return session, gset


def _check_genome_id_attr(attr: GenomeAttr) -> InstrumentedAttribute:
	"""Check that Genome ID attribute is valid, and convert from string argument.
	"""
	if isinstance(attr, str) and attr in Genome.ID_ATTRS:
		return getattr(Genome, attr)

	elif isinstance(attr, InstrumentedAttribute):
		for allowed_name in Genome.ID_ATTRS:
			allowed = getattr(Genome, allowed_name)
			if attr is allowed:
				return attr

	raise ValueError('Genome ID attribute must be one of the following: ' + ', '.join(Genome.ID_ATTRS))


def _get_genome_id(genome: Union[Genome, AnnotatedGenome], attr: InstrumentedAttribute):
	"""Get value of ID attribute for genome."""
	if isinstance(genome, AnnotatedGenome):
		genome = genome.genome
	return attr.__get__(genome, Genome)


def _check_genomes_have_ids(genomeset: ReferenceGenomeSet, id_attr: InstrumentedAttribute):
	"""Check all genomes in ReferenceGenomeSet have values for the given ID attribute or raise a ``RuntimeError``."""
	c = genomeset.genomes \
		.join(AnnotatedGenome.genome) \
		.filter(id_attr == None) \
		.count()

	if c > 0:
		raise RuntimeError(f'{c} genomes missing value for ID attribute {id_attr.key}')


def _map_ids_to_genomes(genomeset: ReferenceGenomeSet, id_attr: Union[str, InstrumentedAttribute]) -> Dict[AnnotatedGenome, Any]:
	"""Get dict mapping ID values to AnnotatedGenome."""
	q = genomeset.genomes.join(AnnotatedGenome.genome).add_columns(id_attr)
	return {id_: g for g, id_ in q}


def genomes_by_id(genomeset: ReferenceGenomeSet, id_attr: GenomeAttr, ids: Sequence, strict: bool = True) -> List[Optional[AnnotatedGenome]]:
	"""Match a :class:`ReferenceGenomeSet`'s genomes to a set of ID values.

	This is primarily used to match genomes to signatures based on the ID values stored in a
	signature file. It is expected that the signature file may contain signatures for more genomes
	than are present in the genome set, see also :func:`.genomes_by_id_subset` for that condition.

	Parameters
	----------
	genomeset
	id_attr
		ID attribute of :class:`gambit.db.models.Genome` to use for lookup. Can be used as the
		attribute itself (e.g. ``Genome.refseq_acc``) or just the name (``'refsec_acc'``).
		See :data:`.GENOME_IDS` for the set of allowed values.
	ids
		Sequence of ID values (strings or integers, matching type of attribute).
	strict
		Raise an exception if a matching genome cannot be found for any ID value.

	Returns
	-------
	List[Optional[AnnotatedGenome]]
		List of genomes of same length as ``ids``. If ``strict=False`` and a genome cannot be found
		for a given ID the list will contain ``None`` at the corresponding position.

	Raises
	------
	KeyError
		If ``strict=True`` and any ID value cannot be found.
	"""
	id_attr = _check_genome_id_attr(id_attr)
	_check_genomes_have_ids(genomeset, id_attr)
	d = _map_ids_to_genomes(genomeset, id_attr)
	if strict:
		return [d[id_] for id_ in ids]
	else:
		return [d.get(id_) for id_ in ids]


def genomes_by_id_subset(genomeset: ReferenceGenomeSet,
                         id_attr: GenomeAttr,
                         ids: Sequence,
                         ) -> Tuple[List[AnnotatedGenome], List[int]]:
	"""Match a :class:`ReferenceGenomeSet`'s genomes to a set of ID values, allowing missing genomes.

	This calls :func:`.genomes_by_id` with ``strict=False`` and filters any ``None`` values from the
	output. The filtered list is returned along with the indices of all values in ``ids`` which were
	not filtered out. The indices can be used to load only those signatures which have a matched
	genome from a signature file.

	Note that it is not checked that every genome in ``genomeset`` is matched by an ID. Check the
	size of the returned lists for this.

	Parameters
	----------
	genomeset
	id_attr
		ID attribute of :class:`gambit.db.models.Genome` to use for lookup. Can be used as the
		attribute itself (e.g. ``Genome.refseq_acc``) or just the name (``'refsec_acc'``).
		See :data:`.GENOME_IDS` for the set of allowed values.
	ids
		Sequence of ID values (strings or integers, matching type of attribute).

	Returns
	-------
	Tuple[List[AnnotatedGenome], List[int]]
	"""
	genomes = genomes_by_id(genomeset, id_attr, ids, strict=False)
	genomes_out = []
	idxs_out = []

	for i, g in enumerate(genomes):
		if g is not None:
			genomes_out.append(g)
			idxs_out.append(i)

	return genomes_out, idxs_out


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

	@classmethod
	def locate_files(cls, path: FilePath) -> Tuple[Path, Path]:
		"""Locate an SQLite genome database file and HDF5 signatures file in a directory.

		Files are located by extension, ``.gdb`` or ``.db`` for SQLite file and ``.gs`` or ``.h5``
		for signatures file. Does not look in subdirectories.

		Parameters
		----------
		path
			Path to directory to look within.

		Returns
		-------
			Paths to genomes database file and signatures file.

		Raises
		------
		FileNotFoundError
			If the path does not exist.
		NotADirectoryError
			If the given path does not point to a directory.
		DatabaseLoadError
			If files could not be located or if multiple files with the same extension exist in the
			directory.
		"""
		path = Path(path)

		def check_single_match(matches, desc: str):
			n = len(matches)
			if n != 1:
				raise DatabaseLoadError(
					f'{"Multiple" if n else "No"} {desc} files found in directory {path}',
					directory=path,
				)

		genomes_matches = {f for f in path.iterdir() if f.suffix in ('.gdb', '.db')}
		check_single_match(genomes_matches, 'genome database (.gdb or .db)')
		genomes_file = genomes_matches.pop()

		signatures_matches = {f for f in path.iterdir() if f.suffix in ('.gs', '.h5')}
		check_single_match(signatures_matches, 'signature (.gs or .h5)')
		signatures_file = signatures_matches.pop()

		return genomes_file, signatures_file

	@classmethod
	def load(cls, genomes_file: FilePath, signatures_file: FilePath) -> 'ReferenceDatabase':
		"""Load complete database given paths to SQLite genomes database file and HDF5 signatures file."""
		session, gset = load_genomeset(genomes_file)
		sigs = load_signatures(signatures_file)
		return cls(gset, sigs)

	@classmethod
	def load_from_dir(cls, path: FilePath) -> 'ReferenceDatabase':
		"""
		Load complete database given directory containing SQLite genomes database file and HDF5
		signatures file.

		See :func:`.locate_db_files` for how these files are located within the directory.

		Raises
		------
		RuntimeError
			If files cannot be located in directory.
		"""
		genomes_file, signatures_file = cls.locate_files(path)
		return cls.load(genomes_file, signatures_file)
