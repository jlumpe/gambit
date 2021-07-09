"""SQLAlchemy models for storing reference genomes and taxonomy information."""

from typing import Sequence, Union, Dict, List, Any, Optional, Tuple, Iterable, Collection

import sqlalchemy as sa
from sqlalchemy import Column, Integer, String, Boolean, Float
from sqlalchemy import ForeignKey, UniqueConstraint
from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm.attributes import InstrumentedAttribute

from .sqla import JsonString


__all__ = [
	'Genome',
	'ReferenceGenomeSet',
	'AnnotatedGenome',
	'Taxon',
]


# Naming convention for constraints and indices, used by SQLAlchemy when creating schema.
# Important for Alembic migration scripts, see https://alembic.sqlalchemy.org/en/latest/naming.html
NAMING_CONVENTION = {
  "ix": "ix_%(column_0_label)s",
  "uq": "uq_%(table_name)s_%(column_0_name)s",
  "ck": "ck_%(table_name)s_%(constraint_name)s",
  "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
  "pk": "pk_%(table_name)s",
}

# SqlAlchemy metadata object and declarative base
metadata = sa.MetaData(naming_convention=NAMING_CONVENTION)
Base = declarative_base(metadata=metadata)


class Genome(Base):
	"""Base model for a reference genome that can be compared to query.

	Corresponds to a single assembly (one or more contigs, but at least partially assembled) from
	what should be a single sequencing run. The same organism or strain may have several genome
	entries for it. Typically this will correspond directly to a record in Genbank (assembly
	database).

	The data on this model should primarily pertain to the sample and sequencing run itself. It
	would be updated if for example a better assembly was produced from the original raw data,
	however more advanced interpretation such as taxonomy assignments belong on an attached
	:class:`AnnotatedGenome` object.

	Attributes
	----------
	id : int
		Integer column (primary key).
	key : str
		String column (unique). Unique "external id" used to reference the genome from outside the
		SQL database, e.g. from a file containing K-mer signatures.
	description : Optional[str]
		String column (optional). Short one-line description. Recommended to be unique but this is
		not enforced.
	ncbi_db : Optional[str]
		String column (optional). If the genome corresponds to a record downloaded from an NCBI
		database this column should be the database name (e.g. ``'assembly'``) and ``ncbi_id``
		should be the entry's UID. Unique along with ``ncbi_id``.
	ncbi_id : Optional[int]
		Integer column (optional). See previous.
	genbank_acc : Optional[str]
		String column (optional, unique). GenBank accession number for this genome, if any.
	refseq_acc : Optional[str]
		String column (optional, unique). RefSeq accession number for this genome, if any.
	extra : Optional[dict]
		JSON column (optional). Additional arbitrary metadata.
	annotations : Collection[.AnnotatedGenome]
		One-to-many relationship to :class:`.AnnotatedGenome`.
	"""

	__tablename__ = 'genomes'

	@declared_attr
	def __table_args__(cls):
		return (
			UniqueConstraint('ncbi_db', 'ncbi_id'),
		)

	id = Column(Integer(), primary_key=True)
	key = Column(String(), unique=True, nullable=False)
	description = Column(String(), nullable=False)
	ncbi_db = Column(String())
	ncbi_id = Column(Integer())
	genbank_acc = Column(String(), unique=True)
	refseq_acc = Column(String(), unique=True)
	extra = deferred(Column(JsonString()))

	annotations = relationship('AnnotatedGenome', lazy=True, cascade='all, delete-orphan')

	def __repr__(self):
		return f'<{type(self).__name__}:{self.id} {self.key!r}>'


class ReferenceGenomeSet(Base):
	"""
	A collection of reference genomes along with additional annotations and data. A full GAMBIT
	database which can be used for queries consists of a genome set plus a set of k-mer signatures
	for those genomes (stored separately).

	Membership of :class:`.Genome`s in the set is determined by the presence of an associated
	:class:`.AnnotatedGenomes` object, which also holds additional annotation data for the genome.
	The genome set also includes a set of associated :class:`.Taxon` entries, which form a taxonomy
	tree under which all its genomes are categorized.

	This schema technically allows for multiple genome sets within the same database (which can
	share :class:`.Genome`\\ s but with different annotations), but the GAMBIT application generally
	expects that genome sets are stored in their own SQLite files.

	Attributes
	----------
	id : int
		Integer primary key.
	key : str
		String column. An "external id"  used to uniquely identify this genome set. Unique along
		with ``version``.
	version : str
		Optional version string, an updated version of a previous genome set should have the same key
		with a later version number.
		Should be in the format defined by `PEP 440 <https://www.python.org/dev/peps/pep-0440/>`_.
	name : str
		String column. Unique name.
	description : Optional[str]
		Text column. Optional description.
	extra : Optional[dict]
		JSON column. Additional arbitrary data.
	genomes : Collection[.AnnotatedGenome]
		Many-to-many relationship with :class:`.AnnotatedGenome`, annotated versions of genomes in
		this set.
	base_genomes : Collection[.Genome]
		Unannotated :class:`Genome`\\ s in this set. Association proxy to the ``genome``
		relationship of members of :attr:`genome`.
	taxa : Collection[.Taxon]
		One-to-many relationship to :class:`.Taxon`. The taxa that form the classification system
		for this genome set.
	"""
	__tablename__ = 'genome_sets'

	@declared_attr
	def __table_args__(cls):
		return (
			UniqueConstraint('key', 'version'),
		)

	id = Column(Integer(), primary_key=True)
	key = Column(String(), index=True, nullable=False)
	version = Column(String())
	name = Column(String(), nullable=False)
	description = Column(String())
	extra = Column(JsonString())

	genomes = relationship('AnnotatedGenome', lazy='dynamic', cascade='all, delete-orphan')
	base_genomes = relationship('Genome', secondary='genome_annotations', lazy='dynamic', viewonly=True)

	def __repr__(self):
		return f'<{type(self).__name__}:{self.id} {self.key!r}:{self.version!r}>'

	def root_taxa(self) -> Collection['Taxon']:
		"""Query for root taxa belonging to the set.

		Returns
		-------
		sqlalchemy.orm.query.Query
		"""
		return self.taxa.filter_by(parent=None)


class AnnotatedGenome(Base):
	"""A genome with additional annotations as part of a genome set.

	This object serves to attach a genome to a :class:`.ReferenceGenomeSet`, and to assign a
	taxonomy classification to that genome. Hybrid attributes mirroring the attributes of the
	connected genome effectively make this behave as an extended ``Genome`` object.

	Attributes
	----------
	genome_id : int
		Integer column, part of composite primary key. ID of :class:`.Genome` the annotations are
		for.
	genome_set_id : int
		Integer column, part of composite primary key. ID of the :class:`.ReferenceGenomeSet` the
		annotations are under.
	organism : str
		String column. Single string describing the organism. May be "Genus species [strain]" but
		could contain more specific information. Intended to be human-readable and shouldn't have
		any semantic meaning for the application (in contrast to the :attr:`taxa` relationship).
	taxon_id : int
		Integer column. ID of the :class:`Taxon` this genome is classified as.
	genome : .Genome
		Many-to-one relationship to :class:`.Genome`.
	genome_set : .ReferenceGenomeSet
		Many-to-one relationship to :class:`.ReferenceGenomeSet`.
	taxon : .Taxon
		Many-to-one relationship to :class:`.Taxon`. The primary taxon this genome is classified as
		under the associated ``ReferenceGenomeSet``. Should be the most specific and "regular"
		(ideally defined on NCBI) taxon this genome belongs to.
	key : str
		Hybrid property connected to attribute on :attr:`genome`.
	description : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	ncbi_db : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	ncbi_id : Optional[int]
		Hybrid property connected to attribute on :attr:`genome`.
	genbank_acc : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	refseq_acc : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	"""
	__tablename__ = 'genome_annotations'

	genome_id = Column(ForeignKey('genomes.id', ondelete='CASCADE'), primary_key=True)
	genome_set_id = Column(ForeignKey('genome_sets.id', ondelete='CASCADE'), primary_key=True)
	taxon_id = Column(ForeignKey('taxa.id', ondelete='SET NULL'), index=True)
	organism = Column(String())

	genome = relationship('Genome', back_populates='annotations')
	genome_set = relationship('ReferenceGenomeSet', back_populates='genomes')
	taxon = relationship('Taxon', backref=backref('genomes', lazy='dynamic'))

	key = hybrid_property(lambda self: self.genome.key)
	description = hybrid_property(lambda self: self.genome.description)
	ncbi_db = hybrid_property(lambda self: self.genome.ncbi_db)
	ncbi_id = hybrid_property(lambda self: self.genome.ncbi_id)
	genbank_acc = hybrid_property(lambda self: self.genome.genbank_acc)
	refseq_acc = hybrid_property(lambda self: self.genome.refseq_acc)

	def __repr__(self):
		return '<{}:{}:{} {!r}/{!r}>'.format(
			type(self).__name__,
			self.genome_set_id,
			self.genome_id,
			self.genome_set.key,
			self.genome.key,
		)


class Taxon(Base):
	"""A taxon used for classifying genomes.

	Taxa are specific to a :class:`.ReferenceGenomeSet` and form a tree/forest structure through the
	:attr:`parent` and :attr`children` relationships.

	Attributes
	----------
	id : int
		Integer column (primary key).
	key : str
		String column (unique). An "external id"  used to uniquely identify this taxon.
	name : str
		String column. Human-readable name for the taxon, typically the standard scientific name.
	rank : Optional[str]
		String column (optional). Taxonomic rank, if any. Species, genus, family, etc.
	description : Optional[str]
		String column (optional). Optional description of taxon.
	distance_threshold : Optional[float]
		Float column (optional). Query genomes within this distance of one of the taxon's reference
		genomes will be classified as that taxon. If NULL the taxon is just used establish the tree
		structure and is not used directly in classification.
	report : Bool
		Boolean column. Whether to report this taxon directly as a match when producing a
		human-readable query result. Some custom taxa might need to be "hidden" from the user,
		in which case the value should be false. The application should then ascend the taxon's
		lineage and choose the first ancestor where this field is true. Defaults to true.
	extra : Optional[dict]
		JSON column (optional). Additional arbitrary data.
	genome_set_id : int
		Integer column. ID of :class:`.ReferenceGenomeSet` the taxon belongs to.
	parent_id : Optional[int]
		Integer column. ID of Taxon that is the direct parent of this one.
	ncbi_id : Optional[int]
		Integer column (optional). ID of the entry in the NCBI taxonomy database this taxon
		corresponds to, if any.
	parent : Optional[.Taxon]
		Many-to-one relationship with :class:`.Taxon`, the parent of this taxon (if any).
	children : Collection[.Taxon]
		One-to-many relationship with :class:`.Taxon`, the children of this taxon.
	genome_set : .ReferenceGenomeSet
		Many-to-one relationship to :class:`.ReferenceGenomeSet`.
	genomes : Collection[.AnnotatedGenome]
		One-to-many relationship with :class:`.AnnotatedGenome`, genomes which are assigned to this
		taxon.
	"""

	__tablename__ = 'taxa'

	id = Column(Integer(), primary_key=True)
	key = Column(String(), unique=True, nullable=False)
	name = Column(String(), index=True, nullable=False)
	rank = Column(String(), index=True)
	description = Column(String())
	distance_threshold = Column(Float())
	report = Column(Boolean(), nullable=False, default=True, server_default=sa.true())

	genome_set_id = Column(ForeignKey('genome_sets.id', ondelete='CASCADE'), nullable=False, index=True)
	parent_id = Column(ForeignKey('taxa.id', ondelete='SET NULL'), index=True)
	ncbi_id = Column(Integer(), index=True)
	extra = deferred(Column(JsonString()))

	genome_set = relationship(
		'ReferenceGenomeSet',
		backref=backref('taxa', lazy='dynamic', cascade='all, delete-orphan')
	)
	parent = relationship('Taxon', remote_side=[id], backref=backref('children', lazy=True))

	def ancestors(self, incself=False) -> Iterable['Taxon']:
		"""Iterate through the taxon's ancestors from bottom to top.

		Parameters
		----------
		incself : bool
			If True start with self, otherwise start with parent.
		"""
		taxon = self if incself else self.parent
		while taxon is not None:
			yield taxon
			taxon = taxon.parent

	def lineage(self) -> List['Taxon']:
		"""Get a sorted list of the taxon's ancestors from top to bottom, including itself."""
		l = list(self.ancestors(incself=True))
		l.reverse()
		return l

	def root(self) -> 'Taxon':
		"""Get the root taxon of this taxon's tree.

		The set of taxa in a :class:`.ReferenceGenomeSet` will generally form
		a forest instead of a single tree, so there can be multiple root taxa.

		Returns self if the taxon has no parent.
		"""
		if self.parent is None:
			return self
		else:
			return self.parent.root()

	def isleaf(self) -> bool:
		"""Check if the taxon is a leaf (has no children)."""
		return not self.children

	def descendants(self, incself=False) -> Iterable['Taxon']:
		"""Iterate through taxa all of the taxon's descendants (pre-order depth-first).

		Parameters
		----------
		incself : bool
			Yield self first.
		"""
		if incself:
			yield self
		for child in self.children:
			yield from child.descendants(incself=True)

	def leaves(self) -> Iterable['Taxon']:
		"""Iterate through all leaves in the taxon's subtree.

		For leaf taxa this will just yield the taxon itself.
		"""
		if self.isleaf():
			yield self
		else:
			for child in self.children:
				yield from child.leaves()

	def print_tree(self, indent='  ', *, _depth=0):
		"""Print the taxon's subtree for debugging.

		Parameters
		---------
		indent : str
			String used to indent each level of descendants.
		"""
		print(indent * _depth + self.name)
		for child in sorted(self.children, key=lambda c: c.name):
			child.print_tree(indent=indent, _depth=_depth + 1)

	def __repr__(self):
		return f'<{type(self).__name__}:{self.id} {self.name!r}>'


#: Attributes of :class:`gambit.db.models.Genome` which serve as unique IDs.
GENOME_ID_ATTRS = {
	'key': Genome.key,
	'genbank_acc': Genome.genbank_acc,
	'refseq_acc': Genome.refseq_acc,
	'ncbi_id': Genome.ncbi_id,
}


# Type alias for argument specifying genome id attribute
GenomeAttr = Union[str, InstrumentedAttribute]


def _check_genome_id_attr(attr: GenomeAttr) -> InstrumentedAttribute:
	"""Check that Genome ID attribute is valid, and convert from string argument.
	"""
	if isinstance(attr, str):
		try:
			return GENOME_ID_ATTRS[attr]
		except KeyError:
			pass

	elif isinstance(attr, InstrumentedAttribute):
		for allowed in GENOME_ID_ATTRS.values():
			if attr is allowed:
				return attr

	raise ValueError('Genome ID attribute must be one of the following: ' + ', '.join(GENOME_ID_ATTRS))


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
