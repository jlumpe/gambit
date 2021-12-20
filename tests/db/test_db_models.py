"""Test gambit.db.models.

Uses the included testdb_210818 database.
"""

import random

import pytest
from sqlalchemy.orm import sessionmaker

from gambit.db import models
from gambit.db import Genome, ReferenceGenomeSet, AnnotatedGenome, Taxon


# Some arbitrary JSON data
JSON_DATA = {
	'int': 1,
	'float': 3.14,
	'string': 'foo',
	'bool': True,
	'null': None,
	'array': [1, 3.14, 'foo', True, None, [1, 2, 3, 4], {'foo': 'bar', 'baz': 1}],
	'object': {
		'foo': 'bar',
		'baz': 1,
		'sub_list': [1, 2, 3],
		'sub_object': {'blah': 123}
	},
}

JSON_DATA2 = JSON_DATA['object']


@pytest.fixture()
def empty_db_session(make_empty_db):
	"""Session factory for empty in-memory database."""
	engine = make_empty_db()
	return sessionmaker(engine)



class TestGenome:
	"""Test Genome model."""

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		session = empty_db_session()

		# Save genome with JSON data
		genome = Genome(
			key='foo',
			description='test genome',
			extra=JSON_DATA,
		)
		session.add(genome)
		session.commit()

		# Reload in fresh session and check value
		session = empty_db_session()
		genome = session.query(Genome).one()
		assert genome.extra == JSON_DATA

		# Assign different data, save
		genome.extra = JSON_DATA2
		session.commit()

		# Check new value
		session = empty_db_session()
		genome = session.query(Genome).one()
		assert genome.extra == JSON_DATA2

		# Assign NULL, save
		genome.extra = None
		session.commit()

		# Check new value
		session = empty_db_session()
		genome = session.query(Genome).one()
		assert genome.extra is None


class TestReferenceGenomeSet:
	"""Test ReferenceGenomeSet model."""

	def test_root_taxa(self, testdb_session):
		session = testdb_session()
		gset = session.query(ReferenceGenomeSet).one()
		assert {taxon.name for taxon in gset.root_taxa()} == {'A1', 'A2', 'A3'}

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		session = empty_db_session()

		# Save genome set with JSON data
		gset = ReferenceGenomeSet(
			key='foo',
			version='1.0',
			name='test',
			extra=JSON_DATA,
		)
		session.add(gset)
		session.commit()

		# Reload in fresh session and check value
		session = empty_db_session()
		gset = session.query(ReferenceGenomeSet).one()
		assert gset.extra == JSON_DATA

		# Assign different data, save
		gset.extra = JSON_DATA2
		session.commit()

		# Check new value
		session = empty_db_session()
		gset = session.query(ReferenceGenomeSet).one()
		assert gset.extra == JSON_DATA2

		# Assign NULL, save
		gset.extra = None
		session.commit()

		# Check new value
		session = empty_db_session()
		gset = session.query(ReferenceGenomeSet).one()
		assert gset.extra is None


class TestAnnotatedGenome:
	"""Test AnnotatedGEnome model."""

	def test_hybrid_props(self, testdb_session):
		session = testdb_session()

		hybrid_attrs = [
			'key',
			'description',
			'ncbi_db',
			'ncbi_id',
			'genbank_acc',
			'refseq_acc',
		]

		for annotated in session.query(AnnotatedGenome):
			for attr in hybrid_attrs:
				assert getattr(annotated, attr) == getattr(annotated.genome, attr)


class TestTaxon:
	"""Test Taxon model."""

	def test_tree(self, testdb_session):
		"""Test tree structure."""
		session = testdb_session()
		gset = session.query(ReferenceGenomeSet).one()
		roots = gset.root_taxa()

		for taxon in gset.taxa:
			root = taxon.root()
			assert root in roots
			assert taxon.isleaf() == (len(taxon.children) == 0)

			# Test parent/child relationships match
			for child in taxon.children:
				assert child.parent == taxon

			# Test ancestors() and lineage() methods
			ancestors = list(taxon.ancestors(incself=True))
			assert ancestors[0] is taxon
			assert ancestors[-1] is root
			assert list(taxon.ancestors()) == ancestors[1:]
			assert list(reversed(taxon.lineage())) == ancestors
			assert taxon.depth() == len(ancestors) - 1

			for i in range(len(ancestors) - 1):
				assert ancestors[i].parent is ancestors[i + 1]

			# Check traversal methods
			descendants_set = {t for t in gset.taxa if taxon in t.ancestors()}
			self.check_traversal(taxon.descendants(False), False, descendants_set)
			self.check_traversal(taxon.descendants(True), True, descendants_set)

			subtree_set = descendants_set | {taxon}
			self.check_traversal(taxon.traverse(False), False, subtree_set)
			self.check_traversal(taxon.traverse(True), True, subtree_set)

			# Check leaves
			assert set(taxon.leaves()) == {d for d in subtree_set if d.isleaf()}

	def check_traversal(self, iterator, postorder, expected):
		seen = set()

		for taxon in iterator:
			if postorder:
				assert all(child in seen for child in taxon.children)
			elif taxon.parent is not None and taxon.parent in expected:
				assert taxon.parent in seen

			seen.add(taxon)

		assert seen == expected

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		session = empty_db_session()

		# Save taxon with JSON data
		gset = ReferenceGenomeSet(
			key='test',
			version='1.0',
			name='test genome set',
		)
		taxon = Taxon(
			genome_set=gset,
			key='test',
			name='test taxon',
			extra=JSON_DATA,
		)
		session.add(gset)
		session.add(taxon)
		session.commit()

		# Reload in fresh session and check value
		session = empty_db_session()
		taxon = session.query(Taxon).one()
		assert taxon.extra == JSON_DATA

		# Assign different data, save
		taxon.extra = JSON_DATA2
		session.commit()

		# Check new value
		session = empty_db_session()
		taxon = session.query(Taxon).one()
		assert taxon.extra == JSON_DATA2

		# Assign NULL, save
		taxon.extra = None
		session.commit()

		# Check new value
		session = empty_db_session()
		taxon = session.query(Taxon).one()
		assert taxon.extra is None


class TestGenomeIDMapping:
	"""Test mapping genomes to ID values."""

	@pytest.fixture()
	def session(self, make_empty_db):
		"""In-memory database containing genomes which have values for all ID attributes."""
		engine = make_empty_db()
		Session = sessionmaker(engine)
		session = Session()

		gset = ReferenceGenomeSet(
			key='test_gset',
			version='1.0',
			name='Test genome set',
		)
		session.add(gset)

		roottaxon = Taxon(
			key='root',
			name='root',
			genome_set=gset,
		)
		session.add(roottaxon)

		for i in range(20):
			g = Genome(
				key=f'test/genome_{i}',
				description=f'Test genome {i}',
				ncbi_db='assembly',
				ncbi_id=i,
				genbank_acc=f'GCA_{i:09d}.1',
				refseq_acc=f'GCF_{i:09d}.1',
			)
			ag = AnnotatedGenome(
				genome_set=gset,
				genome=g,
				taxon=roottaxon,
			)
			session.add(ag)

		session.commit()
		return session

	def test__genome_id_attr(self):
		"""Test _check_genome_id_attr() function."""

		for key, attr in models.GENOME_ID_ATTRS.items():
			assert models._check_genome_id_attr(key) is attr
			assert models._check_genome_id_attr(attr) is attr

		for arg in ['description', Genome.description, AnnotatedGenome.key]:
			with pytest.raises(ValueError):
				models._check_genome_id_attr(arg)

	def test__get_genome_id(self, session):
		"""Test _get_genome_id() function."""

		for genome in session.query(Genome):
			for key, attr in models.GENOME_ID_ATTRS.items():
				assert models._get_genome_id(genome, attr) == getattr(genome, key)

	def test_genomes_by_id(self, session):
		"""Test genomes_by_id() and genomes_by_id_subset() functions."""
		random.seed(0)

		gset = session.query(ReferenceGenomeSet).one()

		# Desired result with only valid IDs
		genomes = list(gset.genomes)
		random.shuffle(genomes)

		# Desired result with some invalid IDs mixed in
		genomes_missing = list(genomes)
		genomes_missing.insert(3, None)
		genomes_missing.insert(8, None)

		for key, attr in models.GENOME_ID_ATTRS.items():
			# Check with all valid IDs
			ids = [models._get_genome_id(g, attr) for g in genomes]
			assert models.genomes_by_id(gset, attr, ids) == genomes
			assert models.genomes_by_id(gset, key, ids) == genomes

			# Add IDs which cannot be mapped to genomes
			ids_missing = [
				models._get_genome_id(g, attr) if g is not None else '!INVALID!'
				for g in genomes_missing
			]

			# Raises KeyError with strict=True (default)
			with pytest.raises(KeyError):
				models.genomes_by_id(gset, attr, ids_missing)

			# strict=False returns list with None values
			assert models.genomes_by_id(gset, attr, ids_missing, strict=False) == genomes_missing

			# Test genomes_by_id_subset
			genomes_sub, idxs_sub = models.genomes_by_id_subset(gset, attr, ids_missing)
			assert genomes_sub == genomes
			assert [genomes_missing[i] for i in idxs_sub] == genomes

			# Incomplete set of IDs which does not encompass all genomes
			ids_incomplete = ids[:-1]

			# TODO


def test_reportable_taxon():
	"Test reportable_taxon() function."
	assert models.reportable_taxon(None) is None

	t = Taxon(report=True)
	assert models.reportable_taxon(t) is t

	t = Taxon(report=False)
	assert models.reportable_taxon(t) is None

	t1 = Taxon(report=False)
	t2 = Taxon(report=True)
	t1.parent = t2
	assert models.reportable_taxon(t1) is t2
