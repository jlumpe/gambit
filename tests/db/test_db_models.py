"""Test gambit.db.models.

Uses the included testdb_210818 database.
"""

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


def check_json_col(empty_db_session, instance, col: str):
	"""Check storing/loading data in JSON column."""
	model = type(instance)

	# Save instance with JSON data
	session = empty_db_session()
	setattr(instance, col, JSON_DATA)
	session.add(instance)
	session.commit()

	# Reload in fresh session and check value
	session = empty_db_session()
	instance = session.query(model).one()
	assert getattr(instance, col) == JSON_DATA

	# Assign different data, save
	setattr(instance, col, JSON_DATA2)
	session.commit()

	# Check new value
	session = empty_db_session()
	instance = session.query(model).one()
	assert getattr(instance, col) == JSON_DATA2

	# Assign NULL, save
	setattr(instance, col, None)
	session.commit()

	# Check new value
	session = empty_db_session()
	instance = session.query(model).one()
	assert getattr(instance, col) is None


class TestGenome:
	"""Test Genome model."""

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		genome = Genome(
			key='foo',
			description='test genome',
		)
		check_json_col(empty_db_session, genome, 'extra')


class TestReferenceGenomeSet:
	"""Test ReferenceGenomeSet model."""

	def test_root_taxa(self, testdb):
		session = testdb.Session()
		gset = session.query(ReferenceGenomeSet).one()
		assert {taxon.name for taxon in gset.root_taxa()} == {'A1', 'A2', 'A3'}

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		gset = ReferenceGenomeSet(
			key='foo',
			version='1.0',
			name='test',
			extra=JSON_DATA,
		)
		check_json_col(empty_db_session, gset, 'extra')


class TestAnnotatedGenome:
	"""Test AnnotatedGEnome model."""

	def test_hybrid_props(self, testdb):
		session = testdb.Session()

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

	def test_tree(self, testdb):
		"""Test tree structure."""
		session = testdb.Session()
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

	def test_genome_membership(self, testdb):
		"""Test the subtree_genomes() and has_genome() methods."""
		session = testdb.Session()

		for taxon in session.query(Taxon):
			sg = set()

			for genome in session.query(AnnotatedGenome):
				if genome.key.split('/')[-1].startswith(taxon.name):
					sg.add(genome)
					assert taxon.has_genome(genome)
				else:
					assert not taxon.has_genome(genome)

			assert set(taxon.subtree_genomes()) == sg

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		session = empty_db_session()

		# Save taxon with JSON data
		gset = ReferenceGenomeSet(
			key='test',
			version='1.0',
			name='test genome set',
		)
		session.add(gset)
		session.commit()

		taxon = Taxon(
			genome_set_id=gset.id,
			key='test',
			name='test taxon',
		)
		check_json_col(empty_db_session, taxon, 'extra')

	def taxon_by_name(self, session, name):
		return session.query(Taxon).filter_by(name=name).one()

	def check_common_ancestry(self, session, names, expected_names):
		taxa = [self.taxon_by_name(session, name) for name in names]
		ca = Taxon.common_ancestors(taxa)
		lca = Taxon.lca(taxa)

		assert [taxon.name for taxon in ca] == expected_names
		if ca:
			assert lca is ca[-1]
		else:
			assert lca is None

	def test_common_ancestry(self, testdb):
		"""Test the common_ancestors() and lca() methods."""

		session = testdb.Session()

		self.check_common_ancestry(session, [], [])

		self.check_common_ancestry(session, ['A1'], ['A1'])
		self.check_common_ancestry(session, ['A1_B1'], ['A1', 'A1_B1'])
		self.check_common_ancestry(session, ['A1_B1_C1'], ['A1', 'A1_B1', 'A1_B1_C1'])

		self.check_common_ancestry(session, ['A1_B1', 'A1_B2'], ['A1'])
		self.check_common_ancestry(session, ['A1_B1_C1', 'A1_B1_C2'], ['A1', 'A1_B1'])

		self.check_common_ancestry(session, ['A1', 'A1_B1'], ['A1'])
		self.check_common_ancestry(session, ['A1_B1', 'A1_B1_C1'], ['A1', 'A1_B1'])

		self.check_common_ancestry(session, ['A1', 'A2'], [])
		self.check_common_ancestry(session, ['A1_B1', 'A1_B2', 'A2_B1'], [])


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
