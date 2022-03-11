"""Test the gambit.query.classify module."""

import pytest

from gambit.classify import matching_taxon, find_matches, consensus_taxon, GenomeMatch
from gambit.db import Taxon, AnnotatedGenome
from gambit.test import make_lineage


def test_matching_taxon():
	taxa = make_lineage([.1, .2, None, .3])

	assert matching_taxon(taxa[0], .05) == taxa[0]
	assert matching_taxon(taxa[0], .1) == taxa[0]
	assert matching_taxon(taxa[0], .15) == taxa[1]
	assert matching_taxon(taxa[0], .25) == taxa[3]
	assert matching_taxon(taxa[0], .35) is None


def test_find_matches():
	pass  # TODO


def test_consensus_taxon(testdb):
	session = testdb.Session()
	get_taxon = lambda name: session.query(Taxon).filter_by(name=name).one()

	A1 = get_taxon('A1')
	A1_B1 = get_taxon('A1_B1')
	A1_B1_C1 = get_taxon('A1_B1_C1')
	A1_B1_C2 = get_taxon('A1_B1_C2')
	A1_B2 = get_taxon('A1_B2')
	A2 = get_taxon('A2')

	# Empty
	assert consensus_taxon([]) == (None, set())

	# Single
	assert consensus_taxon([A1]) == (A1, set())

	# In single lineage
	assert consensus_taxon([A1, A1_B1]) == (A1_B1, set())
	assert consensus_taxon([A1, A1_B1_C1]) == (A1_B1_C1, set())
	assert consensus_taxon([A1, A1_B1, A1_B1_C1]) == (A1_B1_C1, set())

	# Split, with common ancestor
	assert consensus_taxon([A1_B1, A1_B2]) == (A1, {A1_B1, A1_B2})
	assert consensus_taxon([A1_B1_C1, A1_B2]) == (A1, {A1_B1_C1, A1_B2})
	assert consensus_taxon([A1_B1, A1_B1_C1, A1_B2]) == (A1, {A1_B1, A1_B1_C1, A1_B2})
	assert consensus_taxon([A1, A1_B1_C1, A1_B1_C2]) == (A1_B1, {A1_B1_C1, A1_B1_C2})

	# Split, no common ancestor
	assert consensus_taxon([A1, A2]) == (None, {A1, A2})
	assert consensus_taxon([A1_B1, A2]) == (None, {A1_B1, A2})
	assert consensus_taxon([A1, A1_B1, A2]) == (None, {A1, A1_B1, A2})


class TestGenomeMatch:
	"""Test the GenomeMatch class."""

	def test_matched_taxon_default(self):
		"""Test default for the matched_taxon attribute."""

		t1, t2 = make_lineage([.2, .5])
		g = AnnotatedGenome(taxon=t1)

		assert GenomeMatch(g, .1) == GenomeMatch(g, .1, t1)
		assert GenomeMatch(g, .3) == GenomeMatch(g, .3, t2)
		assert GenomeMatch(g, .6) == GenomeMatch(g, .6, None)

	def test_next_taxon(self):
		"""Test the next_taxon() method."""

		taxa = make_lineage([.2, .4, None, .6, None])
		g = AnnotatedGenome(taxon=taxa[0])
		matches = [GenomeMatch(g, d) for d in [.1, .4, .5, .7]]

		assert matches[0].matched_taxon == taxa[0]
		assert matches[0].next_taxon() is None
		assert matches[1].matched_taxon == taxa[1]
		assert matches[1].next_taxon() == taxa[0]
		assert matches[2].matched_taxon == taxa[3]
		assert matches[2].next_taxon() == taxa[1]
		assert matches[3].matched_taxon is None
		assert matches[3].next_taxon() == taxa[3]
