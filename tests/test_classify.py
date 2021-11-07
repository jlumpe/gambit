"""Test the gambit.query.classify module."""

import pytest

from gambit.classify import matching_taxon, find_matches, consensus_taxon
from gambit.db import Taxon


def test_matching_taxon():
	taxa = []
	for d in [.3, None, .2, .1]:
		taxa.insert(0, Taxon(
			distance_threshold=d,
			parent=taxa[0] if taxa else None,
		))

	assert matching_taxon(taxa[0], .05) == taxa[0]
	assert matching_taxon(taxa[0], .1) == taxa[0]
	assert matching_taxon(taxa[0], .15) == taxa[1]
	assert matching_taxon(taxa[0], .25) == taxa[3]
	assert matching_taxon(taxa[0], .35) is None


def test_find_matches():
	pass  # TODO


def test_consensus_taxon(testdb_session):
	session = testdb_session()
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
