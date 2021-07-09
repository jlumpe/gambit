"""Test the gambit.query.classify module."""

import pytest

from gambit.query.classify import matching_taxon, find_matches, consensus_taxon, reportable_taxon
from gambit.db.models import Taxon


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


def test_consensus_taxon():
	pass  # TODO


def test_reportable_taxon():
	pass  # TODO
