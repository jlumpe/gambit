"""Tests for the "gambit tree" command."""

from io import StringIO

import pytest
from Bio import Phylo

from gambit.metric import jaccarddist_pairwise
from gambit.cli.test import invoke_cli
from gambit.cluster import hclust, check_tree_matches_linkage
from gambit.cli import common

@pytest.fixture()
def expected_dmat(testdb):
	sigs = testdb.query_signatures
	return jaccarddist_pairwise(sigs)

@pytest.fixture()
def expected_linkage(expected_dmat):
	return hclust(expected_dmat)

@pytest.mark.parametrize('from_sigs', [False, True])
def test_tree_command(from_sigs, expected_linkage, testdb):
	seqfiles = [str(f.path) for f in testdb.get_query_files()]

	args = ['tree']
	if from_sigs:
		args += ['-s', testdb.paths.query_signatures]
	else:
		kspec = testdb.kmerspec
		args += ['-k', kspec.k, '--prefix', kspec.prefix_str]
		args += seqfiles

	result = invoke_cli(args)
	result_buf = StringIO(result.stdout)
	result_tree = Phylo.read(result_buf, 'newick')

	expected_labels = list(map(common.get_file_id, seqfiles))

	check_tree_matches_linkage(result_tree, expected_linkage, expected_labels)

