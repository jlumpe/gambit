"""Benchmarks for gambit.search module."""

import pytest
import numpy as np

from gambit.kmers import KmerSpec
from gambit.search import calc_signature
from gambit.test import random_seq


@pytest.fixture(scope='module', params=[10**4, 10**6])
def seq(request):
	np.random.seed(0)
	return random_seq(request.param)


@pytest.fixture(scope='module', params=[8, 11, 14])
def kspec_k(request):
	return request.param


@pytest.fixture(scope='module', params=[3, 5, 7])
def kspec_prefix(request):
	return 'ATGACCT'[:request.param]


@pytest.fixture()
def kspec(kspec_k, kspec_prefix):
	return KmerSpec(kspec_k, kspec_prefix)


def benchmark_calc_signature(seq, kspec, benchmark):
	benchmark(calc_signature, kspec, seq)
