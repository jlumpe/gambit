"""Benchmarks for gambit.search module."""

import pytest
import numpy as np

from gambit.kmers import KmerSpec
from gambit.sigs.calc import calc_signature, ArrayAccumulator, SetAccumulator
from gambit.test import random_seq


@pytest.fixture(scope='module', params=[10**4, 10**6])
def seq(request):
	np.random.seed(0)
	return random_seq(request.param)


@pytest.fixture(scope='module', params=[8, 11, 14])
def k(request):
	return request.param


@pytest.fixture(scope='module', params=[3, 5, 7])
def prefix_len(request):
	return request.param


@pytest.fixture()
def kspec(k, prefix_len):
	prefix ='ATGACCT'[:prefix_len]
	return KmerSpec(k, prefix)


@pytest.fixture(
	scope='module',
	params=[pytest.param(ArrayAccumulator, id='array'), pytest.param(SetAccumulator, id='set')],
)
def accumulator(request):
	return request.param


def benchmark_calc_signature(seq, kspec, benchmark, accumulator):
	acc = accumulator(kspec.k)
	benchmark(calc_signature, kspec, seq, accumulator=acc)
