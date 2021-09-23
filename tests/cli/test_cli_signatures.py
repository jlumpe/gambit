"""Tests for the "signatures" command group."""

import json

import pytest
import numpy as np

from gambit.cli.test import invoke_cli
import gambit.io.json as gjson


class TestInfoCommand:

	@pytest.fixture()
	def base_args(self, testdb_files):
		return ['signatures', 'info', str(testdb_files['ref_signatures'])]

	def test_standard(self, base_args, testdb_signatures):
		result = invoke_cli(base_args)
		assert result.exit_code == 0

		# TODO: check

	def test_json(self, base_args, testdb_signatures):
		result = invoke_cli([*base_args, '-j'])
		assert result.exit_code == 0

		data = json.loads(result.stdout)
		assert data['count'] == len(testdb_signatures)
		assert data['kmerspec'] == gjson.to_json(testdb_signatures.kmerspec)
		assert data['metadata'] == gjson.to_json(testdb_signatures.meta)

	def test_ids(self, base_args, testdb_signatures):
		result = invoke_cli([*base_args, '-i'])
		assert result.exit_code == 0

		assert np.array_equal(result.stdout.splitlines(), testdb_signatures.ids)
