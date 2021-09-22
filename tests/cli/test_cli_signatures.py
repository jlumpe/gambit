"""Tests for the "signatures" command group."""

import json

import numpy as np
from click.testing import CliRunner

from gambit.cli import cli
import gambit.io.json as gjson


class TestInfoCommand:
	def test_standard(self):
		pass  # TODO

	def test_json(self, testdb_signatures, testdb_files):
		args = ['signatures', 'info', str(testdb_files['ref_signatures']), '-j']

		runner = CliRunner()
		result = runner.invoke(cli, args, catch_exceptions=False)
		assert result.exit_code == 0

		data = json.loads(result.stdout)
		assert data['count'] == len(testdb_signatures)
		assert data['kmerspec'] == gjson.to_json(testdb_signatures.kmerspec)
		assert data['metadata'] == gjson.to_json(testdb_signatures.meta)

	def test_ids(self, testdb_signatures, testdb_files):
		args = ['signatures', 'info', str(testdb_files['ref_signatures']), '-i']

		runner = CliRunner()
		result = runner.invoke(cli, args, catch_exceptions=False)
		assert result.exit_code == 0

		assert np.array_equal(result.stdout.splitlines(), testdb_signatures.ids)
