from typing import Sequence

import click


def print_table(rows: Sequence[Sequence], colsep: str = ' ', left: str = '', right: str = ''):
	"""Print a basic table."""

	echo = lambda s: click.echo(s, nl=False)

	rows = [list(map(str, row)) for row in rows]
	ncol = max(map(len, rows))

	widths = [0] * ncol
	for row in rows:
		for i, val in enumerate(row):
			widths[i] = max(widths[i], len(val))

	for row in rows:
		echo(left)

		for i, val in enumerate(row):
			echo(val.ljust(widths[i]))

			if i < ncol - 1:
				echo(colsep)

		echo(right)
		echo('\n')
