"""Distance matrices and basic clustering/trees."""

from typing import Union, Optional, Sequence, TextIO, Tuple, List
import csv

import numpy as np

from gambit.util.io import FilePath, maybe_open
from gambit.util.misc import zip_strict


def dump_dmat_csv(file: Union[FilePath, TextIO],
                  dmat: np.ndarray,
                  row_ids: Sequence,
                  col_ids: Sequence,
                  corner: Optional[str] = None,
                  fmt: str = '0.4f',
                  ):
	"""Write distance matrix to file in CSV format."""

	with maybe_open(file, 'w', newline='') as fobj:
		writer = csv.writer(fobj)
		writer.writerow([corner or '', *map(str, col_ids)])
		for row_id, values in zip_strict(row_ids, dmat):
			values_str = (format(d, fmt) for d in values)
			writer.writerow([str(row_id), *values_str])


def load_dmat_csv(file: Union[FilePath, TextIO]) -> Tuple[np.ndarray, List[str], List[str]]:
	"""Load distance matrix from CSV file.

	Returns
	-------
	Tuple
		``(matrix, row_ids, col_ids)`` tuple.
	"""

	with maybe_open(file, newline='') as fobj:
		reader = csv.reader(fobj)
		col_ids = next(reader)[1:]
		nc = len(col_ids)

		rows = []
		row_ids = []

		for rid, *row_str in reader:
			row_ids.append(rid)
			rows.append(np.fromiter(map(float, row_str), np.float32, count=nc))

		values = np.stack(rows)
		return values, row_ids, col_ids
