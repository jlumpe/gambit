"""Distance matrices and basic clustering/trees."""

from typing import Union, Optional, Sequence, TextIO
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
