from typing import NewType

import numpy as np

#: Type for k-mer signatures (k-mer sets in sparse coordinate format)
KmerSignature = NewType('KmerSignature', np.ndarray)
# TODO - use nptyping package to specify dimensions and data type?
