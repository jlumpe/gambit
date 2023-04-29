#!/usr/bin/env python3

# Transpose output table for better display in docs

import sys
import pandas as pd

INT_DT = pd.Int64Dtype()


def format_float(x):
    return '' if pd.isnull(x) else format(x, '.4f')


# Parse args
if len(sys.argv) != 3:
    print('Usage:', sys.argv[0], 'IN_FILE', 'OUT_FILE')
    sys.exit(1)

infile, outfile = sys.argv[1:]


table = pd.read_csv(infile)

# Format columns
for colname, values in table.items():
    # Convert integer columns to proper dtype
    if colname.endswith('.ncbi_id'):
        table[colname] = values.astype(INT_DT)
        continue

    # Format floats as strings with proper precision
    if values.dtype.kind == 'f':
        table[colname] = values.map(format_float)

# Write transposed version
table_t = table.transpose()
table_t.to_csv(outfile, header=False, na_rep='')
