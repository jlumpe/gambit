"""Efficient formats for storing k-mer signatures in memory and in the file system."""

from .base import SignaturesMeta
from .array import SignatureArray, sigarray_eq
