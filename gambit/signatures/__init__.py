"""Efficient formats for storing k-mer signatures in memory and in the file system."""

from .base import SignaturesMeta, sigarray_eq
from .array import SignatureArray
