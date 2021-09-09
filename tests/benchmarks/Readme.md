# Benchmarks

This directory contains pytest files which run benchmarks using the
`pytest-benchmark` library. The file names deliberately omit the `test_`
prefix so that they are not discovered as part of a normal run. To run the
benchmarks use `pytest tests/benchmarks/*`.

Note that discovery has been altered to allow functions/classes prefixed with
`benchmark_`/`Benchmark`to be recognized as tests.
