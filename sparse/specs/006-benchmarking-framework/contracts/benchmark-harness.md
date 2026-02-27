# Contract: Benchmark Harness Functions

**Purpose**: Defines the public API for running benchmarks, managing baselines, and generating reports.

## Function: `run_component_benchmarks`

- **Input**: `Criterion` instance, `BenchmarkConfig`, reference to a `Benchmarkable` solver
- **Output**: Criterion benchmark groups registered and executed
- **Semantics**: For each phase in config, creates a `BenchmarkGroup` named `ssids/{phase}`. For each matrix matching the filter, registers a benchmark parameterized by matrix name. Skips missing matrices and unimplemented phases with `eprintln!` warnings.

## Function: `run_e2e_benchmarks`

- **Input**: `Criterion` instance, `BenchmarkConfig`, reference to a `Benchmarkable` solver
- **Output**: Criterion benchmark group for roundtrip benchmarks
- **Semantics**: Creates a single `BenchmarkGroup` named `ssids/roundtrip`. For each matrix matching the filter, benchmarks the full analyze → factor → solve pipeline.

## Function: `read_peak_rss_kb`

- **Input**: None
- **Output**: `Option<u64>` — peak RSS in kilobytes
- **Semantics**: Reads `/proc/self/status`, parses `VmHWM` line. Returns `None` on non-Linux platforms or if the file cannot be read.

## Function: `collect_results`

- **Input**: Criterion output directory path, `BenchmarkConfig`
- **Output**: `BenchmarkSuiteResult`
- **Semantics**: Parses `raw.csv` and `estimates.json` files from `target/criterion/` for benchmarks matching the config. Assembles into structured results with matrix metadata.

## Function: `save_baseline`

- **Input**: `BenchmarkSuiteResult`, baseline name, output directory
- **Output**: JSON file written to `{output_dir}/{name}.json`
- **Semantics**: Serializes suite results with timestamp and name. Overwrites existing baseline with same name.

## Function: `load_baseline`

- **Input**: Baseline name, directory
- **Output**: `Result<Baseline>` — deserialized baseline or error if not found
- **Semantics**: Reads and deserializes JSON baseline file.

## Function: `detect_regressions`

- **Input**: Current `BenchmarkSuiteResult`, `Baseline`, threshold percentage (default: 5.0)
- **Output**: `RegressionReport`
- **Semantics**: For each matching (matrix, phase) pair, compares mean times. Flags regressions (current > baseline by more than threshold%), improvements (current < baseline by more than threshold%), and unchanged entries.

## Function: `export_csv`

- **Input**: `BenchmarkSuiteResult`, output path
- **Output**: CSV file with one row per (matrix, phase) result
- **Semantics**: Columns: matrix_name, phase, mean_ns, std_dev_ns, median_ns, iterations, matrix_size, matrix_nnz, throughput_nnz_per_sec.

## Function: `export_json`

- **Input**: `BenchmarkSuiteResult`, output path
- **Output**: JSON file with full suite result
- **Semantics**: Serializes the complete `BenchmarkSuiteResult` via serde.

## Function: `generate_markdown_table`

- **Input**: `BenchmarkSuiteResult` or `RegressionReport`
- **Output**: `String` containing a Markdown table
- **Semantics**: Formats results into aligned Markdown table with columns for matrix, phase, mean time, std dev, and optionally throughput or regression percentage.
