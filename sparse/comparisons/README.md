# External Solver Comparisons

This directory contains tools for benchmarking rivrs-sparse against external
reference solvers. Each solver has a standalone Fortran or C driver that runs
as a subprocess — no FFI linking is required. This keeps the rivrs-sparse crate
free of build-time dependencies on external solver libraries.

## Supported Solvers

| Solver | License | Driver | Status |
|--------|---------|--------|--------|
| SPRAL SSIDS | BSD-3 | `drivers/spral_benchmark.f90` | Working |
| MUMPS | Public domain | `drivers/mumps_benchmark.f90` | Working |

## Directory Structure

```
comparisons/
├── README.md                          # This file
├── drivers/                           # Fortran/C driver source and build scripts
│   ├── build_spral.sh                 # Builds libspral.a from SPRAL source
│   ├── build_mumps.sh                 # Builds MUMPS driver (requires libmumps-seq-dev)
│   ├── build_spral_dense_factor.sh    # Builds dense factor comparison driver
│   ├── spral_benchmark.f90            # SPRAL SSIDS benchmark driver (subprocess)
│   ├── mumps_benchmark.f90            # MUMPS benchmark driver (subprocess)
│   ├── spral_full_solve.f90           # SPRAL full solve driver
│   ├── spral_match_order.f90          # SPRAL ordering comparison driver
│   └── spral_dense_factor.cpp         # SPRAL dense APTP comparison
└── src/
    ├── common.rs                      # Shared types, subprocess protocol, formatters
    ├── spral_benchmark.rs             # SPRAL orchestration binary
    └── mumps_benchmark.rs             # MUMPS orchestration binary
```

## Prerequisites

### SPRAL

1. Ensure SPRAL source is available at `/opt/references/spral/` (or adjust paths
   in the build script).

2. Build the SPRAL library and driver binaries:

```sh
comparisons/drivers/build_spral.sh
```

This produces `/tmp/spral_benchmark` (and other driver binaries).

### MUMPS

MUMPS is available as a system package (sequential version with dummy MPI).

1. Install the library:

```sh
# Debian/Ubuntu
apt-get install libmumps-seq-dev
```

2. Build the MUMPS driver:

```sh
comparisons/drivers/build_mumps.sh
```

This produces `/tmp/mumps_benchmark`.

**Ordering**: MUMPS auto-selects its ordering by default (ICNTL(7)=7). Override
via the `MUMPS_ORDERING` environment variable: `auto`, `metis`, `amd`, `scotch`, `pord`.

### METIS

The SPRAL drivers require METIS. The `metis-sys` crate vendored by rivrs-sparse
provides a static `libmetis.a`:

```sh
METIS_LIB=$(find target -name "libmetis.a" | head -1)
```

## Running Comparisons

### SPRAL Benchmark Suite

```sh
# SPRAL-only on CI subset
cargo run --bin spral-comparison --release -- --ci-only

# Side-by-side with rivrs
cargo run --bin spral-comparison --release -- --ci-only --rivrs

# Control thread count (sets OMP_NUM_THREADS for SPRAL, Par::rayon for rivrs)
cargo run --bin spral-comparison --release -- --ci-only --threads 4

# Compare against a previously collected rivrs baseline
cargo run --bin spral-comparison --release -- --ci-only \
  --compare target/benchmarks/baselines/baseline-latest.json

# Filter by category
cargo run --bin spral-comparison --release -- --category hard-indefinite
```

JSON output: `target/benchmarks/spral/`

### MUMPS Benchmark Suite

```sh
# MUMPS-only on CI subset
cargo run --bin mumps-comparison --release -- --ci-only

# Side-by-side with rivrs
cargo run --bin mumps-comparison --release -- --ci-only --rivrs

# Control rivrs thread count (MUMPS sequential is single-threaded)
cargo run --bin mumps-comparison --release -- --ci-only --rivrs --threads 4

# Compare against a previously collected rivrs baseline
cargo run --bin mumps-comparison --release -- --ci-only \
  --compare target/benchmarks/baselines/baseline-latest.json
```

JSON output: `target/benchmarks/mumps/`

### Common CLI Options

All comparison binaries support:

| Flag | Description |
|------|-------------|
| `--ci-only` | Run on CI subset only (10 small matrices) |
| `--rivrs` | Also run rivrs solver for side-by-side comparison |
| `--threads N` | Set thread count for rivrs (and SPRAL via OMP_NUM_THREADS) |
| `--category CAT` | Filter by category: `positive-definite`, `easy-indefinite`, `hard-indefinite` |
| `--compare FILE` | Compare solver timings against a rivrs baseline JSON file |

## Adding a New Solver

To add a comparison against a new solver:

1. Write a Fortran/C driver in `drivers/` that:
   - Reads a matrix from stdin or a file (lower-triangle COO or CSC, 1-indexed)
   - Runs analyze → factor → solve
   - Prints structured results between sentinel markers (see existing drivers)
2. Add a build script in `drivers/`
3. Add a Rust binary in `src/` using `common.rs` shared infrastructure
4. Register the binary in `Cargo.toml`

The subprocess protocol uses sentinel-delimited key-value output, making it easy
to add new solvers without modifying the Rust crate itself.

### Matrix Input Formats

- **CSC** (SPRAL): `n nnz` header, then column pointers (1-indexed), then `row val` pairs
- **COO** (MUMPS): `n nnz` header, then `row col val` lines (1-indexed, lower triangle)

Use `common::format_spral_input()` or `common::format_lower_coo_text()` from the Rust side.

## Results

### MUMPS (sequential)

```bash
comparisons/drivers/build_mumps.sh
cargo run --bin mumps-comparison --release
```

```
=== positive-definite ===
Matrix                                     n        nnz   ana_s   fac_s   slv_s   bwd_err
------------------------------------------------------------------------------------------
AMD/G3_circuit                       1585478    4623152   3.814   6.361   0.703   1.8e-19
Schenk_AFE/af_0_k101                  503625    9027150   0.399   5.055   0.408   3.6e-19
Schenk_AFE/af_shell7                  504855    9046865   0.439   6.482   0.397   3.7e-19
GHS_psdef/apache2                     715176    2766523   2.825  10.024   0.636   2.4e-19
GHS_psdef/bmwcra_1                    148770    5396386   0.593   3.826   0.219   8.8e-19
Oberwolfach/boneS01                   127224    3421188   1.381   2.238   0.172   9.2e-19
Rothberg/cfd2                         123440    1605669   1.111   2.031   0.116   7.9e-19
GHS_psdef/crankseg_1                   52804    5333507   0.321   2.556   0.103   2.5e-17
GHS_psdef/crankseg_2                   63838    7106348   0.378   2.555   0.240   1.1e-17
GHS_psdef/inline_1                    503712   18660027   2.317  10.221   0.517   5.3e-19
GHS_psdef/ldoor                       952203   23737339   1.190  10.793   0.610   3.5e-19
ND/nd12k                               36000    7128473   2.609   9.645   0.233   2.0e-18
ND/nd6k                                18000    3457658   1.188   3.063   0.090   2.6e-18
Um/offshore                           259789    2251231   2.768   4.918   0.302   6.0e-19
DNVS/ship_003                         121728    4103881   1.287   4.102   0.198   6.7e-19
DNVS/shipsec1                         140874    3977139   1.217   2.924   0.177   6.9e-19
DNVS/shipsec5                         179860    5146478   1.603   3.956   0.204   5.9e-19
DNVS/shipsec8                         114919    3384159   1.132   2.547   0.122   7.7e-19
DNVS/thread                            29736    2249892   0.301   2.944   0.092   2.3e-18

=== easy-indefinite ===
Matrix                                     n        nnz   ana_s   fac_s   slv_s   bwd_err
------------------------------------------------------------------------------------------
BenElechi/BenElechi1                  245874    6698185   0.240   3.269   0.213   6.3e-19
Koutsovasilis/F2                       71505    2682895   0.505   1.972   0.038   1.3e-18
PARSEC/H2O                             67024    1141880   1.441  22.720   0.471   1.7e-18
PARSEC/Si10H16                         17077     446500   0.423   2.982   0.078   2.6e-16
PARSEC/Si5H12                          19896     379247   0.449   3.129   0.094   3.9e-17
PARSEC/SiNa                             5743     102265   0.018   0.551   0.021   1.3e-17
Boeing/bcsstk39                        46772    1068033   0.114   0.281   0.071   2.1e-18
GHS_indef/bloweybq                     10001      39996   0.007   0.004   0.001   1.7e-18
GHS_indef/copter2                      55476     407714   0.372   0.520   0.023   4.0e-15
Boeing/crystk02                        13965     491274   0.049   0.173   0.009   2.8e-18
Boeing/crystk03                        24696     887937   0.081   0.535   0.017   2.1e-18
GHS_indef/dawson5                      51537     531157   0.429   0.146   0.014   1.1e-15
GHS_indef/dixmaanl                     60000     179999   0.046   0.022   0.004   4.7e-18
Oberwolfach/filter3D                  106437    1406808   1.240   1.276   0.097   8.4e-19
Oberwolfach/gas_sensor                 66917     885141   0.520   1.116   0.091   9.0e-19
GHS_indef/helm3d01                     32226     230335   0.210   0.272   0.065   6.2e-15
GHS_indef/linverse                     11999      53988   0.009   0.005   0.001   4.1e-18
INPRO/msdoor                          415863   10328399   0.595   2.413   0.179   6.4e-19
ND/nd3k                                 9000    1644345   0.068   1.275   0.074   4.5e-18
Boeing/pwtk                           217918    5926171   0.342   3.954   0.145   4.3e-19
Cunningham/qa8fk                       66127     863353   0.452   1.848   0.093   1.5e-18
Oberwolfach/rail_79841                 79841     316881   0.226   0.116   0.015   5.5e-19
GHS_indef/sparsine                     50000     799494   1.199  19.749   0.281   4.7e-14
GHS_indef/spmsrtls                     29995     129971   0.023   0.013   0.004   1.6e-17
Oberwolfach/t2dal                       4257      20861   0.006   0.004   0.001   1.9e-18
Oberwolfach/t3dh                       79171    2215638   1.274   2.939   0.118   9.0e-19
Cote/vibrobox                          12328     177578   0.146   0.065   0.004   2.9e-18

=== hard-indefinite ===
Matrix                                     n        nnz   ana_s   fac_s   slv_s   bwd_err
------------------------------------------------------------------------------------------
TSOPF/TSOPF_FS_b162_c1                 10798     305732   0.137   0.044   0.002   3.8e-16
TSOPF/TSOPF_FS_b39_c7                  28216     368599   0.161   0.031   0.006   4.6e-16
GHS_indef/aug3dcqp                     35543      77829   0.026   0.067   0.073   6.8e-19
GHS_indef/blockqp1                     60012     340022   0.368   0.047   0.012   2.2e-14
GHS_indef/bratu3d                      27792      88627   0.022   0.379   0.064   6.5e-15
GHS_indef/c-71                         76638     468096   0.511   1.698   0.115   1.3e-18
Schenk_IBMNA/c-big                    345241    1343126   1.494   6.324   0.234   1.9e-17
GHS_indef/cont-201                     80595     239596   0.054   0.189   0.032   7.0e-17
GHS_indef/cont-300                    180895     539396   0.128   1.388   0.110   9.5e-17
GHS_indef/cvxqp3                       17500      69981   0.209   0.907   0.065   8.8e-14
GHS_indef/d_pretok                    182730     885416   0.984   1.339   0.118   1.2e-18
GHS_indef/mario001                     38434     114643   0.122   0.033   0.008   4.1e-18
GHS_indef/ncvxqp1                      12111      40537   0.086   0.435   0.059   3.9e-15
GHS_indef/ncvxqp3                      75000     274982   3.020   8.344   0.183   5.1e-12
GHS_indef/ncvxqp5                      62500     237483   0.805   3.398   0.262   1.9e-14
GHS_indef/ncvxqp7                      87500     312481   2.916  18.192   0.211   2.3e-13
GHS_indef/stokes128                    49666     295938   0.157   0.162   0.018   2.1e-17
GHS_indef/turon_m                     189924     912345   1.194   1.562   0.135   8.4e-18

64/65 completed successfully
```

The failure on astro-ph indicates that it is numerically singular, which is handled correctly by SPRAL + rivrs, but not by MUMPS.

### SPRAL (sequential)


```bash
comparisons/drivers/build_spral.sh
cargo run --bin spral-comparison --release -- --rivrs --threads 1
```

```
=== Comparison: SPRAL vs rivrs (threads=1) ===
Matrix                                     n  spral_fac  rivrs_fac   ratio   spral_be   rivrs_be  sl_st  sl_nd
------------------------------------------------------------------------------------------------------------------
BenElechi/BenElechi1                  245874      1.262      1.468    1.16    6.2e-19    5.9e-19    778   4762
Koutsovasilis/F2                       71505      0.463      0.385    0.83    9.7e-19    9.4e-19    241   1065
PARSEC/H2O                             67024     28.317     26.581    0.94    1.5e-18    1.3e-18      0      0
PARSEC/Si10H16                         17077      2.112      1.637    0.78    4.5e-17    4.5e-17      0      0
PARSEC/Si5H12                          19896      3.156      3.409    1.08    1.8e-17    2.6e-18      0      0
PARSEC/SiNa                             5743      0.217      0.159    0.73    7.5e-18    4.2e-18      0      0
Newman/astro-ph                        16706      0.701      0.361    0.51     8.9e-9    8.0e-18     53    229
Boeing/bcsstk39                        46772      0.132      0.130    0.99    2.6e-18    2.3e-18    126   1172
GHS_indef/bloweybq                     10001      0.002      0.004    1.51    1.6e-18    1.7e-18      1    293
GHS_indef/copter2                      55476      0.467      0.370    0.79    1.2e-15    2.6e-16    207   1092
Boeing/crystk02                        13965      0.087      0.075    0.86    1.9e-18    1.8e-18     62    184
Boeing/crystk03                        24696      0.208      0.178    0.86    1.4e-18    1.4e-18     45    102
GHS_indef/dawson5                      51537      0.122      0.112    0.92    3.8e-16    1.5e-16    103   1337
GHS_indef/dixmaanl                     60000      0.014      0.025    1.86    3.0e-18    5.8e-18      1   1605
Oberwolfach/filter3D                  106437      0.349      0.345    0.99    6.2e-19    6.2e-19    356   2491
Oberwolfach/gas_sensor                 66917      0.562      0.488    0.87    8.0e-19    7.0e-19    335   1727
GHS_indef/helm3d01                     32226      0.174      0.148    0.85    6.3e-16    3.6e-16    117    852
GHS_indef/linverse                     11999      0.003      0.005    1.64    7.4e-18    6.9e-18      1    327
INPRO/msdoor                          415863      0.937      1.134    1.21    6.7e-19    6.5e-19    680  10268
ND/nd3k                                 9000      0.738      0.682    0.92    5.4e-18    2.8e-18      0      0
Boeing/pwtk                           217918      0.938      0.959    1.02    4.4e-19    4.4e-19    727   4404
Cunningham/qa8fk                       66127      0.589      0.566    0.96    8.1e-19    7.8e-19    326   1528
Oberwolfach/rail_79841                 79841      0.037      0.067    1.83    5.9e-19    6.3e-19     12   2188
GHS_indef/sparsine                     50000     32.116     27.338    0.85    8.4e-15    3.0e-15      0      0
GHS_indef/spmsrtls                     29995      0.007      0.013    1.85    1.0e-17    2.0e-17      1    825
Oberwolfach/t2dal                       4257      0.002      0.004    1.54    1.8e-18    2.1e-18      1    110
Oberwolfach/t3dh                       79171      1.583      1.633    1.03    8.8e-19    7.4e-19    346   1298
Cote/vibrobox                          12328      0.050      0.045    0.91    3.1e-18    2.7e-18     52    287
TSOPF/TSOPF_FS_b162_c1                 10798      0.027      0.035    1.28    3.6e-17    5.4e-17     14    790
TSOPF/TSOPF_FS_b39_c7                  28216      0.025      0.034    1.38    1.2e-16    1.6e-16      1   2123
GHS_indef/aug3dcqp                     35543      0.069      0.065    0.95    8.4e-19    1.2e-18     36   1495
GHS_indef/blockqp1                     60012      0.033      0.049    1.49    9.0e-14    9.3e-14      1  19991
GHS_indef/bratu3d                      27792      0.193      0.153    0.79    7.3e-18    1.3e-18    101    750
GHS_indef/c-71                         76638      2.495      3.623    1.45    8.3e-19    1.0e-18    243   2497
Schenk_IBMNA/c-big                    345241      8.459     11.675    1.38    1.3e-17    5.1e-18    904  18712
GHS_indef/cont-201                     80595      0.076      0.082    1.08    2.3e-18    1.9e-18     51   2248
GHS_indef/cont-300                    180895      0.206      0.206    1.00    1.9e-18    1.1e-18    113   5039
GHS_indef/cvxqp3                       17500      0.284      0.138    0.48    1.5e-13    2.0e-14     64    332
GHS_indef/d_pretok                    182730      0.363      0.336    0.92    1.0e-18    9.3e-19    238   4779
GHS_indef/mario001                     38434      0.017      0.032    1.92    1.8e-18    2.0e-18      1   1061
GHS_indef/ncvxqp1                      12111      0.100      0.064    0.64    5.9e-16    1.5e-16     33    243
GHS_indef/ncvxqp3                      75000      2.331      1.412    0.61    1.5e-13    3.3e-14    251   1442
GHS_indef/ncvxqp5                      62500      0.901      0.547    0.61    1.2e-15    2.4e-16    219   1342
GHS_indef/ncvxqp7                      87500      3.244      2.228    0.69    8.7e-14    4.9e-14    277   1651
GHS_indef/stokes128                    49666      0.071      0.071    1.01    9.8e-19    2.2e-18     46   1452
GHS_indef/turon_m                     189924      0.297      0.309    1.04    1.1e-18    2.7e-18    214   4970
AMD/G3_circuit                       1585478      2.306      3.691    1.60    2.7e-19    3.0e-19   1582  48225
Schenk_AFE/af_0_k101                  503625      2.132      2.128    1.00    3.5e-19    3.3e-19   1172  11569
Schenk_AFE/af_shell7                  504855      1.907      1.953    1.02    3.2e-19    3.3e-19   1209  11639
GHS_psdef/apache2                     715176      4.967      5.003    1.01    1.7e-19    2.0e-19   1736  20965
GHS_psdef/bmwcra_1                    148770      1.683      1.535    0.91    6.4e-19    6.0e-19    368    962
Oberwolfach/boneS01                   127224      1.146      1.060    0.93    6.9e-19    6.7e-19    451   3446
Rothberg/cfd2                         123440      0.904      0.850    0.94    5.7e-19    5.5e-19    539   3359
GHS_psdef/crankseg_1                   52804      0.901      0.820    0.91    1.7e-17    2.6e-17     56    145
GHS_psdef/crankseg_2                   63838      1.273      1.138    0.89    3.6e-18    1.0e-17     58    143
GHS_psdef/inline_1                    503712      4.089      4.450    1.09    6.1e-19    4.2e-19   2059   6888
GHS_psdef/ldoor                       952203      2.797      3.375    1.21    4.4e-19    3.8e-19   1654  23174
ND/nd12k                               36000     12.973     10.865    0.84    2.9e-18    1.6e-18      0      0
ND/nd6k                                18000      3.052      2.546    0.83    3.9e-18    2.2e-18      0      0
Um/offshore                           259789      2.347      2.111    0.90    9.0e-19    7.2e-19   1065   5591
DNVS/ship_003                         121728      2.337      2.003    0.86    5.3e-19    5.3e-19    461   1584
DNVS/shipsec1                         140874      1.086      1.037    0.96    7.7e-19    7.4e-19    679   3480
DNVS/shipsec5                         179860      1.905      1.541    0.81    6.3e-19    5.9e-19    859   5009
DNVS/shipsec8                         114919      1.316      0.985    0.75    9.3e-19    8.2e-19    496   2586
DNVS/thread                            29736      3.559      1.492    0.42    1.7e-18    1.4e-18     21     78

65/65 completed successfully
JSON written to: target/benchmarks/spral/spral-benchmark-1772178613.json
```

### SPRAL (`--threads 8`)

```
=== Comparison: SPRAL vs rivrs (threads=8) ===
Matrix                                     n  spral_fac  rivrs_fac   ratio   spral_be   rivrs_be  sl_st  sl_nd
------------------------------------------------------------------------------------------------------------------
BenElechi/BenElechi1                  245874      0.825      0.801    0.97    6.2e-19    5.8e-19    778   4762
Koutsovasilis/F2                       71505      0.287      0.220    0.77    9.8e-19    9.3e-19    241   1065
PARSEC/H2O                             67024     13.067     10.245    0.78    1.6e-18    1.1e-18      0      0
PARSEC/Si10H16                         17077      1.009      0.998    0.99    3.0e-17    5.3e-17      0      0
PARSEC/Si5H12                          19896      1.495      1.920    1.28    8.5e-18    2.2e-18      0      0
PARSEC/SiNa                             5743      0.109      0.097    0.89    1.3e-17    4.0e-18      0      0
Newman/astro-ph                        16706      0.497      0.262    0.53     1.3e-6    1.8e-18     53    229
Boeing/bcsstk39                        46772      0.102      0.100    0.98    2.6e-18    2.3e-18    126   1172
GHS_indef/bloweybq                     10001      0.002      0.004    1.45    1.6e-18    1.7e-18      1    293
GHS_indef/copter2                      55476      0.284      0.204    0.72    9.6e-16    2.4e-16    207   1092
Boeing/crystk02                        13965      0.060      0.044    0.73    2.0e-18    1.8e-18     62    184
Boeing/crystk03                        24696      0.137      0.085    0.62    1.5e-18    1.3e-18     45    102
GHS_indef/dawson5                      51537      0.098      0.087    0.89    2.3e-16    1.5e-16    103   1337
GHS_indef/dixmaanl                     60000      0.014      0.025    1.74    3.0e-18    5.8e-18      1   1605
Oberwolfach/filter3D                  106437      0.256      0.221    0.86    6.3e-19    6.2e-19    356   2491
Oberwolfach/gas_sensor                 66917      0.347      0.258    0.75    8.1e-19    7.0e-19    335   1727
GHS_indef/helm3d01                     32226      0.116      0.091    0.78    7.7e-16    3.8e-16    117    852
GHS_indef/linverse                     11999      0.003      0.005    1.64    7.4e-18    6.9e-18      1    327
INPRO/msdoor                          415863      0.790      0.861    1.09    6.8e-19    6.4e-19    680  10268
ND/nd3k                                 9000      0.359      0.402    1.12    5.5e-18    2.8e-18      0      0
Boeing/pwtk                           217918      0.658      0.562    0.85    4.3e-19    4.3e-19    727   4404
Cunningham/qa8fk                       66127      0.326      0.345    1.06    8.0e-19    7.8e-19    326   1528
Oberwolfach/rail_79841                 79841      0.039      0.068    1.76    5.9e-19    6.3e-19     12   2188
GHS_indef/sparsine                     50000     14.219     11.549    0.81    6.8e-15    3.2e-15      0      0
GHS_indef/spmsrtls                     29995      0.007      0.013    1.82    1.0e-17    2.0e-17      1    825
Oberwolfach/t2dal                       4257      0.003      0.004    1.48    1.8e-18    2.1e-18      1    110
Oberwolfach/t3dh                       79171      0.763      0.789    1.03    8.9e-19    7.7e-19    346   1298
Cote/vibrobox                          12328      0.039      0.030    0.77    3.0e-18    2.6e-18     52    287
TSOPF/TSOPF_FS_b162_c1                 10798      0.031      0.031    0.99    3.9e-17    5.4e-17     14    790
TSOPF/TSOPF_FS_b39_c7                  28216      0.025      0.034    1.35    1.2e-16    1.6e-16      1   2123
GHS_indef/aug3dcqp                     35543      0.055      0.047    0.85    8.4e-19    1.2e-18     36   1495
GHS_indef/blockqp1                     60012      0.033      0.050    1.54    9.0e-14    9.3e-14      1  19991
GHS_indef/bratu3d                      27792      0.122      0.090    0.74    1.2e-17    1.2e-18    101    750
GHS_indef/c-71                         76638      1.402      2.478    1.77    8.3e-19    9.7e-19    243   2497
Schenk_IBMNA/c-big                    345241      4.502      7.858    1.75    1.3e-17    5.2e-18    904  18712
GHS_indef/cont-201                     80595      0.069      0.072    1.04    2.5e-18    1.7e-18     51   2248
GHS_indef/cont-300                    180895      0.174      0.168    0.96    1.8e-18    8.8e-19    113   5039
GHS_indef/cvxqp3                       17500      0.160      0.085    0.53    8.1e-14    3.5e-14     64    332
GHS_indef/d_pretok                    182730      0.273      0.233    0.85    1.0e-18    9.1e-19    238   4779
GHS_indef/mario001                     38434      0.017      0.032    1.84    1.8e-18    2.0e-18      1   1061
GHS_indef/ncvxqp1                      12111      0.067      0.044    0.65    2.9e-16    9.0e-17     33    243
GHS_indef/ncvxqp3                      75000      1.293      0.781    0.60    1.4e-13    3.4e-14    251   1442
GHS_indef/ncvxqp5                      62500      0.521      0.305    0.58    8.0e-16    2.1e-16    219   1342
GHS_indef/ncvxqp7                      87500      1.886      1.197    0.63    9.2e-14    1.8e-14    277   1651
GHS_indef/stokes128                    49666      0.063      0.063    1.00    9.8e-19    2.3e-18     46   1452
GHS_indef/turon_m                     189924      0.251      0.229    0.91    1.1e-18    1.9e-18    214   4970
AMD/G3_circuit                       1585478      1.730      2.876    1.66    2.6e-19    2.2e-19   1582  48225
Schenk_AFE/af_0_k101                  503625      1.418      1.406    0.99    3.5e-19    3.3e-19   1172  11569
Schenk_AFE/af_shell7                  504855      1.351      1.218    0.90    3.2e-19    3.3e-19   1209  11639
GHS_psdef/apache2                     715176      2.841      2.625    0.92    1.7e-19    2.0e-19   1736  20965
GHS_psdef/bmwcra_1                    148770      0.997      0.668    0.67    6.3e-19    6.0e-19    368    962
Oberwolfach/boneS01                   127224      0.698      0.599    0.86    6.9e-19    6.7e-19    451   3446
Rothberg/cfd2                         123440      0.535      0.455    0.85    5.7e-19    5.5e-19    539   3359
GHS_psdef/crankseg_1                   52804      0.522      0.394    0.76    1.8e-17    2.5e-17     56    145
GHS_psdef/crankseg_2                   63838      0.694      0.581    0.84    3.3e-19    6.1e-18     58    143
GHS_psdef/inline_1                    503712      2.579      2.067    0.80    5.9e-19    4.3e-19   2059   6888
GHS_psdef/ldoor                       952203      2.133      2.426    1.14    4.5e-19    3.7e-19   1654  23174
ND/nd12k                               36000      5.767      4.789    0.83    3.0e-18    1.6e-18      0      0
ND/nd6k                                18000      1.505      1.335    0.89    4.1e-18    2.2e-18      0      0
Um/offshore                           259789      1.347      1.220    0.91    9.0e-19    6.9e-19   1065   5591
DNVS/ship_003                         121728      1.200      1.046    0.87    5.3e-19    5.2e-19    461   1584
DNVS/shipsec1                         140874      0.627      0.518    0.83    7.6e-19    6.6e-19    679   3480
DNVS/shipsec5                         179860      1.046      0.905    0.87    6.3e-19    5.5e-19    859   5009
DNVS/shipsec8                         114919      0.745      0.538    0.72    9.2e-19    7.3e-19    496   2586
DNVS/thread                            29736      1.808      0.900    0.50    1.7e-18    1.4e-18     21     78

65/65 completed successfully
```

### Sequential comparison

```bash
python comparisons/analyze_results.py target/benchmarks/spral/spral-benchmark-1772178613.json target/benchmarks/mumps/mumps-benchmark-1772176331.json
```

========================================================================
Multi-Solver Benchmark Analysis
========================================================================
Platform: linux x86_64
Solvers: MUMPS, SPRAL
  MUMPS: 64 matrices, version MUMPS-seq
  SPRAL: 65 matrices, version SPRAL-git

**Multi-solver factor time comparison (seconds):**

| Matrix                              |        n |      MUMPS |      SPRAL |      rivrs |
|-------------------------------------|----------|------------|------------|------------|
| Oberwolfach/t2dal                   |     4257 |      0.004 |      0.002 |      0.004 |
| PARSEC/SiNa                         |     5743 |      0.551 |      0.217 |      0.159 |
| ND/nd3k                             |     9000 |      1.275 |      0.738 |      0.682 |
| GHS_indef/bloweybq                  |    10001 |      0.004 |      0.002 |      0.004 |
| TSOPF/TSOPF_FS_b162_c1              |    10798 |      0.044 |      0.027 |      0.035 |
| GHS_indef/linverse                  |    11999 |      0.005 |      0.003 |      0.005 |
| GHS_indef/ncvxqp1                   |    12111 |      0.435 |      0.100 |      0.064 |
| Cote/vibrobox                       |    12328 |      0.065 |      0.050 |      0.045 |
| Boeing/crystk02                     |    13965 |      0.173 |      0.087 |      0.075 |
| Newman/astro-ph                     |    16706 |       FAIL |      0.701 |      0.361 |
| PARSEC/Si10H16                      |    17077 |      2.982 |      2.112 |      1.637 |
| GHS_indef/cvxqp3                    |    17500 |      0.907 |      0.284 |      0.138 |
| ND/nd6k                             |    18000 |      3.063 |      3.052 |      2.546 |
| PARSEC/Si5H12                       |    19896 |      3.129 |      3.156 |      3.409 |
| Boeing/crystk03                     |    24696 |      0.535 |      0.208 |      0.178 |
| GHS_indef/bratu3d                   |    27792 |      0.379 |      0.193 |      0.153 |
| TSOPF/TSOPF_FS_b39_c7               |    28216 |      0.031 |      0.025 |      0.034 |
| DNVS/thread                         |    29736 |      2.944 |      3.559 |      1.492 |
| GHS_indef/spmsrtls                  |    29995 |      0.013 |      0.007 |      0.013 |
| GHS_indef/helm3d01                  |    32226 |      0.272 |      0.174 |      0.148 |
| GHS_indef/aug3dcqp                  |    35543 |      0.067 |      0.069 |      0.065 |
| ND/nd12k                            |    36000 |      9.645 |     12.973 |     10.865 |
| GHS_indef/mario001                  |    38434 |      0.033 |      0.017 |      0.032 |
| Boeing/bcsstk39                     |    46772 |      0.281 |      0.132 |      0.130 |
| GHS_indef/stokes128                 |    49666 |      0.162 |      0.071 |      0.071 |
| GHS_indef/sparsine                  |    50000 |     19.749 |     32.116 |     27.338 |
| GHS_indef/dawson5                   |    51537 |      0.146 |      0.122 |      0.112 |
| GHS_psdef/crankseg_1                |    52804 |      2.556 |      0.901 |      0.820 |
| GHS_indef/copter2                   |    55476 |      0.520 |      0.467 |      0.370 |
| GHS_indef/dixmaanl                  |    60000 |      0.022 |      0.014 |      0.025 |
| GHS_indef/blockqp1                  |    60012 |      0.047 |      0.033 |      0.049 |
| GHS_indef/ncvxqp5                   |    62500 |      3.398 |      0.901 |      0.547 |
| GHS_psdef/crankseg_2                |    63838 |      2.555 |      1.273 |      1.138 |
| Cunningham/qa8fk                    |    66127 |      1.848 |      0.589 |      0.566 |
| Oberwolfach/gas_sensor              |    66917 |      1.116 |      0.562 |      0.488 |
| PARSEC/H2O                          |    67024 |     22.720 |     28.317 |     26.581 |
| Koutsovasilis/F2                    |    71505 |      1.972 |      0.463 |      0.385 |
| GHS_indef/ncvxqp3                   |    75000 |      8.344 |      2.331 |      1.412 |
| GHS_indef/c-71                      |    76638 |      1.698 |      2.495 |      3.623 |
| Oberwolfach/t3dh                    |    79171 |      2.939 |      1.583 |      1.633 |
| Oberwolfach/rail_79841              |    79841 |      0.116 |      0.037 |      0.067 |
| GHS_indef/cont-201                  |    80595 |      0.189 |      0.076 |      0.082 |
| GHS_indef/ncvxqp7                   |    87500 |     18.192 |      3.244 |      2.228 |
| Oberwolfach/filter3D                |   106437 |      1.276 |      0.349 |      0.345 |
| DNVS/shipsec8                       |   114919 |      2.547 |      1.316 |      0.985 |
| DNVS/ship_003                       |   121728 |      4.102 |      2.337 |      2.003 |
| Rothberg/cfd2                       |   123440 |      2.031 |      0.904 |      0.850 |
| Oberwolfach/boneS01                 |   127224 |      2.238 |      1.146 |      1.060 |
| DNVS/shipsec1                       |   140874 |      2.924 |      1.086 |      1.037 |
| GHS_psdef/bmwcra_1                  |   148770 |      3.826 |      1.683 |      1.535 |
| DNVS/shipsec5                       |   179860 |      3.956 |      1.905 |      1.541 |
| GHS_indef/cont-300                  |   180895 |      1.388 |      0.206 |      0.206 |
| GHS_indef/d_pretok                  |   182730 |      1.339 |      0.363 |      0.336 |
| GHS_indef/turon_m                   |   189924 |      1.562 |      0.297 |      0.309 |
| Boeing/pwtk                         |   217918 |      3.954 |      0.938 |      0.959 |
| BenElechi/BenElechi1                |   245874 |      3.269 |      1.262 |      1.468 |
| Um/offshore                         |   259789 |      4.918 |      2.347 |      2.111 |
| Schenk_IBMNA/c-big                  |   345241 |      6.324 |      8.459 |     11.675 |
| INPRO/msdoor                        |   415863 |      2.413 |      0.937 |      1.134 |
| Schenk_AFE/af_0_k101                |   503625 |      5.055 |      2.132 |      2.128 |
| GHS_psdef/inline_1                  |   503712 |     10.221 |      4.089 |      4.450 |
| Schenk_AFE/af_shell7                |   504855 |      6.482 |      1.907 |      1.953 |
| GHS_psdef/apache2                   |   715176 |     10.024 |      4.967 |      5.003 |
| GHS_psdef/ldoor                     |   952203 |     10.793 |      2.797 |      3.375 |
| AMD/G3_circuit                      |  1585478 |      6.361 |      2.306 |      3.691 |

**Factor time ratio vs rivrs (rivrs_fac / solver_fac):**

| Solver     |   Median |     Mean |   Time-wtd |  Count |
|------------|----------|----------|------------|--------|
| MUMPS      |    0.46x |    0.60x |      0.65x |     64 |
| SPRAL      |    0.94x |    1.03x |      0.94x |     65 |

**SPRAL (threads=1) — 65 matrices:**

| Statistic | Value |
|-----------|-------|
| Median ratio | 0.94x |
| Mean ratio | 1.03x |
| 25th percentile | 0.85x |
| 75th percentile | 1.09x |
| Time-weighted ratio | 0.94x |
| Rivrs faster (<0.90x) | 24/65 (36%) |
| Comparable (0.90-1.10x) | 25/65 (38%) |
| SPRAL faster (>1.10x) | 16/65 (24%) |


### Parallel comparison (`--threads 8`)

```bash
python comparisons/analyze_results.py target/benchmarks/spral/spral-benchmark-1772124268.json
```

========================================================================
SPRAL vs rivrs — Benchmark Analysis
========================================================================
Platform: linux x86_64
rivrs: d1994d0  SPRAL: SPRAL-git
Runs: T8 (8 threads, 65 matrices)

**T8 (threads=8) — 65 matrices:**

| Statistic | Value |
|-----------|-------|
| Median ratio | 0.89x |
| Mean ratio | 0.99x |
| 25th percentile | 0.77x |
| 75th percentile | 1.06x |
| Time-weighted ratio | 0.92x |
| Rivrs faster (<0.90x) | 35/65 (53%) |
| Comparable (0.90-1.10x) | 15/65 (23%) |
| SPRAL faster (>1.10x) | 15/65 (23%) |

**By category:**

| Category | Count | Median | Mean | Time-wtd | Min | Max |
|----------|------:|-------:|-----:|---------:|----:|----:|
| easy-indefinite | 28 | 0.93x | 1.04x | 0.84x | 0.53x | 1.82x |
| hard-indefinite | 18 | 0.94x | 1.03x | 1.26x | 0.53x | 1.84x |
| positive-definite | 19 | 0.86x | 0.88x | 0.90x | 0.50x | 1.66x |

**By application domain:**

| Domain | Count | Median | Mean | Time-wtd |
|--------|------:|-------:|-----:|---------:|
| graph/network | 1 | 0.53x | 0.53x | 0.53x |
| structural | 22 | 0.85x | 0.84x | 0.84x |
| cfd | 3 | 0.85x | 0.86x | 0.82x |
| pde | 11 | 0.89x | 1.02x | 0.86x |
| acoustics | 2 | 0.92x | 0.92x | 1.03x |
| quantum chemistry | 4 | 0.94x | 0.99x | 0.85x |
| model reduction | 6 | 0.95x | 1.12x | 0.92x |
| optimization | 13 | 1.04x | 1.17x | 1.29x |
| power systems | 2 | 1.17x | 1.17x | 1.15x |
| circuit | 1 | 1.66x | 1.66x | 1.66x |

**Phase breakdown (median, seconds):**

| Phase | SPRAL median | rivrs median | Ratio |
|-------|------------:|-----------:|------:|
| Analyse | 0.436 | 0.557 | 1.28x |
| Factor | 0.497 | 0.345 | 0.69x |
| Solve | 0.020 | 0.020 | 1.01x |
| Total | 0.981 | 1.038 | 1.06x |
