/speckit.specify implement Phase 8.2 in ssids-plan.md.  Treat any code as notional and descriptive of *what* functionality is generally needed, not *how* to implement it. The actual spec should be based on the current codebase status and the high-level progress through the plan - mark anything that's unclear or conflicting for further clarification. Make sure to note algorithm references and the locations of the corresponding Markdown files in /workspace/rivrs-linalg/references/ssids/*.md

Can you review PR #22, looking for code quality, consistency with the plan outlined in docs/ssids-plan.md, Rust best practices, duplication of code or functionality, and code clarity? Also be sure to evaluate test comprehensiveness following the core principle that correctness is the critical goal of the implementation - do we have sufficient test coverage to have high confidence in the current implementation?

%%

Can you do a comprehensive review of the examples?  If we're coming up to release, which of these do we really need to keep, and which can be removed or condensed to a minimal set?  Then, let's also add a README.md file to this directory including an explanation of what each example does and how to run it.

%%%

Can you update the README.md file?  At a minimum it should have the basic Rust API and CLI examples for how to run the basic tests, the CI subset, the full SuiteSparse set, and benchmarking (vs SPRAL, and parallel).  Otherwise let's keep the README.md concise since the project is still private and a work in progress.

%%%

Let's revisit the "CI subset", which actually contains some fairly large matrices and hence needs to be ignored, which is counter to the idea of a "CI subset" in the first place. Here is the current performance report - can we select a set of ~10 matrices that should run quickly in CI with up to 4 threads and spans the easy-indefinite, hard-indefinite, and positive-definite cases?  All of the data is in test-data/suitesparse, so I think that swapping these should just be a matter of copying new matrices over and updating metadata.json

=== Comparison: SPRAL vs rivrs (threads=1) ===
Matrix                                     n  spral_fac  rivrs_fac   ratio   spral_be   rivrs_be
----------------------------------------------------------------------------------------------------
BenElechi/BenElechi1                  245874      1.053      1.349    1.28    5.7e-16    5.5e-19
Koutsovasilis/F2                       71505      0.380      0.438    1.15    6.1e-16    9.7e-19
PARSEC/H2O                             67024     27.760     37.070    1.34    2.5e-15    1.3e-18
PARSEC/Si10H16                         17077      2.078      2.305    1.11    1.1e-13    5.0e-17
PARSEC/Si5H12                          19896      3.114      5.198    1.67    1.0e-13    2.6e-18
PARSEC/SiNa                             5743      0.203      0.236    1.17    1.2e-14    4.0e-18
Newman/astro-ph                        16706      0.693      0.431    0.62     1.3e-7    5.9e-18
Boeing/bcsstk39                        46772      0.131      0.149    1.14    6.9e-16    2.2e-18
GHS_indef/bloweybq                     10001      0.002      0.006    2.44    8.3e-17    1.9e-18
GHS_indef/copter2                      55476      0.460      0.444    0.97    2.5e-12    2.9e-16
Boeing/crystk02                        13965      0.086      0.106    1.24    7.6e-16    1.7e-18
Boeing/crystk03                        24696      0.209      0.224    1.07    7.2e-16    1.3e-18
GHS_indef/dawson5                      51537      0.117      0.128    1.09    9.1e-13    9.9e-17
GHS_indef/dixmaanl                     60000      0.014      0.043    3.19    9.1e-15    3.1e-18
Oberwolfach/filter3D                  106437      0.345      0.424    1.23    6.3e-16    6.0e-19
Oberwolfach/gas_sensor                 66917      0.556      0.608    1.09    6.9e-16    6.7e-19
GHS_indef/helm3d01                     32226      0.172      0.183    1.07    3.2e-13    2.8e-16
GHS_indef/linverse                     11999      0.003      0.010    3.21    6.5e-15    9.5e-18
INPRO/msdoor                          415863      0.921      0.982    1.07    9.6e-16    6.4e-19
ND/nd3k                                 9000      0.738      0.780    1.06    2.7e-15    2.8e-18
Boeing/pwtk                           217918      0.928      1.064    1.15    3.9e-16    4.1e-19
Cunningham/qa8fk                       66127      0.554      0.688    1.24    1.4e-15    7.5e-19
Oberwolfach/rail_79841                 79841      0.036      0.206    5.68    6.4e-16    5.7e-19
GHS_indef/sparsine                     50000     31.599     35.627    1.13    1.6e-12    2.9e-15
GHS_indef/spmsrtls                     29995      0.007      0.023    3.32    9.4e-15    1.4e-17
Oberwolfach/t2dal                       4257      0.003      0.009    3.58    2.3e-16    2.1e-18
Oberwolfach/t3dh                       79171      1.551      1.825    1.18    1.3e-15    7.3e-19
Cote/vibrobox                          12328      0.051      0.055    1.10    2.0e-16    2.9e-18
TSOPF/TSOPF_FS_b162_c1                 10798      0.028      0.041    1.48    7.1e-15    3.7e-17
TSOPF/TSOPF_FS_b39_c7                  28216      0.025      0.039    1.57    1.8e-14    1.3e-16
GHS_indef/aug3dcqp                     35543      0.069      0.135    1.95    1.8e-16    6.3e-19
GHS_indef/blockqp1                     60012      0.031      0.057    1.82    6.0e-13    9.2e-14
GHS_indef/bratu3d                      27792      0.177      0.171    0.97    6.3e-14    1.5e-18
GHS_indef/c-71                         76638      2.387     37.563   15.74    1.2e-16    1.7e-18
Schenk_IBMNA/c-big                    345241      8.392     90.549   10.79    2.2e-16    1.1e-17
GHS_indef/cont-201                     80595      0.075      0.094    1.25    5.9e-14    1.8e-18
GHS_indef/cont-300                    180895      0.201      0.228    1.14    8.4e-14    1.0e-18
GHS_indef/cvxqp3                       17500      0.256      0.160    0.63    1.7e-10    1.0e-14
GHS_indef/d_pretok                    182730      0.340      0.370    1.09    5.1e-16    8.9e-19
GHS_indef/mario001                     38434      0.016      0.091    5.57    2.0e-15    2.1e-18
GHS_indef/ncvxqp1                      12111      0.098      0.074    0.75    5.6e-14    1.1e-16
GHS_indef/ncvxqp3                      75000      2.270      1.530    0.67    4.4e-10    3.6e-14
GHS_indef/ncvxqp5                      62500      0.904      0.661    0.73    2.4e-12    2.5e-16
GHS_indef/ncvxqp7                      87500      3.196      2.341    0.73    5.4e-10    8.5e-14
GHS_indef/stokes128                    49666      0.069      0.089    1.30    1.5e-15    1.1e-18
GHS_indef/turon_m                     189924      0.290      0.345    1.19    4.5e-16    1.1e-17
AMD/G3_circuit                       1585478      2.252      3.331    1.48    3.9e-15    2.8e-19
Schenk_AFE/af_0_k101                  503625      2.085      2.247    1.08    1.1e-15    3.2e-19
Schenk_AFE/af_shell7                  504855      1.885      2.082    1.10    6.0e-16    2.9e-19
GHS_psdef/apache2                     715176      4.905      5.049    1.03    1.4e-15    2.1e-19
GHS_psdef/bmwcra_1                    148770      1.661      1.903    1.15    1.0e-15    5.9e-19
Oberwolfach/boneS01                   127224      1.141      1.245    1.09    9.5e-16    6.5e-19
Rothberg/cfd2                         123440      0.887      1.018    1.15    1.0e-15    5.0e-19
GHS_psdef/crankseg_1                   52804      0.894      1.024    1.14    1.9e-15    2.5e-17
GHS_psdef/crankseg_2                   63838      1.258      1.440    1.14    4.7e-16    1.6e-17
GHS_psdef/inline_1                    503712      4.152      4.545    1.09    7.8e-16    4.3e-19
GHS_psdef/ldoor                       952203      2.707      2.846    1.05    9.1e-16    3.6e-19
ND/nd12k                               36000     12.900     13.713    1.06    3.6e-15    1.6e-18
ND/nd6k                                18000      3.024      2.979    0.98    3.4e-15    2.1e-18
Um/offshore                           259789      2.315      2.580    1.11    7.5e-16    7.2e-19
DNVS/ship_003                         121728      2.301      2.415    1.05    5.9e-16    5.2e-19
DNVS/shipsec1                         140874      1.070      1.250    1.17    5.3e-16    7.5e-19
DNVS/shipsec5                         179860      1.877      1.881    1.00    7.9e-16    5.8e-19
DNVS/shipsec8                         114919      1.289      1.175    0.91    4.3e-16    8.1e-19
DNVS/thread                            29736      3.441      2.228    0.65    2.0e-15    1.6e-18
