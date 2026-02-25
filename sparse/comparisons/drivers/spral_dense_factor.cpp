// Standalone driver to factor a dense frontal matrix using SPRAL's APTP kernel.
//
// Reads a dense m x m symmetric matrix (and k = num_fully_summed) from stdin,
// calls SPRAL's ldlt_app_factor, and reports factorization statistics.
//
// Input format (stdin):
//   m k
//   a[0,0] a[0,1] ... a[0,m-1]
//   a[1,0] a[1,1] ... a[1,m-1]
//   ...
//   a[m-1,0] a[m-1,1] ... a[m-1,m-1]
//
// Output (stdout):
//   SPRAL_DENSE_FACTOR
//   m: <m>
//   k: <k>
//   nelim: <nelim>
//   block_size: <block_size>
//   threshold: <u>
//   PERM (pivot permutation, 0-indexed):
//   perm[0] perm[1] ... perm[k-1]
//   D (2*k entries: d[2*i] = D^{-1}_{ii}, d[2*i+1] = off-diag):
//   d[0] d[1] d[2] d[3] ...
//   PIVOT_SUMMARY:
//   num_1x1: <count>
//   num_2x2: <count>
//   num_delayed: <count>
//
// Build (after running comparisons/drivers/build_spral.sh):
//   See comparisons/drivers/build_spral_dense_factor.sh
//
// Usage:
//   ./spral_dense_factor < /tmp/frontal_matrix.txt

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>

// SPRAL headers
#include "ssids/cpu/cpu_iface.hxx"
#include "ssids/cpu/BuddyAllocator.hxx"
#include "ssids/cpu/Workspace.hxx"
#include "ssids/cpu/kernels/ldlt_app.hxx"

using namespace spral::ssids::cpu;

int main(int argc, char* argv[]) {
    // Parse optional arguments
    double u = 0.01;        // threshold (matching our default)
    int block_size = 256;   // outer block size
    bool verbose = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-u") == 0 && i+1 < argc) {
            u = atof(argv[++i]);
        } else if (strcmp(argv[i], "-b") == 0 && i+1 < argc) {
            block_size = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-v") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            fprintf(stderr, "Usage: %s [-u threshold] [-b block_size] [-v] < input.txt\n", argv[0]);
            fprintf(stderr, "  -u threshold   APTP pivot threshold (default: 0.01)\n");
            fprintf(stderr, "  -b block_size  Outer block size (default: 256)\n");
            fprintf(stderr, "  -v             Verbose output\n");
            return 0;
        }
    }

    // Read matrix dimensions
    int m, k;
    if (scanf("%d %d", &m, &k) != 2) {
        fprintf(stderr, "ERROR: Failed to read m and k\n");
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Reading %d x %d matrix (k=%d)...\n", m, m, k);
    }

    // Allocate matrix in column-major format (SPRAL convention)
    // lda = m (no padding for simplicity)
    int lda = m;
    std::vector<double> a(lda * m, 0.0);

    // Read matrix row by row, store in column-major
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double val;
            if (scanf("%lf", &val) != 1) {
                fprintf(stderr, "ERROR: Failed to read a[%d,%d]\n", i, j);
                return 1;
            }
            a[j * lda + i] = val;  // column-major: a[col * lda + row]
        }
    }

    if (verbose) {
        fprintf(stderr, "Matrix read successfully.\n");
        // Print some diagnostics
        double max_val = 0.0;
        double max_diag = 0.0;
        int zero_diag = 0;
        for (int i = 0; i < m; i++) {
            double d = fabs(a[i * lda + i]);
            if (d == 0.0) zero_diag++;
            if (d > max_diag) max_diag = d;
            for (int j = 0; j < m; j++) {
                double v = fabs(a[j * lda + i]);
                if (v > max_val) max_val = v;
            }
        }
        fprintf(stderr, "  max |a_ij| = %.6e\n", max_val);
        fprintf(stderr, "  max |a_ii| = %.6e\n", max_diag);
        fprintf(stderr, "  zero diags = %d\n", zero_diag);
    }

    // Initialize pivot permutation as identity
    std::vector<int> perm(k);
    for (int i = 0; i < k; i++) {
        perm[i] = i;
    }

    // Allocate D array (2 * k entries)
    std::vector<double> d(2 * k, 0.0);

    // Allocate contribution/update matrix (m-k) x (m-k)
    // SPRAL stores the Schur complement update here
    int upd_size = m - k;
    std::vector<double> upd(upd_size > 0 ? upd_size * upd_size : 1, 0.0);

    // Set up SPRAL options
    struct cpu_factor_options options;
    options.print_level = 0;
    options.action = true;
    options.small = 1e-20;
    options.u = u;
    options.multiplier = 1.1;
    options.small_subtree_threshold = 4 * 1024 * 1024;
    options.cpu_block_size = block_size;
    options.pivot_method = PivotMethod::app_block;
    options.failed_pivot_method = FailedPivotMethod::tpp;

    // Set up allocator and workspace
    // BuddyAllocator needs enough memory for the backup copy
    size_t alloc_size = static_cast<size_t>(m) * k + 1024 * 1024;  // generous
    BuddyAllocator<double, std::allocator<double>> alloc(alloc_size);

    // One workspace per "thread" (we run single-threaded)
    size_t work_sz = 2 * static_cast<size_t>(m) * m * sizeof(double) + 1024;
    std::vector<Workspace> work;
    work.emplace_back(work_sz);

    if (verbose) {
        fprintf(stderr, "Calling ldlt_app_factor(m=%d, n=%d, lda=%d, u=%.4f, blk=%d)...\n",
                m, k, lda, u, block_size);
    }

    // Call SPRAL's APTP kernel
    int nelim = ldlt_app_factor(
        m,              // total rows
        k,              // fully-summed cols
        perm.data(),    // pivot permutation
        a.data(),       // dense matrix (column-major)
        lda,            // leading dimension
        d.data(),       // D output
        0.0,            // beta (no prior update)
        upd.data(),     // update matrix
        upd_size > 0 ? upd_size : 1,  // ldupd
        options,        // factor options
        work,           // workspace
        alloc           // allocator
    );

    if (verbose) {
        fprintf(stderr, "Done. nelim = %d\n", nelim);
    }

    // Count pivot types
    int num_1x1 = 0, num_2x2 = 0;
    {
        int j = 0;
        while (j < nelim) {
            if (j + 1 < nelim && !std::isfinite(d[2*j + 2])) {
                // 2x2 pivot at (j, j+1): d[2*(j+1)] = infinity marker
                num_2x2++;
                j += 2;
            } else {
                // 1x1 pivot
                num_1x1++;
                j += 1;
            }
        }
    }
    int num_delayed = k - nelim;

    // Output results
    printf("SPRAL_DENSE_FACTOR\n");
    printf("m: %d\n", m);
    printf("k: %d\n", k);
    printf("nelim: %d\n", nelim);
    printf("block_size: %d\n", block_size);
    printf("threshold: %.6e\n", u);

    printf("PERM:\n");
    for (int i = 0; i < k; i++) {
        printf("%d", perm[i]);
        if (i + 1 < k) printf(" ");
    }
    printf("\n");

    printf("D:\n");
    for (int i = 0; i < 2 * nelim; i++) {
        printf("%.17e", d[i]);
        if (i + 1 < 2 * nelim) printf(" ");
    }
    printf("\n");

    printf("PIVOT_SUMMARY:\n");
    printf("num_1x1: %d\n", num_1x1);
    printf("num_2x2: %d\n", num_2x2);
    printf("num_delayed: %d\n", num_delayed);

    // If verbose, also dump the factored L matrix (first few cols)
    if (verbose) {
        int show = std::min(nelim, 10);
        fprintf(stderr, "\nFactored L (first %d cols, lower triangle):\n", show);
        for (int i = 0; i < std::min(m, 20); i++) {
            for (int j = 0; j < show; j++) {
                fprintf(stderr, "%12.4e ", a[j * lda + i]);
            }
            fprintf(stderr, "\n");
        }

        fprintf(stderr, "\nD entries (first %d):\n", std::min(2*nelim, 20));
        for (int i = 0; i < std::min(2*nelim, 20); i += 2) {
            fprintf(stderr, "  d[%d] = %.6e, d[%d] = %.6e\n",
                    i, d[i], i+1, d[i+1]);
        }
    }

    return (nelim >= 0) ? 0 : 1;
}
