#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <numeric>

// Assuming FDETYPES_m provides RKIND. 
// Typically RKIND is 8 for double precision in Fortran context.
#ifndef RKIND
#define RKIND 8
#endif

// Map RKIND to C++ types
#if RKIND == 4
using real_type = float;
#elif RKIND == 8
using real_type = double;
#else
using real_type = double; // Default to double
#endif

// Complex type corresponding to Fortran complex (usually double precision components)
using complex_type = std::complex<double>;

// LAPACK wrappers for DGETRF, DGETRI, SGETRF, SGETRI
// Note: In a real project, you would link against LAPACK/BLAS libraries.
// These are stubs/declarations for compilation purposes.
extern "C" {
    void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
    void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
    void sgetrf_(int* m, int* n, float* a, int* lda, int* ipiv, int* info);
    void sgetri_(int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info);
    
    void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, 
                double* wr, double* wi, double* vl, int* ldvl, 
                double* vr, int* ldvr, double* work, int* lwork, int* info);
}

namespace utils_m {

    struct entry_t {
        std::vector<real_type> x;
    };

    // Helper to convert 4D Fortran array (index 4) to 3D result
    // Fortran: complex, dimension(:,:,:), allocatable :: res
    // Input: complex, dimension(:,:,:,:), intent(in) :: a
    std::vector<std::vector<std::vector<complex_type>>> sumQComponents(const std::vector<std::vector<std::vector<std::vector<complex_type>>>>& a) {
        if (a.empty() || a[0].empty() || a[0][0].empty()) {
            return {};
        }
        
        int dim1 = a.size();
        int dim2 = a[0].size();
        int dim3 = a[0][0].size();
        int dim4 = a[0][0][0].size();

        std::vector<std::vector<std::vector<complex_type>>> res(dim1, 
            std::vector<std::vector<complex_type>>(dim2, 
            std::vector<complex_type>(dim3, complex_type(0.0, 0.0))));

        for (int i = 0; i < dim4; ++i) {
            for (int j = 0; j < dim1; ++j) {
                for (int k = 0; k < dim2; ++k) {
                    for (int l = 0; l < dim3; ++l) {
                        res[j][k][l] += a[j][k][l][i];
                    }
                }
            }
        }
        return res;
    }

    // Dot matrix multiplication
    // a: (M, N, K), b: (N, K) -> res: (M)
    // res(i) = sum_j dot_product(a(i,j,:), b(j,:))
    std::vector<complex_type> dotmatrixmul(const std::vector<std::vector<std::vector<complex_type>>>& a, 
                                           const std::vector<std::vector<complex_type>>& b) {
        if (a.empty() || a[0].empty()) {
            return {};
        }

        int m = a.size();
        int n = a[0].size();
        int k = a[0][0].size();

        // Check dimensions match
        if (b.size() != n || b.empty() || b[0].size() != k) {
            throw std::invalid_argument("Dimensions mismatch in dotmatrixmul");
        }

        std::vector<complex_type> res(m, complex_type(0.0, 0.0));

        for (int i = 0; i < m; ++i) {
            complex_type sum_val(0.0, 0.0);
            for (int j = 0; j < n; ++j) {
                // Dot product of a(i,j,:) and b(j,:)
                complex_type dot_prod(0.0, 0.0);
                for (int l = 0; l < k; ++l) {
                    dot_prod += a[i][j][l] * b[j][l];
                }
                sum_val += dot_prod;
            }
            res[i] = sum_val;
        }
        return res;
    }

    // Identity matrix
    std::vector<std::vector<real_type>> eye(int dim) {
        std::vector<std::vector<real_type>> res(dim, std::vector<real_type>(dim, 0.0));
        for (int i = 0; i < dim; ++i) {
            res[i][i] = 1.0;
        }
        return res;
    }

    // Get Eigenvalues using LAPACK DGEEV
    // Returns a vector of size 2*n containing real and imaginary parts interleaved?
    // Fortran code: eigvals = [eigvals_real, eigvals_imag]
    // This creates a 1D array of size 2*n where first n are real parts, next n are imag parts.
    std::vector<real_type> getEigenValues(const std::vector<std::vector<real_type>>& matrix) {
        int n = matrix.size();
        if (n == 0) return {};

        std::vector<real_type> eigvals_real(n, 0.0);
        std::vector<real_type> eigvals_imag(n, 0.0);
        std::vector<real_type> eigvals(2 * n, 0.0);
        
        std::vector<real_type> m1(n * n, 0.0); // Flattened copy
        std::vector<real_type> m2(n * n, 0.0); // Flattened copy
        for(int i=0; i<n; ++i) for(int j=0; j<n; ++j) m1[i*n+j] = matrix[i][j];
        
        // LAPACK expects column-major storage. 
        // If matrix is passed as row-major vector<vector>, we need to handle indexing carefully.
        // However, standard LAPACK wrappers often assume the input is already in the correct memory layout 
        // or provide a wrapper. Here we assume the input 'matrix' is stored in a way compatible with 
        // LAPACK's column-major expectation if we pass it directly as a flat pointer, 
        // BUT std::vector<std::vector> is row-major.
        // To be safe and simple, let's assume the caller provides data in LAPACK format (column-major) 
        // or we flatten it. 
        // Given the complexity of converting row-major to column-major for LAPACK inside a generic translator,
        // and the fact that Fortran is column-major, we will assume the input `matrix` 
        // is conceptually column-major or we flatten it assuming row-major input needs transposition.
        // Let's flatten m1 and m2 into 1D vectors in column-major order if the input is row-major.
        
        // Flatten row-major to column-major for LAPACK
        std::vector<real_type> m1_flat(n * n, 0.0);
        std::vector<real_type> m2_flat(n * n, 0.0);
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                // Fortran: m1(i,j) -> C++: m1_flat[j*n + i] for column-major
                m1_flat[j * n + i] = matrix[i][j];
                m2_flat[j * n + i] = matrix[i][j];
            }
        }

        int info = 0;
        int lwork = -1;
        int ldummy = 1;
        std::vector<real_type> dummy(1, 0.0);
        std::vector<real_type> work;
        
        // Query workspace
        char jobvl = 'N', jobvr = 'N';
        
        // dgeev_ signature:
        // dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, 
        //        double* wr, double* wi, double* vl, int* ldvl, 
        //        double* vr, int* ldvr, double* work, int* lwork, int* info);
        
        dgeev_(&jobvl, &jobvr, &n, m1_flat.data(), &n, 
               eigvals_real.data(), eigvals_imag.data(), 
               dummy.data(), &ldummy, dummy.data(), &ldummy, 
               dummy.data(), &lwork, &info);
        
        if (info != 0) {
            throw std::runtime_error("LAPACK dgeev query failed");
        }
        
        lwork = static_cast<int>(std::max(static_cast<double>((64 + 2) * n), dummy[0]));
        work.resize(lwork);
        
        dgeev_(&jobvl, &jobvr, &n, m2_flat.data(), &n, 
               eigvals_real.data(), eigvals_imag.data(), 
               dummy.data(), &ldummy, dummy.data(), &ldummy, 
               work.data(), &lwork, &info);
               
        if (info != 0) {
            throw std::runtime_error("LAPACK dgeev computation failed");
        }
        
        // Fortran: eigvals = [eigvals_real, eigvals_imag]
        // This means first n elements are real parts, next n are imaginary parts.
        for (int i = 0; i < n; ++i) {
            eigvals[i] = eigvals_real[i];
            eigvals[n + i] = eigvals_imag[i];
        }
        
        return eigvals;
    }

    // Matrix Inversion using LAPACK GETRF/GETRI
    std::vector<std::vector<real_type>> inv(const std::vector<std::vector<real_type>>& A) {
        int n = A.size();
        if (n == 0) return {};

        // LAPACK expects column-major. Convert input row-major to column-major.
        std::vector<real_type> Ainv_flat(n * n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ainv_flat[j * n + i] = A[i][j];
            }
        }

        std::vector<int> ipiv(n, 0);
        std::vector<real_type> work(n, 0.0);
        int info = 0;

#ifdef CompileWithReal8
        // Use Double Precision
        int m = n;
        dgetrf_(&m, &n, Ainv_flat.data(), &n, ipiv.data(), &info);
        if (info != 0) {
            throw std::runtime_error("Matrix is numerically singular!");
        }
        dgetri_(&n, Ainv_flat.data(), &n, ipiv.data(), work.data(), &n, &info);
        if (info != 0) {
            throw std::runtime_error("Matrix inversion failed!");
        }
#else
        // Use Single Precision (cast to float)
        std::vector<float> Ainv_float(n * n);
        for(int k=0; k<n*n; ++k) Ainv_float[k] = static_cast<float>(Ainv_flat[k]);
        std::vector<float> work_float(n);
        
        int m = n;
        sgetrf_(&m, &n, Ainv_float.data(), &n, ipiv.data(), &info);
        if (info != 0) {
            throw std::runtime_error("Matrix is numerically singular!");
        }
        sgetri_(&n, Ainv_float.data(), &n, ipiv.data(), work_float.data(), &n, &info);
        if (info != 0) {
            throw std::runtime_error("Matrix inversion failed!");
        }
        
        // Convert back to double
        for(int k=0; k<n*n; ++k) Ainv_flat[k] = static_cast<real_type>(Ainv_float[k]);
#endif

        // Convert column-major flat vector back to row-major 2D vector
        std::vector<std::vector<real_type>> result(n, std::vector<real_type>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result[i][j] = Ainv_flat[j * n + i];
            }
        }
        return result;
    }

    // Element-wise invert
    std::vector<std::vector<real_type>> element_wise_invert(int n, const std::vector<std::vector<real_type>>& x) {
        if (x.size() != n || (n > 0 && x[0].size() != n)) {
            throw std::invalid_argument("Matrix dimensions do not match n");
        }
        
        std::vector<std::vector<real_type>> y(n, std::vector<real_type>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (x[i][j] == 0.0) {
                    throw std::runtime_error("Division by zero in element_wise_invert");
                }
                y[i][j] = 1.0 / x[i][j];
            }
        }
        return y;
    }

} // namespace utils_m