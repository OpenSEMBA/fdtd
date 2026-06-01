#ifndef UTILS_M_H
#define UTILS_M_H

#include <complex>
#include <vector>

namespace utils_m {

using real_type = double;
using complex_type = std::complex<double>;

std::vector<std::vector<std::vector<complex_type>>> sumQComponents(
    const std::vector<std::vector<std::vector<std::vector<complex_type>>>>& a);

std::vector<complex_type> dotmatrixmul(const std::vector<std::vector<std::vector<complex_type>>>& a,
                                       const std::vector<std::vector<complex_type>>& b);

std::vector<real_type> getEigenValues(const std::vector<std::vector<real_type>>& matrix);

std::vector<std::vector<real_type>> inv(const std::vector<std::vector<real_type>>& A);

std::vector<std::vector<real_type>> element_wise_invert(int n,
                                                        const std::vector<std::vector<real_type>>& x);

inline std::vector<real_type> matmul_vec(const std::vector<std::vector<real_type>>& m,
                                         const std::vector<real_type>& v) {
    const int n = static_cast<int>(m.size());
    std::vector<real_type> res(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[static_cast<size_t>(i)] += m[static_cast<size_t>(i)][static_cast<size_t>(j)] *
                                        v[static_cast<size_t>(j)];
        }
    }
    return res;
}

inline std::vector<std::vector<real_type>> matmul(const std::vector<std::vector<real_type>>& a,
                                                  const std::vector<std::vector<real_type>>& b) {
    const int n = static_cast<int>(a.size());
    std::vector<std::vector<real_type>> res(n, std::vector<real_type>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                res[static_cast<size_t>(i)][static_cast<size_t>(j)] +=
                    a[static_cast<size_t>(i)][static_cast<size_t>(k)] *
                    b[static_cast<size_t>(k)][static_cast<size_t>(j)];
            }
        }
    }
    return res;
}

inline double dot_product(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    const int n = static_cast<int>(std::min(a.size(), b.size()));
    for (int i = 0; i < n; ++i) {
        s += a[static_cast<size_t>(i)] * b[static_cast<size_t>(i)];
    }
    return s;
}

} // namespace utils_m

#endif
