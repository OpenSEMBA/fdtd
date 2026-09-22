#ifndef CUDA_TEST_UTILS_H
#define CUDA_TEST_UTILS_H

/*
 * Helpers for CUDA GoogleTest contracts. Only compiled when
 * SEMBA_FDTD_ENABLE_CUDA=ON (see test/CMakeLists.txt + CompileWithCUDA).
 */
#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <vector>

#include "fdtd_cuda.h"

namespace cuda_test {

inline constexpr fdtd_dims3 kTiny = {8, 8, 8};

/* Fortran column-major 0-based linear index (same as fdtd_idx3 on device). */
inline size_t idx3(int i, int j, int k, int nx, int ny)
{
    return (size_t)i + (size_t)j * (size_t)nx + (size_t)k * (size_t)nx * (size_t)ny;
}

inline size_t nelem(fdtd_dims3 d)
{
    return (size_t)d.nx * (size_t)d.ny * (size_t)d.nz;
}

inline fdtd_ibox full_ibox(fdtd_dims3 d)
{
    return fdtd_ibox{0, d.nx - 1, 0, d.ny - 1, 0, d.nz - 1};
}

inline fdtd_ibox empty_ibox()
{
    return fdtd_ibox{2, 1, 0, 0, 0, 0}; /* is > ie → no-op */
}

/* Create ctx, or nullptr if no usable CUDA device. */
inline fdtd_cuda_ctx *TryCreateCuda()
{
    fdtd_cuda_ctx *ctx = fdtd_cuda_create();
    if (!ctx || !fdtd_cuda_ok(ctx)) {
        if (ctx) fdtd_cuda_destroy(ctx);
        return nullptr;
    }
    return ctx;
}

/* Call from TEST body only — GTEST_SKIP must not run inside a helper. */
#define CUDA_REQUIRE_CTX(ctx_var)                                          \
    do {                                                                   \
        (ctx_var) = ::cuda_test::TryCreateCuda();                          \
        if (!(ctx_var))                                                    \
            GTEST_SKIP() << "No usable CUDA device (fdtd_cuda_ok=0)";      \
    } while (0)

inline std::vector<fdtd_real> zeros(size_t n)
{
    return std::vector<fdtd_real>(n, (fdtd_real)0);
}

inline void fill_pattern(std::vector<fdtd_real> &a, fdtd_dims3 d, fdtd_real base)
{
    a.assign(nelem(d), (fdtd_real)0);
    for (int k = 0; k < d.nz; ++k)
        for (int j = 0; j < d.ny; ++j)
            for (int i = 0; i < d.nx; ++i)
                a[idx3(i, j, k, d.nx, d.ny)] =
                    base + (fdtd_real)(i + 10 * j + 100 * k);
}

inline bool exact_eq(const std::vector<fdtd_real> &a, const std::vector<fdtd_real> &b)
{
    if (a.size() != b.size()) return false;
    for (size_t n = 0; n < a.size(); ++n) {
        if (a[n] != b[n]) return false;
    }
    return true;
}

inline bool near_eq(const std::vector<fdtd_real> &a, const std::vector<fdtd_real> &b,
                    fdtd_real atol = (fdtd_real)1e-5)
{
    if (a.size() != b.size()) return false;
    for (size_t n = 0; n < a.size(); ++n) {
        if (std::fabs(a[n] - b[n]) > atol) return false;
    }
    return true;
}

struct SixFields {
    std::vector<fdtd_real> Ex, Ey, Ez, Hx, Hy, Hz;
};

inline SixFields make_six(fdtd_dims3 d, fdtd_real base)
{
    SixFields f;
    fill_pattern(f.Ex, d, base + 0);
    fill_pattern(f.Ey, d, base + 1000);
    fill_pattern(f.Ez, d, base + 2000);
    fill_pattern(f.Hx, d, base + 3000);
    fill_pattern(f.Hy, d, base + 4000);
    fill_pattern(f.Hz, d, base + 5000);
    return f;
}

inline SixFields make_six_zero(fdtd_dims3 d)
{
    SixFields f;
    f.Ex = zeros(nelem(d));
    f.Ey = zeros(nelem(d));
    f.Ez = zeros(nelem(d));
    f.Hx = zeros(nelem(d));
    f.Hy = zeros(nelem(d));
    f.Hz = zeros(nelem(d));
    return f;
}

inline int upload_six(fdtd_cuda_ctx *ctx, const SixFields &f)
{
    return fdtd_cuda_upload_fields(ctx, f.Ex.data(), f.Ey.data(), f.Ez.data(),
                                   f.Hx.data(), f.Hy.data(), f.Hz.data());
}

inline int download_six(fdtd_cuda_ctx *ctx, SixFields &f)
{
    return fdtd_cuda_download_fields(ctx, f.Ex.data(), f.Ey.data(), f.Ez.data(),
                                     f.Hx.data(), f.Hy.data(), f.Hz.data());
}

inline int upload_six_box(fdtd_cuda_ctx *ctx, const SixFields &f, const fdtd_ibox &b)
{
    return fdtd_cuda_upload_fields_box(ctx, f.Ex.data(), f.Ey.data(), f.Ez.data(),
                                       f.Hx.data(), f.Hy.data(), f.Hz.data(),
                                       &b, &b, &b, &b, &b, &b);
}

inline int download_six_box(fdtd_cuda_ctx *ctx, SixFields &f, const fdtd_ibox &b)
{
    return fdtd_cuda_download_fields_box(ctx, f.Ex.data(), f.Ey.data(), f.Ez.data(),
                                         f.Hx.data(), f.Hy.data(), f.Hz.data(),
                                         &b, &b, &b, &b, &b, &b);
}

inline bool in_box(int i, int j, int k, const fdtd_ibox &b)
{
    return i >= b.is && i <= b.ie && j >= b.js && j <= b.je && k >= b.ks && k <= b.ke;
}

/* Host mirror of k_advance_ex2 (uniform dims). */
inline void host_advance_ex(std::vector<fdtd_real> &Ex, const std::vector<fdtd_real> &Hy,
                            const std::vector<fdtd_real> &Hz, const std::vector<fdtd_media> &Mi,
                            const std::vector<fdtd_real> &g1, const std::vector<fdtd_real> &g2,
                            const std::vector<fdtd_real> &Idyh, const std::vector<fdtd_real> &Idzh,
                            fdtd_dims3 d, fdtd_ibox s)
{
    for (int k = s.ks; k <= s.ke; ++k)
        for (int j = s.js; j <= s.je; ++j)
            for (int i = s.is; i <= s.ie; ++i) {
                size_t id = idx3(i, j, k, d.nx, d.ny);
                int medio = (int)Mi[id];
                fdtd_real curl =
                    (Hz[idx3(i, j, k, d.nx, d.ny)] - Hz[idx3(i, j - 1, k, d.nx, d.ny)]) * Idyh[j] -
                    (Hy[idx3(i, j, k, d.nx, d.ny)] - Hy[idx3(i, j, k - 1, d.nx, d.ny)]) * Idzh[k];
                Ex[id] = g1[medio] * Ex[id] + g2[medio] * curl;
            }
}

/* Host mirror of k_advance_hx2 (uniform dims). */
inline void host_advance_hx(std::vector<fdtd_real> &Hx, const std::vector<fdtd_real> &Ey,
                            const std::vector<fdtd_real> &Ez, const std::vector<fdtd_media> &Mi,
                            const std::vector<fdtd_real> &gm1, const std::vector<fdtd_real> &gm2,
                            const std::vector<fdtd_real> &Idze, const std::vector<fdtd_real> &Idye,
                            fdtd_dims3 d, fdtd_ibox s)
{
    for (int k = s.ks; k <= s.ke; ++k)
        for (int j = s.js; j <= s.je; ++j)
            for (int i = s.is; i <= s.ie; ++i) {
                size_t id = idx3(i, j, k, d.nx, d.ny);
                int medio = (int)Mi[id];
                fdtd_real curl =
                    (Ey[idx3(i, j, k + 1, d.nx, d.ny)] - Ey[idx3(i, j, k, d.nx, d.ny)]) * Idze[k] -
                    (Ez[idx3(i, j + 1, k, d.nx, d.ny)] - Ez[idx3(i, j, k, d.nx, d.ny)]) * Idye[j];
                Hx[id] = gm1[medio] * Hx[id] + gm2[medio] * curl;
            }
}

/*
 * Host mirror of k_cpml for one absolute cell (use_fixed_medio path).
 * Field/Ha use remapped indices; Psi uses (i-psi_xi,...).
 */
inline void host_cpml_cell(fdtd_real &Field, fdtd_real &Psi, fdtd_real Ha0, fdtd_real Ha1,
                           fdtd_real Pb, fdtd_real Pc, fdtd_real G, int h_sign)
{
    fdtd_real dH = Ha0 - Ha1;
    Psi = Pb * Psi + dH * Pc;
    Field = Field + (fdtd_real)h_sign * G * Psi;
}

} // namespace cuda_test

#endif /* CUDA_TEST_UTILS_H */
