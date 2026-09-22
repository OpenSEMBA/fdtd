#ifndef CUDA_YEE_TESTS_H
#define CUDA_YEE_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/*
 * Documents fdtd_cuda_advance_ex: g1/g2, media index, Idyh/Idzh, Hy/Hz
 * neighbour offsets (j-1, k-1). ey/ez/hy/hz share the same pattern.
 */
TEST(cuda, yee_advance_ex_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_media(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_coeffs(ctx, 2), 1);
    ASSERT_EQ(fdtd_cuda_alloc_metrics(ctx, d.nx, d.ny, d.nz, d.nx, d.ny, d.nz), 1);

    std::vector<fdtd_media> mEx(nelem(d), (fdtd_media)0);
    std::vector<fdtd_media> mZero(nelem(d), (fdtd_media)0);
    ASSERT_EQ(fdtd_cuda_upload_media(ctx, mEx.data(), mZero.data(), mZero.data(),
                                     mZero.data(), mZero.data(), mZero.data()),
              1);

    std::vector<fdtd_real> g1 = {(fdtd_real)0.9, (fdtd_real)0.0};
    std::vector<fdtd_real> g2 = {(fdtd_real)0.1, (fdtd_real)0.0};
    std::vector<fdtd_real> gm1 = {(fdtd_real)1.0, (fdtd_real)0.0};
    std::vector<fdtd_real> gm2 = {(fdtd_real)1.0, (fdtd_real)0.0};
    ASSERT_EQ(fdtd_cuda_upload_coeffs(ctx, g1.data(), g2.data(), gm1.data(), gm2.data(), 2), 1);

    std::vector<fdtd_real> Idxh(d.nx, (fdtd_real)1), Idyh(d.ny, (fdtd_real)2), Idzh(d.nz, (fdtd_real)3);
    std::vector<fdtd_real> Idxe(d.nx, (fdtd_real)1), Idye(d.ny, (fdtd_real)1), Idze(d.nz, (fdtd_real)1);
    ASSERT_EQ(fdtd_cuda_upload_metrics(ctx, Idxh.data(), Idyh.data(), Idzh.data(),
                                       Idxe.data(), Idye.data(), Idze.data()),
              1);

    SixFields f = make_six_zero(d);
    /* Non-trivial H neighbours so curl ≠ 0 inside the sweep. */
    for (int k = 0; k < d.nz; ++k)
        for (int j = 0; j < d.ny; ++j)
            for (int i = 0; i < d.nx; ++i) {
                f.Hy[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)(0.01 * (i + j + k));
                f.Hz[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)(0.02 * (i + 2 * j + 3 * k));
                f.Ex[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)(0.5);
            }

    SixFields host = f;
    fdtd_ibox sweep = {1, 6, 1, 6, 1, 6};
    host_advance_ex(host.Ex, host.Hy, host.Hz, mEx, g1, g2, Idyh, Idzh, d, sweep);

    ASSERT_EQ(upload_six(ctx, f), 1);
    ASSERT_EQ(fdtd_cuda_advance_ex(ctx, sweep), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);

    EXPECT_TRUE(near_eq(out.Ex, host.Ex, (fdtd_real)1e-5));

    fdtd_cuda_destroy(ctx);
}

/*
 * Documents fdtd_cuda_advance_hx: gm1/gm2, Idze/Idye, Ey(k+1)/Ez(j+1) offsets.
 */
TEST(cuda, yee_advance_hx_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_media(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_coeffs(ctx, 2), 1);
    ASSERT_EQ(fdtd_cuda_alloc_metrics(ctx, d.nx, d.ny, d.nz, d.nx, d.ny, d.nz), 1);

    std::vector<fdtd_media> mHx(nelem(d), (fdtd_media)0);
    std::vector<fdtd_media> mZero(nelem(d), (fdtd_media)0);
    ASSERT_EQ(fdtd_cuda_upload_media(ctx, mZero.data(), mZero.data(), mZero.data(),
                                     mHx.data(), mZero.data(), mZero.data()),
              1);

    std::vector<fdtd_real> g1 = {(fdtd_real)1.0, (fdtd_real)0.0};
    std::vector<fdtd_real> g2 = {(fdtd_real)1.0, (fdtd_real)0.0};
    std::vector<fdtd_real> gm1 = {(fdtd_real)0.8, (fdtd_real)0.0};
    std::vector<fdtd_real> gm2 = {(fdtd_real)0.2, (fdtd_real)0.0};
    ASSERT_EQ(fdtd_cuda_upload_coeffs(ctx, g1.data(), g2.data(), gm1.data(), gm2.data(), 2), 1);

    std::vector<fdtd_real> Idxh(d.nx, (fdtd_real)1), Idyh(d.ny, (fdtd_real)1), Idzh(d.nz, (fdtd_real)1);
    std::vector<fdtd_real> Idxe(d.nx, (fdtd_real)1), Idye(d.ny, (fdtd_real)4), Idze(d.nz, (fdtd_real)5);
    ASSERT_EQ(fdtd_cuda_upload_metrics(ctx, Idxh.data(), Idyh.data(), Idzh.data(),
                                       Idxe.data(), Idye.data(), Idze.data()),
              1);

    SixFields f = make_six_zero(d);
    for (int k = 0; k < d.nz; ++k)
        for (int j = 0; j < d.ny; ++j)
            for (int i = 0; i < d.nx; ++i) {
                f.Ey[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)(0.03 * (i + j + k));
                f.Ez[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)(0.04 * (i + 2 * j));
                f.Hx[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)(0.25);
            }

    SixFields host = f;
    fdtd_ibox sweep = {1, 5, 1, 5, 1, 5}; /* leave k+1 / j+1 in-bounds */
    host_advance_hx(host.Hx, host.Ey, host.Ez, mHx, gm1, gm2, Idze, Idye, d, sweep);

    ASSERT_EQ(upload_six(ctx, f), 1);
    ASSERT_EQ(fdtd_cuda_advance_hx(ctx, sweep), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);

    EXPECT_TRUE(near_eq(out.Hx, host.Hx, (fdtd_real)1e-5));

    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_YEE_TESTS_H */
