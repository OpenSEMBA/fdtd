#ifndef CUDA_PROBE_TESTS_H
#define CUDA_PROBE_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/* Invariant: gather of one registered cell returns the uploaded value. */
TEST(cuda, point_probe_gather_single)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields src = make_six_zero(d);
    src.Ex[idx3(3, 4, 2, d.nx, d.ny)] = (fdtd_real)42.5;
    ASSERT_EQ(upload_six(ctx, src), 1);

    int comp[1] = {0};
    int i[1] = {3}, j[1] = {4}, k[1] = {2};
    ASSERT_EQ(fdtd_cuda_set_point_probes(ctx, 1, comp, i, j, k), 1);
    fdtd_real out[1] = {0};
    ASSERT_EQ(fdtd_cuda_gather_point_probes(ctx, out), 1);
    EXPECT_EQ(out[0], (fdtd_real)42.5);

    fdtd_cuda_destroy(ctx);
}

/* Invariant: mixed components; out[] order matches registration order. */
TEST(cuda, point_probe_gather_multi_components)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields src = make_six_zero(d);
    src.Ex[idx3(1, 1, 1, d.nx, d.ny)] = (fdtd_real)1.0;
    src.Ey[idx3(2, 3, 4, d.nx, d.ny)] = (fdtd_real)2.0;
    src.Hz[idx3(5, 5, 5, d.nx, d.ny)] = (fdtd_real)3.0;
    ASSERT_EQ(upload_six(ctx, src), 1);

    int comp[3] = {0, 1, 5};
    int i[3] = {1, 2, 5}, j[3] = {1, 3, 5}, k[3] = {1, 4, 5};
    ASSERT_EQ(fdtd_cuda_set_point_probes(ctx, 3, comp, i, j, k), 1);
    fdtd_real out[3] = {0, 0, 0};
    ASSERT_EQ(fdtd_cuda_gather_point_probes(ctx, out), 1);
    EXPECT_EQ(out[0], (fdtd_real)1.0);
    EXPECT_EQ(out[1], (fdtd_real)2.0);
    EXPECT_EQ(out[2], (fdtd_real)3.0);

    fdtd_cuda_destroy(ctx);
}

/* Invariant: column-major layout for probe index. */
TEST(cuda, point_probe_gather_layout)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields src = make_six_zero(d);
    const int i0 = 2, j0 = 6, k0 = 3;
    src.Hx[idx3(i0, j0, k0, d.nx, d.ny)] = (fdtd_real)-7.25;
    ASSERT_EQ(upload_six(ctx, src), 1);

    int comp[1] = {3};
    int i[1] = {i0}, j[1] = {j0}, k[1] = {k0};
    ASSERT_EQ(fdtd_cuda_set_point_probes(ctx, 1, comp, i, j, k), 1);
    fdtd_real out[1] = {0};
    ASSERT_EQ(fdtd_cuda_gather_point_probes(ctx, out), 1);
    EXPECT_EQ(out[0], (fdtd_real)-7.25);

    fdtd_cuda_destroy(ctx);
}

/* Invariant: n==0 is a successful no-op. */
TEST(cuda, point_probe_gather_empty_noop)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_set_point_probes(ctx, 0, nullptr, nullptr, nullptr, nullptr), 1);
    ASSERT_EQ(fdtd_cuda_gather_point_probes(ctx, nullptr), 1);
    fdtd_cuda_destroy(ctx);
}

/* Invariant: out-of-range index rejected at set time. */
TEST(cuda, point_probe_gather_oob_rejected)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    int comp[1] = {0};
    int i[1] = {99}, j[1] = {0}, k[1] = {0};
    EXPECT_EQ(fdtd_cuda_set_point_probes(ctx, 1, comp, i, j, k), 0);

    int badc[1] = {9};
    i[0] = 0;
    EXPECT_EQ(fdtd_cuda_set_point_probes(ctx, 1, badc, i, j, k), 0);

    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_PROBE_TESTS_H */
