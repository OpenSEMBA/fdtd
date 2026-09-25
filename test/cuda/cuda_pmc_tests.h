#ifndef CUDA_PMC_TESTS_H
#define CUDA_PMC_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/*
 * Host MinusCloneMagneticPMC writes ghost = -adjacent interior.
 * Host CloneMagneticPeriodic writes ghost = opposite interior.
 * field_comp 3=Hx. wall_axis 2 is z. Origins here are 0.
 */

static fdtd_clone_job hx_z_job(int ghost, int source, int sign)
{
    fdtd_clone_job job{};
    job.field_comp = 3;
    job.wall_axis = 2;
    job.ghost = ghost;
    job.source = source;
    job.sign = sign;
    job.a0 = 0;
    job.a1 = kTiny.nx - 1;
    job.b0 = 0;
    job.b1 = kTiny.ny - 1;
    return job;
}

TEST(cuda, pmc_negates_adjacent_plane)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    SixFields f = make_six_zero(d);
    for (int j = 0; j < d.ny; ++j)
        for (int i = 0; i < d.nx; ++i) {
            f.Hx[idx3(i, j, 1, d.nx, d.ny)] = (fdtd_real)2;
            f.Hx[idx3(i, j, 0, d.nx, d.ny)] = (fdtd_real)9;
            f.Hx[idx3(i, j, 3, d.nx, d.ny)] = (fdtd_real)4;
        }
    ASSERT_EQ(upload_six(ctx, f), 1);
    fdtd_clone_job job = hx_z_job(0, 1, -1);
    ASSERT_EQ(fdtd_cuda_set_clone_jobs(ctx, &job, 1), 1);
    ASSERT_EQ(fdtd_cuda_advance_clones(ctx), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Hx[idx3(2, 3, 0, d.nx, d.ny)], (fdtd_real)-2, 1e-5);
    EXPECT_NEAR(out.Hx[idx3(2, 3, 1, d.nx, d.ny)], (fdtd_real)2, 1e-5);
    EXPECT_NEAR(out.Hx[idx3(2, 3, 3, d.nx, d.ny)], (fdtd_real)4, 1e-5);
    EXPECT_NEAR(out.Hy[idx3(2, 3, 0, d.nx, d.ny)], (fdtd_real)0, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, periodic_copies_opposite_plane)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    SixFields f = make_six_zero(d);
    for (int j = 0; j < d.ny; ++j)
        for (int i = 0; i < d.nx; ++i) {
            f.Hx[idx3(i, j, 0, d.nx, d.ny)] = (fdtd_real)1;
            f.Hx[idx3(i, j, 7, d.nx, d.ny)] = (fdtd_real)4;
            f.Hx[idx3(i, j, 2, d.nx, d.ny)] = (fdtd_real)-6;
        }
    ASSERT_EQ(upload_six(ctx, f), 1);
    fdtd_clone_job job = hx_z_job(0, 7, +1);
    ASSERT_EQ(fdtd_cuda_set_clone_jobs(ctx, &job, 1), 1);
    ASSERT_EQ(fdtd_cuda_advance_clones(ctx), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Hx[idx3(1, 4, 0, d.nx, d.ny)], (fdtd_real)4, 1e-5);
    EXPECT_NEAR(out.Hx[idx3(1, 4, 7, d.nx, d.ny)], (fdtd_real)4, 1e-5);
    EXPECT_NEAR(out.Hx[idx3(1, 4, 2, d.nx, d.ny)], (fdtd_real)-6, 1e-5);
    fdtd_cuda_destroy(ctx);
}

#endif
