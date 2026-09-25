#ifndef CUDA_NODAL_TESTS_H
#define CUDA_NODAL_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/*
 * Host AdvanceNodalE soft Ez:
 *   Ez -= G2(medio) * Idxh(i) * Idyh(j) * amplitude * evolucion(t)
 * evolucion is linear between evol(nprev) and evol(nprev+1), 0 outside [0, numus].
 * PEC media (skip_e != 0) are left unchanged. Hard sources assign amplitude*wave.
 * Initial-value sources apply evol(0) only at step 0.
 */

static int upload_nodal_base(fdtd_cuda_ctx *ctx, fdtd_dims3 d, const std::vector<fdtd_media> &mez,
                             const std::vector<fdtd_real> &g2, const std::vector<fdtd_real> &idxh,
                             const std::vector<fdtd_real> &idyh, const std::vector<int> &skip_e)
{
    std::vector<fdtd_media> m0(nelem(d), (fdtd_media)1);
    std::vector<fdtd_real> z(8, (fdtd_real)1);
    std::vector<fdtd_real> gm(2, (fdtd_real)1);
    std::vector<int> skip_h(skip_e.size(), 0);
    if (fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d) != 1) return 0;
    if (fdtd_cuda_alloc_media(ctx, d, d, d, d, d, d) != 1) return 0;
    if (fdtd_cuda_alloc_coeffs(ctx, 1) != 1) return 0;
    if (fdtd_cuda_alloc_metrics(ctx, 8, 8, 8, 8, 8, 8) != 1) return 0;
    if (fdtd_cuda_upload_media(ctx, m0.data(), m0.data(), mez.data(), m0.data(), m0.data(),
                               m0.data()) != 1)
        return 0;
    std::vector<fdtd_real> g1(2, (fdtd_real)1);
    if (fdtd_cuda_upload_coeffs(ctx, g1.data(), g2.data(), gm.data(), gm.data(), 1) != 1) return 0;
    if (fdtd_cuda_upload_metrics(ctx, idxh.data(), idyh.data(), z.data(), z.data(), z.data(),
                                 z.data()) != 1)
        return 0;
    (void)skip_h;
    return 1;
}

TEST(cuda, nodal_soft_interp_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    std::vector<fdtd_media> mez(nelem(d), (fdtd_media)1);
    std::vector<fdtd_real> g2(2, (fdtd_real)0);
    g2[1] = (fdtd_real)0.25;
    std::vector<fdtd_real> idxh(8, (fdtd_real)0), idyh(8, (fdtd_real)0);
    idxh[2] = (fdtd_real)2;
    idyh[3] = (fdtd_real)3;
    std::vector<int> skip_e = {1, 0};
    ASSERT_EQ(upload_nodal_base(ctx, d, mez, g2, idxh, idyh, skip_e), 1);

    const int i = 2, j = 3, k = 4;
    SixFields f = make_six_zero(d);
    f.Ez[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)1;
    f.Ez[idx3(i + 1, j, k, d.nx, d.ny)] = (fdtd_real)8;
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_real samples[3] = {(fdtd_real)0, (fdtd_real)1, (fdtd_real)3};
    fdtd_nodal_job job{};
    job.field_comp = 2;
    job.hard = 0;
    job.initial_only = 0;
    job.xi = job.xe = i;
    job.yi = job.ye = j;
    job.zi = job.ze = k;
    job.e_xi = job.e_yi = job.e_zi = 0;
    job.evol_off = 0;
    job.numus = 2;
    job.amplitude = (fdtd_real)4;
    job.deltaevol = (fdtd_real)1;
    std::vector<int> skip_h = {0, 0};
    ASSERT_EQ(fdtd_cuda_upload_nodal(ctx, &job, 1, samples, 3, skip_e.data(), skip_h.data(), 2), 1);
    EXPECT_EQ(fdtd_cuda_nodal_ready(ctx), 1);

    const fdtd_real time = (fdtd_real)1.5;
    const fdtd_real wave = (samples[2] - samples[1]) / job.deltaevol * (time - (fdtd_real)1) + samples[1];
    const fdtd_real expect = (fdtd_real)1 - g2[1] * idxh[i] * idyh[j] * job.amplitude * wave;

    ASSERT_EQ(fdtd_cuda_advance_nodal_e(ctx, time, 1), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(i, j, k, d.nx, d.ny)], expect, 1e-5);
    EXPECT_NEAR(out.Ez[idx3(i + 1, j, k, d.nx, d.ny)], (fdtd_real)8, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, nodal_hard_assigns_and_skips_pec)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    std::vector<fdtd_media> mez(nelem(d), (fdtd_media)1);
    const int j = 3, k = 4;
    const int i_pec = 2, i_vac = 3;
    mez[idx3(i_pec, j, k, d.nx, d.ny)] = (fdtd_media)0;
    std::vector<fdtd_real> g2(2, (fdtd_real)1);
    std::vector<fdtd_real> ones(8, (fdtd_real)1);
    std::vector<int> skip_e = {1, 0};
    ASSERT_EQ(upload_nodal_base(ctx, d, mez, g2, ones, ones, skip_e), 1);

    SixFields f = make_six_zero(d);
    f.Ez[idx3(i_pec, j, k, d.nx, d.ny)] = (fdtd_real)4;
    f.Ez[idx3(i_vac, j, k, d.nx, d.ny)] = (fdtd_real)9;
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_real samples[2] = {(fdtd_real)7, (fdtd_real)9};
    fdtd_nodal_job job{};
    job.field_comp = 2;
    job.hard = 1;
    job.xi = i_pec;
    job.xe = i_vac;
    job.yi = job.ye = j;
    job.zi = job.ze = k;
    job.numus = 1;
    job.amplitude = (fdtd_real)5;
    job.deltaevol = (fdtd_real)1;
    std::vector<int> skip_h = {0, 0};
    ASSERT_EQ(fdtd_cuda_upload_nodal(ctx, &job, 1, samples, 2, skip_e.data(), skip_h.data(), 2), 1);

    const fdtd_real time = (fdtd_real)0;
    const fdtd_real wave = samples[0];
    ASSERT_EQ(fdtd_cuda_advance_nodal_e(ctx, time, 0), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(i_pec, j, k, d.nx, d.ny)], (fdtd_real)4, 1e-5);
    EXPECT_NEAR(out.Ez[idx3(i_vac, j, k, d.nx, d.ny)], job.amplitude * wave, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, nodal_initial_value_only_step0)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    std::vector<fdtd_media> mez(nelem(d), (fdtd_media)1);
    std::vector<fdtd_real> g2(2, (fdtd_real)1);
    std::vector<fdtd_real> ones(8, (fdtd_real)1);
    std::vector<int> skip_e = {1, 0};
    ASSERT_EQ(upload_nodal_base(ctx, d, mez, g2, ones, ones, skip_e), 1);

    const int i = 2, j = 2, k = 2;
    SixFields f = make_six_zero(d);
    f.Ez[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)1;
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_real samples[1] = {(fdtd_real)2.5};
    fdtd_nodal_job job{};
    job.field_comp = 2;
    job.hard = 1;
    job.initial_only = 1;
    job.xi = job.xe = i;
    job.yi = job.ye = j;
    job.zi = job.ze = k;
    job.numus = 0;
    job.amplitude = (fdtd_real)4;
    job.deltaevol = (fdtd_real)1;
    std::vector<int> skip_h = {0, 0};
    ASSERT_EQ(fdtd_cuda_upload_nodal(ctx, &job, 1, samples, 1, skip_e.data(), skip_h.data(), 2), 1);

    ASSERT_EQ(fdtd_cuda_advance_nodal_e(ctx, (fdtd_real)0, 0), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(i, j, k, d.nx, d.ny)], job.amplitude * samples[0], 1e-5);

    ASSERT_EQ(fdtd_cuda_advance_nodal_e(ctx, (fdtd_real)1, 1), 1);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(i, j, k, d.nx, d.ny)], job.amplitude * samples[0], 1e-5);
    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_NODAL_TESTS_H */
