#ifndef CUDA_MUR_TESTS_H
#define CUDA_MUR_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/*
 * Host first-order Mur (AdvanceMagneticMUR, left Hx):
 *   Hx(i, j_g, k) = Past(i, j_g+1, k) + CAB1(medio) * (Hx(i, j_g+1, k) - Past(i, j_g, k))
 * then Past slab is overwritten from current Hx.
 */
TEST(cuda, mur_apply_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_media(ctx, d, d, d, d, d, d), 1);

    std::vector<fdtd_media> mZero(nelem(d), (fdtd_media)0);
    ASSERT_EQ(fdtd_cuda_upload_media(ctx, mZero.data(), mZero.data(), mZero.data(),
                                     mZero.data(), mZero.data(), mZero.data()),
              1);

    const int n_cab = 2;
    std::vector<fdtd_real> cab(n_cab, (fdtd_real)0);
    cab[0] = (fdtd_real)0.5;
    ASSERT_EQ(fdtd_cuda_upload_mur_cab(ctx, 0, n_cab, cab.data()), 1);

    /* Ghost plane j=1, interior neighbour j=2. Past slab j in [1,2]. */
    const int jg = 1, jn = 2;
    const int i0 = 1, i1 = 4, k0 = 1, k1 = 4;
    const int p_xi = 1, p_yi = 1, p_zi = 1;
    const int nx_p = i1 - p_xi + 1;
    const int ny_p = jn - p_yi + 1;
    const int nz_p = k1 - p_zi + 1;
    const int n_past = nx_p * ny_p * nz_p;

    fdtd_mur_job job{};
    job.field_comp = 3;
    job.media_comp = 3;
    job.cab_which = 0;
    job.wall_axis = 1;
    job.neigh_sign = 1;
    job.plane_abs = jg;
    job.a0 = i0;
    job.a1 = i1;
    job.b0 = k0;
    job.b1 = k1;
    job.e_xi = 0;
    job.e_yi = 0;
    job.e_zi = 0;
    job.p_xi = p_xi;
    job.p_yi = p_yi;
    job.p_zi = p_zi;
    job.nx_p = nx_p;
    job.ny_p = ny_p;
    job.nz_p = nz_p;
    job.store_xi = i0;
    job.store_xe = i1;
    job.store_yi = jg;
    job.store_ye = jn;
    job.store_zi = k0;
    job.store_ze = k1;

    ASSERT_EQ(fdtd_cuda_set_mur_jobs(ctx, &job, 1), 1);
    EXPECT_EQ(fdtd_cuda_mur_ready(ctx), 1);

    std::vector<fdtd_real> past(n_past, (fdtd_real)0);
    const fdtd_real past_g = (fdtd_real)0.1;
    const fdtd_real past_n = (fdtd_real)0.3;
    const fdtd_real hx_n = (fdtd_real)2.0;
    for (int k = k0; k <= k1; ++k)
        for (int i = i0; i <= i1; ++i) {
            past[idx3(i - p_xi, jg - p_yi, k - p_zi, nx_p, ny_p)] = past_g;
            past[idx3(i - p_xi, jn - p_yi, k - p_zi, nx_p, ny_p)] = past_n;
        }
    ASSERT_EQ(fdtd_cuda_upload_mur_past(ctx, 0, past.data(), n_past), 1);

    SixFields f = make_six_zero(d);
    for (int k = k0; k <= k1; ++k)
        for (int i = i0; i <= i1; ++i) {
            f.Hx[idx3(i, jg, k, d.nx, d.ny)] = (fdtd_real)9.0;
            f.Hx[idx3(i, jn, k, d.nx, d.ny)] = hx_n;
        }
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_real host_ghost = past_n + cab[0] * (hx_n - past_g);

    ASSERT_EQ(fdtd_cuda_advance_mur(ctx), 1);

    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    std::vector<fdtd_real> past_out(n_past, (fdtd_real)0);
    ASSERT_EQ(fdtd_cuda_download_mur_past(ctx, 0, past_out.data(), n_past), 1);

    for (int k = k0; k <= k1; ++k)
        for (int i = i0; i <= i1; ++i) {
            EXPECT_NEAR(out.Hx[idx3(i, jg, k, d.nx, d.ny)], host_ghost, 1e-5);
            EXPECT_NEAR(out.Hx[idx3(i, jn, k, d.nx, d.ny)], hx_n, 1e-5);
            EXPECT_NEAR(past_out[idx3(i - p_xi, jg - p_yi, k - p_zi, nx_p, ny_p)], host_ghost,
                        1e-5);
            EXPECT_NEAR(past_out[idx3(i - p_xi, jn - p_yi, k - p_zi, nx_p, ny_p)], hx_n, 1e-5);
        }

    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_MUR_TESTS_H */
