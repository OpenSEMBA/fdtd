#ifndef CUDA_CPML_TESTS_H
#define CUDA_CPML_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/*
 * Documents fdtd_cpml_job + fdtd_cuda_cpml_apply for one absolute cell:
 *   Psi = P_b * Psi + (Ha[i,j,k] - Ha[i,j-1,k]) * P_c   (h_diff_axis=1)
 *   Ex += h_sign * G2[medio] * Psi
 * with use_fixed_medio so media lookup is not the focus.
 *
 * Job fields exercised: field_comp, psi_slot, h_comp_a, h_diff_axis,
 * free_axis, p_base, p_b_which/p_c_which, e_xi/psi_xi origins, h_sign.
 */
TEST(cuda, cpml_apply_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_media(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_coeffs(ctx, 2), 1);

    std::vector<fdtd_media> mZero(nelem(d), (fdtd_media)0);
    ASSERT_EQ(fdtd_cuda_upload_media(ctx, mZero.data(), mZero.data(), mZero.data(),
                                     mZero.data(), mZero.data(), mZero.data()),
              1);

    std::vector<fdtd_real> g1 = {(fdtd_real)1.0, (fdtd_real)0.0};
    std::vector<fdtd_real> g2 = {(fdtd_real)0.5, (fdtd_real)0.0};
    std::vector<fdtd_real> gm1 = {(fdtd_real)1.0, (fdtd_real)0.0};
    std::vector<fdtd_real> gm2 = {(fdtd_real)1.0, (fdtd_real)0.0};
    ASSERT_EQ(fdtd_cuda_upload_coeffs(ctx, g1.data(), g2.data(), gm1.data(), gm2.data(), 2), 1);

    /* Absolute region is a single cell at (i,j,k)=(3,4,2); field origin 0. */
    const int ai = 3, aj = 4, ak = 2;
    const int psi_n = 4;
    ASSERT_EQ(fdtd_cuda_alloc_psi(ctx, 0, psi_n * psi_n * psi_n), 1);

    std::vector<fdtd_real> psi(psi_n * psi_n * psi_n, (fdtd_real)0);
    /* psi absolute origin = (2,3,1) so remapped (1,1,1) for cell (3,4,2) */
    const int psi_xi = 2, psi_yi = 3, psi_zi = 1;
    psi[idx3(ai - psi_xi, aj - psi_yi, ak - psi_zi, psi_n, psi_n)] = (fdtd_real)0.1;
    ASSERT_EQ(fdtd_cuda_upload_psi(ctx, 0, psi.data(), (int)psi.size()), 1);

    /* P_be_y = which 1, P_ce_y = which 7 (see fdtd_cuda.h comment). */
    const int n1d = 16;
    std::vector<fdtd_real> Pb(n1d, (fdtd_real)0), Pc(n1d, (fdtd_real)0);
    const int p_base = 0;
    Pb[aj - p_base] = (fdtd_real)0.7;
    Pc[aj - p_base] = (fdtd_real)0.3;
    ASSERT_EQ(fdtd_cuda_alloc_cpml_1d(ctx, 1, n1d, Pb.data()), 1);
    ASSERT_EQ(fdtd_cuda_alloc_cpml_1d(ctx, 7, n1d, Pc.data()), 1);

    SixFields f = make_six_zero(d);
    f.Ex[idx3(ai, aj, ak, d.nx, d.ny)] = (fdtd_real)1.0;
    f.Hy[idx3(ai, aj, ak, d.nx, d.ny)] = (fdtd_real)2.0;     /* Ha at cell */
    f.Hy[idx3(ai, aj - 1, ak, d.nx, d.ny)] = (fdtd_real)0.5; /* neighbour j-1 */
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_cpml_job job{};
    job.field_comp = 0; /* Ex */
    job.psi_slot = 0;
    job.media_comp = 0;
    job.h_comp_a = 4; /* Hy */
    job.h_comp_b = 4;
    job.h_diff_axis = 1; /* y */
    job.p_b_which = 1;
    job.p_c_which = 7;
    job.p_base = p_base;
    job.free_axis = 1; /* y feeds P_b/P_c */
    job.xi = ai;
    job.xe = ai;
    job.yi = aj;
    job.ye = aj;
    job.zi = ak;
    job.ze = ak;
    job.e_xi = 0;
    job.e_yi = 0;
    job.e_zi = 0;
    job.psi_xi = psi_xi;
    job.psi_yi = psi_yi;
    job.psi_zi = psi_zi;
    job.nx_psi = psi_n;
    job.ny_psi = psi_n;
    job.nz_psi = psi_n;
    job.h_sign = 1;
    job.use_fixed_medio = 1;
    job.medio_fixed = 0;

    fdtd_real host_Ex = f.Ex[idx3(ai, aj, ak, d.nx, d.ny)];
    fdtd_real host_Psi = psi[idx3(ai - psi_xi, aj - psi_yi, ak - psi_zi, psi_n, psi_n)];
    host_cpml_cell(host_Ex, host_Psi, f.Hy[idx3(ai, aj, ak, d.nx, d.ny)],
                   f.Hy[idx3(ai, aj - 1, ak, d.nx, d.ny)], Pb[aj - p_base], Pc[aj - p_base],
                   g2[0], job.h_sign);

    ASSERT_EQ(fdtd_cuda_cpml_apply(ctx, &job), 1);

    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    std::vector<fdtd_real> psi_out(psi.size(), (fdtd_real)0);
    ASSERT_EQ(fdtd_cuda_download_psi(ctx, 0, psi_out.data(), (int)psi_out.size()), 1);

    EXPECT_NEAR(out.Ex[idx3(ai, aj, ak, d.nx, d.ny)], host_Ex, 1e-5);
    EXPECT_NEAR(psi_out[idx3(ai - psi_xi, aj - psi_yi, ak - psi_zi, psi_n, psi_n)], host_Psi,
                1e-5);

    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_CPML_TESTS_H */
