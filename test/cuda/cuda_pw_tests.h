#ifndef CUDA_PW_TESTS_H
#define CUDA_PW_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

#include <vector>

using namespace cuda_test;

namespace {

/* Host mirror of device evolucion + Incid for one mode, one wave. */
fdtd_real host_evolucion(fdtd_real t, fdtd_real d, fdtd_real cluz, fdtd_real delta,
                         int numus, const fdtd_real *evol)
{
    long long nprev = (long long)((t - d / cluz) / delta);
    if (nprev + 1 <= (long long)numus && nprev > 0) {
        fdtd_real e0 = evol[nprev];
        fdtd_real e1 = evol[nprev + 1];
        return (e1 - e0) / delta * ((t - d / cluz) - (fdtd_real)nprev * delta) + e0;
    }
    return (fdtd_real)0;
}

fdtd_real host_incid(fdtd_real time, int i, int j, int k, fdtd_real cluz, fdtd_real px, fdtd_real py,
                     fdtd_real pz, fdtd_real d0, fdtd_real fpw, fdtd_real delta, int numus,
                     const fdtd_real *evol, const fdtd_real *phys_x, const fdtd_real *phys_y,
                     const fdtd_real *phys_z, int x0, int y0, int z0)
{
    fdtd_real xf = phys_x[i - x0];
    fdtd_real yf = phys_y[j - y0];
    fdtd_real zf = phys_z[k - z0];
    fdtd_real dd = xf * px + yf * py + zf * pz - d0;
    return fpw * host_evolucion(time, dd, cluz, delta, numus, evol);
}

} // namespace

/* Tiny E-face: Ez back plane, one cell; device matches host G2*Incid*Id. */
TEST(cuda, planewave_e_face_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_coeffs(ctx, 4), 1);
    ASSERT_EQ(fdtd_cuda_alloc_metrics(ctx, 8, 8, 8, 8, 8, 8), 1);

    std::vector<fdtd_real> g1(5, 0), g2(5, 0), gm1(5, 0), gm2(5, 0);
    g2[1] = (fdtd_real)0.25;
    gm2[1] = (fdtd_real)0.5;
    ASSERT_EQ(fdtd_cuda_upload_coeffs(ctx, g1.data(), g2.data(), gm1.data(), gm2.data(), 4), 1);

    std::vector<fdtd_real> Idxh(8, (fdtd_real)2.0), Idyh(8, 1), Idzh(8, 1);
    std::vector<fdtd_real> Idxe(8, 1), Idye(8, 1), Idze(8, 1);
    ASSERT_EQ(fdtd_cuda_upload_metrics(ctx, Idxh.data(), Idyh.data(), Idzh.data(), Idxe.data(),
                                       Idye.data(), Idze.data()),
              1);

    /* PhysCoor for iHy (field index 5 → CUDA 4): uniform coords. */
    std::vector<fdtd_real> px(16), py(16), pz(16);
    for (int t = 0; t < 16; ++t) {
        px[t] = (fdtd_real)t;
        py[t] = (fdtd_real)t;
        pz[t] = (fdtd_real)t;
    }
    const int base = 0;
    ASSERT_EQ(fdtd_cuda_upload_planewave_phys(ctx, 4 /*Hy*/, 0, base, 16, px.data()), 1);
    ASSERT_EQ(fdtd_cuda_upload_planewave_phys(ctx, 4, 1, base, 16, py.data()), 1);
    ASSERT_EQ(fdtd_cuda_upload_planewave_phys(ctx, 4, 2, base, 16, pz.data()), 1);

    const int numus = 8;
    const fdtd_real delta = (fdtd_real)1.0;
    const fdtd_real cluz = (fdtd_real)1.0;
    std::vector<fdtd_real> evol(numus + 1);
    for (int s = 0; s <= numus; ++s) evol[s] = (fdtd_real)(0.1 * s);

    int num_modes[1] = {1};
    int numus_a[1] = {numus};
    fdtd_real delta_a[1] = {delta};
    fdtd_real p_x[1] = {(fdtd_real)1.0}, p_y[1] = {0}, p_z[1] = {0}, d0[1] = {0};
    /* fpw[wave*6*modes + (nfield-1)*modes + mode]; nfield=iHy=5 → index 4 */
    std::vector<fdtd_real> fpw(6, (fdtd_real)0);
    fpw[4] = (fdtd_real)1.0;

    fdtd_pw_face face{};
    face.field_comp = 2; /* Ez */
    face.incid_nfield = 5; /* iHy */
    face.wave = 0;
    face.free_mode = 0;
    face.fixed_abs = 3;
    face.a0 = 2;
    face.a1 = 2;
    face.b0 = 4;
    face.b1 = 4;
    face.incid_di = -1;
    face.incid_dj = 0;
    face.incid_dk = 0;
    face.field_xi = 0;
    face.field_yi = 0;
    face.field_zi = 0;
    face.id_axis = 0;
    face.use_e_metric = 0;
    face.sign = -1;

    fdtd_pw_face dummy{};
    ASSERT_EQ(fdtd_cuda_upload_planewave(ctx, 1, 1, numus, cluz, num_modes, numus_a, delta_a,
                                         evol.data(), p_x, p_y, p_z, d0, fpw.data(), &face, 1,
                                         &dummy, 0),
              1);

    SixFields src = make_six_zero(d);
    src.Ez[idx3(3, 2, 4, d.nx, d.ny)] = (fdtd_real)10.0;
    ASSERT_EQ(upload_six(ctx, src), 1);

    const fdtd_real time = (fdtd_real)3.5;
    int still = 0;
    ASSERT_EQ(fdtd_cuda_advance_planewave_e(ctx, time, &still), 1);

    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);

    fdtd_real inc = host_incid(time, 3 - 1, 2, 4, cluz, p_x[0], p_y[0], p_z[0], d0[0], fpw[4], delta,
                               numus, evol.data(), px.data(), py.data(), pz.data(), base, base, base);
    fdtd_real expect = (fdtd_real)10.0 + (fdtd_real)(-1) * g2[1] * inc * Idxh[3];
    EXPECT_NEAR(out.Ez[idx3(3, 2, 4, d.nx, d.ny)], expect, (fdtd_real)1e-5);
    EXPECT_EQ(still, 1);

    fdtd_cuda_destroy(ctx);
}

/* Tiny H-face: Hz back, one cell vs host Gm2*Incid*Id. */
TEST(cuda, planewave_h_face_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_coeffs(ctx, 4), 1);
    ASSERT_EQ(fdtd_cuda_alloc_metrics(ctx, 8, 8, 8, 8, 8, 8), 1);

    std::vector<fdtd_real> g1(5, 0), g2(5, 0), gm1(5, 0), gm2(5, 0);
    g2[1] = (fdtd_real)0.25;
    gm2[1] = (fdtd_real)0.5;
    ASSERT_EQ(fdtd_cuda_upload_coeffs(ctx, g1.data(), g2.data(), gm1.data(), gm2.data(), 4), 1);

    std::vector<fdtd_real> Idxh(8, 1), Idyh(8, 1), Idzh(8, 1);
    std::vector<fdtd_real> Idxe(8, (fdtd_real)3.0), Idye(8, 1), Idze(8, 1);
    ASSERT_EQ(fdtd_cuda_upload_metrics(ctx, Idxh.data(), Idyh.data(), Idzh.data(), Idxe.data(),
                                       Idye.data(), Idze.data()),
              1);

    std::vector<fdtd_real> px(16), py(16), pz(16);
    for (int t = 0; t < 16; ++t) {
        px[t] = (fdtd_real)t;
        py[t] = (fdtd_real)t;
        pz[t] = (fdtd_real)t;
    }
    ASSERT_EQ(fdtd_cuda_upload_planewave_phys(ctx, 1 /*Ey*/, 0, 0, 16, px.data()), 1);
    ASSERT_EQ(fdtd_cuda_upload_planewave_phys(ctx, 1, 1, 0, 16, py.data()), 1);
    ASSERT_EQ(fdtd_cuda_upload_planewave_phys(ctx, 1, 2, 0, 16, pz.data()), 1);

    const int numus = 8;
    const fdtd_real delta = (fdtd_real)1.0;
    const fdtd_real cluz = (fdtd_real)1.0;
    std::vector<fdtd_real> evol(numus + 1);
    for (int s = 0; s <= numus; ++s) evol[s] = (fdtd_real)(0.2 * s);

    int num_modes[1] = {1};
    int numus_a[1] = {numus};
    fdtd_real delta_a[1] = {delta};
    fdtd_real p_x[1] = {0}, p_y[1] = {(fdtd_real)1.0}, p_z[1] = {0}, d0[1] = {0};
    std::vector<fdtd_real> fpw(6, (fdtd_real)0);
    fpw[1] = (fdtd_real)1.0; /* iEy */

    fdtd_pw_face face{};
    face.field_comp = 5; /* Hz */
    face.incid_nfield = 2; /* iEy */
    face.wave = 0;
    face.free_mode = 0;
    face.fixed_abs = 1;
    face.a0 = 3;
    face.a1 = 3;
    face.b0 = 2;
    face.b1 = 2;
    face.incid_di = +1;
    face.incid_dj = 0;
    face.incid_dk = 0;
    face.field_xi = 0;
    face.field_yi = 0;
    face.field_zi = 0;
    face.id_axis = 0;
    face.use_e_metric = 1;
    face.sign = +1;

    fdtd_pw_face dummy{};
    ASSERT_EQ(fdtd_cuda_upload_planewave(ctx, 1, 1, numus, cluz, num_modes, numus_a, delta_a,
                                         evol.data(), p_x, p_y, p_z, d0, fpw.data(), &dummy, 0, &face,
                                         1),
              1);

    SixFields src = make_six_zero(d);
    src.Hz[idx3(1, 3, 2, d.nx, d.ny)] = (fdtd_real)5.0;
    ASSERT_EQ(upload_six(ctx, src), 1);

    const fdtd_real time = (fdtd_real)4.0;
    int still = 0;
    ASSERT_EQ(fdtd_cuda_advance_planewave_h(ctx, time, &still), 1);

    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);

    fdtd_real inc = host_incid(time, 1 + 1, 3, 2, cluz, p_x[0], p_y[0], p_z[0], d0[0], fpw[1], delta,
                               numus, evol.data(), px.data(), py.data(), pz.data(), 0, 0, 0);
    fdtd_real expect = (fdtd_real)5.0 + (fdtd_real)(+1) * gm2[1] * inc * Idxe[1];
    EXPECT_NEAR(out.Hz[idx3(1, 3, 2, d.nx, d.ny)], expect, (fdtd_real)1e-5);

    fdtd_cuda_destroy(ctx);
}

/* No faces / empty upload → advance is a successful no-op. */
TEST(cuda, planewave_disabled_noop)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    ASSERT_EQ(fdtd_cuda_alloc_coeffs(ctx, 4), 1);
    ASSERT_EQ(fdtd_cuda_alloc_metrics(ctx, 8, 8, 8, 8, 8, 8), 1);

    std::vector<fdtd_real> g1(5, 0), g2(5, (fdtd_real)1), gm1(5, 0), gm2(5, (fdtd_real)1);
    ASSERT_EQ(fdtd_cuda_upload_coeffs(ctx, g1.data(), g2.data(), gm1.data(), gm2.data(), 4), 1);
    std::vector<fdtd_real> ones(8, 1);
    ASSERT_EQ(fdtd_cuda_upload_metrics(ctx, ones.data(), ones.data(), ones.data(), ones.data(),
                                       ones.data(), ones.data()),
              1);

    ASSERT_EQ(fdtd_cuda_upload_planewave(ctx, 0, 0, 0, 1.0f, nullptr, nullptr, nullptr, nullptr,
                                         nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0,
                                         nullptr, 0),
              1);
    EXPECT_EQ(fdtd_cuda_planewave_ready(ctx), 1);

    SixFields src = make_six(d, (fdtd_real)1.0);
    ASSERT_EQ(upload_six(ctx, src), 1);
    int still = 99;
    ASSERT_EQ(fdtd_cuda_advance_planewave_e(ctx, (fdtd_real)1.0, &still), 1);
    ASSERT_EQ(fdtd_cuda_advance_planewave_h(ctx, (fdtd_real)1.0, &still), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_TRUE(exact_eq(src.Ex, out.Ex));
    EXPECT_TRUE(exact_eq(src.Hy, out.Hy));

    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_PW_TESTS_H */
