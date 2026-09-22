#ifndef CUDA_SYNC_TESTS_H
#define CUDA_SYNC_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/* Invariant: bulk H2D then D2H is an identity on all six fields. */
TEST(cuda, fields_bulk_roundtrip)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields src = make_six(d, (fdtd_real)1);
    ASSERT_EQ(upload_six(ctx, src), 1);

    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);

    EXPECT_TRUE(exact_eq(src.Ex, out.Ex));
    EXPECT_TRUE(exact_eq(src.Ey, out.Ey));
    EXPECT_TRUE(exact_eq(src.Ez, out.Ez));
    EXPECT_TRUE(exact_eq(src.Hx, out.Hx));
    EXPECT_TRUE(exact_eq(src.Hy, out.Hy));
    EXPECT_TRUE(exact_eq(src.Hz, out.Hz));

    fdtd_cuda_destroy(ctx);
}

/*
 * Invariant: boxed H2D of a pattern onto a zero device, then full D2H —
 * cells inside ibox match the pattern; outside stay 0.
 */
TEST(cuda, fields_box_roundtrip_preserves_outside)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields zero = make_six_zero(d);
    ASSERT_EQ(upload_six(ctx, zero), 1);

    SixFields pat = make_six(d, (fdtd_real)7);
    fdtd_ibox box = {2, 4, 2, 5, 1, 3};
    ASSERT_EQ(upload_six_box(ctx, pat, box), 1);

    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);

    for (int k = 0; k < d.nz; ++k)
        for (int j = 0; j < d.ny; ++j)
            for (int i = 0; i < d.nx; ++i) {
                size_t id = idx3(i, j, k, d.nx, d.ny);
                if (in_box(i, j, k, box)) {
                    EXPECT_EQ(out.Ex[id], pat.Ex[id]) << "i=" << i << " j=" << j << " k=" << k;
                } else {
                    EXPECT_EQ(out.Ex[id], (fdtd_real)0) << "i=" << i << " j=" << j << " k=" << k;
                }
            }

    fdtd_cuda_destroy(ctx);
}

/*
 * Invariant: boxed D2H only overwrites the host region inside ibox;
 * host cells outside keep their prefilled sentinel.
 */
TEST(cuda, fields_box_download_only_touches_box)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields device_pat = make_six(d, (fdtd_real)3);
    ASSERT_EQ(upload_six(ctx, device_pat), 1);

    SixFields host = make_six(d, (fdtd_real)9000); /* sentinel B */
    SixFields host_before = host;
    fdtd_ibox box = {1, 3, 2, 4, 2, 6};
    ASSERT_EQ(download_six_box(ctx, host, box), 1);

    for (int k = 0; k < d.nz; ++k)
        for (int j = 0; j < d.ny; ++j)
            for (int i = 0; i < d.nx; ++i) {
                size_t id = idx3(i, j, k, d.nx, d.ny);
                if (in_box(i, j, k, box)) {
                    EXPECT_EQ(host.Ex[id], device_pat.Ex[id]);
                } else {
                    EXPECT_EQ(host.Ex[id], host_before.Ex[id]);
                }
            }

    fdtd_cuda_destroy(ctx);
}

/* Invariant: full-domain ibox path matches bulk upload/download exactly. */
TEST(cuda, fields_full_ibox_matches_bulk)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields src = make_six(d, (fdtd_real)11);
    fdtd_ibox full = full_ibox(d);

    ASSERT_EQ(upload_six(ctx, src), 1);
    SixFields via_bulk = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, via_bulk), 1);

    SixFields zero = make_six_zero(d);
    ASSERT_EQ(upload_six(ctx, zero), 1);
    ASSERT_EQ(upload_six_box(ctx, src, full), 1);
    SixFields via_box = make_six_zero(d);
    ASSERT_EQ(download_six_box(ctx, via_box, full), 1);

    EXPECT_TRUE(exact_eq(via_bulk.Ex, via_box.Ex));
    EXPECT_TRUE(exact_eq(via_bulk.Ey, via_box.Ey));
    EXPECT_TRUE(exact_eq(via_bulk.Ez, via_box.Ez));
    EXPECT_TRUE(exact_eq(via_bulk.Hx, via_box.Hx));
    EXPECT_TRUE(exact_eq(via_bulk.Hy, via_box.Hy));
    EXPECT_TRUE(exact_eq(via_bulk.Hz, via_box.Hz));

    fdtd_cuda_destroy(ctx);
}

/* Invariant: inverted/empty ibox is a safe no-op (ok return, no clobber). */
TEST(cuda, fields_empty_ibox_noop)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields device_pat = make_six(d, (fdtd_real)5);
    ASSERT_EQ(upload_six(ctx, device_pat), 1);

    SixFields host = make_six(d, (fdtd_real)42);
    SixFields host_before = host;
    fdtd_ibox empty = empty_ibox();
    ASSERT_EQ(download_six_box(ctx, host, empty), 1);
    EXPECT_TRUE(exact_eq(host.Ex, host_before.Ex));
    EXPECT_TRUE(exact_eq(host.Hz, host_before.Hz));

    SixFields upload_try = make_six(d, (fdtd_real)77);
    ASSERT_EQ(upload_six_box(ctx, upload_try, empty), 1);
    SixFields after = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, after), 1);
    EXPECT_TRUE(exact_eq(after.Ex, device_pat.Ex));

    fdtd_cuda_destroy(ctx);
}

/*
 * Invariant: host linear index is Fortran column-major
 *   id = i + j*nx + k*nx*ny
 * after a bulk roundtrip of sparse markers.
 */
TEST(cuda, fortran_column_major_layout)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);

    SixFields src = make_six_zero(d);
    const int i0 = 3, j0 = 5, k0 = 2;
    src.Ex[idx3(i0, j0, k0, d.nx, d.ny)] = (fdtd_real)1234.5;
    src.Hz[idx3(1, 2, 7, d.nx, d.ny)] = (fdtd_real)-99.25;

    ASSERT_EQ(upload_six(ctx, src), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);

    EXPECT_EQ(out.Ex[idx3(i0, j0, k0, d.nx, d.ny)], (fdtd_real)1234.5);
    EXPECT_EQ(out.Hz[idx3(1, 2, 7, d.nx, d.ny)], (fdtd_real)-99.25);

    int hits = 0;
    for (size_t n = 0; n < out.Ex.size(); ++n)
        if (out.Ex[n] != (fdtd_real)0) ++hits;
    EXPECT_EQ(hits, 1);

    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_SYNC_TESTS_H */
