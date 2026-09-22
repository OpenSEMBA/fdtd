#ifndef CUDA_TESTS_H
#define CUDA_TESTS_H

/*
 * CUDA GoogleTest suite — compiled and linked ONLY when
 * SEMBA_FDTD_ENABLE_CUDA=ON (test/CMakeLists.txt + CompileWithCUDA).
 * Not part of basic CPU or MPI CPU builds.
 */
#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include <gtest/gtest.h>

#include "fdtd_cuda.h"

#include "cuda_sync_tests.h"
#include "cuda_yee_tests.h"
#include "cuda_cpml_tests.h"

TEST(cuda, context_create_destroy)
{
    fdtd_cuda_ctx *ctx = fdtd_cuda_create();
    ASSERT_NE(ctx, nullptr);
    if (!fdtd_cuda_ok(ctx)) {
        fdtd_cuda_destroy(ctx);
        GTEST_SKIP() << "No usable CUDA device (fdtd_cuda_ok=0)";
    }
    EXPECT_EQ(fdtd_cuda_ok(ctx), 1);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, alloc_tiny_fields)
{
    fdtd_cuda_ctx *ctx = fdtd_cuda_create();
    ASSERT_NE(ctx, nullptr);
    if (!fdtd_cuda_ok(ctx)) {
        fdtd_cuda_destroy(ctx);
        GTEST_SKIP() << "No usable CUDA device (fdtd_cuda_ok=0)";
    }

    fdtd_dims3 d = {8, 8, 8};
    EXPECT_EQ(fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d), 1);
    EXPECT_EQ(fdtd_cuda_alloc_media(ctx, d, d, d, d, d, d), 1);
    EXPECT_EQ(fdtd_cuda_alloc_coeffs(ctx, 4), 1);
    EXPECT_EQ(fdtd_cuda_alloc_metrics(ctx, 8, 8, 8, 8, 8, 8), 1);

    fdtd_cuda_destroy(ctx);
}

#endif /* CUDA_TESTS_H */
