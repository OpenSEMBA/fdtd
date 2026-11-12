#include <gtest/gtest.h>

extern "C" int test_lorentz_res2_s11();
extern "C" int test_lorentz_plasma_s11();

// The Lorentz slab simulations need a single-rank run; the system tests use the
// same guard. S11 references come from dgtd and scikit-rf (PEC-backed slab, S21 = 0).
#ifndef CompileWithMPI
TEST(lorentz, res2_s11_matches_dgtd_and_scikit_rf)   { EXPECT_EQ(0, test_lorentz_res2_s11()); }
TEST(lorentz, plasma_s11_matches_dgtd_and_scikit_rf) { EXPECT_EQ(0, test_lorentz_plasma_s11()); }
#endif
