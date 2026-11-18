#include <gtest/gtest.h>

extern "C" int test_lorentz_ade_single_node();

// Unit test of the electric dispersive subroutines (InitEDispersives /
// AdvanceEDispersiveE) on a one-node mesh. The full-slab S11 checks against
// the dgtd and scikit-rf vectors live in test/pyWrapper/test_lorentz.py.
TEST(lorentz, ade_single_node) { EXPECT_EQ(0, test_lorentz_ade_single_node()); }
