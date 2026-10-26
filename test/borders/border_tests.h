#include <gtest/gtest.h>

extern "C" int test_periodic_x_ghost_planes();
extern "C" int test_periodic_y_ghost_planes();
extern "C" int test_periodic_z_ghost_planes();

TEST(periodic_borders, x_ghost_planes) {
    EXPECT_EQ(0, test_periodic_x_ghost_planes());
}

TEST(periodic_borders, y_ghost_planes) {
    EXPECT_EQ(0, test_periodic_y_ghost_planes());
}

TEST(periodic_borders, z_ghost_planes) {
    EXPECT_EQ(0, test_periodic_z_ghost_planes());
}
