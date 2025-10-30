#include <gtest/gtest.h>

extern "C" int test_init_vtk();

TEST(vtkModule, test_initialize_vtk)    { EXPECT_EQ(0, test_init_vtk()); }


