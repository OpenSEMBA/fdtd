#ifdef CompileWithNewOutputModule
#include <gtest/gtest.h>

extern "C" int test_create_h5_file();
extern "C" int test_write_1d_dataset();
extern "C" int test_write_2d_dataset();
extern "C" int test_write_3d_dataset();
extern "C" int test_xdmf_file_creation();
extern "C" int test_xdmf_file_with_h5data();

TEST(xdmfapi, test_create_h5) { EXPECT_EQ(0, test_create_h5_file()); }
TEST(xdmfapi, test_write_1d) { EXPECT_EQ(0, test_write_1d_dataset()); }
TEST(xdmfapi, test_write_2d) { EXPECT_EQ(0, test_write_2d_dataset()); }
TEST(xdmfapi, test_write_3d) { EXPECT_EQ(0, test_write_3d_dataset()); }
TEST(xdmfapi, test_xdmf_file) { EXPECT_EQ(0, test_xdmf_file_creation()); }
TEST(xdmfapi, test_xdmf_file_with_h5) { EXPECT_EQ(0, test_xdmf_file_with_h5data()); }

#endif