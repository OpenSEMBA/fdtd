#ifdef CompileWithNewOutputModule
#include <gtest/gtest.h>

extern "C" int test_vtkapi_points_allocation();
extern "C" int test_vtkapi_point_scalar();
extern "C" int test_vtkapi_point_vector();
extern "C" int test_vtkapi_cell_scalar();
extern "C" int test_vtkapi_cell_vector();
extern "C" int test_vtkapi_vts_file_creation();
extern "C" int test_vtkapi_vtu_file_creation();
extern "C" int test_vtkapi_vtu_cell_data();
extern "C" int test_vtkapi_vts_content();
extern "C" int test_vtkapi_vtu_content();

TEST(vtkapi, test_points_allocation) {EXPECT_EQ(0, test_vtkapi_points_allocation());}
TEST(vtkapi, test_point_scalar) {EXPECT_EQ(0, test_vtkapi_point_scalar());}
TEST(vtkapi, test_point_vector) {EXPECT_EQ(0, test_vtkapi_point_vector());}
TEST(vtkapi, test_cell_scalar) {EXPECT_EQ(0, test_vtkapi_cell_scalar());}
TEST(vtkapi, test_cell_vector) {EXPECT_EQ(0, test_vtkapi_cell_vector());}
TEST(vtkapi, test_vts_file_creation) {EXPECT_EQ(0, test_vtkapi_vts_file_creation());}
TEST(vtkapi, test_vtu_file_creation) {EXPECT_EQ(0, test_vtkapi_vtu_file_creation());}
TEST(vtkapi, test_vtu_cell_data) {EXPECT_EQ(0, test_vtkapi_vtu_cell_data());}
TEST(vtkapi, test_vts_content) {EXPECT_EQ(0, test_vtkapi_vts_content());}
TEST(vtkapi, test_vtu_content) {EXPECT_EQ(0, test_vtkapi_vtu_content());}


#endif
