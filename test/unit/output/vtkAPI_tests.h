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

// Allocates the structured-grid point array with the expected shape.
TEST(vtkapi, test_points_allocation) {EXPECT_EQ(0, test_vtkapi_points_allocation());}
// Stores scalar point data in a structured grid.
TEST(vtkapi, test_point_scalar) {EXPECT_EQ(0, test_vtkapi_point_scalar());}
// Stores vector point data in a structured grid.
TEST(vtkapi, test_point_vector) {EXPECT_EQ(0, test_vtkapi_point_vector());}
// Stores scalar cell data in a structured grid.
TEST(vtkapi, test_cell_scalar) {EXPECT_EQ(0, test_vtkapi_cell_scalar());}
// Stores vector cell data in a structured grid.
TEST(vtkapi, test_cell_vector) {EXPECT_EQ(0, test_vtkapi_cell_vector());}
// Creates a readable VTS structured-grid file.
TEST(vtkapi, test_vts_file_creation) {EXPECT_EQ(0, test_vtkapi_vts_file_creation());}
// Creates a readable VTU unstructured-grid file.
TEST(vtkapi, test_vtu_file_creation) {EXPECT_EQ(0, test_vtkapi_vtu_file_creation());}
// Stores scalar and vector data associated with unstructured-grid cells.
TEST(vtkapi, test_vtu_cell_data) {EXPECT_EQ(0, test_vtkapi_vtu_cell_data());}
// Writes scalar and vector point-data declarations to a VTS file.
TEST(vtkapi, test_vts_content) {EXPECT_EQ(0, test_vtkapi_vts_content());}
// Writes points, cells, and point/cell data declarations to a VTU file.
TEST(vtkapi, test_vtu_content) {EXPECT_EQ(0, test_vtkapi_vtu_content());}


#endif
