#include <gtest/gtest.h>


extern "C" int test_geometry_side_position();
extern "C" int test_geometry_coord_position();
extern "C" int test_geometry_triangle_on_face();
extern "C" int test_geometry_triangle_normal();
extern "C" int test_geometry_triangle_edges();
extern "C" int test_geometry_triangle_cell();
extern "C" int test_geometry_elements_in_cell();
extern "C" int test_geometry_path();
extern "C" int test_geometry_map_sides();
extern "C" int test_geometry_vertex_vertex_contour();
extern "C" int test_geometry_vertex_side_contour();
extern "C" int test_geometry_side_vertex_contour();
extern "C" int test_geometry_side_side_contour();
extern "C" int test_geometry_side_side_contour_2();
extern "C" int test_geometry_side_side_contour_3();
extern "C" int test_geometry_areas();
extern "C" int test_fhash_coords();
extern "C" int test_fhash_array();
extern "C" int test_fhash_add_triangle();
extern "C" int test_fhash_cellmap_set_get();
extern "C" int test_conformal_edges_open();
extern "C" int test_conformal_faces_open();
extern "C" int test_conformal_edges_closed();
extern "C" int test_conformal_faces_closed();
extern "C" int test_conformal_edge_next_cell();


TEST(conformal, geometry_coord_position)   { EXPECT_EQ(0, test_geometry_coord_position()); }
TEST(conformal, geometry_side_position)    { EXPECT_EQ(0, test_geometry_side_position()); }
TEST(conformal, geometry_triangle_on_face) { EXPECT_EQ(0, test_geometry_triangle_on_face()); }
TEST(conformal, geometry_triangle_normal)  { EXPECT_EQ(0, test_geometry_triangle_normal()); }
TEST(conformal, geometry_triangle_edges)   { EXPECT_EQ(0, test_geometry_triangle_edges()); }
TEST(conformal, geometry_triangle_cells)   { EXPECT_EQ(0, test_geometry_triangle_cell()); }
TEST(conformal, geometry_elements_in_cell) { EXPECT_EQ(0, test_geometry_elements_in_cell()); }
TEST(conformal, geometry_path)             { EXPECT_EQ(0, test_geometry_path()); }
TEST(conformal, geometry_map_sides)        { EXPECT_EQ(0, test_geometry_map_sides()); }
TEST(conformal, geometry_vv_contour)       { EXPECT_EQ(0, test_geometry_vertex_vertex_contour()); }
TEST(conformal, geometry_vs_contour)       { EXPECT_EQ(0, test_geometry_vertex_side_contour()); }
TEST(conformal, geometry_sv_contour)       { EXPECT_EQ(0, test_geometry_side_vertex_contour()); }
TEST(conformal, geometry_ss_contour)       { EXPECT_EQ(0, test_geometry_side_side_contour()); }
TEST(conformal, geometry_ss_contour2)      { EXPECT_EQ(0, test_geometry_side_side_contour_2()); }
TEST(conformal, geometry_ss_contour3)      { EXPECT_EQ(0, test_geometry_side_side_contour_3()); }
TEST(conformal, geometry_areas)            { EXPECT_EQ(0, test_geometry_areas()); }

TEST(conformal, fhash_coords)              { EXPECT_EQ(0, test_fhash_coords()); }
TEST(conformal, fhash_array)               { EXPECT_EQ(0, test_fhash_array()); }
TEST(conformal, fhash_add_triangle)        { EXPECT_EQ(0, test_fhash_add_triangle()); }
TEST(conformal, fhash_cellmap_set_get)     { EXPECT_EQ(0, test_fhash_cellmap_set_get()); }

TEST(conformal, conformal_edges_open)        { EXPECT_EQ(0, test_conformal_edges_open()); }
TEST(conformal, conformal_faces_open)        { EXPECT_EQ(0, test_conformal_faces_open()); }
TEST(conformal, conformal_edges_closed)        { EXPECT_EQ(0, test_conformal_edges_closed()); }
TEST(conformal, conformal_faces_closed)        { EXPECT_EQ(0, test_conformal_faces_closed()); }
TEST(conformal, conformal_edge_next_cell)        { EXPECT_EQ(0, test_conformal_edge_next_cell()); }
