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
extern "C" int test_cell_map_coords();
extern "C" int test_cell_map_array();
extern "C" int test_cell_map_add_triangle();
extern "C" int test_cell_map_cellmap_set_get();
extern "C" int test_conformal_filling_off_face_triangle_x();
extern "C" int test_conformal_filling_off_face_triangle_y();
extern "C" int test_conformal_filling_off_face_triangle_z();
extern "C" int test_conformal_filling_open();
extern "C" int test_conformal_filling_closed();
extern "C" int test_conformal_edge_next_cell();
extern "C" int test_conformal_pec_media();
extern "C" int test_conformal_pec_media_raytracing();
extern "C" int test_conformal_pec_corner();
extern "C" int test_conformal_filling_closed_corner();
extern "C" int test_conformal_filling_block_and_corner();

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

TEST(conformal, cell_map_coords)              { EXPECT_EQ(0, test_cell_map_coords()); }
TEST(conformal, cell_map_array)               { EXPECT_EQ(0, test_cell_map_array()); }
TEST(conformal, cell_map_add_triangle)        { EXPECT_EQ(0, test_cell_map_add_triangle()); }
TEST(conformal, cell_map_cellmap_set_get)     { EXPECT_EQ(0, test_cell_map_cellmap_set_get()); }

TEST(conformal, conformal_filling_off_face_triangle_x)     { EXPECT_EQ(0, test_conformal_filling_off_face_triangle_x()); }
TEST(conformal, conformal_filling_off_face_triangle_y)     { EXPECT_EQ(0, test_conformal_filling_off_face_triangle_y()); }
TEST(conformal, conformal_filling_off_face_triangle_z)     { EXPECT_EQ(0, test_conformal_filling_off_face_triangle_z()); }
TEST(conformal, conformal_filling_open)                        { EXPECT_EQ(0, test_conformal_filling_open()); }
TEST(conformal, conformal_filling_closed)                      { EXPECT_EQ(0, test_conformal_filling_closed()); }
TEST(conformal, conformal_edge_next_cell)                      { EXPECT_EQ(0, test_conformal_edge_next_cell()); }
TEST(conformal, conformal_filling_closed_corner)               { EXPECT_EQ(0, test_conformal_filling_closed_corner()); }
TEST(conformal, conformal_filling_block_and_corner)            { EXPECT_EQ(0, test_conformal_filling_block_and_corner()); }

TEST(conformal, conformal_pec_media)        { EXPECT_EQ(0, test_conformal_pec_media()); }
TEST(conformal, conformal_pec_corner)        { EXPECT_EQ(0, test_conformal_pec_corner()); }
TEST(conformal, conformal_pec_media_raytracing)        { EXPECT_EQ(0, test_conformal_pec_media_raytracing()); }
