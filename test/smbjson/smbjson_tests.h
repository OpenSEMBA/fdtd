#include <gtest/gtest.h>

extern "C" int test_idchildtable_fhash();
extern "C" int test_idchildtable();

extern "C" int test_cells();
extern "C" int test_mesh_add_get();
extern "C" int test_mesh_add_get_long_list();
extern "C" int test_mesh_node_to_pixel();
extern "C" int test_mesh_polyline_to_linel();

extern "C" int test_parser_ctor();
extern "C" int test_parser_read_mesh();
extern "C" int test_parser_read_conformal_volume();

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
extern "C" int test_fhash_coords();
extern "C" int test_fhash_array();
extern "C" int test_fhash_add_triangle();
extern "C" int test_fhash_cellmap_set_get();
extern "C" int test_read_planewave();
extern "C" int test_read_sgbc();
extern "C" int test_read_dielectricslab();
extern "C" int test_read_thinslot();
extern "C" int test_read_holland1981();
extern "C" int test_read_towelhanger();
extern "C" int test_read_currentinjection();
extern "C" int test_read_sphere();
extern "C" int test_read_airplane();
extern "C" int test_read_mtln();
extern "C" int test_read_holland1981();
extern "C" int test_read_connectedwires();
extern "C" int test_read_shieldedpair();
extern "C" int test_read_large_airplane_mtln();

TEST(smbjson, idchildtable_fhash)     {EXPECT_EQ(0, test_idchildtable_fhash()); }
TEST(smbjson, idchildtable_add_get)   {EXPECT_EQ(0, test_idchildtable()); }

TEST(smbjson, mesh_cells)                { EXPECT_EQ(0, test_cells()); }
TEST(smbjson, mesh_add_get)              { EXPECT_EQ(0, test_mesh_add_get()); }
TEST(smbjson, mesh_add_get_long_list)    { EXPECT_EQ(0, test_mesh_add_get_long_list()); }
TEST(smbjson, mesh_node_to_pixel)        { EXPECT_EQ(0, test_mesh_node_to_pixel()); }
TEST(smbjson, mesh_polyline_to_linel)    { EXPECT_EQ(0, test_mesh_polyline_to_linel()); }


TEST(smbjson, geometry_coord_position)   { EXPECT_EQ(0, test_geometry_coord_position()); }
TEST(smbjson, geometry_side_position)    { EXPECT_EQ(0, test_geometry_side_position()); }
TEST(smbjson, geometry_triangle_on_face) { EXPECT_EQ(0, test_geometry_triangle_on_face()); }
TEST(smbjson, geometry_triangle_normal)  { EXPECT_EQ(0, test_geometry_triangle_normal()); }
TEST(smbjson, geometry_triangle_edges)   { EXPECT_EQ(0, test_geometry_triangle_edges()); }
TEST(smbjson, geometry_triangle_cells)   { EXPECT_EQ(0, test_geometry_triangle_cell()); }
TEST(smbjson, geometry_elements_in_cell) { EXPECT_EQ(0, test_geometry_elements_in_cell()); }
TEST(smbjson, geometry_path)             { EXPECT_EQ(0, test_geometry_path()); }
TEST(smbjson, geometry_map_sides)        { EXPECT_EQ(0, test_geometry_map_sides()); }
TEST(smbjson, geometry_vv_contour)       { EXPECT_EQ(0, test_geometry_vertex_vertex_contour()); }
TEST(smbjson, geometry_vs_contour)       { EXPECT_EQ(0, test_geometry_vertex_side_contour()); }
TEST(smbjson, geometry_sv_contour)       { EXPECT_EQ(0, test_geometry_side_vertex_contour()); }
TEST(smbjson, geometry_ss_contour)       { EXPECT_EQ(0, test_geometry_side_side_contour()); }
TEST(smbjson, geometry_ss_contour2)      { EXPECT_EQ(0, test_geometry_side_side_contour_2()); }
TEST(smbjson, fhash_coords)              { EXPECT_EQ(0, test_fhash_coords()); }
TEST(smbjson, fhash_array)               { EXPECT_EQ(0, test_fhash_array()); }
TEST(smbjson, fhash_add_triangle)        { EXPECT_EQ(0, test_fhash_add_triangle()); }
TEST(smbjson, fhash_cellmap_set_get)     { EXPECT_EQ(0, test_fhash_cellmap_set_get()); }

TEST(smbjson, parser_ctor)               { EXPECT_EQ(0, test_parser_ctor()); }
TEST(smbjson, parser_read_mesh)          { EXPECT_EQ(0, test_parser_read_mesh()); }
TEST(smbjson, parser_read_conf_volume)   { EXPECT_EQ(0, test_parser_read_conformal_volume()); }
TEST(smbjson, read_planewave)            { EXPECT_EQ(0, test_read_planewave()); }
TEST(smbjson, read_dielectricslab)       { EXPECT_EQ(0, test_read_dielectricslab()); }
TEST(smbjson, read_thinslot)             { EXPECT_EQ(0, test_read_thinslot()); }
TEST(smbjson, read_sgbc)                 { EXPECT_EQ(0, test_read_sgbc()); }
TEST(smbjson, read_sphere)               { EXPECT_EQ(0, test_read_sphere()); }
TEST(smbjson, read_airplane)             { EXPECT_EQ(0, test_read_airplane()); }


#ifdef CompileWithMTLN
TEST(smbjson, read_towelhanger)          { EXPECT_EQ(0, test_read_towelhanger()); }
TEST(smbjson, read_holland1981)          { EXPECT_EQ(0, test_read_holland1981()); }
TEST(smbjson, read_connectedwires)       { EXPECT_EQ(0, test_read_connectedwires()); }
TEST(smbjson, read_currentinjection)     { EXPECT_EQ(0, test_read_currentinjection()); }
// TEST(smbjson, read_shieldedpair)         { EXPECT_EQ(0, test_read_shieldedpair()); }
TEST(smbjson, read_mtln)                 { EXPECT_EQ(0, test_read_mtln()); }
TEST(smbjson, read_large_airplane_mtln)  { EXPECT_EQ(0, test_read_large_airplane_mtln()); }
#endif
