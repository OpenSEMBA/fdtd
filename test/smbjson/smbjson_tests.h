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

extern "C" int test_read_planewave();
extern "C" int test_read_sgbc();
extern "C" int test_read_dielectricslab();
extern "C" int test_read_holland1981();
extern "C" int test_read_towelhanger();
extern "C" int test_read_connectedwires();
extern "C" int test_read_currentinjection();
extern "C" int test_read_shieldedpair();
extern "C" int test_read_mtln();
extern "C" int test_read_sphere();
extern "C" int test_read_airplane();
extern "C" int test_read_large_airplane_mtln();

TEST(smbjson, idchildtable_fhash)     {EXPECT_EQ(0, test_idchildtable_fhash()); }
TEST(smbjson, idchildtable_add_get)   {EXPECT_EQ(0, test_idchildtable()); }

TEST(smbjson, mesh_cells)                { EXPECT_EQ(0, test_cells()); }
TEST(smbjson, mesh_add_get)              { EXPECT_EQ(0, test_mesh_add_get()); }
TEST(smbjson, mesh_add_get_long_list)    { EXPECT_EQ(0, test_mesh_add_get_long_list()); }
TEST(smbjson, mesh_node_to_pixel)        { EXPECT_EQ(0, test_mesh_node_to_pixel()); }
TEST(smbjson, mesh_polyline_to_linel)    { EXPECT_EQ(0, test_mesh_polyline_to_linel()); }
   
TEST(smbjson, parser_ctor)               { EXPECT_EQ(0, test_parser_ctor()); }
TEST(smbjson, parser_read_mesh)          { EXPECT_EQ(0, test_parser_read_mesh()); }
TEST(smbjson, read_planewave)            { EXPECT_EQ(0, test_read_planewave()); }
TEST(smbjson, read_dielectricslab)       { EXPECT_EQ(0, test_read_dielectricslab()); }
TEST(smbjson, read_sgbc)                 { EXPECT_EQ(0, test_read_sgbc()); }
TEST(smbjson, read_holland1981)          { EXPECT_EQ(0, test_read_holland1981()); }
TEST(smbjson, read_towelhanger)          { EXPECT_EQ(0, test_read_towelhanger()); }
TEST(smbjson, read_connectedwires)       { EXPECT_EQ(0, test_read_connectedwires()); }
TEST(smbjson, read_currentinjection)     { EXPECT_EQ(0, test_read_currentinjection()); }
// TEST(smbjson, read_shieldedpair)         { EXPECT_EQ(0, test_read_shieldedpair()); }
TEST(smbjson, read_mtln)                 { EXPECT_EQ(0, test_read_mtln()); }
TEST(smbjson, read_sphere)               { EXPECT_EQ(0, test_read_sphere()); }
TEST(smbjson, read_airplane)             { EXPECT_EQ(0, test_read_airplane()); }
TEST(smbjson, read_large_airplane_mtln)  { EXPECT_EQ(0, test_read_large_airplane_mtln()); }