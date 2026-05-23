#ifndef TEST_SMBJSON_PARSER_MESH_H
#define TEST_SMBJSON_PARSER_MESH_H

#include <gtest/gtest.h>
#include "test_smbjson_helpers.h"
#include "smbjson_m.h"

TEST(smbjson_cpp, parser_read_mesh) {
    smbjson::parser_t parser(testDataPath("mtln.fdtd.json"));
    ASSERT_TRUE(parser.isInitialized);
    Mesh::mesh_t mesh = parser.readMesh();

    bool found = false;
    Mesh::coordinate_t c59 = mesh.getCoordinate(59, found);
    ASSERT_TRUE(found);
    EXPECT_DOUBLE_EQ(c59.position[0], 10.0);
    EXPECT_DOUBLE_EQ(c59.position[1], 0.0);
    EXPECT_DOUBLE_EQ(c59.position[2], 1.0);

    Mesh::coordinate_t c64 = mesh.getCoordinate(64, found);
    ASSERT_TRUE(found);
    EXPECT_DOUBLE_EQ(c64.position[0], 10.0);
    EXPECT_DOUBLE_EQ(c64.position[1], 0.0);
    EXPECT_DOUBLE_EQ(c64.position[2], 1.0);

    Mesh::coordinate_t c61 = mesh.getCoordinate(61, found);
    ASSERT_TRUE(found);
    EXPECT_DOUBLE_EQ(c61.position[0], 10.0);
    EXPECT_DOUBLE_EQ(c61.position[1], 0.0);
    EXPECT_DOUBLE_EQ(c61.position[2], 1.0);
}

TEST(smbjson_cpp, parser_read_conformal_volume) {
    smbjson::parser_t parser(testDataPath("conformal.fdtd.json"));
    ASSERT_TRUE(parser.isInitialized);
    Mesh::mesh_t mesh = parser.readMesh();

    bool found = false;
    auto conformal_regions = mesh.getConformalRegions({5});
    ASSERT_EQ(conformal_regions.size(), 1u);
    EXPECT_EQ(conformal_regions[0].triangles.size(), 24u);

    auto cell_regions = mesh.getCellRegions({5});
    ASSERT_EQ(cell_regions.size(), 1u);
    ASSERT_EQ(cell_regions[0].intervals.size(), 1u);
}

#endif // TEST_SMBJSON_PARSER_MESH_H
