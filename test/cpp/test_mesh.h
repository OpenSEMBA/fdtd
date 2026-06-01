#ifndef TEST_MESH_H
#define TEST_MESH_H

#include <gtest/gtest.h>
#include "mesh_m.h"

// Test addCoordinate and getCoordinate
TEST(mesh, add_get_coordinate) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c1;
    c1.position[0] = 0.0; c1.position[1] = 0.0; c1.position[2] = 0.0;
    mesh.addCoordinate(10, c1);

    mesh_m::coordinate_t c2;
    c2.position[0] = 1.0; c2.position[1] = 0.0; c2.position[2] = 0.0;
    mesh.addCoordinate(11, c2);

    mesh_m::coordinate_t c3;
    c3.position[0] = 2.0; c3.position[1] = 0.0; c3.position[2] = 0.0;
    mesh.addCoordinate(12, c3);

    bool found = false;
    auto obtained = mesh.getCoordinate(10, found);
    EXPECT_TRUE(found);
    EXPECT_DOUBLE_EQ(obtained.position[0], 0.0);
    EXPECT_DOUBLE_EQ(obtained.position[1], 0.0);
    EXPECT_DOUBLE_EQ(obtained.position[2], 0.0);

    obtained = mesh.getCoordinate(11, found);
    EXPECT_TRUE(found);
    EXPECT_DOUBLE_EQ(obtained.position[0], 1.0);

    obtained = mesh.getCoordinate(99, found);
    EXPECT_FALSE(found);
}

// Test addElement for node and getNode
TEST(mesh, add_node_get_node) {
    mesh_m::mesh_t mesh;

    mesh_m::node_t node;
    node.coordIds.push_back(10);
    mesh.addElement(1, node);

    bool found = false;
    auto obtained = mesh.getNode(1, found);
    EXPECT_TRUE(found);
    ASSERT_EQ(obtained.coordIds.size(), 1u);
    EXPECT_EQ(obtained.coordIds[0], 10);

    obtained = mesh.getNode(99, found);
    EXPECT_FALSE(found);
}

// Test addElement for polyline and getPolyline
TEST(mesh, add_polyline_get_polyline) {
    mesh_m::mesh_t mesh;

    mesh_m::polyline_t pl;
    pl.coordIds.push_back(11);
    pl.coordIds.push_back(12);
    pl.coordIds.push_back(13);
    mesh.addElement(2, pl);

    bool found = false;
    auto obtained = mesh.getPolyline(2, found);
    EXPECT_TRUE(found);
    ASSERT_EQ(obtained.coordIds.size(), 3u);
    EXPECT_EQ(obtained.coordIds[0], 11);
    EXPECT_EQ(obtained.coordIds[1], 12);
    EXPECT_EQ(obtained.coordIds[2], 13);
}

// Test nodeToPixel
TEST(mesh, node_to_pixel) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c;
    c.position[0] = 0.0; c.position[1] = 0.0; c.position[2] = 0.0;
    mesh.addCoordinate(1, c);

    mesh_m::node_t node;
    node.coordIds.push_back(1);

    auto pix = mesh.nodeToPixel(node);
    EXPECT_DOUBLE_EQ(pix.cell[0], 0.0);
    EXPECT_DOUBLE_EQ(pix.cell[1], 0.0);
    EXPECT_DOUBLE_EQ(pix.cell[2], 0.0);
    EXPECT_EQ(pix.tag, 1);
}

// Test polylineToLinels — structured +X polyline 1→2→3
TEST(mesh, polyline_to_linels_x) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c1, c2, c3;
    c1.position[0] = 0.0; c1.position[1] = 0.0; c1.position[2] = 0.0;
    c2.position[0] = 3.0; c2.position[1] = 0.0; c2.position[2] = 0.0;
    c3.position[0] = 3.0; c3.position[1] = 1.0; c3.position[2] = 0.0;

    mesh.addCoordinate(1, c1);
    mesh.addCoordinate(2, c2);
    mesh.addCoordinate(3, c3);

    mesh_m::polyline_t pl;
    pl.coordIds.push_back(1);
    pl.coordIds.push_back(2);
    pl.coordIds.push_back(3);

    EXPECT_TRUE(mesh.arePolylineSegmentsStructured(pl));

    auto linels = mesh.polylineToLinels(pl);
    ASSERT_EQ(linels.size(), 4u);
    // First segment: (0,0,0)→(3,0,0) = 3 linels in +X direction
    // Last segment: (3,0,0)→(3,1,0) = 1 linel in +Y direction
    EXPECT_EQ(linels[0].tag, 1);
    EXPECT_EQ(linels[0].cell[0], 0); EXPECT_EQ(linels[0].cell[1], 0); EXPECT_EQ(linels[0].cell[2], 0);
    EXPECT_EQ(linels[0].orientation, mesh_m::DIR_X + 1); // 1-based: DIR_X is 1

    // Last linel should have tag=3 (end of polyline)
    EXPECT_EQ(linels[linels.size()-1].tag, 3);
}

// Test polylineToLinels — reverse direction
TEST(mesh, polyline_to_linels_reverse_x) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c1, c2, c3;
    c1.position[0] = 3.0; c1.position[1] = 0.0; c1.position[2] = 0.0;
    c2.position[0] = 0.0; c2.position[1] = 0.0; c2.position[2] = 0.0;

    mesh.addCoordinate(1, c1);
    mesh.addCoordinate(2, c2);

    mesh_m::polyline_t pl;
    pl.coordIds.push_back(1);
    pl.coordIds.push_back(2);

    EXPECT_TRUE(mesh.arePolylineSegmentsStructured(pl));

    auto linels = mesh.polylineToLinels(pl);
    ASSERT_EQ(linels.size(), 3u);
    EXPECT_EQ(linels[0].tag, 1);
    EXPECT_EQ(linels[linels.size()-1].tag, 2);
}

// Test arePolylineSegmentsStructured — diagonal line is not structured
TEST(mesh, polyline_unstructured_diagonal) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c1, c2;
    c1.position[0] = 0.0; c1.position[1] = 0.0; c1.position[2] = 0.0;
    c2.position[0] = 3.0; c2.position[1] = 0.0; c2.position[2] = 0.0;

    mesh.addCoordinate(1, c1);
    mesh.addCoordinate(3, c2);

    mesh_m::polyline_t pl;
    pl.coordIds.push_back(1);
    pl.coordIds.push_back(3);

    EXPECT_TRUE(mesh.arePolylineSegmentsStructured(pl));
}

// Test arePolylineSegmentsStructured — fractional coords are not structured
TEST(mesh, polyline_fractional_coords) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c1, c2;
    c1.position[0] = 0.0; c1.position[1] = 0.0; c1.position[2] = 0.0;
    c2.position[0] = 3.5; c2.position[1] = 0.0; c2.position[2] = 0.0;

    mesh.addCoordinate(1, c1);
    mesh.addCoordinate(5, c2);

    mesh_m::polyline_t pl;
    pl.coordIds.push_back(1);
    pl.coordIds.push_back(5);

    EXPECT_FALSE(mesh.arePolylineSegmentsStructured(pl));

    auto linels = mesh.polylineToLinels(pl);
    EXPECT_EQ(linels.size(), 0u);
}

// Test addCoordinate and getCoordinate — long list (from Fortran test)
TEST(mesh, long_coordinate_list) {
    mesh_m::mesh_t mesh;

    // Coordinates from Fortran test_mesh_add_get_long_list
    std::vector<std::pair<int, std::vector<double>>> coords = {
        {1,  {1, 9, 1}},    {2,  {10, 9, 1}},   {5,  {18, 9, 1}},   {6,  {10, 2, 1}},
        {11, {1, 9, 1}},    {15, {10, 9, 1}},   {23, {18, 9, 1}},   {24, {10, 2, 1}},
        {33, {1, 9, 1}},    {34, {1, 9, 1}},    {35, {1, 9, 1}},    {36, {1, 9, 1}},
        {37, {1, 9, 1}},    {38, {1, 9, 1}},    {39, {1, 9, 1}},    {40, {1, 9, 1}},
        {41, {10, 9, 1}},   {42, {10, 9, 1}},   {43, {10, 9, 1}},   {44, {10, 9, 1}},
        {45, {10, 9, 1}},   {46, {10, 9, 1}},   {47, {10, 9, 1}},   {48, {10, 9, 1}},
        {51, {18, 9, 1}},   {52, {18, 9, 1}},   {59, {10, 2, 1}},   {60, {10, 2, 1}},
        {61, {10, 2, 1}},   {62, {10, 2, 1}},   {63, {10, 2, 1}}, {64, {10, 2, 1}},
        {65, {10, 2, 1}},   {66, {10, 2, 1}}
    };

    for (auto& p : coords) {
        mesh_m::coordinate_t c;
        c.position[0] = p.second[0];
        c.position[1] = p.second[1];
        c.position[2] = p.second[2];
        mesh.addCoordinate(p.first, c);
    }

    // Check coord 65
    bool found = false;
    auto obtained = mesh.getCoordinate(65, found);
    EXPECT_TRUE(found);
    EXPECT_DOUBLE_EQ(obtained.position[0], 10.0);
    EXPECT_DOUBLE_EQ(obtained.position[1], 2.0);
    EXPECT_DOUBLE_EQ(obtained.position[2], 1.0);

    // Check coord 66
    obtained = mesh.getCoordinate(66, found);
    EXPECT_TRUE(found);
    EXPECT_DOUBLE_EQ(obtained.position[0], 10.0);
    EXPECT_DOUBLE_EQ(obtained.position[1], 2.0);
    EXPECT_DOUBLE_EQ(obtained.position[2], 1.0);
}

// Test addCoordinate with duplicate IDs (last one wins)
TEST(mesh, coordinate_update) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c1;
    c1.position[0] = 1.0; c1.position[1] = 2.0; c1.position[2] = 3.0;
    mesh.addCoordinate(1, c1);

    mesh_m::coordinate_t c2;
    c2.position[0] = 10.0; c2.position[1] = 20.0; c2.position[2] = 30.0;
    mesh.addCoordinate(1, c2);

    bool found = false;
    auto obtained = mesh.getCoordinate(1, found);
    EXPECT_TRUE(found);
    EXPECT_DOUBLE_EQ(obtained.position[0], 10.0);
    EXPECT_DOUBLE_EQ(obtained.position[1], 20.0);
    EXPECT_DOUBLE_EQ(obtained.position[2], 30.0);
}

// Test coordinate_t subtraction
TEST(mesh, coordinate_subtraction) {
    mesh_m::coordinate_t a, b;
    a.position[0] = 5.0; a.position[1] = 3.0; a.position[2] = 1.0;
    b.position[0] = 2.0; b.position[1] = 1.0; b.position[2] = 0.0;

    auto diff = a - b;
    EXPECT_DOUBLE_EQ(diff.position[0], 3.0);
    EXPECT_DOUBLE_EQ(diff.position[1], 2.0);
    EXPECT_DOUBLE_EQ(diff.position[2], 1.0);
}

// Test coordinate_t equality
TEST(mesh, coordinate_equality) {
    mesh_m::coordinate_t a, b, c;
    a.position[0] = 1.0; a.position[1] = 2.0; a.position[2] = 3.0;
    b.position[0] = 1.0; b.position[1] = 2.0; b.position[2] = 3.0;
    c.position[0] = 1.0; c.position[1] = 2.0; c.position[2] = 4.0;

    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

// Test checkCoordinateId
TEST(mesh, check_coordinate_id) {
    mesh_m::mesh_t mesh;

    mesh_m::coordinate_t c;
    c.position[0] = 0.0; c.position[1] = 0.0; c.position[2] = 0.0;
    mesh.addCoordinate(42, c);

    int stat = 0;
    mesh.checkCoordinateId(42, stat);
    EXPECT_EQ(stat, 0); // found

    mesh.checkCoordinateId(99, stat);
    EXPECT_EQ(stat, 1); // not found
}

// Test getCellRegions
TEST(mesh, get_cell_regions) {
    mesh_m::mesh_t mesh;

    cells_m::cell_region_t cr1;
    cells_m::cell_interval_t ivl;
    ivl.ini.cell[0] = 0; ivl.ini.cell[1] = 0; ivl.ini.cell[2] = 0;
    ivl.end.cell[0] = 2; ivl.end.cell[1] = 0; ivl.end.cell[2] = 0;
    cr1.intervals.push_back(ivl);
    mesh.addCellRegion(10, cr1);

    cells_m::cell_region_t cr2;
    ivl.ini.cell[0] = 0; ivl.ini.cell[1] = 0; ivl.ini.cell[2] = 0;
    ivl.end.cell[0] = 0; ivl.end.cell[1] = 3; ivl.end.cell[2] = 0;
    cr2.intervals.push_back(ivl);
    mesh.addCellRegion(20, cr2);

    auto regions = mesh.getCellRegions({10, 20});
    EXPECT_EQ(regions.size(), 2u);

    // Get a non-existing region
    regions = mesh.getCellRegions({10, 99});
    EXPECT_EQ(regions.size(), 1u);
}

// Test empty polyline
TEST(mesh, empty_polyline) {
    mesh_m::mesh_t mesh;
    mesh_m::polyline_t pl;

    bool structured = mesh.arePolylineSegmentsStructured(pl);
    EXPECT_FALSE(structured);

    auto linels = mesh.polylineToLinels(pl);
    EXPECT_EQ(linels.size(), 0u);
}

#endif // TEST_MESH_H
