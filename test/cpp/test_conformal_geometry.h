#ifndef TEST_CONFORMAL_GEOMETRY_H
#define TEST_CONFORMAL_GEOMETRY_H

#include <gtest/gtest.h>
#include <vector>

#include "conformal_types.h"
#include "geometry_m.h"
#include "cell_map_m.h"
#include "test_conformal_helpers.h"

using namespace conformal_types_m;
using conformal_test::expectCellEqual;
using conformal_test::expectPositionsEqual;
using conformal_test::expectSideEndpointPositions;
using conformal_test::makeCoord;
using conformal_test::makeSide;
using conformal_test::makeTriangle;

TEST(conformal, geometry_coord_position) {
    coord_t c = makeCoord({0.0, 0.0, 0.0}, 1);
    EXPECT_TRUE(c.isOnVertex());
    EXPECT_FALSE(c.isOnAnyEdge());
    EXPECT_FALSE(c.isOnAnyFace());

    c = makeCoord({1.0, 0.0, 0.0}, 1);
    EXPECT_TRUE(c.isOnVertex());
    EXPECT_FALSE(c.isOnAnyEdge());
    EXPECT_FALSE(c.isOnAnyFace());

    c = makeCoord({0.0, 1.0, 0.0}, 1);
    EXPECT_TRUE(c.isOnVertex());
    EXPECT_FALSE(c.isOnAnyEdge());
    EXPECT_FALSE(c.isOnAnyFace());

    c = makeCoord({0.0, 0.0, 1.0}, 1);
    EXPECT_TRUE(c.isOnVertex());
    EXPECT_FALSE(c.isOnAnyEdge());
    EXPECT_FALSE(c.isOnAnyFace());

    c = makeCoord({0.5, 0.0, 0.0}, 1);
    EXPECT_FALSE(c.isOnVertex());
    EXPECT_TRUE(c.isOnEdge(EDGE_X));
    EXPECT_FALSE(c.isOnAnyFace());

    c = makeCoord({0.0, 0.5, 0.0}, 1);
    EXPECT_FALSE(c.isOnVertex());
    EXPECT_TRUE(c.isOnEdge(EDGE_Y));
    EXPECT_FALSE(c.isOnAnyFace());

    c = makeCoord({0.0, 0.0, 0.5}, 1);
    EXPECT_FALSE(c.isOnVertex());
    EXPECT_TRUE(c.isOnEdge(EDGE_Z));
    EXPECT_FALSE(c.isOnAnyFace());

    c = makeCoord({0.5, 0.5, 0.0}, 1);
    EXPECT_FALSE(c.isOnVertex());
    EXPECT_FALSE(c.isOnAnyEdge());
    EXPECT_TRUE(c.isOnFace(FACE_Z));

    c = makeCoord({0.5, 0.0, 0.5}, 1);
    EXPECT_FALSE(c.isOnVertex());
    EXPECT_FALSE(c.isOnAnyEdge());
    EXPECT_TRUE(c.isOnFace(FACE_Y));

    c = makeCoord({0.0, 0.5, 0.5}, 1);
    EXPECT_FALSE(c.isOnVertex());
    EXPECT_FALSE(c.isOnAnyEdge());
    EXPECT_TRUE(c.isOnFace(FACE_X));
}

TEST(conformal, geometry_side_position) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 1.0}, 3);
    const coord_t c4 = makeCoord({0.0, 1.0, 1.0}, 4);
    const coord_t c5 = makeCoord({0.0, 1.0, 0.0}, 5);

    side_t side = makeSide(c1, c2);
    EXPECT_TRUE(side.isOnEdge(EDGE_X));
    EXPECT_FALSE(side.isOnAnyFace());

    side = makeSide(c1, c5);
    EXPECT_TRUE(side.isOnEdge(EDGE_Y));
    EXPECT_FALSE(side.isOnAnyFace());

    side = makeSide(c1, c3);
    EXPECT_TRUE(side.isOnEdge(EDGE_Z));
    EXPECT_FALSE(side.isOnAnyFace());

    side = makeSide(c1, c4);
    EXPECT_FALSE(side.isOnAnyEdge());
    EXPECT_TRUE(side.isOnFace(FACE_X));
}

TEST(conformal, geometry_triangle_on_face) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 1.0}, 3);
    const coord_t c4 = makeCoord({0.0, 1.0, 1.0}, 4);
    const coord_t c5 = makeCoord({0.0, 1.0, 0.0}, 5);

    triangle_t t = makeTriangle(c1, c3, c5);
    EXPECT_TRUE(t.isOnAnyFace());
    EXPECT_EQ(t.getFace(), FACE_X);

    t = makeTriangle(c1, c2, c3);
    EXPECT_TRUE(t.isOnAnyFace());
    EXPECT_EQ(t.getFace(), FACE_Y);

    t = makeTriangle(c1, c2, c5);
    EXPECT_TRUE(t.isOnAnyFace());
    EXPECT_EQ(t.getFace(), FACE_Z);

    t = makeTriangle(c1, c2, c4);
    EXPECT_FALSE(t.isOnAnyFace());
    EXPECT_EQ(t.getFace(), NOT_ON_FACE);
}

TEST(conformal, geometry_triangle_normal) {
    triangle_t t;

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 0.0, 1.0};
    t.vertices[2].position = {0.0, 1.0, 0.0};
    EXPECT_EQ(t.getFace(), FACE_X);

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 0.0, 1.0};
    t.vertices[2].position = {1.0, 0.0, 0.0};
    EXPECT_EQ(t.getFace(), FACE_Y);

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 1.0, 0.0};
    t.vertices[2].position = {1.0, 0.0, 0.0};
    EXPECT_EQ(t.getFace(), FACE_Z);

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 1.0, 0.0};
    t.vertices[2].position = {1.0, 0.0, 1.0};
    EXPECT_EQ(t.getFace(), NOT_ON_FACE);
}

TEST(conformal, geometry_triangle_edges) {
    triangle_t t;

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 0.0, 1.0};
    t.vertices[2].position = {0.0, 1.0, 0.0};
    std::vector<side_t> sides = t.getSides();
    EXPECT_EQ(sides[0].getEdge(), EDGE_Z);
    EXPECT_EQ(sides[1].getEdge(), NOT_ON_EDGE);
    EXPECT_EQ(sides[2].getEdge(), EDGE_Y);

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 0.0, 1.0};
    t.vertices[2].position = {1.0, 0.0, 0.0};
    sides = t.getSides();
    EXPECT_EQ(sides[0].getEdge(), EDGE_Z);
    EXPECT_EQ(sides[1].getEdge(), NOT_ON_EDGE);
    EXPECT_EQ(sides[2].getEdge(), EDGE_X);

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 1.0, 0.0};
    t.vertices[2].position = {1.0, 0.0, 0.0};
    sides = t.getSides();
    EXPECT_EQ(sides[0].getEdge(), EDGE_Y);
    EXPECT_EQ(sides[1].getEdge(), NOT_ON_EDGE);
    EXPECT_EQ(sides[2].getEdge(), EDGE_X);

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 1.0, 0.0};
    t.vertices[2].position = {1.0, 0.0, 1.0};
    sides = t.getSides();
    EXPECT_EQ(sides[0].getEdge(), EDGE_Y);
    EXPECT_EQ(sides[1].getEdge(), NOT_ON_EDGE);
    EXPECT_EQ(sides[2].getEdge(), NOT_ON_EDGE);
}

TEST(conformal, geometry_triangle_cell) {
    triangle_t t;

    t.vertices[0].position = {0.0, 0.0, 0.0};
    t.vertices[1].position = {0.0, 0.0, 1.0};
    t.vertices[2].position = {0.0, 1.0, 0.0};
    expectCellEqual(t.getCell(), {0, 0, 0});

    t.vertices[0].position = {1.0, 0.0, 0.0};
    t.vertices[1].position = {1.0, 0.0, 1.0};
    t.vertices[2].position = {1.0, 1.0, 0.0};
    expectCellEqual(t.getCell(), {1, 0, 0});

    t.vertices[0].position = {1.0, 0.0, 1.0};
    t.vertices[1].position = {1.0, 0.0, 2.0};
    t.vertices[2].position = {1.0, 1.0, 1.0};
    expectCellEqual(t.getCell(), {1, 0, 1});

    t.vertices[0].position = {0.0, 0.0, 1.0};
    t.vertices[1].position = {1.0, 0.0, 1.0};
    t.vertices[2].position = {1.0, 1.0, 1.0};
    expectCellEqual(t.getCell(), {0, 0, 1});
}

TEST(conformal, geometry_elements_in_cell) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 1.0}, 3);
    const coord_t c4 = makeCoord({0.0, 1.0, 1.0}, 4);
    const coord_t c5 = makeCoord({0.0, 1.0, 0.0}, 5);

    const std::vector<triangle_t> triangles = {
        makeTriangle(c1, c5, c3),
        makeTriangle(c1, c2, c5),
        makeTriangle(c1, c2, c4),
    };

    cell_map_m::triangle_map_t tri_map;
    cell_map_m::side_map_t side_map;
    cell_map_m::buildMapOfTrisOnFaces(tri_map, triangles);

    const std::vector<int> cell = triangles[0].getCell();
    const std::vector<triangle_t> tris = tri_map.getTrianglesInCell(cell);
    EXPECT_EQ(tris.size(), 2u);

    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    const std::vector<side_t> sides = side_map.getSidesInCell(cell);
    EXPECT_EQ(sides.size(), 2u);
}

namespace {

std::vector<triangle_t> makeEightTriangleSet() {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.75, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.75, 0.25, 0.0}, 3);
    const coord_t c4 = makeCoord({0.25, 0.75, 0.0}, 4);
    const coord_t c5 = makeCoord({0.0, 0.75, 0.0}, 5);
    const coord_t c6 = makeCoord({0.0, 0.0, 0.75}, 6);

    return {
        makeTriangle(c1, c2, c6),
        makeTriangle(c1, c6, c5),
        makeTriangle(c1, c3, c2),
        makeTriangle(c1, c4, c3),
        makeTriangle(c1, c5, c4),
        makeTriangle(c2, c3, c6),
        makeTriangle(c3, c4, c6),
        makeTriangle(c4, c5, c6),
    };
}

} // namespace

TEST(conformal, geometry_map_sides) {
    const std::vector<triangle_t> triangles = makeEightTriangleSet();
    const std::vector<int> cell = triangles[0].getCell();

    cell_map_m::triangle_map_t tri_map;
    cell_map_m::side_map_t side_map;
    cell_map_m::buildMapOfTrisOnFaces(tri_map, triangles);

    const std::vector<triangle_t> tris = tri_map.getTrianglesInCell(cell);
    EXPECT_EQ(tris.size(), 5u);

    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    const std::vector<side_t> sides = side_map.getSidesInCell(cell);
    EXPECT_EQ(sides.size(), 5u);

    EXPECT_EQ(geometry_m::getSidesOnFace(sides, FACE_X).size(), 1u);
    EXPECT_EQ(geometry_m::getSidesOnFace(sides, FACE_Y).size(), 1u);
    EXPECT_EQ(geometry_m::getSidesOnFace(sides, FACE_Z).size(), 3u);

    EXPECT_EQ(geometry_m::getSidesOnEdge(sides, EDGE_X).size(), 0u);
    EXPECT_EQ(geometry_m::getSidesOnEdge(sides, EDGE_Y).size(), 0u);
    EXPECT_EQ(geometry_m::getSidesOnEdge(sides, EDGE_Z).size(), 0u);
}

TEST(conformal, geometry_path) {
    const std::vector<triangle_t> triangles = makeEightTriangleSet();
    const std::vector<int> cell = triangles[0].getCell();

    cell_map_m::side_map_t side_map;
    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    const std::vector<side_t> sides = side_map.getSidesInCell(cell);

    const std::vector<side_t> sides_on_face = geometry_m::getSidesOnFace(sides, FACE_Z);
    const std::vector<side_t> path = geometry_m::getPathOnFace(sides_on_face);

    EXPECT_EQ(path[0].init.id, 2);
    EXPECT_EQ(path[0].end.id, 3);
    EXPECT_EQ(path[1].init.id, 3);
    EXPECT_EQ(path[1].end.id, 4);
    EXPECT_EQ(path[2].init.id, 4);
    EXPECT_EQ(path[2].end.id, 5);
}

namespace {

void expectContourClosureOnZFace(const std::vector<side_t>& contour,
                                 const coord_t& c1,
                                 const coord_t& c2,
                                 const coord_t& c5) {
    EXPECT_EQ(contour.size(), 5u);
    expectSideEndpointPositions(contour[3], c5, c1);
    expectSideEndpointPositions(contour[4], c1, c2);
}

std::vector<side_t> buildZFaceContour(const std::vector<triangle_t>& triangles) {
    cell_map_m::side_map_t side_map;
    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    const std::vector<side_t> sides = side_map.getSidesInCell(triangles[0].getCell());
    const std::vector<side_t> sides_on_face = geometry_m::getSidesOnFace(sides, FACE_Z);
    const std::vector<side_t> path = geometry_m::getPathOnFace(sides_on_face);
    return geometry_m::buildSidesContour(path);
}

} // namespace

TEST(conformal, geometry_vertex_vertex_contour) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.75, 0.25, 0.0}, 3);
    const coord_t c4 = makeCoord({0.25, 0.75, 0.0}, 4);
    const coord_t c5 = makeCoord({0.0, 1.0, 0.0}, 5);
    const coord_t c6 = makeCoord({0.0, 0.0, 1.0}, 6);

    const std::vector<triangle_t> triangles = {
        makeTriangle(c1, c2, c6),
        makeTriangle(c1, c6, c5),
        makeTriangle(c1, c3, c2),
        makeTriangle(c1, c4, c3),
        makeTriangle(c1, c5, c4),
        makeTriangle(c2, c3, c6),
        makeTriangle(c3, c4, c6),
        makeTriangle(c4, c5, c6),
    };

    const std::vector<side_t> contour = buildZFaceContour(triangles);
    expectContourClosureOnZFace(contour, c1, c2, c5);
}

TEST(conformal, geometry_vertex_side_contour) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.75, 0.25, 0.0}, 3);
    const coord_t c4 = makeCoord({0.25, 0.75, 0.0}, 4);
    const coord_t c5 = makeCoord({0.0, 0.75, 0.0}, 5);
    const coord_t c6 = makeCoord({0.0, 0.0, 1.00}, 6);

    const std::vector<triangle_t> triangles = {
        makeTriangle(c1, c2, c6),
        makeTriangle(c1, c6, c5),
        makeTriangle(c1, c3, c2),
        makeTriangle(c1, c4, c3),
        makeTriangle(c1, c5, c4),
        makeTriangle(c2, c3, c6),
        makeTriangle(c3, c4, c6),
        makeTriangle(c4, c5, c6),
    };

    const std::vector<side_t> contour = buildZFaceContour(triangles);
    expectContourClosureOnZFace(contour, c1, c2, c5);
}

TEST(conformal, geometry_side_vertex_contour) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.75, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.75, 0.25, 0.0}, 3);
    const coord_t c4 = makeCoord({0.25, 0.75, 0.0}, 4);
    const coord_t c5 = makeCoord({0.0, 1.0, 0.0}, 5);
    const coord_t c6 = makeCoord({0.0, 0.0, 1.00}, 6);

    const std::vector<triangle_t> triangles = {
        makeTriangle(c1, c2, c6),
        makeTriangle(c1, c6, c5),
        makeTriangle(c1, c3, c2),
        makeTriangle(c1, c4, c3),
        makeTriangle(c1, c5, c4),
        makeTriangle(c2, c3, c6),
        makeTriangle(c3, c4, c6),
        makeTriangle(c4, c5, c6),
    };

    const std::vector<side_t> contour = buildZFaceContour(triangles);
    expectContourClosureOnZFace(contour, c1, c2, c5);
}

TEST(conformal, geometry_side_side_contour) {
    const std::vector<triangle_t> triangles = makeEightTriangleSet();
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.75, 0.0, 0.0}, 2);
    const coord_t c5 = makeCoord({0.0, 0.75, 0.0}, 5);

    const std::vector<side_t> contour = buildZFaceContour(triangles);
    expectContourClosureOnZFace(contour, c1, c2, c5);
}

TEST(conformal, geometry_side_side_contour_2) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.5, 0.0}, 3);
    const coord_t c4 = makeCoord({1.0, 0.5, 0.0}, 4);
    const coord_t c5 = makeCoord({0.25, 0.25, 0.0}, 5);
    const coord_t c6 = makeCoord({0.0, 0.0, 1.0}, 6);
    const coord_t c7 = makeCoord({1.0, 0.0, 1.0}, 7);

    const std::vector<triangle_t> triangles = {
        makeTriangle(c1, c7, c6),
        makeTriangle(c1, c2, c7),
        makeTriangle(c1, c3, c6),
        makeTriangle(c2, c4, c7),
        makeTriangle(c1, c3, c5),
        makeTriangle(c1, c5, c2),
        makeTriangle(c2, c5, c4),
        makeTriangle(c5, c3, c6),
        makeTriangle(c4, c5, c7),
        makeTriangle(c5, c6, c7),
    };

    cell_map_m::side_map_t side_map;
    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    const std::vector<side_t> sides = side_map.getSidesInCell(triangles[0].getCell());
    EXPECT_EQ(sides.size(), 3u);

    const std::vector<side_t> sides_on_face = geometry_m::getSidesOnFace(sides, FACE_Z);
    const std::vector<side_t> path = geometry_m::getPathOnFace(sides_on_face);
    const std::vector<side_t> contour = geometry_m::buildSidesContour(path);

    EXPECT_EQ(contour.size(), 5u);
    expectSideEndpointPositions(contour[2], c3, c1);
    expectSideEndpointPositions(contour[3], c1, c2);
    expectSideEndpointPositions(contour[4], c2, c4);
}

TEST(conformal, geometry_side_side_contour_3) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.25}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.25}, 2);
    const coord_t c3 = makeCoord({0.0, 1.0, 0.25}, 3);
    const coord_t c4 = makeCoord({0.0, 0.0, 0.0}, 4);
    const coord_t c5 = makeCoord({1.0, 0.0, 0.0}, 5);

    const std::vector<triangle_t> triangles = {makeTriangle(c1, c2, c3)};

    cell_map_m::side_map_t side_map;
    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    const std::vector<side_t> sides = side_map.getSidesInCell(triangles[0].getCell());
    EXPECT_EQ(sides.size(), 2u);

    const std::vector<side_t> sides_on_face = geometry_m::getSidesOnFace(sides, FACE_Y);
    const std::vector<side_t> path = geometry_m::getPathOnFace(sides_on_face);
    const std::vector<side_t> contour = geometry_m::buildSidesContour(path);

    EXPECT_EQ(contour.size(), 4u);
    expectSideEndpointPositions(contour[0], c1, c2);
    expectSideEndpointPositions(contour[1], c2, c5);
    expectSideEndpointPositions(contour[2], c5, c4);
    expectSideEndpointPositions(contour[3], c4, c1);
}

TEST(conformal, geometry_areas) {
    const coord_t c1 = makeCoord({1.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.0, 1.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 1.0}, 3);

    std::vector<triangle_t> triangles = {makeTriangle(c1, c2, c3)};
    const std::vector<int> cell = triangles[0].getCell();

    cell_map_m::side_map_t side_map;
    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    std::vector<side_t> sides = side_map.getSidesInCell(cell);
    EXPECT_EQ(sides.size(), 3u);

    std::vector<side_t> sides_on_face = geometry_m::getSidesOnFace(sides, FACE_Z);
    std::vector<side_t> path = geometry_m::getPathOnFace(sides_on_face);
    std::vector<side_t> contour = geometry_m::buildSidesContour(path);
    EXPECT_EQ(geometry_m::contourArea(contour), 0.5);

    side_map.unset(cell_map_m::key(cell));
    triangles[0].vertices[0].position = {0.5, 0.0, 0.0};
    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    sides = side_map.getSidesInCell(cell);
    EXPECT_EQ(sides.size(), 3u);

    sides_on_face = geometry_m::getSidesOnFace(sides, FACE_Y);
    path = geometry_m::getPathOnFace(sides_on_face);
    contour = geometry_m::buildSidesContour(path);
    EXPECT_EQ(geometry_m::contourArea(contour), 0.25);

    side_map.unset(cell_map_m::key(cell));
    triangles[0].vertices[1].position = {0.0, 0.25, 0.0};
    cell_map_m::buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, triangles);
    sides = side_map.getSidesInCell(cell);
    EXPECT_EQ(sides.size(), 3u);

    sides_on_face = geometry_m::getSidesOnFace(sides, FACE_X);
    path = geometry_m::getPathOnFace(sides_on_face);
    contour = geometry_m::buildSidesContour(path);
    EXPECT_EQ(geometry_m::contourArea(contour), 0.125);
}

#endif // TEST_CONFORMAL_GEOMETRY_H
