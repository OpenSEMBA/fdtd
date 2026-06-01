#ifndef TEST_CONFORMAL_CELL_MAP_H
#define TEST_CONFORMAL_CELL_MAP_H

#include <any>
#include <gtest/gtest.h>
#include <vector>

#include "cell_map_m.h"
#include "conformal_types.h"
#include "fhash_m.h"
#include "test_conformal_helpers.h"

using namespace conformal_types_m;
using conformal_test::expectCellEqual;
using conformal_test::expectPositionsEqual;
using conformal_test::makeCoord;
using conformal_test::makeTriangle;

TEST(conformal, cell_map_coords) {
    fhash_m::fhash_tbl_t tbl;

    triangle_t t_set;
    t_set.vertices[0].position = {0.0, 0.0, 0.0};
    t_set.vertices[1].position = {0.0, 1.0, 0.0};
    t_set.vertices[2].position = {1.0, 0.0, 1.0};

    const std::vector<int> cell_set = t_set.getCell();
    tbl.set(cell_map_m::key(cell_set), t_set);

    std::any t_alloc;
    tbl.get_raw(cell_map_m::key(cell_set), t_alloc);
    const auto t_get = std::any_cast<triangle_t>(t_alloc);

    expectPositionsEqual(t_get.vertices[0], t_set.vertices[0]);
    expectPositionsEqual(t_get.vertices[1], t_set.vertices[1]);
    expectPositionsEqual(t_get.vertices[2], t_set.vertices[2]);
}

TEST(conformal, cell_map_array) {
    cell_map_m::triangle_map_t tbl;

    triangle_t t1;
    t1.vertices[0].position = {0.0, 0.0, 0.0};
    t1.vertices[1].position = {0.0, 1.0, 0.0};
    t1.vertices[2].position = {1.0, 0.0, 1.0};

    cell_map_m::element_set_t tri_list;
    tri_list.triangles = {t1};
    tbl.set(cell_map_m::key(t1.getCell()), tri_list);

    t1.vertices[2].position = {1.0, 0.0, 0.0};
    const std::vector<int> cell_set = t1.getCell();
    ASSERT_TRUE(tbl.hasKey(cell_set));

    std::any t_alloc;
    tbl.get_raw(cell_map_m::key(cell_set), t_alloc);
    auto tri_list_aux = std::any_cast<cell_map_m::element_set_t>(t_alloc);
    ASSERT_EQ(tri_list_aux.triangles.size(), 1u);

    tri_list_aux.triangles.push_back(t1);
    tbl.set(cell_map_m::key(cell_set), tri_list_aux);

    ASSERT_TRUE(tbl.hasKey(cell_set));
    tbl.get_raw(cell_map_m::key(cell_set), t_alloc);
    const auto t_alloc_set = std::any_cast<cell_map_m::element_set_t>(t_alloc);
    ASSERT_EQ(t_alloc_set.triangles.size(), 2u);

    expectPositionsEqual(t_alloc_set.triangles[0].vertices[0], makeCoord({0, 0, 0}, 1));
    expectPositionsEqual(t_alloc_set.triangles[0].vertices[1], makeCoord({0, 1, 0}, 2));
    expectPositionsEqual(t_alloc_set.triangles[0].vertices[2], makeCoord({1, 0, 1}, 3));
    expectPositionsEqual(t_alloc_set.triangles[1].vertices[0], makeCoord({0, 0, 0}, 1));
    expectPositionsEqual(t_alloc_set.triangles[1].vertices[1], makeCoord({0, 1, 0}, 2));
    expectPositionsEqual(t_alloc_set.triangles[1].vertices[2], makeCoord({1, 0, 0}, 3));
}

TEST(conformal, cell_map_add_triangle) {
    cell_map_m::triangle_map_t tbl;
    tbl.keys.clear();

    triangle_t t1;
    t1.vertices[0].position = {0.0, 0.0, 0.0};
    t1.vertices[1].position = {0.0, 1.0, 0.0};
    t1.vertices[2].position = {1.0, 0.0, 1.0};
    std::vector<int> cell_set = t1.getCell();
    tbl.addTriangle(t1);

    ASSERT_TRUE(tbl.hasKey(cell_set));
    std::any t_alloc;
    tbl.get_raw(cell_map_m::key(cell_set), t_alloc);
    const auto elems = std::any_cast<cell_map_m::element_set_t>(t_alloc);
    ASSERT_EQ(elems.triangles.size(), 1u);
    expectPositionsEqual(elems.triangles[0].vertices[0], makeCoord({0, 0, 0}, 1));
    expectPositionsEqual(elems.triangles[0].vertices[1], makeCoord({0, 1, 0}, 2));
    expectPositionsEqual(elems.triangles[0].vertices[2], makeCoord({1, 0, 1}, 3));

    t1.vertices[2].position = {1.0, 0.0, 0.0};
    cell_set = t1.getCell();
    tbl.addTriangle(t1);

    tbl.get_raw(cell_map_m::key(cell_set), t_alloc);
    const auto elems2 = std::any_cast<cell_map_m::element_set_t>(t_alloc);
    ASSERT_EQ(elems2.triangles.size(), 2u);
    expectPositionsEqual(elems2.triangles[0].vertices[2], makeCoord({1, 0, 1}, 3));
    expectPositionsEqual(elems2.triangles[1].vertices[2], makeCoord({1, 0, 0}, 3));

    ASSERT_EQ(tbl.keys.size(), 1u);
    expectCellEqual(tbl.keys[0].cell, {0, 0, 0});
}

TEST(conformal, cell_map_cellmap_set_get) {
    const coord_t c1 = makeCoord({0, 0, 0}, 1);
    const coord_t c2 = makeCoord({1, 0, 0}, 2);
    const coord_t c3 = makeCoord({0, 0, 1}, 3);
    const coord_t c4 = makeCoord({0, 1, 1}, 4);
    const coord_t c5 = makeCoord({0, 1, 0}, 5);

    const triangle_t t1 = makeTriangle(c1, c5, c3);
    const triangle_t t2 = makeTriangle(c1, c2, c5);
    const triangle_t t3 = makeTriangle(c1, c2, c4);

    cell_map_m::triangle_map_t tri_map;
    std::vector<triangle_t> tri_set = {t1};
    cell_map_m::buildMapOfTrisOnFaces(tri_map, tri_set);

    std::vector<int> cell = {static_cast<int>(std::floor(c1.position[0])),
                             static_cast<int>(std::floor(c1.position[1])),
                             static_cast<int>(std::floor(c1.position[2]))};
    std::vector<triangle_t> tri_get = tri_map.getTrianglesInCell(cell);
    ASSERT_EQ(tri_get.size(), 1u);
    EXPECT_EQ(tri_get[0].vertices[0].id, tri_set[0].vertices[0].id);
    EXPECT_EQ(tri_get[0].vertices[1].id, tri_set[0].vertices[1].id);
    EXPECT_EQ(tri_get[0].vertices[2].id, tri_set[0].vertices[2].id);

    tri_set = {t1, t2};
    tri_map.unset(cell_map_m::key(t1.getCell()));
    cell_map_m::buildMapOfTrisOnFaces(tri_map, tri_set);

    cell = {static_cast<int>(std::floor(c1.position[0])),
            static_cast<int>(std::floor(c1.position[1])),
            static_cast<int>(std::floor(c1.position[2]))};
    tri_get = tri_map.getTrianglesInCell(cell);
    ASSERT_EQ(tri_get.size(), 2u);
    EXPECT_EQ(tri_get[0].vertices[0].id, tri_set[0].vertices[0].id);
    EXPECT_EQ(tri_get[0].vertices[1].id, tri_set[0].vertices[1].id);
    EXPECT_EQ(tri_get[0].vertices[2].id, tri_set[0].vertices[2].id);
    EXPECT_EQ(tri_get[1].vertices[0].id, tri_set[1].vertices[0].id);
    EXPECT_EQ(tri_get[1].vertices[1].id, tri_set[1].vertices[1].id);
    EXPECT_EQ(tri_get[1].vertices[2].id, tri_set[1].vertices[2].id);
}

#endif // TEST_CONFORMAL_CELL_MAP_H
