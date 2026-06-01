#ifndef TEST_CONFORMAL_FILLING_H
#define TEST_CONFORMAL_FILLING_H

#include <gtest/gtest.h>
#include <vector>

#include "conformal_m.h"
#include "test_conformal_helpers.h"

using namespace conformal_types_m;
using conformal_m::ConformalMedia_t;
using conformal_test::makeCoord;
using conformal_test::makeTriangle;

namespace {

constexpr double kTol = 0.01;

conformal_m::ConformalPECRegions_t makeRegion(const std::vector<triangle_t>& tris) {
    conformal_m::ConformalPECRegions_t regions;
    cell_map_m::ConformalPECElements_t volume;
    volume.triangles = tris;
    regions.volumes.push_back(volume);
    return regions;
}

void expectRatioNear(double actual, double expected) {
    EXPECT_NEAR(actual, expected, kTol);
}

void checkOffFaceTriangleXPlus(const ConformalMedia_t& cM) {
    ASSERT_EQ(cM.edge_media.size(), 2u);
    expectRatioNear(cM.edge_media[0].ratio, 0.75);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.75);
    expectRatioNear(cM.edge_media[1].ratio, 0.0);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 2u);
    expectRatioNear(cM.edge_media[1].edges[0].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[1].ratio, 0.0);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.75);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.75);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.75);
}

void checkOffFaceTriangleXMinus(const ConformalMedia_t& cM) {
    ASSERT_EQ(cM.edge_media.size(), 2u);
    expectRatioNear(cM.edge_media[0].ratio, 0.25);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.25);
    expectRatioNear(cM.edge_media[1].ratio, 0.0);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 2u);
    expectRatioNear(cM.edge_media[1].edges[0].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[1].ratio, 0.0);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.25);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.25);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.25);
}

} // namespace

TEST(conformal, filling_off_face_triangle_x) {
    const coord_t c1 = makeCoord({0.75, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.75, 0.0, 1.0}, 2);
    const coord_t c3 = makeCoord({0.75, 1.0, 0.0}, 3);

    auto regions = makeRegion({makeTriangle(c1, c2, c3)});
    ConformalMedia_t cM = conformal_m::buildConformalMedia(regions)[0];
    checkOffFaceTriangleXPlus(cM);

    regions.volumes[0].triangles = {makeTriangle(c1, c3, c2)};
    cM = conformal_m::buildConformalMedia(regions)[0];
    checkOffFaceTriangleXMinus(cM);
}

TEST(conformal, filling_off_face_triangle_y) {
    const coord_t c1 = makeCoord({0.0, 0.75, 0.0}, 1);
    const coord_t c2 = makeCoord({0.0, 0.75, 1.0}, 2);
    const coord_t c3 = makeCoord({1.0, 0.75, 0.0}, 3);

    auto regions = makeRegion({makeTriangle(c1, c2, c3)});
    ConformalMedia_t cM = conformal_m::buildConformalMedia(regions)[0];

    ASSERT_EQ(cM.edge_media.size(), 2u);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.25);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 2u);
    expectRatioNear(cM.edge_media[1].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[0].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[1].ratio, 0.0);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.25);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.25);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.25);

    regions.volumes[0].triangles = {makeTriangle(c1, c3, c2)};
    cM = conformal_m::buildConformalMedia(regions)[0];

    ASSERT_EQ(cM.edge_media.size(), 2u);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.75);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 2u);
    expectRatioNear(cM.edge_media[1].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[0].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[1].ratio, 0.0);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.75);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.75);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.75);
}

TEST(conformal, filling_off_face_triangle_z) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.75}, 1);
    const coord_t c2 = makeCoord({0.0, 1.0, 0.75}, 2);
    const coord_t c3 = makeCoord({1.0, 0.0, 0.75}, 3);

    auto regions = makeRegion({makeTriangle(c1, c2, c3)});
    ConformalMedia_t cM = conformal_m::buildConformalMedia(regions)[0];

    ASSERT_EQ(cM.edge_media.size(), 2u);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.75);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.75);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 2u);
    expectRatioNear(cM.edge_media[1].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[0].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[1].ratio, 0.0);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.75);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.75);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.75);

    regions.volumes[0].triangles = {makeTriangle(c1, c3, c2)};
    cM = conformal_m::buildConformalMedia(regions)[0];

    ASSERT_EQ(cM.edge_media.size(), 2u);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.25);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.25);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 2u);
    expectRatioNear(cM.edge_media[1].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[0].ratio, 0.0);
    expectRatioNear(cM.edge_media[1].edges[1].ratio, 0.0);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.25);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.25);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.25);
}

TEST(conformal, filling_open) {
    const coord_t c1 = makeCoord({0.6, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.0, 0.6, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 0.6}, 3);

    const ConformalMedia_t cM =
        conformal_m::buildConformalMedia(makeRegion({makeTriangle(c1, c2, c3)}))[0];

    ASSERT_EQ(cM.edge_media.size(), 1u);
    expectRatioNear(cM.edge_media[0].ratio, 0.4);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.4);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.4);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.4);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.82);
    ASSERT_EQ(cM.face_media[0].faces.size(), 3u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.82);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.82);
    expectRatioNear(cM.face_media[0].faces[2].ratio, 0.82);
    expectRatioNear(cM.time_step_scale_factor, 0.9055);
}

TEST(conformal, filling_closed) {
    const coord_t c1 = makeCoord({0.6, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.0, 0.6, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 0.6}, 3);
    const coord_t c4 = makeCoord({0.0, 0.0, 0.0}, 4);

    const std::vector<triangle_t> tris = {
        makeTriangle(c1, c2, c3),
        makeTriangle(c2, c4, c3),
        makeTriangle(c1, c3, c4),
        makeTriangle(c1, c4, c2),
    };

    const ConformalMedia_t cM = conformal_m::buildConformalMedia(makeRegion(tris))[0];

    ASSERT_EQ(cM.edge_media.size(), 1u);
    expectRatioNear(cM.edge_media[0].ratio, 0.4);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 3u);
    expectRatioNear(cM.edge_media[0].edges[0].ratio, 0.4);
    expectRatioNear(cM.edge_media[0].edges[1].ratio, 0.4);
    expectRatioNear(cM.edge_media[0].edges[2].ratio, 0.4);
    ASSERT_EQ(cM.face_media.size(), 1u);
    expectRatioNear(cM.face_media[0].ratio, 0.82);
    ASSERT_EQ(cM.face_media[0].faces.size(), 3u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.82);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.82);
    expectRatioNear(cM.face_media[0].faces[2].ratio, 0.82);
    expectRatioNear(cM.time_step_scale_factor, 0.9055);
}

TEST(conformal, filling_edge_next_cell) {
    const coord_t c1 = makeCoord({0.0, 0.25, 0.0}, 1);
    const coord_t c2 = makeCoord({0.0, 1.0, 1.0}, 2);
    const coord_t c3 = makeCoord({0.25, 1.0, 1.0}, 3);
    const coord_t c4 = makeCoord({1.0, 1.0, 0.5}, 4);
    const coord_t c5 = makeCoord({1.0, 1.0, 0.0}, 5);
    const coord_t c6 = makeCoord({0.0, 1.0, 0.0}, 6);

    const std::vector<triangle_t> tris = {
        makeTriangle(c1, c3, c2),
        makeTriangle(c1, c4, c3),
        makeTriangle(c1, c5, c4),
        makeTriangle(c1, c2, c6),
        makeTriangle(c1, c6, c5),
        makeTriangle(c6, c2, c3),
        makeTriangle(c6, c3, c4),
        makeTriangle(c6, c4, c5),
    };

    const ConformalMedia_t cM = conformal_m::buildConformalMedia(makeRegion(tris))[0];

    ASSERT_EQ(cM.edge_media.size(), 4u);
    expectRatioNear(cM.edge_media[0].ratio, 0.25);
    expectRatioNear(cM.edge_media[1].ratio, 0.0);
    expectRatioNear(cM.edge_media[2].ratio, 0.75);
    expectRatioNear(cM.edge_media[3].ratio, 0.5);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 1u);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 2u);
    ASSERT_EQ(cM.edge_media[2].edges.size(), 1u);
    ASSERT_EQ(cM.edge_media[3].edges.size(), 1u);
    ASSERT_EQ(cM.face_media.size(), 2u);
    expectRatioNear(cM.face_media[0].ratio, 0.625);
    expectRatioNear(cM.face_media[1].ratio, 0.1875);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    ASSERT_EQ(cM.face_media[1].faces.size(), 1u);
}

TEST(conformal, filling_closed_corner) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({1.0, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 1.0}, 3);
    const coord_t c4 = makeCoord({1.0, 0.0, 1.0}, 4);
    const coord_t c5 = makeCoord({0.0, 1.0, 0.0}, 5);
    const coord_t c6 = makeCoord({0.0, 1.0, 1.0}, 6);

    const std::vector<triangle_t> tris = {
        makeTriangle(c1, c2, c3),
        makeTriangle(c2, c4, c3),
        makeTriangle(c1, c3, c6),
        makeTriangle(c1, c6, c5),
        makeTriangle(c3, c4, c6),
        makeTriangle(c1, c5, c2),
        makeTriangle(c2, c5, c4),
        makeTriangle(c5, c6, c4),
    };

    const ConformalMedia_t cM = conformal_m::buildConformalMedia(makeRegion(tris))[0];

    ASSERT_EQ(cM.edge_media.size(), 1u);
    expectRatioNear(cM.edge_media[0].ratio, 0.0);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 7u);
    for (std::size_t i = 0; i < 7; ++i) {
        expectRatioNear(cM.edge_media[0].edges[i].ratio, 0.0);
    }
    ASSERT_EQ(cM.face_media.size(), 2u);
    expectRatioNear(cM.face_media[0].ratio, 0.0);
    ASSERT_EQ(cM.face_media[0].faces.size(), 2u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.0);
    expectRatioNear(cM.face_media[0].faces[1].ratio, 0.0);
    expectRatioNear(cM.face_media[1].ratio, 0.5);
    ASSERT_EQ(cM.face_media[1].faces.size(), 2u);
    // Fortran uses 0.5 tolerance for per-face ratios in this group (see test_filling.F90)
    EXPECT_NEAR(cM.face_media[1].faces[0].ratio, 0.0, 0.5);
    EXPECT_NEAR(cM.face_media[1].faces[1].ratio, 0.0, 0.5);
}

TEST(conformal, filling_block_and_corner) {
    const coord_t c1 = makeCoord({0.0, 0.0, 0.0}, 1);
    const coord_t c2 = makeCoord({0.5, 0.0, 0.0}, 2);
    const coord_t c3 = makeCoord({0.0, 0.0, 1.0}, 3);
    const coord_t c4 = makeCoord({0.5, 0.0, 1.0}, 4);
    const coord_t c5 = makeCoord({0.0, 1.0, 1.0}, 5);
    const coord_t c6 = makeCoord({0.5, 1.0, 1.0}, 6);
    const coord_t c7 = makeCoord({0.0, 1.0, 0.0}, 7);
    const coord_t c8 = makeCoord({0.5, 1.0, 0.0}, 8);
    const coord_t c9 = makeCoord({0.0, 1.5, 1.0}, 9);
    const coord_t c10 = makeCoord({0.0, 1.5, 0.0}, 10);

    const std::vector<triangle_t> tris = {
        makeTriangle(c1, c2, c3),
        makeTriangle(c2, c4, c3),
        makeTriangle(c1, c3, c5),
        makeTriangle(c1, c5, c7),
        makeTriangle(c3, c4, c5),
        makeTriangle(c4, c6, c5),
        makeTriangle(c1, c7, c2),
        makeTriangle(c2, c7, c8),
        makeTriangle(c2, c8, c6),
        makeTriangle(c2, c6, c4),
        makeTriangle(c5, c6, c9),
        makeTriangle(c7, c10, c8),
        makeTriangle(c5, c9, c10),
        makeTriangle(c5, c10, c7),
        makeTriangle(c9, c6, c10),
        makeTriangle(c8, c10, c6),
    };

    const ConformalMedia_t cM = conformal_m::buildConformalMedia(makeRegion(tris))[0];

    ASSERT_EQ(cM.edge_media.size(), 2u);
    expectRatioNear(cM.edge_media[0].ratio, 0.0);
    ASSERT_EQ(cM.edge_media[0].edges.size(), 4u);
    for (std::size_t i = 0; i < 4; ++i) {
        expectRatioNear(cM.edge_media[0].edges[i].ratio, 0.0);
    }
    expectRatioNear(cM.edge_media[1].ratio, 0.5);
    ASSERT_EQ(cM.edge_media[1].edges.size(), 6u);
    for (std::size_t i = 0; i < 6; ++i) {
        expectRatioNear(cM.edge_media[1].edges[i].ratio, 0.5);
    }
    ASSERT_EQ(cM.face_media.size(), 3u);
    expectRatioNear(cM.face_media[0].ratio, 0.0);
    ASSERT_EQ(cM.face_media[0].faces.size(), 1u);
    expectRatioNear(cM.face_media[0].faces[0].ratio, 0.0);
    expectRatioNear(cM.face_media[1].ratio, 0.5);
    ASSERT_EQ(cM.face_media[1].faces.size(), 5u);
    for (std::size_t i = 0; i < 5; ++i) {
        expectRatioNear(cM.face_media[1].faces[i].ratio, 0.5);
    }
    expectRatioNear(cM.face_media[2].ratio, 0.875);
    ASSERT_EQ(cM.face_media[2].faces.size(), 2u);
    expectRatioNear(cM.face_media[2].faces[0].ratio, 0.875);
    expectRatioNear(cM.face_media[2].faces[1].ratio, 0.875);
}

#endif // TEST_CONFORMAL_FILLING_H
