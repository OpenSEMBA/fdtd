#ifndef TEST_CONFORMAL_HELPERS_H
#define TEST_CONFORMAL_HELPERS_H

#include <gtest/gtest.h>
#include <initializer_list>
#include <vector>

#include "conformal_types.h"

namespace conformal_test {

using conformal_types_m::coord_t;
using conformal_types_m::side_t;
using conformal_types_m::triangle_t;

inline coord_t makeCoord(std::initializer_list<double> xyz, int id) {
    return coord_t(std::vector<double>(xyz.begin(), xyz.end()), id);
}

inline side_t makeSide(const coord_t& init, const coord_t& end) {
    side_t side;
    side.init = init;
    side.end = end;
    return side;
}

inline triangle_t makeTriangle(const coord_t& c1, const coord_t& c2, const coord_t& c3) {
    triangle_t tri;
    tri.vertices = {c1, c2, c3};
    return tri;
}

inline bool positionsEqual(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == 3u && b.size() == 3u &&
           a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

inline void expectPositionsEqual(const coord_t& a, const coord_t& b) {
    EXPECT_TRUE(positionsEqual(a.position, b.position));
}

inline void expectSideEndpointPositions(const side_t& side,
                                        const coord_t& init,
                                        const coord_t& end) {
    expectPositionsEqual(side.init, init);
    expectPositionsEqual(side.end, end);
}

inline void expectCellEqual(const std::vector<int>& actual,
                            const std::vector<int>& expected) {
    ASSERT_EQ(actual.size(), 3u);
    ASSERT_EQ(expected.size(), 3u);
    EXPECT_EQ(actual[0], expected[0]);
    EXPECT_EQ(actual[1], expected[1]);
    EXPECT_EQ(actual[2], expected[2]);
}

} // namespace conformal_test

#endif // TEST_CONFORMAL_HELPERS_H
