#ifndef TEST_CELLS_H
#define TEST_CELLS_H

#include <gtest/gtest.h>
#include "cells_m.h"

// Test cell_interval_t for linel in +X direction
TEST(cells, linel_positive_x) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 1; interval.ini.cell[1] = 1; interval.ini.cell[2] = 1;
    interval.end.cell[0] = 3; interval.end.cell[1] = 1; interval.end.cell[2] = 1;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_LINEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_X);
    EXPECT_EQ(interval.getSize(), 2);
}

// Test cell_interval_t for linel in -X direction
TEST(cells, linel_negative_x) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 2; interval.ini.cell[1] = 1; interval.ini.cell[2] = 1;
    interval.end.cell[0] = 0; interval.end.cell[1] = 1; interval.end.cell[2] = 1;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_LINEL);
    EXPECT_EQ(interval.getOrientation(), -cells_m::DIR_X);
    EXPECT_EQ(interval.getSize(), 2);
}

// Test cell_interval_t for linel in +Y direction
TEST(cells, linel_positive_y) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 11; interval.ini.cell[1] = 16; interval.ini.cell[2] = 11;
    interval.end.cell[0] = 11; interval.end.cell[1] = 21; interval.end.cell[2] = 11;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_LINEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_Y);
    EXPECT_EQ(interval.getSize(), 5);
}

// Test cell_interval_t for linel in +Z direction
TEST(cells, linel_positive_z) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 0; interval.ini.cell[1] = 0; interval.ini.cell[2] = 2;
    interval.end.cell[0] = 0; interval.end.cell[1] = 0; interval.end.cell[2] = 5;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_LINEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_Z);
    EXPECT_EQ(interval.getSize(), 3);
}

// Test cell_interval_t for linel in -Z direction
TEST(cells, linel_negative_z) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 0; interval.ini.cell[1] = 0; interval.ini.cell[2] = 5;
    interval.end.cell[0] = 0; interval.end.cell[1] = 0; interval.end.cell[2] = 2;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_LINEL);
    EXPECT_EQ(interval.getOrientation(), -cells_m::DIR_Z);
    EXPECT_EQ(interval.getSize(), 3);
}

// Test cell_interval_t for surfel +Z (2x2)
TEST(cells, surfel_positive_z) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 3; interval.ini.cell[1] = 3; interval.ini.cell[2] = 3;
    interval.end.cell[0] = 5; interval.end.cell[1] = 5; interval.end.cell[2] = 3;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_SURFEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_Z);
    EXPECT_EQ(interval.getSize(), 4);
}

// Test cell_interval_t for surfel -X (3x1)
TEST(cells, surfel_negative_x) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 5; interval.ini.cell[1] = 5; interval.ini.cell[2] = 5;
    interval.end.cell[0] = 5; interval.end.cell[1] = 4; interval.end.cell[2] = 2;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_SURFEL);
    EXPECT_EQ(interval.getOrientation(), -cells_m::DIR_X);
    EXPECT_EQ(interval.getSize(), 3);
}

// Test cell_interval_t for surfel -Y (XZ face, varying Y)
TEST(cells, surfel_negative_y) {
    cells_m::cell_interval_t interval;
    // X and Z vary, Y is fixed -> surface normal is along Y
    interval.ini.cell[0] = 1; interval.ini.cell[1] = 1; interval.ini.cell[2] = 1;
    interval.end.cell[0] = 3; interval.end.cell[1] = 1; interval.end.cell[2] = 3;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_SURFEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_Y);
    EXPECT_EQ(interval.getSize(), 4);
}

// Test cell_interval_t for voxel (3D)
TEST(cells, voxel) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 0; interval.ini.cell[1] = 0; interval.ini.cell[2] = 0;
    interval.end.cell[0] = 2; interval.end.cell[1] = 3; interval.end.cell[2] = 4;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_VOXEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_NULL);
    EXPECT_EQ(interval.getSize(), 2 * 3 * 4);
}

// Test cell_interval_t for single cell (pixel)
TEST(cells, pixel) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 5; interval.ini.cell[1] = 5; interval.ini.cell[2] = 5;
    interval.end.cell[0] = 5; interval.end.cell[1] = 5; interval.end.cell[2] = 5;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_PIXEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_NULL);
    EXPECT_EQ(interval.getSize(), 0);
}

// Test operator== for pixel_t
TEST(cells, pixel_equality) {
    cells_m::pixel_t a;
    a.cell[0] = 1; a.cell[1] = 2; a.cell[2] = 3; a.tag = 10;
    cells_m::pixel_t b;
    b.cell[0] = 1; b.cell[1] = 2; b.cell[2] = 3; b.tag = 10;
    cells_m::pixel_t c;
    c.cell[0] = 1; c.cell[1] = 2; c.cell[2] = 3; c.tag = 20;
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

// Test operator== for linel_t
TEST(cells, linel_equality) {
    cells_m::linel_t a;
    a.cell[0] = 1; a.cell[1] = 2; a.cell[2] = 3; a.tag = 10; a.orientation = cells_m::DIR_X;
    cells_m::linel_t b;
    b.cell[0] = 1; b.cell[1] = 2; b.cell[2] = 3; b.tag = 10; b.orientation = cells_m::DIR_X;
    cells_m::linel_t c;
    c.cell[0] = 1; c.cell[1] = 2; c.cell[2] = 3; c.tag = 10; c.orientation = cells_m::DIR_Y;
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

// Test cell_region_t::getIntervalsOfType
TEST(cells, region_intervals_of_type) {
    cells_m::cell_region_t region;
    // Add a linel interval
    cells_m::cell_interval_t linel;
    linel.ini.cell[0] = 0; linel.ini.cell[1] = 0; linel.ini.cell[2] = 0;
    linel.end.cell[0] = 2; linel.end.cell[1] = 0; linel.end.cell[2] = 0;
    region.intervals.push_back(linel);
    // Add a surfel interval
    cells_m::cell_interval_t surfel;
    surfel.ini.cell[0] = 0; surfel.ini.cell[1] = 0; surfel.ini.cell[2] = 0;
    surfel.end.cell[0] = 2; surfel.end.cell[1] = 2; surfel.end.cell[2] = 0;
    region.intervals.push_back(surfel);
    // Add another linel interval
    cells_m::cell_interval_t linel2;
    linel2.ini.cell[0] = 0; linel2.ini.cell[1] = 0; linel2.ini.cell[2] = 0;
    linel2.end.cell[0] = 0; linel2.end.cell[1] = 3; linel2.end.cell[2] = 0;
    region.intervals.push_back(linel2);

    auto linels = region.getIntervalsOfType(cells_m::CELL_TYPE_LINEL);
    EXPECT_EQ(linels.size(), 2u);
    auto surfels = region.getIntervalsOfType(cells_m::CELL_TYPE_SURFEL);
    EXPECT_EQ(surfels.size(), 1u);
    auto pixels = region.getIntervalsOfType(cells_m::CELL_TYPE_PIXEL);
    EXPECT_EQ(pixels.size(), 0u);
}

// Test cell_region_t::toPixels
TEST(cells, region_to_pixels) {
    cells_m::cell_region_t region;
    // Add two pixel intervals
    for (int x = 0; x < 3; x++) {
        cells_m::cell_interval_t pixel;
        pixel.ini.cell[0] = x; pixel.ini.cell[1] = 0; pixel.ini.cell[2] = 0;
        pixel.end.cell[0] = x; pixel.end.cell[1] = 0; pixel.end.cell[2] = 0;
        region.intervals.push_back(pixel);
    }
    // Add a linel (should not appear in toPixels result)
    cells_m::cell_interval_t linel;
    linel.ini.cell[0] = 0; linel.ini.cell[1] = 2; linel.ini.cell[2] = 0;
    linel.end.cell[0] = 0; linel.end.cell[1] = 5; linel.end.cell[2] = 0;
    region.intervals.push_back(linel);

    auto pixels = region.toPixels();
    EXPECT_EQ(pixels.size(), 3u);
    for (int x = 0; x < 3; x++) {
        EXPECT_EQ(pixels[x].cell[0], x);
        EXPECT_EQ(pixels[x].cell[1], 0);
        EXPECT_EQ(pixels[x].cell[2], 0);
    }
}

// Test linel of size 1
TEST(cells, linel_single_step) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 5; interval.ini.cell[1] = 5; interval.ini.cell[2] = 5;
    interval.end.cell[0] = 6; interval.end.cell[1] = 5; interval.end.cell[2] = 5;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_LINEL);
    EXPECT_EQ(interval.getOrientation(), cells_m::DIR_X);
    EXPECT_EQ(interval.getSize(), 1);
}

// Test surfel 1x1
TEST(cells, surfel_single_face) {
    cells_m::cell_interval_t interval;
    interval.ini.cell[0] = 2; interval.ini.cell[1] = 3; interval.ini.cell[2] = 4;
    interval.end.cell[0] = 3; interval.end.cell[1] = 4; interval.end.cell[2] = 4;
    EXPECT_EQ(interval.getType(), cells_m::CELL_TYPE_SURFEL);
    EXPECT_EQ(interval.getSize(), 1);
}

#endif // TEST_CELLS_H
