#ifndef WIRES_TESTS_H
#define WIRES_TESTS_H

#include <gtest/gtest.h>

extern "C" {
    int test_getHwires_returns_associated();
    int test_destroyWires_deallocates_wires_data();
    int test_destroyWires_keeps_non_thinwire_media();
    int test_destroyWireMedia_thinwire();
    int test_destroyWireMedia_non_thinwire();
    int test_evolucion_out_of_range_low();
    int test_evolucion_out_of_range_high();
    int test_evolucion_at_zero();
    int test_evolucion_midpoint_interp();
    int test_evolucion_exact_sample();
}

TEST(wires, getHwires_returns_associated) {
    EXPECT_EQ(0, test_getHwires_returns_associated());
}

TEST(wires, destroyWires_deallocates_wires_data) {
    EXPECT_EQ(0, test_destroyWires_deallocates_wires_data());
}

TEST(wires, destroyWires_keeps_non_thinwire_media) {
    EXPECT_EQ(0, test_destroyWires_keeps_non_thinwire_media());
}

TEST(wires, destroyWireMedia_thinwire) {
    EXPECT_EQ(0, test_destroyWireMedia_thinwire());
}

TEST(wires, destroyWireMedia_non_thinwire) {
    EXPECT_EQ(0, test_destroyWireMedia_non_thinwire());
}

TEST(wires, evolucion_out_of_range_low) {
    EXPECT_EQ(0, test_evolucion_out_of_range_low());
}

TEST(wires, evolucion_out_of_range_high) {
    EXPECT_EQ(0, test_evolucion_out_of_range_high());
}

TEST(wires, evolucion_at_zero) {
    EXPECT_EQ(0, test_evolucion_at_zero());
}

TEST(wires, evolucion_midpoint_interp) {
    EXPECT_EQ(0, test_evolucion_midpoint_interp());
}

TEST(wires, evolucion_exact_sample) {
    EXPECT_EQ(0, test_evolucion_exact_sample());
}

#endif // WIRES_TESTS_H
