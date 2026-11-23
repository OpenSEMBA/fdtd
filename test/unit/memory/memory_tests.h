#ifndef MEMORY_TESTS_H
#define MEMORY_TESTS_H

#include <gtest/gtest.h>

extern "C" {
    int test_memory_amount_format();
    int test_solver_memory_estimate();
    int test_system_memory_probe();
    int test_memory_failure_message();
    int test_cpml_estimate();
    int test_holland_wire_estimate();
    int test_conformal_estimate();
#ifndef CompileWithStochastic
    int test_sgbc_estimate();
#endif
}

TEST(memory, amount_format) {
    EXPECT_EQ(0, test_memory_amount_format());
}

TEST(memory, solver_estimate) {
    EXPECT_EQ(0, test_solver_memory_estimate());
}

TEST(memory, system_probe) {
    EXPECT_EQ(0, test_system_memory_probe());
}

TEST(memory, failure_message) {
    EXPECT_EQ(0, test_memory_failure_message());
}

TEST(memory, cpml_estimate) {
    EXPECT_EQ(0, test_cpml_estimate());
}

TEST(memory, holland_wire_estimate) {
    EXPECT_EQ(0, test_holland_wire_estimate());
}

TEST(memory, conformal_estimate) {
    EXPECT_EQ(0, test_conformal_estimate());
}

#ifndef CompileWithStochastic
TEST(memory, sgbc_estimate) {
    EXPECT_EQ(0, test_sgbc_estimate());
}
#endif

#endif // MEMORY_TESTS_H
