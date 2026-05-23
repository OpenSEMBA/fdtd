#ifndef TEST_MTLN_SPICE_H
#define TEST_MTLN_SPICE_H

#include <gtest/gtest.h>
#include <cstdlib>
#include <string>

#include "test_mtln_helpers.h"

namespace {

bool spiceEnvironmentReady() {
    const char* semba = std::getenv("SEMBA");
    return semba != nullptr && std::string(semba).size() > 0;
}

} // namespace

// ngspice circuit tests require the full circuit_m C++ port (src_mtln/circuit.cpp).

TEST(mtln, codemodels) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set (codemodel path)";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

TEST(mtln, spice_tran) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

TEST(mtln, spice_tran_2) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

TEST(mtln, spice_multiple) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

TEST(mtln, spice_current_source) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

TEST(mtln, spice_dc) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

TEST(mtln, spice_read_message) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

TEST(mtln, spice_mod_times) {
    if (!spiceEnvironmentReady()) {
        GTEST_SKIP() << "SEMBA environment variable not set";
    }
    GTEST_SKIP() << "circuit_m ngspice interface not yet ported in src_cpp";
}

#endif
