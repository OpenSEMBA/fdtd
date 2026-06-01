#ifndef TEST_SMBJSON_READ_MTLN_H
#define TEST_SMBJSON_READ_MTLN_H

#ifdef CompileWithMTLN

#include "test_smbjson_helpers.h"
#include "expected/test_read_connectedWires.h"
#include "expected/test_read_shieldedPair.h"
#include "expected/test_read_mtln.h"
#include "expected/test_read_towelHanger.h"
#include "expected/test_read_holland1981_unshielded.h"
#include "expected/test_read_unshielded_multiwires_multipolar_expansion.h"
#include "expected/test_read_currentInjection.h"

TEST(smbjson_cpp, read_connectedwires) {
    auto expected = expectedReadConnectedWires();
    auto actual = parseFile(testDataPath("connectedWires.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_shieldedpair) {
    auto expected = expectedReadShieldedPair();
    auto actual = parseFile(testDataPath("shieldedPair.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_mtln) {
    auto expected = expectedReadMtln();
    auto actual = parseFile(testDataPath("mtln.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_towelhanger) {
    auto expected = expectedReadTowelHanger();
    auto actual = parseFile(testDataPath("towelHanger.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_holland1981_unshielded) {
    auto expected = expectedReadHolland1981Unshielded();
    auto actual = parseFile(testDataPath("holland1981_unshielded.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_unshielded_multiwires_multipolar_expansion) {
    auto expected = expectedReadUnshieldedMultiwiresMultipolarExpansion();
    auto actual = parseFile(testDataPath("unshielded_multiwires_multipolar_expansion.fdtd.json"));
    expect_eq(expected, actual, true);
}

TEST(smbjson_cpp, read_currentinjection_mtln) {
    auto expected = expectedReadCurrentInjection();
    auto actual = parseFile(testDataPath("currentInjection.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_large_airplane_mtln) {
    EXPECT_NO_THROW({
        smbjson::parser_t parser(testDataPath("large_airplane_mtln.fdtd.json"));
        parser.readProblemDescription();
    });
}

// Endpoint voltage generators (e.g. paul start node) must not become TLM wireGenerators;
// excitation is applied via Spice termination only (Fortran IsGeneratorOnWire interior-only).
TEST(smbjson_cpp, read_paul_8_6_square_no_endpoint_wire_generators) {
    auto actual = parseFile(testDataPathCases("paul/paul_8_6_square.fdtd.json"));
    ASSERT_NE(actual.mtln, nullptr);
    EXPECT_TRUE(actual.mtln->wireGenerators.empty());
    ASSERT_FALSE(actual.mtln->networks.empty());
    bool has_excitation_on_network = false;
    for (const auto& net : actual.mtln->networks) {
        for (const auto& conn : net.connections) {
            for (const auto& node : conn.nodes) {
                if (!node.termination.source.path_to_excitation.empty()) {
                    has_excitation_on_network = true;
                }
            }
        }
    }
    EXPECT_TRUE(has_excitation_on_network);
}

#endif // CompileWithMTLN

#endif // TEST_SMBJSON_READ_MTLN_H
