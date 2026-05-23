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

#endif // CompileWithMTLN

#endif // TEST_SMBJSON_READ_MTLN_H
