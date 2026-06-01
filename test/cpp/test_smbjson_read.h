#ifndef TEST_SMBJSON_READ_H
#define TEST_SMBJSON_READ_H

#include "test_smbjson_helpers.h"
#include "expected/test_read_planewave.h"
#include "expected/test_read_dielectricSlab.h"
#include "expected/test_read_thinSlot.h"
#include "expected/test_read_sgbc.h"
#include "expected/test_read_sphere.h"
#include "expected/test_read_airplane.h"
#include "expected/test_read_lumped_fixture.h"
#include "expected/test_read_currentInjection.h"
#include "expected/test_read_holland1981.h"
#include "expected/test_read_nodal_source_resistance.h"

TEST(smbjson_cpp, read_planewave) {
    auto expected = expectedReadPlanewave();
    auto actual = parseFile(testDataPath("planewave.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_dielectricslab) {
    auto expected = expectedReadDielectricSlab();
    auto actual = parseFile(testDataPath("dielectric_slab.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_thinslot) {
    auto expected = expectedReadThinSlot();
    auto actual = parseFile(testDataPath("thinSlot.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_sgbc) {
    auto expected = expectedReadSgbc();
    auto actual = parseFile(testDataPath("sgbc.fdtd.json"));
    expect_eq(expected, actual);
}

TEST(smbjson_cpp, read_sphere) {
    auto expected = expectedReadSphere();
    auto actual = parseFile(testDataPath("sphere.fdtd.json"));
    expect_eq(expected, actual, true);
}

TEST(smbjson_cpp, read_airplane) {
    auto expected = expectedReadAirplane();
    auto actual = parseFile(testDataPath("airplane.fdtd.json"));
    expect_eq(expected, actual, true);
}

TEST(smbjson_cpp, read_lumped_fixture) {
    auto expected = expectedReadLumpedFixture();
    auto actual = parseFile(testDataPath("lumped_fixture.fdtd.json"));
    expect_eq(expected, actual);
}

#ifndef CompileWithMTLN
TEST(smbjson_cpp, read_currentinjection) {
    auto expected = expectedReadCurrentInjection();
    auto actual = parseFile(testDataPath("currentInjection.fdtd.json"));
    expect_eq(expected, actual);
}
#endif

#ifndef CompileWithMTLN
TEST(smbjson_cpp, read_holland1981) {
    auto expected = expectedReadHolland1981();
    auto actual = parseFile(testDataPath("holland1981.fdtd.json"));
    expect_eq(expected, actual);
}
#endif

TEST(smbjson_cpp, read_nodal_source_resistance_per_meter) {
    smbjson::parser_t parser(nodalSourceResistanceJsonPath());
    auto problem = parser.readProblemDescription();
    expectCableResistancePerMeter(problem, EXPECTED_CABLE_RESISTANCE_PER_METER);
}

TEST(smbjson_cpp, read_nodal_source_total_resistance) {
    smbjson::parser_t parser(nodalSourceTotalResistanceJsonPath());
    auto problem = parser.readProblemDescription();
    expectTotalResistanceOverride(problem, EXPECTED_CABLE_RESISTANCE_PER_METER);
}

#endif // TEST_SMBJSON_READ_H
