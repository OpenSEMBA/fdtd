#ifndef TEST_SMBJSON_HELPERS_H
#define TEST_SMBJSON_HELPERS_H

#include <gtest/gtest.h>
#include "smbjson_m.h"
#include "NFDETypes_extension_m.h"

using namespace NFDETypes_m;

// Path constant matching Fortran tests
static const std::string TEST_DATA_DIR = "testData/";
static const std::string INPUT_EXAMPLES = "input_examples/";

inline std::string testDataPath(const std::string& filename) {
    return TEST_DATA_DIR + INPUT_EXAMPLES + filename;
}

// Parse a JSON file and return the Parseador_t
inline Parseador_t parseFile(const std::string& jsonFile) {
    smbjson::parser_t parser(jsonFile);
    return parser.readProblemDescription();
}

// Build expected Parseador_t with defaults, then caller sets specific fields
inline Parseador_t buildExpected() {
    Parseador_t expected;
    NFDETypes_extension_m::initializeProblemDescription(expected);
    return expected;
}

// Deep comparison with GTest failure reporting
inline void expect_eq(const Parseador_t& expected, const Parseador_t& actual,
                      bool ignoreRegions = false) {
    if (expected.switches != actual.switches)
        ADD_FAILURE() << "switches do not match";
    if (!(*expected.general == *actual.general))
        ADD_FAILURE() << "general do not match";
    if (!(*expected.matriz == *actual.matriz))
        ADD_FAILURE() << "media matrix do not match";
    if (!(*expected.despl == *actual.despl))
        ADD_FAILURE() << "grid do not match";
    if (!(*expected.front == *actual.front))
        ADD_FAILURE() << "boundary do not match";
    if (!(*expected.Mats == *actual.Mats))
        ADD_FAILURE() << "materials do not match";

    if (!ignoreRegions) {
        if (!(*expected.pecRegs == *actual.pecRegs))
            ADD_FAILURE() << "pec regions do not match";
        if (!(*expected.pmcRegs == *actual.pmcRegs))
            ADD_FAILURE() << "pmc regions do not match";
        if (!(*expected.DielRegs == *actual.DielRegs))
            ADD_FAILURE() << "dielectric regions do not match";
        if (!(*expected.LossyThinSurfs == *actual.LossyThinSurfs))
            ADD_FAILURE() << "lossy thin surfaces do not match";
    }

    if (!(*expected.plnSrc == *actual.plnSrc))
        ADD_FAILURE() << "planewave sources do not match";
    if (!(*expected.nodSrc == *actual.nodSrc))
        ADD_FAILURE() << "nodal sources do not match";
    if (!(*expected.oldSONDA == *actual.oldSONDA))
        ADD_FAILURE() << "old probes do not match";
    if (!(*expected.Sonda == *actual.Sonda))
        ADD_FAILURE() << "new probes do not match";
    if (!(*expected.BloquePrb == *actual.BloquePrb))
        ADD_FAILURE() << "block probes do not match";
    if (!(*expected.VolPrb == *actual.VolPrb))
        ADD_FAILURE() << "vol probes do not match";
    if (!(*expected.tWires == *actual.tWires))
        ADD_FAILURE() << "thin wires do not match";
    if (!(*expected.tSlots == *actual.tSlots))
        ADD_FAILURE() << "thin slots do not match";
}

#endif // TEST_SMBJSON_HELPERS_H
