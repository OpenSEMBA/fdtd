#ifndef TEST_SMBJSON_HELPERS_H
#define TEST_SMBJSON_HELPERS_H

#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include "smbjson_m.h"
#include "NFDETypes_extension_m.h"

using namespace NFDETypes_m;

#ifdef CompileWithMTLN
using mtln_types_m::cable_t;
using mtln_types_m::shielded_multiwire_t;
using mtln_types_m::unshielded_multiwire_t;

inline void initializeCablePULParameters(cable_t* cable, int dim = 1) {
    if (auto* sh = dynamic_cast<shielded_multiwire_t*>(cable)) {
        sh->inductance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
        sh->capacitance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
        sh->resistance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
        sh->conductance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
    } else if (auto* un = dynamic_cast<unshielded_multiwire_t*>(cable)) {
        un->cell_inductance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
        un->cell_capacitance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
        un->resistance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
        un->conductance_per_meter.assign(dim, std::vector<double>(dim, 0.0));
        un->multipolar_expansion.clear();
    }
}
#endif

static const std::string TEST_DATA_DIR = "testData/";
static const std::string INPUT_EXAMPLES = "input_examples/";

inline std::string testDataPath(const std::string& filename) {
    return TEST_DATA_DIR + INPUT_EXAMPLES + filename;
}

inline std::string testDataPathCases(const std::string& filename) {
    return TEST_DATA_DIR + "cases/" + filename;
}

inline Parseador_t parseFile(const std::string& jsonFile) {
    smbjson::parser_t parser(jsonFile);
    return parser.readProblemDescription();
}

inline Parseador_t buildExpected() {
    Parseador_t expected;
    NFDETypes_extension_m::initializeProblemDescription(expected);
    return expected;
}

inline bool expect_near(double a, double b, double tol) {
    return std::abs(a - b) < tol;
}

inline void expect_eq(const Parseador_t& expected, const Parseador_t& actual,
                      bool ignoreRegions = false) {
    if (expected.switches != actual.switches)
        ADD_FAILURE() << "Expected and read switches do not match";
    if (!(*expected.general == *actual.general))
        ADD_FAILURE() << "Expected and read \"general\" do not match";
    if (!(*expected.matriz == *actual.matriz))
        ADD_FAILURE() << "Expected and read \"media matrix\" do not match";
    if (!(*expected.despl == *actual.despl))
        ADD_FAILURE() << "Expected and read \"grid\" do not match";
    if (!(*expected.front == *actual.front))
        ADD_FAILURE() << "Expected and read \"boundary\" do not match";
    if (!(*expected.Mats == *actual.Mats))
        ADD_FAILURE() << "Expected and read \"materials\" do not match";

    if (!ignoreRegions) {
        if (!(*expected.pecRegs == *actual.pecRegs))
            ADD_FAILURE() << "Expected and read \"pec regions\" do not match";
        if (!(*expected.pmcRegs == *actual.pmcRegs))
            ADD_FAILURE() << "Expected and read \"pmc regions\" do not match";
        if (!(*expected.DielRegs == *actual.DielRegs))
            ADD_FAILURE() << "Expected and read \"dielectric regions\" do not match";
        if (!(*expected.LossyThinSurfs == *actual.LossyThinSurfs))
            ADD_FAILURE() << "Expected and read \"lossy thin surfaces\" do not match";
    }

    if (!(*expected.plnSrc == *actual.plnSrc))
        ADD_FAILURE() << "Expected and read \"planewave sources\" do not match";
    if (!(*expected.nodSrc == *actual.nodSrc))
        ADD_FAILURE() << "Expected and read \"nodal sources\" do not match";
    if (!(*expected.oldSONDA == *actual.oldSONDA))
        ADD_FAILURE() << "Expected and read \"old probes\" do not match";
    if (!(*expected.Sonda == *actual.Sonda))
        ADD_FAILURE() << "Expected and read \"new probes\" do not match";
    if (!(*expected.BloquePrb == *actual.BloquePrb))
        ADD_FAILURE() << "Expected and read \"block probes\" do not match";
    if (!(*expected.VolPrb == *actual.VolPrb))
        ADD_FAILURE() << "Expected and read \"vol probes\" do not match";
    if (!(*expected.tWires == *actual.tWires))
        ADD_FAILURE() << "Expected and read \"thin wires\" do not match";
    if (!(*expected.tSlots == *actual.tSlots))
        ADD_FAILURE() << "Expected and read \"thin slots\" do not match";
#ifdef CompileWithMTLN
    if (expected.mtln && actual.mtln) {
        for (size_t i = 0; i < expected.mtln->cables.size() && i < actual.mtln->cables.size(); ++i) {
            if (!expected.mtln->cables[i].ptr || !actual.mtln->cables[i].ptr)
                continue;
            auto& es = expected.mtln->cables[i].ptr->segments;
            const auto& as = actual.mtln->cables[i].ptr->segments;
            if (es.size() != as.size())
                continue;
            for (size_t j = 0; j < es.size(); ++j) {
                es[j].dualBox = as[j].dualBox;
                es[j].d1 = as[j].d1;
                es[j].d2 = as[j].d2;
            }
        }
        if (!(*expected.mtln == *actual.mtln))
            ADD_FAILURE() << "Expected and read \"mtln\" do not match";
    } else if (expected.mtln || actual.mtln) {
        ADD_FAILURE() << "Expected and read \"mtln\" presence do not match";
    }
#endif
}

#endif // TEST_SMBJSON_HELPERS_H
