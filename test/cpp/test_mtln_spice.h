#ifndef TEST_MTLN_SPICE_H
#define TEST_MTLN_SPICE_H

#include <gtest/gtest.h>
#include <vector>

#ifndef _WIN32
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "circuit_m.h"
#include "test_mtln_helpers.h"

namespace {

using circuit_m::circuit_t;
using circuit_m::string_t;

std::vector<string_t> dcNodeNames() {
    return {
        string_t("node1", 5),
        string_t("node2", 5),
        string_t("node3", 5),
        string_t("v-sweep", 7),
    };
}

std::vector<string_t> tranNodeNames() {
    return {
        string_t("in", 2),
        string_t("int", 3),
        string_t("out", 3),
        string_t("time", 4),
    };
}

std::vector<string_t> multipleNodeNames() {
    return {
        string_t("n1_in", 5),
        string_t("n1_int", 6),
        string_t("n1_out", 6),
        string_t("time", 4),
        string_t("n2_in", 5),
        string_t("n2_int", 6),
        string_t("n2_out", 6),
    };
}

std::vector<string_t> codemodelNodeNames() {
    return {
        string_t("wire1_1_initial_R", 17),
        string_t("wire1_1_initial", 15),
        string_t("wire1_1_initial_S", 17),
        string_t("wire1_2_initial", 15),
        string_t("wire1_1_end", 11),
        string_t("wire1_wire1_inter", 17),
        string_t("wire1_2_end", 11),
    };
}

void expectDcVoltages(const circuit_t& circuit, int& error_cnt) {
    const std::vector<double> expected = {24.0, 9.7469741675197206, 15.0, 24.0};
    if (circuit.nodes.values.size() != expected.size()) {
        error_cnt += 1;
    }
    for (size_t i = 0; i < expected.size() && i < circuit.nodes.values.size(); ++i) {
        if (!mtln_test::checkNear(expected[i], circuit.nodes.values[i].voltage, 0.01)) {
            error_cnt += 1;
        }
    }
}

int runTranNetlist(const std::string& netlist_path, const std::vector<double>& expected) {
    int error_cnt = 0;
    circuit_t circuit;
    circuit.time = 0.0;
    circuit.dt = 50e-6;
    const double finalTime = 200e-6;
    circuit.init(tranNodeNames(), {}, netlist_path);
    circuit.setStopTimes(finalTime, circuit.dt);
    while (circuit.time < finalTime) {
        circuit.step();
        circuit.time += circuit.dt;
        if (!mtln_test::checkNearTime(circuit.getTime(), circuit.time, 0.01)) {
            error_cnt += 1;
        }
    }
    if (!mtln_test::checkNear(expected[0], circuit.getNodeVoltage("in"), 0.01)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(expected[1], circuit.getNodeVoltage("int"), 0.01)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(expected[2], circuit.getNodeVoltage("out"), 0.01)) {
        error_cnt += 1;
    }
    return error_cnt;
}

#ifndef _WIN32
bool runInIsolatedProcess(int (*test_fn)()) {
    const pid_t pid = fork();
    if (pid == 0) {
        _exit(test_fn() == 0 ? 0 : 1);
    }
    if (pid < 0) {
        return false;
    }
    int status = 0;
    waitpid(pid, &status, 0);
    return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}
#endif

int spiceTranBody() {
    const std::vector<double> expected = {5.0, 0.092995181699999999, 0.053166680000000001};
    return runTranNetlist(mtln_test::PATH_TO_TEST_DATA + "netlists/netlist_tran.cir", expected);
}

int spiceTran2Body() {
    const std::vector<double> expected = {5.0, 0.0039656539400000001, 0.00069279532199999997};
    return runTranNetlist(mtln_test::PATH_TO_TEST_DATA + "netlists/netlist_tran_2.cir", expected);
}

} // namespace

TEST(mtln, codemodels) {
    circuit_t circuit;
    circuit.init(codemodelNodeNames(),
                 {},
                 mtln_test::PATH_TO_TEST_DATA + "netlists/saturation.cir");
}

TEST(mtln, spice_tran_2) {
#ifndef _WIN32
    EXPECT_TRUE(runInIsolatedProcess(spiceTran2Body));
#else
    EXPECT_EQ(0, spiceTran2Body());
#endif
}

TEST(mtln, spice_tran) {
#ifndef _WIN32
    EXPECT_TRUE(runInIsolatedProcess(spiceTranBody));
#else
    EXPECT_EQ(0, spiceTranBody());
#endif
}

TEST(mtln, spice_multiple) {
    int error_cnt = 0;
    circuit_t circuit;
    circuit.time = 0.0;
    circuit.dt = 50e-6;
    const double finalTime = 200e-6;
    circuit.init(multipleNodeNames(), {}, mtln_test::PATH_TO_TEST_DATA + "netlists/netlist_multiple.cir");
    circuit.setStopTimes(finalTime, circuit.dt);
    while (circuit.time < finalTime) {
        circuit.step();
        circuit.time += circuit.dt;
        if (!mtln_test::checkNearTime(circuit.getTime(), circuit.time, 0.01)) {
            error_cnt += 1;
        }
    }
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, spice_current_source) {
    int error_cnt = 0;
    const double resistance = 10.0;
    circuit_t circuit;
    circuit.time = 0.0;
    circuit.dt = 50e-6;
    const double finalTime = 200e-6;
    circuit.init({string_t("1_initial", 9)},
                 {},
                 mtln_test::PATH_TO_TEST_DATA + "netlists/netlist_current_source.cir");
    circuit.setStopTimes(finalTime, circuit.dt);
    double current = 0.1;
    while (circuit.time < finalTime) {
        circuit.updateNodeCurrent("1_initial", current);
        circuit.step();
        circuit.time += circuit.dt;
        if (!mtln_test::checkNear(circuit.getNodeVoltage("1_initial"), current * resistance, 0.01)) {
            error_cnt += 1;
        }
        current = 2.0 * current;
    }
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, spice_dc) {
    int error_cnt = 0;
    circuit_t circuit;
    circuit.init(dcNodeNames(), {}, mtln_test::PATH_TO_TEST_DATA + "netlists/netlist_dc.cir");
    circuit.step();
    expectDcVoltages(circuit, error_cnt);
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, spice_read_message) {
    int error_cnt = 0;
    circuit_t circuit;
    const std::vector<std::string> input = {
        "* Multiple dc sources",
        "vn1 node1 0 dc 24",
        "vn2 node3 0 dc 15",
        "rn1 node1 node2 10k",
        "rn2 node2 node3 8.1k",
        "rn3 node2 0 4.7k",
        ".dc vn1 24 24 1",
        ".save v(node3) v(node2) v(node1)",
        ".end",
        "NULL",
    };
    circuit.init(dcNodeNames());
    circuit.readInput(input);
    circuit.step();
    expectDcVoltages(circuit, error_cnt);
    EXPECT_EQ(error_cnt, 0);
}

TEST(mtln, spice_mod_times) {
    int error_cnt = 0;
    const std::vector<double> expected = {5.0, 0.092995181699999999, 0.053166680000000001};
    circuit_t circuit;
    circuit.time = 0.0;
    circuit.dt = 50e-6;
    const double finalTime = 200e-6;
    circuit.init(tranNodeNames(), {}, mtln_test::PATH_TO_TEST_DATA + "netlists/netlist_tran.cir");
    circuit.setModStopTimes(circuit.dt);
    while (circuit.time < finalTime) {
        circuit.step();
        circuit.time += circuit.dt;
        if (!mtln_test::checkNearTime(circuit.getTime(), circuit.time, 0.01)) {
            error_cnt += 1;
        }
    }
    if (!mtln_test::checkNear(expected[0], circuit.getNodeVoltage("in"), 0.01)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(expected[1], circuit.getNodeVoltage("int"), 0.01)) {
        error_cnt += 1;
    }
    if (!mtln_test::checkNear(expected[2], circuit.getNodeVoltage("out"), 0.01)) {
        error_cnt += 1;
    }
    EXPECT_EQ(error_cnt, 0);
}

#endif
