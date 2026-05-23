#include "circuit_m.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "ngspice/sharedspice.h"
#include "ngspice_interface_m.h"
#include "Report_m.h"

namespace circuit_m {

namespace {

constexpr char kNullChar = '\0';

void sendCommand(const std::string& cmd) {
    command(cmd.c_str());
}

} // namespace

RKIND source_t::interpolate(RKIND_TIEMPO eval_time, RKIND_TIEMPO dt_step) const {
    const int n = static_cast<int>(time.size());
    if (n == 0) {
        return 0.0;
    }

    const RKIND_TIEMPO t_eval = eval_time - dt_step;
    if (t_eval <= time[0]) {
        return value[0];
    }
    if (t_eval >= time[n - 1]) {
        return value[n - 1];
    }

    int index = 0;
    RKIND_TIEMPO max_timediff = -std::numeric_limits<RKIND_TIEMPO>::infinity();
    for (int i = 0; i < n; ++i) {
        const RKIND_TIEMPO timediff = time[static_cast<size_t>(i)] - t_eval;
        if (timediff <= 0.0 && timediff >= max_timediff) {
            max_timediff = timediff;
            index = i + 1;
        }
    }
    if (index == 0) {
        index = 1;
    }
    if (index >= n) {
        index = n - 1;
    }

    const int i1 = index - 1;
    const int i2 = i1 + 1;
    const RKIND x1 = time[static_cast<size_t>(i1)];
    const RKIND y1 = value[static_cast<size_t>(i1)];
    const RKIND x2 = time[static_cast<size_t>(i2)];
    const RKIND y2 = value[static_cast<size_t>(i2)];

    if (x2 == x1) {
        return y2;
    }
    return (t_eval * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1);
}

source_t setSource(const std::string& source_path) {
    source_t res;
    if (source_path.empty()) {
        res.has_source = false;
        return res;
    }

    std::ifstream in(source_path);
    if (!in) {
        Report_m::WarnErrReport("Cannot open excitation file: " + source_path, true);
        res.has_source = false;
        return res;
    }

    std::vector<RKIND_TIEMPO> times;
    std::vector<RKIND> values;
    RKIND_TIEMPO t = 0.0;
    RKIND v = 0.0;
    while (in >> t >> v) {
        times.push_back(t);
        values.push_back(v);
    }

    if (times.empty()) {
        Report_m::WarnErrReport("Error reading excitation file: " + source_path, true);
        res.has_source = false;
        return res;
    }

    res.has_source = true;
    res.time = std::move(times);
    res.value = std::move(values);
    return res;
}

void circuit_t::printCWD() {
    sendCommand("getcwd");
}

void circuit_t::init(const std::vector<string_t>& names, const std::vector<node_source_t>& sources,
                     const std::string& netlist) {
    start();

    if (!netlist.empty()) {
        loadNetlist(netlist);
    }

    if (names.empty()) {
        throw std::runtime_error("Missing node names");
    }

    nodes.names = names;
    nodes.values.assign(names.size(), VI_t{});
    nodes.sources.assign(names.size(), source_t{});

    if (!sources.empty()) {
        for (size_t i = 0; i < sources.size(); ++i) {
            nodes.sources[i] = setSource(sources[i].path_to_excitation);
            nodes.sources[i].source_type = sources[i].source_type;
        }
    }
}

void circuit_t::loadNetlist(const std::string& netlist) {
    sendCommand("source " + netlist);
    sendCommand("save none");
    sendCommand("set history = 1");
}

void circuit_t::updateCircuitSources(RKIND_TIEMPO eval_time) {
    for (size_t i = 0; i < nodes.sources.size(); ++i) {
        if (!nodes.sources[i].has_source) {
            continue;
        }
        const RKIND interp = nodes.sources[i].interpolate(eval_time, 0.0);
        std::ostringstream oss;
        oss << interp;
        const std::string source_value = oss.str();
        if (nodes.sources[i].source_type == SOURCE_TYPE_VOLTAGE) {
            sendCommand("alter @V" + nodes.names[i].name + "_s[dc] = " + source_value);
        } else if (nodes.sources[i].source_type == SOURCE_TYPE_CURRENT) {
            sendCommand("alter @I" + nodes.names[i].name + "_s[dc] = " + source_value);
        }
    }
}

void circuit_t::run() {
    sendCommand("run ");
}

void circuit_t::resume() {
    sendCommand("resume ");
}

void circuit_t::quit() {
    sendCommand("quit 0");
}

void circuit_t::setStopTimes(RKIND_TIEMPO finalTime, RKIND_TIEMPO dt_in) {
    sendCommand("delete all");
    RKIND_TIEMPO t = 0.0;
    while (t < finalTime) {
        t += dt_in;
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "% .16g", static_cast<double>(t));
        sendCommand(std::string("stop when time = ") + buffer);
    }
}

void circuit_t::setModStopTimes(RKIND_TIEMPO dt_in) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%g", static_cast<float>(dt_in));
    sendCommand(std::string("stop when time mod ") + buffer);
}

void circuit_t::readInput(const std::vector<std::string>& input, bool printInput) {
    if (printInput) {
        for (const auto& line : input) {
            std::cout << line << '\n';
        }
    }

    std::vector<std::string> tmp(input.size());
    std::vector<char*> argv_c(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        tmp[i] = input[i];
        tmp[i].push_back(kNullChar);
        argv_c[i] = tmp[i].data();
    }
    circ(argv_c.data());

    sendCommand("save none");
    sendCommand("set history = 1");
}

void circuit_t::updateNodeCurrent(const std::string& node_name, RKIND current) {
    std::string sCurrent;
    if (node_name.find("initial") != std::string::npos) {
        std::ostringstream oss;
        oss << current;
        sCurrent = oss.str();
    } else if (node_name.find("end") != std::string::npos) {
        std::ostringstream oss;
        oss << -current;
        sCurrent = oss.str();
    } else {
        sCurrent = "0";
    }
    sendCommand("alter @I" + node_name + "[dc] = " + sCurrent);
}

void circuit_t::modifyLineCapacitorValue(const std::string& name, RKIND c) {
    std::ostringstream oss;
    oss << c;
    sendCommand("alter @CL" + name + " = " + oss.str());
}

void circuit_t::updateNodes() {
    if (has_error() != 0) {
        Report_m::WarnErrReport("Ngspice reported a controlled exit while updating nodes.", true);
        return;
    }

    for (size_t i = 0; i < nodes.names.size(); ++i) {
        const std::string name = nodes.names[i].name.substr(0, static_cast<size_t>(nodes.names[i].length));
        void* info_ptr = get_vector_info(const_cast<char*>(name.c_str()));
        if (info_ptr == nullptr) {
            Report_m::WarnErrReport("Ngspice returned null vector info for " + name, true);
            return;
        }

        auto* info = static_cast<pvector_info>(info_ptr);
        if (info->v_realdata == nullptr) {
            Report_m::WarnErrReport("Ngspice returned null vector data for " + name, true);
            return;
        }
        if (info->v_length <= 0) {
            Report_m::WarnErrReport("Ngspice returned empty vector for " + name, true);
            return;
        }

        const double last_val = info->v_realdata[info->v_length - 1];
        if (name != "time") {
            nodes.values[i].voltage = static_cast<RKIND>(last_val);
        } else {
            nodes.values[i].time = static_cast<RKIND_TIEMPO>(last_val);
        }
    }
}

void circuit_t::step() {
    if (has_error() != 0) {
        Report_m::WarnErrReport("Ngspice reported a controlled exit before MTLN step.", true);
        return;
    }

    updateCircuitSources(time);
    if (time == 0.0) {
        run();
    } else {
        resume();
    }

    if (has_error() != 0) {
        Report_m::WarnErrReport("Ngspice reported a controlled exit after run/resume.", true);
        return;
    }

    updateNodes();
}

RKIND circuit_t::getNodeVoltage(const std::string& name) const {
    const int idx = findVoltageIndexByName(nodes.names, name);
    if (idx <= 0 || idx > static_cast<int>(nodes.values.size())) {
        return 0.0;
    }
    return nodes.values[static_cast<size_t>(idx - 1)].voltage;
}

RKIND circuit_t::getNodeCurrent(const std::string& name) const {
    const int idx = findVoltageIndexByName(nodes.names, name);
    if (idx <= 0 || idx > static_cast<int>(nodes.values.size())) {
        return 0.0;
    }
    return nodes.values[static_cast<size_t>(idx - 1)].current;
}

RKIND_TIEMPO circuit_t::getTime() const {
    const int idx = findIndexByName(nodes.names, "time");
    if (idx <= 0 || idx > static_cast<int>(nodes.values.size())) {
        return 0.0;
    }
    return nodes.values[static_cast<size_t>(idx - 1)].time;
}

int findIndexByName(const std::vector<string_t>& names, const std::string& name) {
    for (size_t i = 0; i < names.size(); ++i) {
        const std::string sub = names[i].name.substr(0, static_cast<size_t>(names[i].length));
        if (sub == name) {
            return static_cast<int>(i) + 1;
        }
    }
    return 0;
}

int findVoltageIndexByName(const std::vector<string_t>& names, const std::string& name) {
    const std::string target_v = "V(" + name + ")";
    for (size_t i = 0; i < names.size(); ++i) {
        const std::string sub = names[i].name.substr(0, static_cast<size_t>(names[i].length));
        if (sub == target_v || sub == name) {
            return static_cast<int>(i) + 1;
        }
    }
    return 0;
}

} // namespace circuit_m
