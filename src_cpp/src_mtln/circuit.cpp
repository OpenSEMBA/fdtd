#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <optional>
#include <memory>

// Assuming these headers exist based on the Fortran 'use' statements
// #include "ngspice_interface_m.h"
// #include "mtln_types_m.h"
// #include "Report_m.h"
// #include "FDETYPES_m.h"

// Forward declarations for types used from other modules
struct node_source_t;
struct vectorInfo_t;

// Mapping Fortran kinds to C++ types
using RKIND = double;
using RKIND_TIEMPO = double;
using SINGLE = float;

// Placeholder for external C functions/variables if not defined in included headers
extern "C" {
    void start();
    int has_error();
    void WarnErrReport(const char* msg, bool is_error);
    void command(const char* cmd);
    void circ(const char** argv);
    void* get_vector_info(const char* name);
}

// Placeholder for external functions if not defined in included headers
int findVoltageIndexByName(const std::vector<std::string>& names, const std::string& name);
int findIndexByName(const std::vector<std::string>& names, const std::string& name);
// Note: The Fortran code uses intrinsic-like functions for index finding internally, 
// but calls findVoltageIndexByName and findIndexByName which seem to be defined elsewhere or implied.
// Based on the code provided, findIndexByName and findVoltageIndexByName are defined in this module.
// However, getNodeVoltage calls findVoltageIndexByName. 
// We will implement them below.

namespace circuit_m {

    struct string_t {
        std::string name;
        int length;
    };

    struct source_t {
        bool has_source = false;
        std::vector<RKIND_TIEMPO> time;
        std::vector<RKIND> value;
        int source_type;

        // Method declaration
        RKIND interpolate(RKIND_TIEMPO time, RKIND_TIEMPO dt) const;
    };

    struct VI_t {
        RKIND voltage;
        RKIND current;
        RKIND_TIEMPO time;
    };

    struct nodes_t {
        std::vector<VI_t> values;
        std::vector<source_t> sources;
        std::vector<string_t> names;
    };

    struct circuit_t {
        std::string name;
        RKIND_TIEMPO time = 0.0;
        RKIND_TIEMPO dt = 0.0;
        bool errorFlag = false;
        nodes_t nodes;
        nodes_t saved_nodes;

        void init(const std::vector<string_t>& names, const std::vector<node_source_t>& sources, const std::string& netlist);
        void init(const std::vector<string_t>& names, const std::vector<node_source_t>& sources);
        void init(const std::vector<string_t>& names);
        
        void run();
        void step();
        void quit();
        void readInput(const std::vector<std::string>& input, bool printInput = false);
        void setStopTimes(RKIND_TIEMPO finalTime, RKIND_TIEMPO dt);
        void setModStopTimes(RKIND_TIEMPO dt);
        RKIND getNodeVoltage(const std::string& name) const;
        RKIND getNodeCurrent(const std::string& name) const;
        void updateNodes();
        RKIND_TIEMPO getTime() const;
        void updateNodeCurrent(const std::string& node_name, RKIND current);
        void updateCircuitSources(RKIND_TIEMPO time);
        void modifyLineCapacitorValue(const std::string& name, RKIND c);
        void printCWD();

    private:
        void resume();
        void loadNetlist(const std::string& netlist);
    };

    // Implementation of source_t::interpolate
    RKIND source_t::interpolate(RKIND_TIEMPO time, RKIND_TIEMPO dt) const {
        RKIND res = 0.0;
        int n = static_cast<int>(time.size());
        if (n == 0) {
            return res;
        }

        RKIND_TIEMPO t_eval = time - dt;

        // Clamp to avoid extrapolation and division by zero at source tail.
        if (t_eval <= time[0]) {
            res = value[0];
            return res;
        }
        if (t_eval >= time[n - 1]) {
            res = value[n - 1];
            return res;
        }

        // Calculate timediff = this%time - t_eval
        std::vector<RKIND> timediff(n);
        for (int i = 0; i < n; ++i) {
            timediff[i] = static_cast<RKIND>(time[i] - t_eval);
        }

        // Find index of max element where timediff <= 0
        // Fortran: index = maxloc(timediff, 1, (timediff) <= 0)
        // This finds the first index of the maximum value in the subset where condition is true.
        // Since timediff is decreasing (time is increasing), the max value <= 0 will be the one closest to 0 from below.
        // Actually, maxloc returns the index of the maximum element. 
        // If we filter by <= 0, we want the largest value among those <= 0.
        // Since time is sorted ascending, timediff[i] = time[i] - t_eval.
        // As i increases, time[i] increases, so timediff[i] increases.
        // We want the largest timediff[i] that is <= 0. This corresponds to the largest i such that time[i] <= t_eval.
        
        int index = 0;
        bool found = false;
        for (int i = 0; i < n; ++i) {
            if (timediff[i] <= 0) {
                index = i;
                found = true;
            } else {
                break;
            }
        }
        
        if (!found) index = 1; // Fallback, though logic suggests index should be valid if t_eval < time[n-1]
        if (index >= n) index = n - 1;
        if (index >= n - 1) index = n - 2; // Ensure we have a next element for interpolation

        RKIND x1 = static_cast<RKIND>(time[index]);
        RKIND y1 = value[index];
        RKIND x2 = static_cast<RKIND>(time[index + 1]);
        RKIND y2 = value[index + 1];

        if (x2 == x1) {
            res = y2;
            return res;
        }

        res = (t_eval * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1);
        return res;
    }

    void circuit_t::printCWD() {
        // command('getcwd' // c_null_char)
        // Note: c_null_char is typically '\0'
        std::string cmd = "getcwd";
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::init(const std::vector<string_t>& names, const std::vector<node_source_t>& sources, const std::string& netlist) {
        start();
        if (!netlist.empty()) {
            std::cout << "load netlist" << std::endl;
            loadNetlist(netlist);
        }

        if (names.empty()) {
            // error stop 'Missing node names'
            throw std::runtime_error("Missing node names");
        }

        nodes.names.resize(names.size());
        nodes.values.resize(names.size());
        nodes.sources.resize(names.size());

        for (size_t i = 0; i < names.size(); ++i) {
            nodes.names[i] = names[i];
        }

        for (size_t i = 0; i < sources.size(); ++i) {
            nodes.sources[i] = setSource(sources[i].path_to_excitation);
            nodes.sources[i].source_type = sources[i].source_type;
        }
    }

    // Overloads for optional arguments
    void circuit_t::init(const std::vector<string_t>& names, const std::vector<node_source_t>& sources) {
        init(names, sources, "");
    }

    void circuit_t::init(const std::vector<string_t>& names) {
        init(names, {}, "");
    }

    source_t setSource(const std::string& source_path) {
        source_t res;
        res.has_source = true;

        if (source_path.empty()) {
            res.time = {};
            res.value = {};
            res.has_source = false;
            return res;
        }

        int line_count = 0;
        std::ifstream file(source_path);
        if (!file.is_open()) {
            std::string msg = "Cannot open excitation file: " + source_path;
            WarnErrReport(msg.c_str(), true);
            res.time = {};
            res.value = {};
            res.has_source = false;
            return res;
        }

        double time_val, value_val;
        while (file >> time_val >> value_val) {
            line_count++;
        }
        file.close();

        if (line_count == 0) {
             // Handle empty file case if necessary, though loop didn't execute
             res.time = {};
             res.value = {};
             res.has_source = false;
             return res;
        }

        res.time.resize(line_count);
        res.value.resize(line_count);

        file.open(source_path);
        for (int i = 0; i < line_count; ++i) {
            file >> res.time[i] >> res.value[i];
        }
        file.close();

        return res;
    }

    void circuit_t::loadNetlist(const std::string& netlist) {
        std::string cmd = "source " + netlist;
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::step() {
        if (has_error() != 0) {
            std::string msg = "Ngspice reported a controlled exit before MTLN step.";
            WarnErrReport(msg.c_str(), true);
            return;
        }

        updateCircuitSources(time);
        if (time == 0) {
            run();
        } else {
            resume();
        }

        if (has_error() != 0) {
            std::string msg = "Ngspice reported a controlled exit after run/resume.";
            WarnErrReport(msg.c_str(), true);
            return;
        }

        updateNodes();
    }

    void circuit_t::run() {
        std::string cmd = "run";
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::setStopTimes(RKIND_TIEMPO finalTime, RKIND_TIEMPO dt) {
        RKIND time = 0.0;
        while (time < finalTime) {
            time += dt;
            char charTime[50];
            snprintf(charTime, sizeof(charTime), "%f", static_cast<double>(time));
            std::string cmd = "stop when time = ";
            cmd += charTime;
            cmd += '\0';
            command(cmd.c_str());
        }
    }

    void circuit_t::setModStopTimes(RKIND_TIEMPO dt) {
        char charTime[50];
        snprintf(charTime, sizeof(charTime), "%f", static_cast<float>(dt));
        std::string cmd = "stop when time mod ";
        cmd += charTime;
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::resume() {
        std::string cmd = "resume";
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::quit() {
        std::string cmd = "quit 0";
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::readInput(const std::vector<std::string>& input, bool printInput) {
        if (printInput) {
            for (const auto& s : input) {
                std::cout << s << std::endl;
            }
        }

        std::vector<const char*> argv_c(input.size());
        std::vector<std::string> tmp_strs(input.size());

        for (size_t i = 0; i < input.size(); ++i) {
            tmp_strs[i] = input[i] + '\0';
            argv_c[i] = tmp_strs[i].c_str();
        }

        // circ expects const char**
        circ(argv_c.data());

        // Prevent ngspice from accumulating per-step vector history
        std::string cmd1 = "save none";
        cmd1 += '\0';
        command(cmd1.c_str());

        // Limit command history
        std::string cmd2 = "set history = 1";
        cmd2 += '\0';
        command(cmd2.c_str());
    }

    string_t getName(const void* cName) {
        string_t res;
        res.name = "";
        res.length = 0;
        
        if (cName == nullptr) return res;

        const char* f_output = static_cast<const char*>(cName);
        int i = 0;
        for (i = 0; i < 100; ++i) {
            if (f_output[i] == '\0') break;
            res.name += f_output[i];
        }
        res.length = i;
        return res;
    }

    void circuit_t::updateCircuitSources(RKIND_TIEMPO time) {
        for (size_t i = 0; i < nodes.sources.size(); ++i) {
            if (nodes.sources[i].has_source) {
                if (nodes.sources[i].source_type == 1) { // SOURCE_TYPE_VOLTAGE
                    RKIND interp = nodes.sources[i].interpolate(time, 0.0);
                    char source_value[50];
                    snprintf(source_value, sizeof(source_value), "%f", static_cast<double>(interp));
                    
                    std::string cmd = "alter @V";
                    cmd += nodes.names[i].name;
                    cmd += "_s[dc] = ";
                    cmd += source_value;
                    cmd += '\0';
                    command(cmd.c_str());
                } else if (nodes.sources[i].source_type == 2) { // SOURCE_TYPE_CURRENT
                    RKIND interp = nodes.sources[i].interpolate(time, 0.0);
                    char source_value[50];
                    snprintf(source_value, sizeof(source_value), "%f", static_cast<double>(interp));
                    
                    std::string cmd = "alter @I";
                    cmd += nodes.names[i].name;
                    cmd += "_s[dc] = ";
                    cmd += source_value;
                    cmd += '\0';
                    command(cmd.c_str());
                }
            }
        }
    }

    void circuit_t::modifyLineCapacitorValue(const std::string& name, RKIND c) {
        char sC[50];
        snprintf(sC, sizeof(sC), "%f", static_cast<double>(c));
        
        std::string cmd = "alter @CL";
        cmd += name;
        cmd += " = ";
        cmd += sC;
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::updateNodeCurrent(const std::string& node_name, RKIND current) {
        char sCurrent[50];
        snprintf(sCurrent, sizeof(sCurrent), "%f", static_cast<double>(current));
        
        std::string cmd = "alter @I";
        cmd += node_name;
        cmd += "[dc] = ";
        cmd += sCurrent;
        cmd += '\0';
        command(cmd.c_str());
    }

    void circuit_t::updateNodes() {
        if (has_error() != 0) {
            std::string msg = "Ngspice reported a controlled exit while updating nodes.";
            WarnErrReport(msg.c_str(), true);
            return;
        }

        for (size_t i = 0; i < nodes.names.size(); ++i) {
            std::string name_str = nodes.names[i].name + '\0';
            void* info_ptr = get_vector_info(name_str.c_str());
            
            if (info_ptr == nullptr) {
                std::string msg = "Ngspice returned null vector info for " + nodes.names[i].name;
                WarnErrReport(msg.c_str(), true);
                return;
            }

            vectorInfo_t* info = static_cast<vectorInfo_t*>(info_ptr);
            
            if (info->vRealData == nullptr) {
                std::string msg = "Ngspice returned null vector data for " + nodes.names[i].name;
                WarnErrReport(msg.c_str(), true);
                return;
            }
            
            if (info->vLength <= 0) {
                std::string msg = "Ngspice returned empty vector for " + nodes.names[i].name;
                WarnErrReport(msg.c_str(), true);
                return;
            }

            double* values = static_cast<double*>(info->vRealData);
            
            if (nodes.names[i].name != "time") {
                nodes.values[i].voltage = static_cast<RKIND>(values[info->vLength - 1]);
            } else {
                nodes.values[i].time = static_cast<RKIND_TIEMPO>(values[info->vLength - 1]);
            }
        }
    }

    RKIND circuit_t::getNodeVoltage(const std::string& name) const {
        int idx = findVoltageIndexByName(nodes.names, name);
        if (idx == 0 || idx > static_cast<int>(nodes.values.size())) {
            return 0.0; // Or handle error
        }
        return nodes.values[idx - 1].voltage;
    }

    RKIND circuit_t::getNodeCurrent(const std::string& name) const {
        int idx = findVoltageIndexByName(nodes.names, name);
        if (idx == 0 || idx > static_cast<int>(nodes.values.size())) {
            return 0.0; // Or handle error
        }
        return nodes.values[idx - 1].current;
    }

    RKIND_TIEMPO circuit_t::getTime() const {
        int idx = findIndexByName(nodes.names, "time");
        if (idx == 0 || idx > static_cast<int>(nodes.values.size())) {
            return 0.0; // Or handle error
        }
        return nodes.values[idx - 1].time;
    }

    int findIndexByName(const std::vector<string_t>& names, const std::string& name) {
        for (size_t i = 0; i < names.size(); ++i) {
            if (names[i].name.substr(0, names[i].length) == name) {
                return static_cast<int>(i) + 1;
            }
        }
        return 0;
    }

    int findVoltageIndexByName(const std::vector<string_t>& names, const std::string& name) {
        for (size_t i = 0; i < names.size(); ++i) {
            std::string v_name = "V(" + name + ")";
            if (names[i].name.substr(0, names[i].length) == v_name) {
                return static_cast<int>(i) + 1;
            } else if (names[i].name.substr(0, names[i].length) == name) {
                return static_cast<int>(i) + 1;
            }
        }
        return 0;
    }

}