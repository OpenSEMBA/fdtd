#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>
#include <cassert>

// Forward declarations and includes for external dependencies
// Assuming these headers exist based on the Fortran use statements
#include "ngspice_interface_m.h"
#include "mtln_types_m.h"
#include "Report_m.h"
#include "FDETYPES_m.h"

// Mapping Fortran types to C++ types based on rules
using rkind = double;
using rkind_tiempo = double;
using single = float;

// Helper to simulate c_null_char
constexpr char c_null_char = '\0';

// Helper to simulate c_char type
using c_char = char;

// Helper to simulate c_loc and c_f_pointer behavior roughly
// In a real implementation, these would wrap actual C pointer operations.
// For this translation, we assume standard C++ pointer semantics where applicable,
// but we must keep the structure compatible with the logic.

// Mocking external functions/types that are not defined in the snippet but used
// These would normally come from the included headers.

namespace circuit_m {

    // Derived type: string_t
    struct string_t {
        std::string name;
        int length;

        string_t() : name(""), length(0) {}
    };

    // Derived type: source_t
    class source_t {
    public:
        bool has_source = false;
        std::vector<rkind_tiempo> time;
        std::vector<rkind> value;
        int source_type;

        source_t() : has_source(false), source_type(0) {}

        // Method: interpolate
        rkind interpolate(rkind_tiempo time, rkind dt) const {
            int n = static_cast<int>(time.size());
            if (n == 0) {
                return 0.0_rkind;
            }

            rkind_tiempo t_eval = time - dt;

            // Clamp to avoid extrapolation and division by zero at source tail.
            if (t_eval <= time[0]) {
                return value[0];
            }
            if (t_eval >= time[n - 1]) {
                return value[n - 1];
            }

            // Calculate timediff
            std::vector<rkind_tiempo> timediff(n);
            for (int i = 0; i < n; ++i) {
                timediff[i] = time[i] - t_eval;
            }

            // Find index where timediff <= 0. 
            // Fortran maxloc with mask returns the index of the first element satisfying the condition 
            // if multiple exist? No, maxloc returns the index of the maximum value. 
            // The mask logic in Fortran: maxloc(array, dim, mask) returns index of max value where mask is true.
            // Here mask is (timediff <= 0). We want the index of the maximum value among those <= 0.
            // Since timediff is likely decreasing (time is increasing), the first one <= 0 is the largest negative/zero.
            
            int index = 0;
            bool found = false;
            for (int i = 0; i < n; ++i) {
                if (timediff[i] <= 0) {
                    index = i + 1; // Fortran is 1-based
                    found = true;
                    break; 
                }
            }
            
            // If no element <= 0, maxloc returns 0 in Fortran (or undefined behavior depending on version, usually 0 or 1)
            // The code handles index == 0 by setting to 1.
            if (!found) {
                index = 1;
            }

            // Adjust for 1-based indexing logic in Fortran vs 0-based in C++
            // In Fortran: index is 1-based. 
            // If index >= n, set to n-1.
            // In C++, we need to map this carefully.
            
            // Let's re-evaluate the Fortran logic:
            // index = maxloc(timediff, 1, (timediff) <= 0)
            // If no element satisfies mask, maxloc returns 0.
            // if (index == 0) index = 1;
            // if (index >= n) index = n - 1;
            
            // In C++, let's find the first index i such that time[i] >= t_eval.
            // This corresponds to the interpolation interval.
            // timediff[i] = time[i] - t_eval.
            // We want the largest time[i] such that time[i] <= t_eval? 
            // No, maxloc finds the MAXIMUM value in timediff where timediff <= 0.
            // Since time is sorted ascending, timediff is sorted descending.
            // The first element <= 0 is the largest value <= 0.
            // So we want the first i where time[i] <= t_eval.
            
            // Let's stick to the explicit calculation to be safe.
            int idx_f = 0;
            for (int i = 0; i < n; ++i) {
                if (timediff[i] <= 0) {
                    idx_f = i + 1; // 1-based
                    break;
                }
            }
            if (idx_f == 0) idx_f = 1;
            if (idx_f >= n) idx_f = n - 1;

            // Convert to 0-based for access
            int i1 = idx_f - 1;
            int i2 = i1 + 1;

            // Ensure bounds
            if (i1 < 0) i1 = 0;
            if (i2 >= n) i2 = n - 1;

            rkind x1 = time[i1];
            rkind y1 = value[i1];
            rkind x2 = time[i2];
            rkind y2 = value[i2];

            if (x2 == x1) {
                return y2;
            }

            return (t_eval * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1);
        }
    };

    // Derived type: VI_t
    struct VI_t {
        rkind voltage;
        rkind current;
        rkind_tiempo time;

        VI_t() : voltage(0.0), current(0.0), time(0.0) {}
    };

    // Derived type: nodes_t
    struct nodes_t {
        std::vector<VI_t> values;
        std::vector<source_t> sources;
        std::vector<string_t> names;
    };

    // Forward declaration of external functions used in circuit_t
    void start();
    void command(const char* cmd);
    bool has_error();
    void WarnErrReport(const std::string& msg, bool fatal);
    source_t setSource(const std::string& source_path);
    void circ(std::vector<void*> argv_c);
    void* get_vector_info(const char* name);
    int findVoltageIndexByName(const std::vector<string_t>& names, const std::string& name);
    int findIndexByName(const std::vector<string_t>& names, const std::string& name);

    // External struct definitions assumed from headers
    struct vectorInfo_t {
        double* vRealData;
        int vLength;
    };

    // Derived type: circuit_t
    class circuit_t {
    public:
        std::string name;
        rkind_tiempo time = 0.0;
        rkind_tiempo dt = 0.0;
        bool errorFlag = false;
        nodes_t nodes;
        nodes_t saved_nodes;

        // Methods
        void init(const std::vector<string_t>& names, const std::vector<node_source_t>& sources, const std::string& netlist = "");
        void run();
        void step();
        void quit();
        void readInput(const std::vector<std::string>& input, bool printInput = false);
        void setStopTimes(rkind_tiempo finalTime, rkind_tiempo dt);
        void setModStopTimes(rkind_tiempo dt);
        rkind getNodeVoltage(const std::string& name);
        rkind getNodeCurrent(const std::string& name);
        rkind_tiempo getTime();
        void updateNodes();
        void updateNodeCurrent(const std::string& node_name, rkind current);
        void updateCircuitSources(rkind_tiempo time);
        void modifyLineCapacitorValue(const std::string& name, rkind c);
        void printCWD();

    private:
        void resume();
        void loadNetlist(const std::string& netlist);
    };

    // Implementation of circuit_t methods

    void circuit_t::printCWD() {
        // command('getcwd' // c_null_char)
        // Note: command expects a null-terminated string.
        std::string cmd = "getcwd";
        cmd += c_null_char;
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

        if (!sources.empty()) {
            for (size_t i = 0; i < sources.size(); ++i) {
                nodes.sources[i] = setSource(sources[i].path_to_excitation);
                nodes.sources[i].source_type = sources[i].source_type;
            }
        }
    }

    void circuit_t::loadNetlist(const std::string& netlist) {
        std::string cmd = "source ";
        cmd += netlist;
        cmd += c_null_char;
        command(cmd.c_str());
    }

    void circuit_t::step() {
        if (has_error() != 0) {
            WarnErrReport("Ngspice reported a controlled exit before MTLN step.", true);
            return;
        }

        updateCircuitSources(time);
        if (time == 0) {
            run();
        } else {
            resume();
        }

        if (has_error() != 0) {
            WarnErrReport("Ngspice reported a controlled exit after run/resume.", true);
            return;
        }

        updateNodes();
    }

    void circuit_t::run() {
        std::string cmd = "run";
        cmd += c_null_char;
        command(cmd.c_str());
    }

    void circuit_t::setStopTimes(rkind_tiempo finalTime, rkind_tiempo dt) {
        rkind time = 0.0;
        while (time < finalTime) {
            time += dt;
            std::ostringstream oss;
            oss << time;
            std::string charTime = oss.str();
            std::string cmd = "stop when time = " + charTime;
            cmd += c_null_char;
            command(cmd.c_str());
        }
    }

    void circuit_t::setModStopTimes(rkind_tiempo dt) {
        std::ostringstream oss;
        oss << static_cast<float>(dt);
        std::string charTime = oss.str();
        std::string cmd = "stop when time mod " + charTime;
        cmd += c_null_char;
        command(cmd.c_str());
    }

    void circuit_t::resume() {
        std::string cmd = "resume";
        cmd += c_null_char;
        command(cmd.c_str());
    }

    void circuit_t::quit() {
        std::string cmd = "quit 0";
        cmd += c_null_char;
        command(cmd.c_str());
    }

    void circuit_t::readInput(const std::vector<std::string>& input, bool printInput) {
        if (printInput) {
            for (size_t i = 0; i < input.size(); ++i) {
                std::cout << input[i] << std::endl;
            }
        }

        std::vector<std::string> tmp(input.size());
        std::vector<void*> argv_c(input.size());

        for (size_t i = 0; i < input.size(); ++i) {
            tmp[i] = input[i] + c_null_char;
            argv_c[i] = const_cast<char*>(tmp[i].c_str());
        }

        circ(argv_c);
        
        // Prevent ngspice from accumulating per-step vector history
        std::string cmd1 = "save none";
        cmd1 += c_null_char;
        command(cmd1.c_str());

        // Limit command history
        std::string cmd2 = "set history = 1";
        cmd2 += c_null_char;
        command(cmd2.c_str());
    }

    string_t getName(void* cName) {
        string_t res;
        res.name = "";
        res.length = 0;
        
        // c_f_pointer(cName, f_output, [100])
        // Assuming cName points to a null-terminated char array
        char* f_output = static_cast<char*>(cName);
        if (!f_output) return res;

        int i = 0;
        for (i = 0; i < 100; ++i) {
            if (f_output[i] == c_null_char) break;
            res.name += f_output[i];
        }
        res.length = i;

        return res;
    }

    void circuit_t::updateCircuitSources(rkind_tiempo time) {
        for (size_t i = 0; i < nodes.sources.size(); ++i) {
            if (nodes.sources[i].has_source) {
                if (nodes.sources[i].source_type == SOURCE_TYPE_VOLTAGE) {
                    rkind interp = nodes.sources[i].interpolate(time, 0.0_rkind_tiempo);
                    std::ostringstream oss;
                    oss << interp;
                    std::string source_value = oss.str();
                    
                    std::string cmd = "alter @V" + nodes.names[i].name + "_s[dc] = " + source_value;
                    cmd += c_null_char;
                    command(cmd.c_str());
                } else if (nodes.sources[i].source_type == SOURCE_TYPE_CURRENT) {
                    rkind interp = nodes.sources[i].interpolate(time, 0.0_rkind_tiempo);
                    std::ostringstream oss;
                    oss << interp;
                    std::string source_value = oss.str();
                    
                    std::string cmd = "alter @I" + nodes.names[i].name + "_s[dc] = " + source_value;
                    cmd += c_null_char;
                    command(cmd.c_str());
                }
            }
        }
    }

    void circuit_t::modifyLineCapacitorValue(const std::string& name, rkind c) {
        std::ostringstream oss;
        oss << c;
        std::string sC = oss.str();
        
        std::string cmd = "alter @CL" + name + " = " + sC;
        cmd += c_null_char;
        command(cmd.c_str());
    }

    void circuit_t::updateNodeCurrent(const std::string& node_name, rkind current) {
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
            // Default case if neither found? Fortran code doesn't set sCurrent in this case.
            // This might be a bug in original Fortran or implies sCurrent is uninitialized.
            // We'll assume it's handled by the caller or defaults to 0.
            sCurrent = "0";
        }

        std::string cmd = "alter @I" + node_name + "[dc] = " + sCurrent;
        cmd += c_null_char;
        command(cmd.c_str());
    }

    void circuit_t::updateNodes() {
        if (has_error() != 0) {
            WarnErrReport("Ngspice reported a controlled exit while updating nodes.", true);
            return;
        }

        for (size_t i = 0; i < nodes.names.size(); ++i) {
            std::string name_str = nodes.names[i].name + c_null_char;
            void* info_ptr = get_vector_info(name_str.c_str());
            
            if (info_ptr == nullptr) {
                WarnErrReport("Ngspice returned null vector info for " + nodes.names[i].name, true);
                return;
            }

            vectorInfo_t* info = static_cast<vectorInfo_t*>(info_ptr);
            
            if (info->vRealData == nullptr) {
                WarnErrReport("Ngspice returned null vector data for " + nodes.names[i].name, true);
                return;
            }
            
            if (info->vLength <= 0) {
                WarnErrReport("Ngspice returned empty vector for " + nodes.names[i].name, true);
                return;
            }

            // c_f_pointer(info%vRealData, values, shape=[info%vLength])
            // In C++, we just access the array directly.
            // values is a pointer to double.
            double* values = info->vRealData;
            
            // ubound(values,1) is the last index.
            double last_val = values[info->vLength - 1];

            if (nodes.names[i].name != "time") {
                nodes.values[i].voltage = static_cast<rkind>(last_val);
            } else {
                nodes.values[i].time = static_cast<rkind_tiempo>(last_val);
            }
        }
    }

    rkind circuit_t::getNodeVoltage(const std::string& name) {
        int idx = findVoltageIndexByName(nodes.names, name);
        if (idx == 0 || idx > static_cast<int>(nodes.values.size())) {
            return 0.0; // Error handling or default
        }
        return nodes.values[idx - 1].voltage;
    }

    rkind circuit_t::getNodeCurrent(const std::string& name) {
        int idx = findVoltageIndexByName(nodes.names, name);
        if (idx == 0 || idx > static_cast<int>(nodes.values.size())) {
            return 0.0; // Error handling or default
        }
        return nodes.values[idx - 1].current;
    }

    rkind_tiempo circuit_t::getTime() {
        int idx = findIndexByName(nodes.names, "time");
        if (idx == 0 || idx > static_cast<int>(nodes.values.size())) {
            return 0.0; // Error handling or default
        }
        return nodes.values[idx - 1].time;
    }

    int findIndexByName(const std::vector<string_t>& names, const std::string& name) {
        for (size_t i = 0; i < names.size(); ++i) {
            // Fortran: names(i)%name(1:names(i)%length) == trim(name)
            std::string sub = names[i].name.substr(0, names[i].length);
            if (sub == name) {
                return static_cast<int>(i) + 1; // 1-based index
            }
        }
        return 0;
    }

    int findVoltageIndexByName(const std::vector<string_t>& names, const std::string& name) {
        for (size_t i = 0; i < names.size(); ++i) {
            std::string sub = names[i].name.substr(0, names[i].length);
            std::string target_v = "V(" + name + ")";
            if (sub == target_v) {
                return static_cast<int>(i) + 1;
            } else if (sub == name) {
                return static_cast<int>(i) + 1;
            }
        }
        return 0;
    }

} // namespace circuit_m