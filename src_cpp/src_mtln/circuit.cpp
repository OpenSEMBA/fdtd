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
    struct node_source_t { int type=0; double value=0.0; std::string path_to_excitation; int source_type=0; double resistance=0.0; std::shared_ptr<struct node_source_t> source; std::string name; double v=0.0; double i=0.0; };
    struct string_t { std::string name; std::string type; int length=0; };
    struct vectorInfo_t { std::vector<std::string> names; std::vector<double> values; };
    struct ngspice_interface_t { void init(){} void step(){} double getTime() const { return 0.0; } };
    struct circuit_t { std::string name; double time=0.0, dt=0.0, final_time=0.0; int number_of_nodes=0; std::vector<std::shared_ptr<node_source_t>> nodes; std::vector<std::shared_ptr<string_t>> node_names; ngspice_interface_t ngspice; circuit_t()=default; void init(){} void step(){} void setSource(int,const node_source_t&){} double getNodeVoltage(int) const { return 0.0; } void updateNodeCurrent(int,double){} void setModStopTimes(double,double){} };
    circuit_t circuitCtor(const std::string& name, double dt, double final_time, int n_nodes) { circuit_t c; c.name=name; c.dt=dt; c.final_time=final_time; c.number_of_nodes=n_nodes; return c; }
} // namespace circuit_m
