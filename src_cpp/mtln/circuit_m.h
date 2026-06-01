#ifndef CIRCUIT_M_H
#define CIRCUIT_M_H

#include <string>
#include <vector>

#include "mtln_types.h"

namespace circuit_m {

using mtln_types_m::RKIND;
using mtln_types_m::RKIND_TIEMPO;
using mtln_types_m::SOURCE_TYPE_CURRENT;
using mtln_types_m::SOURCE_TYPE_VOLTAGE;
using mtln_types_m::node_source_t;

struct string_t {
    std::string name;
    int length = 0;

    string_t() = default;
    string_t(const std::string& n, int len) : name(n), length(len) {}
};

struct source_t {
    bool has_source = false;
    std::vector<RKIND_TIEMPO> time;
    std::vector<RKIND> value;
    int source_type = 0;

    RKIND interpolate(RKIND_TIEMPO eval_time, RKIND_TIEMPO dt_step) const;
};

struct VI_t {
    RKIND voltage = 0.0;
    RKIND current = 0.0;
    RKIND_TIEMPO time = 0.0;
};

struct nodes_t {
    std::vector<VI_t> values;
    std::vector<source_t> sources;
    std::vector<string_t> names;
};

class circuit_t {
public:
    std::string name;
    RKIND_TIEMPO time = 0.0;
    RKIND_TIEMPO dt = 0.0;
    bool errorFlag = false;
    nodes_t nodes;
    nodes_t saved_nodes;

    void init(const std::vector<string_t>& names,
              const std::vector<node_source_t>& sources = {},
              const std::string& netlist = "");
    void run();
    void step();
    void quit();
    void readInput(const std::vector<std::string>& input, bool printInput = false);
    void setStopTimes(RKIND_TIEMPO finalTime, RKIND_TIEMPO dt_in);
    void setModStopTimes(RKIND_TIEMPO dt_in);
    RKIND getNodeVoltage(const std::string& name) const;
    RKIND getNodeCurrent(const std::string& name) const;
    RKIND_TIEMPO getTime() const;
    void updateNodeCurrent(const std::string& node_name, RKIND current);
    void modifyLineCapacitorValue(const std::string& name, RKIND c);
    void printCWD();

private:
    void resume();
    void loadNetlist(const std::string& netlist);
    void updateNodes();
    void updateCircuitSources(RKIND_TIEMPO eval_time);
};

source_t setSource(const std::string& source_path);
int findIndexByName(const std::vector<string_t>& names, const std::string& name);
int findVoltageIndexByName(const std::vector<string_t>& names, const std::string& name);

} // namespace circuit_m

#endif
