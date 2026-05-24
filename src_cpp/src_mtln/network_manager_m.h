#ifndef NETWORK_MANAGER_M_H
#define NETWORK_MANAGER_M_H

#include <vector>

#include "circuit_m.h"
#include "mtln_types.h"
#include "network_m.h"

namespace network_manager_m {

using circuit_m::circuit_t;
using mtln_types_m::RKIND;
using mtln_types_m::RKIND_TIEMPO;
using network_m::network_t;

struct network_manager_t {
    std::vector<network_t> networks;
    circuit_t circuit;
    RKIND time = 0.0;
    RKIND dt = 0.0;

    void advanceVoltage();
    void updateCircuitCurrentsFromNetwork();
    void updateNetworkVoltages();
};

network_manager_t network_managerCtor(const std::vector<network_t>& networks,
                                      const std::vector<std::string>& description,
                                      RKIND_TIEMPO final_time,
                                      RKIND_TIEMPO dt);

} // namespace network_manager_m

#endif
