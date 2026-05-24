#ifndef WIRES_MTLN_M_H
#define WIRES_MTLN_M_H

#include <string>

#include "mtln_types.h"

namespace Wire_bundles_mtln_m {

void solveMTLNProblem(const mtln_types_m::mtln_t& mtln_parsed, const std::string& nEntradaRoot);
void reportSimulationEnd(int layoutnumber);

} // namespace Wire_bundles_mtln_m

#endif
