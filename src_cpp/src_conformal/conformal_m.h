#ifndef CONFORMAL_M_H
#define CONFORMAL_M_H

#include <vector>

#include "cell_map_m.h"
#include "nfde_types.h"

namespace conformal_m {

// Conformal geometry regions (triangle_t volumes), matching Fortran NFDETypes layout.
struct ConformalPECRegions_t {
    std::vector<cell_map_m::ConformalPECElements_t> volumes;
};

using NFDETypes_m::ConformalMedia_t;
using NFDETypes_m::edge_t;
using NFDETypes_m::face_t;
using NFDETypes_m::conformal_edge_media_t;
using NFDETypes_m::conformal_face_media_t;

std::vector<ConformalMedia_t> buildConformalMedia(const ConformalPECRegions_t& regions);

} // namespace conformal_m

#endif // CONFORMAL_M_H
