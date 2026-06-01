#include "xdmf_h5.h"

#include <string>
#include <vector>

namespace snapxdmf_m {

void WRITE_XDMFSNAP(int ninstant, const std::string& filename, int minXabs, int maxXabs, int minYabs,
                    int maxYabs, int minZabs, int maxZabs,
                    const std::vector<std::vector<std::vector<std::vector<float>>>>& valor3D) {
#ifdef SEMBA_CPP_ENABLE_HDF5
    const int nx = maxXabs - minXabs + 1;
    const int ny = maxYabs - minYabs + 1;
    const int nz = maxZabs - minZabs + 1;
    const int finalstep = 1;

    std::vector<float> slab(static_cast<size_t>(nx * ny * nz), 0.0f);
    size_t idx = 0;
    for (int k = minZabs; k <= maxZabs; ++k) {
        for (int j = minYabs; j <= maxYabs; ++j) {
            for (int i = minXabs; i <= maxXabs; ++i) {
                const int ii = i - minXabs;
                const int jj = j - minYabs;
                const int kk = k - minZabs;
                if (ii >= 0 && ii < static_cast<int>(valor3D.size()) && jj >= 0 &&
                    jj < static_cast<int>(valor3D[ii].size()) && kk >= 0 &&
                    kk < static_cast<int>(valor3D[ii][jj].size()) &&
                    !valor3D[ii][jj][kk].empty()) {
                    slab[idx] = valor3D[ii][jj][kk][0];
                }
                ++idx;
            }
        }
    }

    xdmf_h5_m::openh5file(filename, finalstep, minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs);
    xdmf_h5_m::writeh5file(filename, slab.data(), nx, ny, nz, 1, static_cast<double>(ninstant),
                           minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs,
                           static_cast<double>(minZabs), static_cast<double>(minYabs),
                           static_cast<double>(minXabs), 1.0, 1.0, 1.0, minZabs, minYabs, minXabs,
                           finalstep, true);
    xdmf_h5_m::closeh5file(finalstep, {static_cast<double>(ninstant)});
#else
    (void)ninstant;
    (void)filename;
    (void)minXabs;
    (void)maxXabs;
    (void)minYabs;
    (void)maxYabs;
    (void)minZabs;
    (void)maxZabs;
    (void)valor3D;
#endif
}

} // namespace snapxdmf_m
