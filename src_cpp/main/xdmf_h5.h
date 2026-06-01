#ifndef XDMF_H5_H
#define XDMF_H5_H

#include <string>
#include <vector>

namespace xdmf_h5_m {

void openh5file(const std::string& filename, int finalstep, int minXabs, int maxXabs,
                int minYabs, int maxYabs, int minZabs, int maxZabs);

void writeh5file(const std::string& filename,
                 const float* valor3d, int nx, int ny, int nz, int indi, double attindi,
                 int minXabs, int maxXabs, int minYabs, int maxYabs, int minZabs, int maxZabs,
                 double linez_minZabs_primero, double liney_minYabs_primero,
                 double linex_minXabs_primero, double dz_minZabs, double dy_minYabs,
                 double dx_minXabs, int minZabs_primero, int minYabs_primero,
                 int minXabs_primero, int finalstep, bool vtkindex);

void closeh5file(int finalstep, const std::vector<double>& att);

} // namespace xdmf_h5_m

#endif
