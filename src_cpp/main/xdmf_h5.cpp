#include "xdmf_h5.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#ifdef SEMBA_CPP_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace xdmf_h5_m {
namespace {

#ifdef SEMBA_CPP_ENABLE_HDF5
hid_t g_file_id = -1;
hid_t g_dset_id = -1;
hid_t g_dspace_id = -1;
hid_t g_slice_id = -1;
hsize_t g_dims[4] = {0, 0, 0, 0};
hsize_t g_slice_dims[4] = {0, 0, 0, 1};
hsize_t g_offset[4] = {0, 0, 0, 0};
std::ofstream g_xdmf;
#endif

std::string trim(const std::string& value) {
    const size_t first = value.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return {};
    const size_t last = value.find_last_not_of(" \t\n\r");
    return value.substr(first, last - first + 1);
}

} // namespace

void openh5file(const std::string& filename, int finalstep, int minXabs, int maxXabs,
                int minYabs, int maxYabs, int minZabs, int maxZabs) {
#ifdef SEMBA_CPP_ENABLE_HDF5
    // Match Fortran/h5py layout: (time, z, y, x) — see xdmf_h5.F90 HDF5 transpose note.
    const hsize_t nx = static_cast<hsize_t>(maxXabs - minXabs + 1);
    const hsize_t ny = static_cast<hsize_t>(maxYabs - minYabs + 1);
    const hsize_t nz = static_cast<hsize_t>(maxZabs - minZabs + 1);
    g_dims[0] = static_cast<hsize_t>(finalstep);
    g_dims[1] = nz;
    g_dims[2] = ny;
    g_dims[3] = nx;
    g_slice_dims[0] = 1;
    g_slice_dims[1] = nz;
    g_slice_dims[2] = ny;
    g_slice_dims[3] = nx;

    const std::string h5_path = trim(filename) + ".h5";
    g_file_id = H5Fcreate(h5_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    g_dspace_id = H5Screate_simple(4, g_dims, nullptr);
    g_slice_id = H5Screate_simple(4, g_slice_dims, nullptr);
    g_dset_id = H5Dcreate2(g_file_id, "data", H5T_NATIVE_FLOAT, g_dspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    g_xdmf.open(trim(filename) + ".xdmf");
    g_xdmf << "<Xdmf>\n<Domain>\n";
    g_xdmf << "<Grid Name=\"GridTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";
#else
    (void)filename;
    (void)finalstep;
    (void)minXabs;
    (void)maxXabs;
    (void)minYabs;
    (void)maxYabs;
    (void)minZabs;
    (void)maxZabs;
#endif
}

void writeh5file(const std::string& filename, const float* valor3d, int nx, int ny, int nz,
                 int indi, double attindi, int minXabs, int maxXabs, int minYabs, int maxYabs,
                 int minZabs, int maxZabs, double linez_minZabs_primero,
                 double liney_minYabs_primero, double linex_minXabs_primero, double dz_minZabs,
                 double dy_minYabs, double dx_minXabs, int minZabs_primero, int minYabs_primero,
                 int minXabs_primero, int finalstep, bool vtkindex) {
#ifdef SEMBA_CPP_ENABLE_HDF5
    (void)filename;
    g_offset[0] = static_cast<hsize_t>(indi - 1);
    g_offset[1] = 0;
    g_offset[2] = 0;
    g_offset[3] = 0;

    H5Sselect_hyperslab(g_dspace_id, H5S_SELECT_SET, g_offset, nullptr, g_slice_dims, nullptr);
    H5Dwrite(g_dset_id, H5T_NATIVE_FLOAT, g_slice_id, g_dspace_id, H5P_DEFAULT, valor3d);

    char time_text[32];
    std::snprintf(time_text, sizeof(time_text), "% .9E", attindi);

    const int dim_x = maxXabs - minXabs + 1;
    const int dim_y = maxYabs - minYabs + 1;
    const int dim_z = maxZabs - minZabs + 1;

    g_xdmf << "<Grid Name=\"IntGrid\" GridType=\"Uniform\"  CollectionType=\"Spatial\">>\n";
    g_xdmf << "<Time Value=\"" << trim(time_text) << "\" />\n";
    g_xdmf << "<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << dim_z << " "
           << dim_y << " " << dim_x << "\">\n</Topology>\n";
    g_xdmf << "<Geometry Type=\"ORIGIN_DXDYDZ\">\n<DataItem Format=\"XML\" Dimensions=\"3\">\n";
    if (vtkindex) {
        g_xdmf << minZabs_primero << " " << minYabs_primero << " " << minXabs_primero << "\n";
    } else {
        g_xdmf << linez_minZabs_primero << " " << liney_minYabs_primero << " "
               << linex_minXabs_primero << "\n";
    }
    g_xdmf << "</DataItem>\n<DataItem Format=\"XML\" Dimensions=\"3\">\n";
    g_xdmf << dz_minZabs << " " << dy_minYabs << " " << dx_minXabs << "\n</DataItem>\n</Geometry>\n";
    g_xdmf << "<Attribute Name=\"IntValues\" Center=\"Node\">\n";
    g_xdmf << "<DataItem ItemType=\"HyperSlab\" Dimensions=\"1 " << dim_z << " " << dim_y << " "
           << dim_x << "\" Format=\"XML\">\n";
    g_xdmf << "<DataItem Dimensions=\"3 4\" Format=\"XML\">\n";
    g_xdmf << g_offset[0] << " 0 0 0\n1 1 1 1\n1 " << dim_z << " " << dim_y << " " << dim_x
           << "\n</DataItem>\n";
    g_xdmf << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\""
           << finalstep << " " << dim_z << " " << dim_y << " " << dim_x << "\">\n";
    g_xdmf << trim(filename) << ".h5:/data\n</DataItem>\n</DataItem>\n</Attribute>\n</Grid>\n";
#else
    (void)filename;
    (void)valor3d;
    (void)nx;
    (void)ny;
    (void)nz;
    (void)indi;
    (void)attindi;
    (void)minXabs;
    (void)maxXabs;
    (void)minYabs;
    (void)maxYabs;
    (void)minZabs;
    (void)maxZabs;
    (void)linez_minZabs_primero;
    (void)liney_minYabs_primero;
    (void)linex_minXabs_primero;
    (void)dz_minZabs;
    (void)dy_minYabs;
    (void)dx_minXabs;
    (void)minZabs_primero;
    (void)minYabs_primero;
    (void)minXabs_primero;
    (void)finalstep;
    (void)vtkindex;
#endif
}

void closeh5file(int finalstep, const std::vector<double>& att) {
#ifdef SEMBA_CPP_ENABLE_HDF5
    H5Dclose(g_dset_id);
    g_dset_id = -1;

    const hsize_t time_dims[1] = {static_cast<hsize_t>(finalstep)};
    hid_t time_space = H5Screate_simple(1, time_dims, nullptr);
    hid_t time_set = H5Dcreate2(g_file_id, "Time", H5T_NATIVE_FLOAT, time_space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    std::vector<float> att_f(att.size());
    for (size_t i = 0; i < att.size(); ++i) {
        att_f[i] = static_cast<float>(att[i]);
    }
    H5Dwrite(time_set, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, att_f.data());
    H5Dclose(time_set);
    H5Sclose(time_space);

    H5Sclose(g_slice_id);
    H5Sclose(g_dspace_id);
    H5Fclose(g_file_id);
    g_slice_id = -1;
    g_dspace_id = -1;
    g_file_id = -1;

    g_xdmf << "</Grid>\n</Domain>\n</Xdmf>\n";
    g_xdmf.close();
#else
    (void)finalstep;
    (void)att;
#endif
}

} // namespace xdmf_h5_m
