#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cstring>

// Assuming these types are defined elsewhere or standard equivalents
// HID_T, HSIZET, etc. usually map to long long or specific HDF5 types.
// For this translation, we assume a basic HDF5 wrapper or direct calls if headers were included.
// Since we are translating to C++ without the full Fortran context, we will simulate the HDF5 calls
// or assume a C++ HDF5 interface exists. However, to strictly follow "translate", we will keep the
// function signatures similar but adapted to C++.

// Note: The original code uses HDF5 Fortran bindings (h5xxx_f). In C++, we would typically use
// the C++ API or C API. For the purpose of this translation, I will assume a hypothetical
// C++ wrapper or simply translate the logic assuming the underlying HDF5 calls are available
// via a namespace or global functions that mimic the Fortran interface for ease of porting,
// OR more likely, replace them with standard C++ file I/O if the logic allows, but here
// the logic is tightly coupled to HDF5.

// To make this compile-ready as a standalone snippet, I will define stubs for the HDF5 calls
// or assume they are provided by a library. Given the instruction to translate, I will
// convert the structure and logic.

// Type definitions based on context
typedef long long HID_T;
typedef long long HSIZET;
typedef double RKIND; // Default real kind
typedef double RKIND_tiempo; // Time real kind
const int BUFSIZE = 256;

// HDF5 Constants (Simplified)
const int H5F_ACC_TRUNC_F = 0;
const int H5T_NATIVE_DOUBLE = 0;
const int H5T_NATIVE_REAL = 1;
const int H5T_NATIVE_LDOUBLE = 2;
const int H5S_SELECT_SET_F = 0;

// Mock HDF5 functions to make the code structurally complete for translation purposes
// In a real scenario, these would be replaced by actual HDF5 C++ API calls.
namespace HDF5 {
    void h5open_f(int& error) { error = 0; }
    void h5fcreate_f(const std::string& name, int flags, HID_T& file_id, int& error) { error = 0; file_id = 1; }
    void h5screate_simple_f(int rank, const std::vector<HSIZET>& dims, HID_T& dspace_id, int& error) { error = 0; dspace_id = 1; }
    void h5dcreate_f(HID_T file_id, const std::string& name, int dtype, HID_T dspace_id, HID_T& dset_id, int& error) { error = 0; dset_id = 1; }
    void h5sselect_hyperslab_f(HID_T dspace_id, int op, const std::vector<HSIZET>& offset, const std::vector<HSIZET>& count, int& error) { error = 0; }
    void h5dwrite_f(HID_T dset_id, int dtype, const std::vector<std::vector<std::vector<std::vector<RKIND>>>>& data, const std::vector<HSIZET>& dims, int& error, HID_T slice_id, HID_T dspace_id) { error = 0; }
    void h5dclose_f(HID_T dset_id, int& error) { error = 0; }
    void h5sclose_f(HID_T dspace_id, int& error) { error = 0; }
    void h5fclose_f(HID_T file_id, int& error) { error = 0; }
    void h5close_f(int& error) { error = 0; }
}

// Helper to trim whitespace
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (std::string::npos == first) return str;
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

// Helper to adjustl (left justify)
std::string adjustl(const std::string& str) {
    return trim(str);
}

namespace xdmf_h5_m {

    // Global variables moved to namespace scope or class members.
    // Since the original module has module-level variables, we keep them in the namespace.
    HID_T file_id;
    HID_T dset_id;
    HID_T dspace_id;
    HID_T slice2D_id;
    std::vector<HSIZET> DATA_dims;
    std::vector<HSIZET> offset;
    std::vector<HSIZET> valor3d_dims;
    
    // File handle for XDMF output
    std::ofstream xdmf_file;

    // Constants
    const int RANK = 4;

    void openh5file(const std::string& filename, int finalstep, int minXabs, int maxXabs, int minYabs, int maxYabs, int minZabs, int maxZabs) {
        int error = 0;
        std::string dsetname = "data";
        int rank = 4;

        DATA_dims.resize(rank);
        valor3d_dims.resize(rank);
        offset.resize(rank);

        DATA_dims[0] = maxXabs - minXabs + 1;
        DATA_dims[1] = maxYabs - minYabs + 1;
        DATA_dims[2] = maxZabs - minZabs + 1;
        DATA_dims[3] = finalstep;

        valor3d_dims[0] = DATA_dims[0];
        valor3d_dims[1] = DATA_dims[1];
        valor3d_dims[2] = DATA_dims[2];
        valor3d_dims[3] = 1;

        // HDF5 Open/Create
        HDF5::h5open_f(error);
        std::string h5filename = trim(adjustl(filename)) + ".h5";
        HDF5::h5fcreate_f(h5filename, H5F_ACC_TRUNC_F, file_id, error);
        HDF5::h5screate_simple_f(rank, DATA_dims, dspace_id, error);
        HDF5::h5screate_simple_f(rank, valor3d_dims, slice2D_id, error);
        
        // Determine data type based on compile flag (simulated)
        // In C++, we usually template or use templates. Here we assume double for simplicity or check a macro.
        // The original code uses #ifdef CompileWithReal8. We'll assume double for now.
        HDF5::h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error);

        // XDMF part
        std::string xdmffilename = trim(adjustl(filename)) + ".xdmf";
        xdmf_file.open(xdmffilename);
        if (xdmf_file.is_open()) {
            xdmf_file << "<Xdmf>" << std::endl;
            xdmf_file << "<Domain>" << std::endl;
            xdmf_file << "<Grid Name=\"GridTime\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
        }
    }

    void writeh5file(const std::string& filename, 
                     const std::vector<std::vector<std::vector<std::vector<RKIND>>>>& valor3d,
                     int indi, RKIND attindi, 
                     int minXabs, int maxXabs, int minYabs, int maxYabs, int minZabs, int maxZabs,
                     RKIND linez_minZabs_primero, RKIND liney_minYabs_primero, RKIND linex_minXabs_primero,
                     RKIND dz_minZabs, RKIND dy_minYabs, RKIND dx_minXabs,
                     int minZabs_primero, int minYabs_primero, int minXabs_primero,
                     int finalstep, bool vtkindex) {
        
        int error = 0;
        std::string dsetname = "data";
        std::string charc;

        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;
        offset[3] = indi - 1;

        HDF5::h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, valor3d_dims, error);
        
        // Write HDF5 data
        // Note: The Fortran call passes valor3d which is 4D. The C++ vector is 4D.
        // The dims passed to write are valor3d_dims.
        HDF5::h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, valor3d, valor3d_dims, error, slice2D_id, dspace_id);

        // XDMF part
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(9) << attindi;
        charc = oss.str();
        
        DATA_dims[0] = maxXabs - minXabs + 1;
        DATA_dims[1] = maxYabs - minYabs + 1;
        DATA_dims[2] = maxZabs - minZabs + 1;

        if (xdmf_file.is_open()) {
            xdmf_file << "<Grid Name=\"IntGrid\" GridType=\"Uniform\"  CollectionType=\"Spatial\">>" << std::endl;
            xdmf_file << "<Time Value=\"" << trim(adjustl(charc)) << "\" />" << std::endl;
            xdmf_file << "<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" 
                      << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << "\">" << std::endl;
            xdmf_file << "</Topology>" << std::endl;
            xdmf_file << "<Geometry Type=\"ORIGIN_DXDYDZ\">" << std::endl;
            xdmf_file << "<DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
            
            if (vtkindex) {
                xdmf_file << minZabs_primero << " " << minYabs_primero << " " << minXabs_primero << std::endl;
            } else {
                xdmf_file << linez_minZabs_primero << " " << liney_minYabs_primero << " " << linex_minXabs_primero << std::endl;
            }
            
            xdmf_file << "</DataItem>" << std::endl;
            xdmf_file << "<DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
            xdmf_file << dz_minZabs << " " << dy_minYabs << " " << dx_minXabs << std::endl;
            xdmf_file << "</DataItem>" << std::endl;
            xdmf_file << "</Geometry>" << std::endl;
            xdmf_file << "<Attribute Name=\"IntValues\" Center=\"Node\">" << std::endl;
            xdmf_file << "<DataItem ItemType=\"HyperSlab\" Dimensions=\"" 
                      << "1 " << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << "\" Format=\"XML\">" << std::endl;
            xdmf_file << "<DataItem Dimensions=\"3 4\" Format=\"XML\">" << std::endl;
            xdmf_file << offset[3] << " " << 0 << " " << 0 << " " << 0 << std::endl;
            xdmf_file << 1 << " " << 1 << " " << 1 << " " << 1 << std::endl;
            xdmf_file << 1 << " " << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << std::endl;
            xdmf_file << "</DataItem>" << std::endl;
            xdmf_file << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" 
                      << finalstep << " " << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << "\">" << std::endl;
            xdmf_file << trim(adjustl(filename)) << ".h5:/" << trim(adjustl(dsetname)) << std::endl;
            xdmf_file << "</DataItem>" << std::endl;
            xdmf_file << "</DataItem>" << std::endl;
            xdmf_file << "</Attribute>" << std::endl;
            xdmf_file << "</Grid>" << std::endl;
        }
    }

    void closeh5file(int finalstep, const std::vector<RKIND>& att) {
        int error = 0;
        int rank = 1;
        std::string dsetname = "Time";

        DATA_dims.clear();
        valor3d_dims.clear();
        offset.clear();

        HDF5::h5dclose_f(dset_id, error);

        dsetname = "Time";
        rank = 1;
        DATA_dims.resize(rank);
        DATA_dims[0] = finalstep;

        HDF5::h5screate_simple_f(rank, DATA_dims, dspace_id, error);
        
        // Write time data
        HDF5::h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error);
        HDF5::h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, att, DATA_dims, error);
        HDF5::h5dclose_f(dset_id, error);

        HDF5::h5sclose_f(slice2D_id, error);
        HDF5::h5sclose_f(dspace_id, error);

        HDF5::h5fclose_f(file_id, error);
        HDF5::h5close_f(error);

        if (xdmf_file.is_open()) {
            xdmf_file << "</Grid>" << std::endl;
            xdmf_file << "</Domain>" << std::endl;
            xdmf_file << "</Xdmf>" << std::endl;
            xdmf_file.close();
        }
    }

    void createh5filefromsinglebin(const std::string& filename, bool vtkindex) {
        int myunit = 0; // Simulated unit
        int fieldob = 0;
        int pasadas = 0;
        int pasadastotales = 0;
        
        std::string fichin;
        std::vector<RKIND> att;
        std::vector<std::vector<std::vector<std::vector<RKIND>>>> valor3d;
        bool SGGObservationiiTimeDomain = false;
        
        int minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs;
        int minZabs_primero, minYabs_primero, minXabs_primero, finalstep, indi, i1, j1, k1;
        RKIND linez_minZabs_primero, liney_minYabs_primero, linex_minXabs_primero;
        RKIND dz_minZabs, dy_minYabs, dx_minXabs;
        std::string dubuf;

        // Process filename
        std::string base_filename = filename;
        size_t pos = base_filename.find(".h5bin");
        if (pos != std::string::npos) {
            base_filename = base_filename.substr(0, pos);
        }
        base_filename = trim(adjustl(base_filename));

        // Open binary file
        std::ifstream bin_file(base_filename + ".h5bin", std::ios::binary);
        if (!bin_file.is_open()) {
            std::cerr << "Error opening file: " << base_filename << ".h5bin" << std::endl;
            return;
        }

        bin_file >> finalstep >> minXabs >> maxXabs >> minYabs >> maxYabs >> minZabs >> maxZabs >> fieldob >> SGGObservationiiTimeDomain >> pasadastotales;

        // Allocate valor3d
        // Fortran: allocate(valor3d(minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1))
        // C++: Use 0-based indexing internally, but map to Fortran indices if needed.
        // Here we allocate a vector of vectors of vectors of vectors.
        int dim1 = maxXabs - minXabs + 1;
        int dim2 = maxYabs - minYabs + 1;
        int dim3 = maxZabs - minZabs + 1;
        int dim4 = 1;
        
        valor3d.resize(dim1, std::vector<std::vector<std::vector<RKIND>>>(dim2, std::vector<std::vector<RKIND>>(dim3, std::vector<RKIND>(dim4, 0.0))));
        
        att.resize(finalstep + 1, 0.0); // 1-based indexing in Fortran, so size is finalstep + 1

        for (pasadas = 1; pasadas <= pasadastotales; ++pasadas) {
            if (SGGObservationiiTimeDomain) {
                if (pasadas == 1) {
                    fichin = trim(adjustl(base_filename)) + "_time";
                } else {
                    std::cerr << "Buggy error in valor3d." << std::endl;
                    return;
                }
            } else {
                if (pasadas == 1) {
                    fichin = trim(adjustl(base_filename)) + "_mod";
                } else if (pasadas == 2) {
                    fichin = trim(adjustl(base_filename)) + "_phase";
                } else {
                    std::cerr << "Buggy error in valor3d." << std::endl;
                    return;
                }
            }

            // Check if we should write
            // Assuming iMEC and iMHC are constants defined elsewhere. 
            // For translation, we assume they are not equal to fieldob or pasadas != 2.
            bool skip_write = false; // Placeholder for ((fieldob == iMEC) || (fieldob == iMHC)) && (pasadas == 2)
            
            if (!skip_write) {
                openh5file(fichin, finalstep, minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs);
            }

            // Initialize valor3d and att
            for (int i = 0; i < dim1; ++i)
                for (int j = 0; j < dim2; ++j)
                    for (int k = 0; k < dim3; ++k)
                        for (int l = 0; l < dim4; ++l)
                            valor3d[i][j][k][l] = 0.0;
            
            for (int i = 0; i <= finalstep; ++i)
                att[i] = 0.0;

            for (indi = 1; indi <= finalstep; ++indi) {
                bin_file >> minZabs_primero >> minYabs_primero >> minXabs_primero;
                bin_file >> linez_minZabs_primero >> liney_minYabs_primero >> linex_minXabs_primero;
                bin_file >> dz_minZabs >> dy_minYabs >> dx_minXabs;
                bin_file >> att[indi];

                std::ostringstream oss;
                oss << " ----> .xdmf file " << att[indi] << "(" << indi << "/" << finalstep << ")";
                dubuf = oss.str();
                std::cout << trim(adjustl(dubuf)) << std::endl;

                for (k1 = minZabs; k1 <= maxZabs; ++k1) {
                    for (j1 = minYabs; j1 <= maxYabs; ++j1) {
                        // Read row i1 from minXabs to maxXabs
                        // valor3d is indexed [i1-minXabs][j1-minYabs][k1-minZabs][0]
                        for (i1 = minXabs; i1 <= maxXabs; ++i1) {
                            bin_file >> valor3d[i1 - minXabs][j1 - minYabs][k1 - minZabs][0];
                        }
                    }
                }

                if (!skip_write) {
                    writeh5file(fichin, valor3d, indi, att[indi], 
                                minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs,
                                linez_minZabs_primero, liney_minYabs_primero, linex_minXabs_primero,
                                dz_minZabs, dy_minYabs, dx_minXabs,
                                minZabs_primero, minYabs_primero, minXabs_primero, finalstep, vtkindex);
                }
            }

            if (!skip_write) {
                closeh5file(finalstep, att);
            }
        }

        bin_file.close();
        valor3d.clear();
        att.clear();
    }

} // namespace xdmf_h5_m