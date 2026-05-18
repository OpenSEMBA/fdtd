#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdint>

// Note: The original code uses HDF5. In C++, you would typically use the HDF5 C++ API or C API.
// For this translation, we assume a wrapper or direct C API calls are available.
// We also assume FDETYPES_m provides RKIND and BUFSIZE.
// Since we cannot include the actual Fortran modules, we define placeholders/constants.

// Placeholder for constants from FDETYPES_m
constexpr int BUFSIZE = 256;
constexpr int RKIND = 4; // Assuming 4-byte real

// Placeholder for HDF5 types and functions if not using the actual HDF5 headers
// In a real scenario, you would #include <hdf5.h>
// For the purpose of this translation, we assume the following are available:
// - HID_T, HSIZE_T
// - H5F_ACC_TRUNC_F, H5T_NATIVE_REAL, H5S_SELECT_SET_F
// - h5open_f, h5fcreate_f, h5screate_simple_f, h5dcreate_f, h5dwrite_f, h5dclose_f, h5sclose_f, h5fclose_f, h5close_f

// If using actual HDF5 C API, the signatures would be different. 
// Here we mimic the Fortran interface style for clarity, but in practice, you'd use the C++ HDF5 bindings or C API.

namespace snapxdmf_m {

    void WRITE_XDMFSNAP(int ninstant, const std::string& filename, int minXabs, int maxXabs, int minYabs, int maxYabs, int minZabs, int maxZabs, const std::vector<std::vector<std::vector<std::vector<float>>>>& valor3D) {
#ifdef CompileWithHDF
        // Note: Fortran array valor3D is 1-based and has shape (minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1:1)
        // In C++, we assume valor3D is passed as a 4D vector with 0-based indexing, but we need to map it correctly.
        // The Fortran array is contiguous in memory in the order of the last index varying fastest? 
        // Actually, Fortran is column-major. The declaration is dimension(minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1:1).
        // So the first index varies slowest, last index varies fastest.
        // In C++, std::vector of vectors is row-major. We need to be careful with data layout.
        // For simplicity, we assume the data is passed in a way that matches the expected memory layout for HDF5, 
        // or we copy it. Given the complexity, we'll assume a flat vector or a specific layout.
        // However, the original code passes valor3D directly to h5dwrite_f. 
        // Let's assume valor3D is a 4D vector where the indices are [i][j][k][l] corresponding to [minXabs:maxXabs][minYabs:maxYabs][minZabs:maxZabs][1:1].
        // To match Fortran's column-major order, we might need to transpose or reorder.
        // But for the sake of this translation, we will pass a pointer to the data, assuming the caller handles the layout.
        // We'll use a flat vector for simplicity in the C++ signature, but the original is 4D.
        // Let's stick to the 4D vector for type safety, but note that HDF5 expects a contiguous block.
        
        // Extract dimensions
        int finalstep = 1;
        int rank = 4;
        
        // DATA_dims for HDF5
        std::vector<hsize_t> DATA_dims(rank);
        DATA_dims[0] = maxXabs - minXabs + 1;
        DATA_dims[1] = maxYabs - minYabs + 1;
        DATA_dims[2] = maxZabs - minZabs + 1;
        DATA_dims[3] = finalstep;
        
        // valor3d_dims for hyperslab
        std::vector<hsize_t> valor3d_dims(rank);
        valor3d_dims[0] = DATA_dims[0];
        valor3d_dims[1] = DATA_dims[1];
        valor3d_dims[2] = DATA_dims[2];
        valor3d_dims[3] = 1;
        
        // offset for hyperslab
        std::vector<hsize_t> offset(rank);
        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;
        offset[3] = 0;
        
        std::string dsetname = "data";
        
        int error = 0;
        hid_t file_id;
        hid_t dset_id;
        hid_t dspace_id;
        hid_t slice2D_id;
        
        // Open HDF5 file
        // h5open_f is called once per program usually, but here it's called every time.
        // We'll call it as in the original.
        // Note: h5open_f is a Fortran interface. In C++, we might use H5open().
        // Assuming a wrapper function h5open_f exists that calls H5open().
        h5open_f(error);
        
        std::string h5filename = filename + ".h5";
        h5fcreate_f(h5filename.c_str(), H5F_ACC_TRUNC_F, file_id, error);
        
        h5screate_simple_f(rank, DATA_dims.data(), dspace_id, error);
        h5screate_simple_f(rank, valor3d_dims.data(), slice2D_id, error);
        
        h5dcreate_f(file_id, dsetname.c_str(), H5T_NATIVE_REAL, dspace_id, dset_id, error);
        
        // Open XDMF file
        std::string xdmffilename = filename + ".xdmf";
        std::ofstream xdmf_file(xdmffilename);
        if (!xdmf_file.is_open()) {
            std::cerr << "Error opening XDMF file: " << xdmffilename << std::endl;
            return;
        }
        
        xdmf_file << "<Xdmf>" << std::endl;
        xdmf_file << "<Domain>" << std::endl;
        xdmf_file << "<Grid Name=\"GridTime\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
        
        int indi = finalstep;
        std::vector<float> att(finalstep);
        att[0] = static_cast<float>(ninstant);
        
        char charc[BUFSIZE];
        snprintf(charc, BUFSIZE, "%d", indi);
        
        // Select hyperslab
        h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset.data(), valor3d_dims.data(), error);
        
        // Write data
        // Note: valor3D is a 4D vector. We need to pass a pointer to the first element.
        // The layout must match what HDF5 expects. Fortran is column-major.
        // If valor3D is stored in row-major (C++ default), we might need to transpose.
        // For now, we assume the data is laid out correctly for HDF5 or we pass a flattened version.
        // To be safe, we can create a flat vector in Fortran order.
        // However, for simplicity, we pass the data directly, assuming the caller ensures correct layout.
        // In practice, you might need to copy and transpose.
        
        // Flatten valor3D to a 1D array in Fortran order
        // Dimensions: (maxXabs-minXabs+1) x (maxYabs-minYabs+1) x (maxZabs-minZabs+1) x 1
        size_t total_size = DATA_dims[0] * DATA_dims[1] * DATA_dims[2] * DATA_dims[3];
        std::vector<float> valor3d_flat(total_size);
        
        // Copy data from 4D vector to 1D flat array in Fortran order
        // Fortran order: i varies slowest, l varies fastest
        // valor3D[i][j][k][l] -> flat index = i + j*dim0 + k*dim0*dim1 + l*dim0*dim1*dim2
        for (int i = 0; i < DATA_dims[0]; ++i) {
            for (int j = 0; j < DATA_dims[1]; ++j) {
                for (int k = 0; k < DATA_dims[2]; ++k) {
                    for (int l = 0; l < DATA_dims[3]; ++l) {
                        // Adjust indices for 0-based C++ vector
                        // The original Fortran array is from minXabs to maxXabs, etc.
                        // So the index in the vector should be (i + minXabs), etc.
                        size_t flat_idx = i + j * DATA_dims[0] + k * DATA_dims[0] * DATA_dims[1] + l * DATA_dims[0] * DATA_dims[1] * DATA_dims[2];
                        valor3d_flat[flat_idx] = valor3D[i + minXabs][j + minYabs][k + minZabs][l];
                    }
                }
            }
        }
        
        h5dwrite_f(dset_id, H5T_NATIVE_REAL, valor3d_flat.data(), valor3d_dims.data(), error, slice2D_id, dspace_id);
        
        // Write XDMF content
        xdmf_file << "<Grid Name=\"IntGrid\" GridType=\"Uniform\"  CollectionType=\"Spatial\">>" << std::endl;
        xdmf_file << "<Time Value=\"" << charc << "\" />" << std::endl;
        xdmf_file << "<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" 
                  << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << "\">" << std::endl;
        xdmf_file << "</Topology>" << std::endl;
        xdmf_file << "<Geometry Type=\"ORIGIN_DXDYDZ\">" << std::endl;
        xdmf_file << "<DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
        xdmf_file << "<DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
        xdmf_file << "1.0 1.0 1.0" << std::endl;
        xdmf_file << "</DataItem>" << std::endl;
        xdmf_file << "</Geometry>" << std::endl;
        xdmf_file << "<Attribute Name=\"IntValues\" Center=\"Node\">" << std::endl;
        xdmf_file << "<DataItem ItemType=\"HyperSlab\" Dimensions=\"" 
                  << 1 << " " << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << "\" Format=\"XML\">" << std::endl;
        xdmf_file << "<DataItem Dimensions=\"3 4\" Format=\"XML\">" << std::endl;
        xdmf_file << offset[3] << " " << 0 << " " << 0 << " " << 0 << std::endl;
        xdmf_file << 1 << " " << 1 << " " << 1 << " " << 1 << std::endl;
        xdmf_file << 1 << " " << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << std::endl;
        xdmf_file << "</DataItem>" << std::endl;
        xdmf_file << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" 
                  << finalstep << " " << DATA_dims[2] << " " << DATA_dims[1] << " " << DATA_dims[0] << "\">" << std::endl;
        xdmf_file << filename << ".h5:/" << dsetname << std::endl;
        xdmf_file << "</DataItem>" << std::endl;
        xdmf_file << "</DataItem>" << std::endl;
        xdmf_file << "</Attribute>" << std::endl;
        xdmf_file << "</Grid>" << std::endl;
        
        // Close dataset and dataspace
        h5dclose_f(dset_id, error);
        h5sclose_f(slice2D_id, error);
        h5sclose_f(dspace_id, error);
        
        // Create time dataset
        dsetname = "Time";
        rank = 1;
        std::vector<hsize_t> DATA_dims_time(rank);
        DATA_dims_time[0] = finalstep;
        
        h5screate_simple_f(rank, DATA_dims_time.data(), dspace_id, error);
        h5dcreate_f(file_id, dsetname.c_str(), H5T_NATIVE_REAL, dspace_id, dset_id, error);
        h5dwrite_f(dset_id, H5T_NATIVE_REAL, att.data(), DATA_dims_time.data(), error);
        
        h5dclose_f(dset_id, error);
        h5sclose_f(dspace_id, error);
        
        xdmf_file << "</Grid>" << std::endl;
        xdmf_file << "</Domain>" << std::endl;
        xdmf_file << "</Xdmf>" << std::endl;
        
        xdmf_file.close();
        
        h5fclose_f(file_id, error);
        h5close_f(error);
        
#endif
    }

} // namespace snapxdmf_m