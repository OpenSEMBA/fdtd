#ifndef TEST_XDMF_H5_H
#define TEST_XDMF_H5_H

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#ifdef SEMBA_CPP_ENABLE_HDF5
#include <hdf5.h>
#endif

#include "xdmf_h5.h"

TEST(XdmfH5, WritesFourDimensionalSlab) {
#ifndef SEMBA_CPP_ENABLE_HDF5
    GTEST_SKIP() << "HDF5 not enabled for this build";
#else
    const std::string stem = "test_xdmf_h5_unit";
    const int nx = 2;
    const int ny = 2;
    const int nz = 2;
    const int nsteps = 3;
    const int minX = 1;
    const int maxX = 2;
    const int minY = 1;
    const int maxY = 2;
    const int minZ = 1;
    const int maxZ = 2;

    xdmf_h5_m::openh5file(stem, nsteps, minX, maxX, minY, maxY, minZ, maxZ);

    std::vector<float> slab(static_cast<size_t>(nx * ny * nz), 0.0f);
    for (int step = 1; step <= nsteps; ++step) {
        const float value = static_cast<float>(step) * 0.25f;
        std::fill(slab.begin(), slab.end(), value);
        xdmf_h5_m::writeh5file(stem, slab.data(), nx, ny, nz, step, static_cast<double>(step),
                               minX, maxX, minY, maxY, minZ, maxZ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                               minZ, minY, minX, nsteps, true);
    }
    xdmf_h5_m::closeh5file(nsteps, {1.0, 2.0, 3.0});

    hid_t file = H5Fopen((stem + ".h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(file, 0);
    hid_t dataset = H5Dopen2(file, "data", H5P_DEFAULT);
    ASSERT_GE(dataset, 0);

    hid_t space = H5Dget_space(dataset);
    hsize_t dims[4] = {0, 0, 0, 0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(space), 4);
    ASSERT_EQ(H5Sget_simple_extent_dims(space, dims, nullptr), 4);
    EXPECT_EQ(dims[0], static_cast<hsize_t>(nsteps));
    EXPECT_EQ(dims[1], static_cast<hsize_t>(nz));
    EXPECT_EQ(dims[2], static_cast<hsize_t>(ny));
    EXPECT_EQ(dims[3], static_cast<hsize_t>(nx));

    std::vector<float> readback(static_cast<size_t>(nx * ny * nz), 0.0f);
    hsize_t offset[4] = {1, 0, 0, 0};
    hsize_t count[4] = {1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny),
                        static_cast<hsize_t>(nx)};
    hid_t memspace = H5Screate_simple(4, count, nullptr);
    H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, space, H5P_DEFAULT, readback.data());

    EXPECT_NEAR(readback.front(), 0.5f, 1e-6f);

    H5Sclose(memspace);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
#endif
}

#endif
