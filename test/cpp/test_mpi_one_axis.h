#pragma once

#include <gtest/gtest.h>

#include <stdexcept>

#ifdef CompileWithMPI
#include <mpi.h>
#endif

#include "semba_fdtd.h"

#ifdef CompileWithMPI
namespace {

bool requireAtLeastTwoMpiRanksForFieldExchange() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized == 0) {
        ADD_FAILURE() << "MPI is not initialized";
        return false;
    }

    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2) {
        return false;
    }
    return true;
}

} // namespace
#endif

TEST(MpiOneAxis, ParsesFortranMpidirFlags) {
    using namespace SEMBA_FDTD_m::SEMBA_FDTD_test;

    EXPECT_EQ(test_mpi_axis_from_flags("-i case.fdtd.json"), 3);
    EXPECT_EQ(test_mpi_axis_from_flags("-i case.fdtd.json -mpidir x"), 1);
    EXPECT_EQ(test_mpi_axis_from_flags("-mpidir y -i case.fdtd.json"), 2);
    EXPECT_EQ(test_mpi_axis_from_flags("-mpidir=z -i case.fdtd.json"), 3);
    EXPECT_EQ(test_mpi_axis_from_flags("-mpidir x -mpidir x"), 1);
    EXPECT_THROW(test_mpi_axis_from_flags("-mpidir x -mpidir y"),
                 std::runtime_error);
}

TEST(MpiOneAxis, SplitsOnlyAlongCanonicalAxis) {
    using namespace SEMBA_FDTD_m::SEMBA_FDTD_test;

    const auto slices = test_mpi_one_axis_slices(16, 4, 0, 0, -1, 1);
    ASSERT_EQ(slices.size(), 4U);
    for (int rank = 0; rank < 4; ++rank) {
        EXPECT_EQ(slices[static_cast<size_t>(rank)].axis, 1);
        EXPECT_EQ(slices[static_cast<size_t>(rank)].com, 4 * rank);
        EXPECT_EQ(slices[static_cast<size_t>(rank)].fin, 4 * (rank + 1));
        EXPECT_EQ(slices[static_cast<size_t>(rank)].allocZI,
                  slices[static_cast<size_t>(rank)].sweepZI - 1);
        EXPECT_EQ(slices[static_cast<size_t>(rank)].allocZE,
                  slices[static_cast<size_t>(rank)].sweepZE + 1);
    }

    EXPECT_TRUE(slices.front().physicalDown);
    EXPECT_FALSE(slices.front().physicalUp);
    EXPECT_FALSE(slices[1].physicalDown);
    EXPECT_FALSE(slices[1].physicalUp);
    EXPECT_FALSE(slices.back().physicalDown);
    EXPECT_TRUE(slices.back().physicalUp);
    EXPECT_EQ(slices[0].sweepZE, 3);
    EXPECT_EQ(slices[1].sweepZE, 7);
    EXPECT_EQ(slices[2].sweepZE, 11);
    EXPECT_EQ(slices[3].sweepZE, 16);
}

TEST(MpiOneAxis, AppliesFortranPmlCpuWeighting) {
    using namespace SEMBA_FDTD_m::SEMBA_FDTD_test;

    const auto slices = test_mpi_one_axis_slices(20, 2, 4, 0, -1, 3);
    ASSERT_EQ(slices.size(), 2U);
    EXPECT_EQ(slices[0].com, 0);
    EXPECT_EQ(slices[0].fin, 8);
    EXPECT_EQ(slices[1].com, 8);
    EXPECT_EQ(slices[1].fin, 20);
    EXPECT_TRUE(slices[0].pmlDown);
    EXPECT_FALSE(slices[0].pmlUp);
    EXPECT_FALSE(slices[1].pmlDown);
    EXPECT_FALSE(slices[1].pmlUp);
}

TEST(MpiOneAxis, SupportsTwoRankForcedCut) {
    using namespace SEMBA_FDTD_m::SEMBA_FDTD_test;

    const auto slices = test_mpi_one_axis_slices(18, 2, 0, 0, 5, 3);
    ASSERT_EQ(slices.size(), 2U);
    EXPECT_EQ(slices[0].fin, 5);
    EXPECT_EQ(slices[1].com, 5);
    EXPECT_THROW(test_mpi_one_axis_slices(18, 3, 0, 0, 5, 3),
                 std::runtime_error);
}

#ifdef CompileWithMPI
TEST(MpiOneAxis, LeavesElectricGhostPlanesUnchangedWithoutExtraInfo) {
    using namespace SEMBA_FDTD_m::SEMBA_FDTD_test;
    if (!requireAtLeastTwoMpiRanksForFieldExchange()) {
        GTEST_SKIP() << "requires at least two MPI ranks";
    }

    EXPECT_EQ(test_mpi_exchange_electric_ghost_planes(1), 0);
    EXPECT_EQ(test_mpi_exchange_electric_ghost_planes(2), 0);
    EXPECT_EQ(test_mpi_exchange_electric_ghost_planes(3), 0);
}

TEST(MpiOneAxis, ExchangesMagneticGhostPlanesOnEveryAxis) {
    using namespace SEMBA_FDTD_m::SEMBA_FDTD_test;
    if (!requireAtLeastTwoMpiRanksForFieldExchange()) {
        GTEST_SKIP() << "requires at least two MPI ranks";
    }

    EXPECT_EQ(test_mpi_exchange_magnetic_ghost_planes(1), 0);
    EXPECT_EQ(test_mpi_exchange_magnetic_ghost_planes(2), 0);
    EXPECT_EQ(test_mpi_exchange_magnetic_ghost_planes(3), 0);
}
#endif
