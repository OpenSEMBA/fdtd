#pragma once

#if defined(CompileWithMPI) && defined(CompileWithMTLN)

#include <mpi.h>

#include <gtest/gtest.h>

#include "mtl_bundle_m.h"

extern MPI_Comm SUBCOMM_MPI;

namespace {

constexpr int kCommSend = 1;
constexpr int kCommRecv = -1;
constexpr int kCommField = 1;
constexpr int kCommV = 2;

void requireTwoMpiRanks() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    ASSERT_NE(initialized, 0);

    int size = 1;
    MPI_Comm_size(SUBCOMM_MPI, &size);
    if (size != 2) {
        GTEST_SKIP() << "requires exactly two MPI ranks";
    }
}

} // namespace

TEST(MtlnMpiBundle, ExchangesVoltageBoundaryValues) {
    requireTwoMpiRanks();

    int rank = 0;
    MPI_Comm_rank(SUBCOMM_MPI, &rank);

    mtl_bundle_m::mtl_bundle_t bundle;
    bundle.v = {
        {10.0 + rank, 12.5 + rank, 15.0 + rank},
        {20.0 + rank, 22.5 + rank, 25.0 + rank},
    };

    mtl_m::communicator_t comm;
    comm.comm_type = kCommV;
    comm.comm_task = (rank == 0) ? kCommSend : kCommRecv;
    comm.delta_rank = (rank == 0) ? 1 : -1;
    comm.v_index = (rank == 0) ? 1 : 0;
    bundle.mpi_comm.comms = {comm};

    bundle.Comm_MPI_V();
    MPI_Barrier(SUBCOMM_MPI);

    if (rank == 1) {
        EXPECT_DOUBLE_EQ(bundle.v[0][0], 12.5);
        EXPECT_DOUBLE_EQ(bundle.v[1][0], 22.5);
    }
}

TEST(MtlnMpiBundle, ExchangesExternalFieldValues) {
    requireTwoMpiRanks();

    int rank = 0;
    MPI_Comm_rank(SUBCOMM_MPI, &rank);

    double field_send = 123.25;
    double field_recv = -7.0;

    mtl_bundle_m::mtl_bundle_t bundle;
    bundle.external_field_segments.resize(2);
    bundle.external_field_segments[0].field = &field_send;
    bundle.external_field_segments[1].field = &field_recv;

    mtl_m::communicator_t comm;
    comm.comm_type = kCommField;
    comm.comm_task = (rank == 0) ? kCommSend : kCommRecv;
    comm.delta_rank = (rank == 0) ? 1 : -1;
    comm.field_index = (rank == 0) ? 0 : 1;
    bundle.mpi_comm.comms = {comm};

    bundle.Comm_MPI_Fields();
    MPI_Barrier(SUBCOMM_MPI);

    if (rank == 1) {
        EXPECT_DOUBLE_EQ(field_recv, 123.25);
    }
}

#endif // defined(CompileWithMPI) && defined(CompileWithMTLN)
