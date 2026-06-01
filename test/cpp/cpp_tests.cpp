#include <gtest/gtest.h>

#ifdef CompileWithMPI
#include <mpi.h>
#endif

#include "test_cells.h"
#include "test_idchildtable.h"
#include "test_mesh.h"
#include "test_idchildtable_fhash.h"

#include "test_smbjson_parser.h"
#include "test_smbjson_parser_mesh.h"
#include "test_smbjson_read.h"
#ifdef CompileWithMTLN
#include "test_smbjson_read_mtln.h"
#endif

#include "test_conformal_geometry.h"
#include "test_conformal_cell_map.h"
#include "test_conformal_filling.h"

#include "test_rotate.h"

#ifndef CompileWithMPI
#include "test_observation.h"
#endif

#include "test_preprocess_geom.h"

#include "test_system_init_solver.h"
#include "test_planewave_evolucion.h"
#include "test_planewave_init.h"
#include "test_planewave_strict.h"
#include "test_planewave_tfsf.h"
#include "test_bordersmur.h"
#include "test_mpi_one_axis.h"
#include "test_pml_boundary.h"
#include "test_planewave_pw_in_box.h"
#include "test_holland_wire.h"
#include "test_bulk_current_probe.h"
#include "test_maloney_nostoch.h"
#include "test_maloney_missing.h"
#include "test_xdmf_h5.h"

#ifdef CompileWithMTLN
#include "test_mtln_mtl.h"
#include "test_mtln_types.h"
#include "test_mtln_math.h"
#include "test_mtln_multipolar.h"
#include "test_mtln_dispersive.h"
#include "test_mtln_spice.h"
#include "test_mtln_preprocess.h"
#include "test_mtln_fhash.h"
#ifdef CompileWithMPI
#include "test_mtln_mpi_bundle.h"
#endif
#endif

int main(int argc, char** argv) {
#ifdef CompileWithMPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    const bool finalize_mpi = (mpi_initialized == 0);
    if (finalize_mpi) {
        MPI_Init(&argc, &argv);
    }
#endif
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef CompileWithMPI
    int mpi_finalized = 0;
    MPI_Finalized(&mpi_finalized);
    if (finalize_mpi && mpi_finalized == 0) {
        MPI_Finalize();
    }
#endif
    return result;
}
