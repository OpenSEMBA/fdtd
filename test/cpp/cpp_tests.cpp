#include <gtest/gtest.h>

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
#include "test_planewave_tfsf.h"
#include "test_bordersmur.h"
#include "test_planewave_pw_in_box.h"

#ifdef CompileWithMTLN
#include "test_mtln_mtl.h"
#include "test_mtln_types.h"
#include "test_mtln_math.h"
#include "test_mtln_multipolar.h"
#include "test_mtln_dispersive.h"
#include "test_mtln_spice.h"
#include "test_mtln_preprocess.h"
#include "test_mtln_fhash.h"
#endif

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
