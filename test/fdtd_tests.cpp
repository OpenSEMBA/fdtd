#include <gtest/gtest.h>

#ifdef CompileWithMPI
#include <mpi.h>
#endif

#ifdef CompileWithMTLN
    #include "mtln/mtln_tests.h"
    //#include "system/system_tests.h"
#endif
#ifdef CompileWithSMBJSON
    #include "smbjson/smbjson_tests.h"
    #include "rotate/rotate_tests.h"
    #ifdef CompileWithNewOutputModule
        #include "unit/output/output_tests.h"
        #include "unit/output/vtkAPI_tests.h"
    #endif
#endif
#if !defined(CompileWithMPI) && !defined(CompileWithNewOutputModule)
    #include "observation/observation_tests.h"
#endif
#include "conformal/conformal_tests.h"
#include "unit/preprocess/preprocess_tests.h"

int main(int argc, char **argv) {
#ifdef CompileWithMPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(&argc, &argv);
    }
#endif
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef CompileWithMPI
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) {
        MPI_Finalize();
    }
#endif
    return result;
}
