#include <gtest/gtest.h>

#ifdef CompileWithMTLN
    #include "mtln/mtln_tests.h"
    //#include "system/system_tests.h"
#endif
#ifdef CompileWithSMBJSON
    #include "smbjson/smbjson_tests.h"
    #include "rotate/rotate_tests.h"
    #include "vtk/vtk_tests.h"
#endif
#ifndef CompileWithMPI
    #include "observation/observation_tests.h"
#endif
#include "conformal/conformal_tests.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
