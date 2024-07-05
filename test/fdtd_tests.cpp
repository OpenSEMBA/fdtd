#include <gtest/gtest.h>

#ifdef CompileWithMTLN
    #include "mtln/mtln_tests.h"
    #include "system/system_tests.h"
#endif

#include "smbjson/smbjson_tests.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
