#include <gtest/gtest.h>

#include "test_cells.h"
#include "test_idchildtable.h"
#include "test_mesh.h"

#include "test_smbjson_parser.h"
#include "test_smbjson_read.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
