#ifndef TEST_SMBJSON_READ_H
#define TEST_SMBJSON_READ_H

#include "test_smbjson_helpers.h"

TEST(smbjson_cpp, read_dielectricslab_loads) {
    smbjson::parser_t parser(testDataPath("dielectric_slab.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

TEST(smbjson_cpp, read_thinSlot_loads) {
    smbjson::parser_t parser(testDataPath("thinSlot.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

TEST(smbjson_cpp, read_lumped_fixture_loads) {
    smbjson::parser_t parser(testDataPath("lumped_fixture.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

#endif // TEST_SMBJSON_READ_H
