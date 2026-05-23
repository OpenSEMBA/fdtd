#ifndef TEST_SMBJSON_PARSER_H
#define TEST_SMBJSON_PARSER_H

#include "test_smbjson_helpers.h"

TEST(smbjson_cpp, parser_ctor) {
    smbjson::parser_t parser(testDataPath("planewave.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

TEST(smbjson_cpp, background_defaults) {
    smbjson::parser_t parser(testDataPath("planewave.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

TEST(smbjson_cpp, background_set) {
    smbjson::parser_t parser(testDataPath("background.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

TEST(smbjson_cpp, read_planewave_loads) {
    smbjson::parser_t parser(testDataPath("planewave.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

#endif // TEST_SMBJSON_PARSER_H
