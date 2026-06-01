#ifndef TEST_SMBJSON_PARSER_H
#define TEST_SMBJSON_PARSER_H

#include "test_smbjson_helpers.h"
#include "expected/test_read_background.h"

TEST(smbjson_cpp, parser_ctor) {
    smbjson::parser_t parser(testDataPath("planewave.fdtd.json"));
    EXPECT_TRUE(parser.isInitialized);
}

TEST(smbjson_cpp, read_background_defaults) {
    smbjson::parser_t parser(testDataPath("planewave.fdtd.json"));
    auto pr = parser.readProblemDescription();
    expectBackgroundDefaults(pr);
}

TEST(smbjson_cpp, read_background_set) {
    smbjson::parser_t parser(testDataPath("background.fdtd.json"));
    auto pr = parser.readProblemDescription();
    expectBackgroundSet(pr);
}

TEST(smbjson_cpp, read_planewave_empty_elementids) {
    smbjson::parser_t parser(testDataPath("planewave_empty_elementids.fdtd.json"));
    EXPECT_THROW(parser.readProblemDescription(), std::runtime_error);
}

#endif // TEST_SMBJSON_PARSER_H
