#ifndef TEST_IDCHILDTABLE_FHASH_H
#define TEST_IDCHILDTABLE_FHASH_H

#include <gtest/gtest.h>
#include <any>
#include <string>
#include "fhash_m.h"

TEST(smbjson_cpp, idchildtable_fhash) {
    fhash_m::fhash_tbl_t tbl;
    constexpr int expectedValue = 10;

    tbl.set(fhash_m::key(std::string("my_key_1")), std::any(expectedValue));
    tbl.set(fhash_m::key(std::string("my_key_2")), std::any(1.0));
    tbl.set(fhash_m::key(std::string("a string value")), std::any(std::string("a string value")));
    tbl.set(fhash_m::key(std::vector<int>{1, 2, 3, 4, 5}), std::any(false));

    int val = 0;
    int stat = 0;
    std::any raw;
    tbl.get_raw(fhash_m::key(std::string("my_key_1")), raw, &stat);
    val = std::any_cast<int>(raw);
    EXPECT_EQ(stat, 0);
    EXPECT_EQ(val, expectedValue);
}

#endif
