#ifndef TEST_MTLN_FHASH_H
#define TEST_MTLN_FHASH_H

#include <gtest/gtest.h>
#include <any>
#include <string>

#include "fhash_m.h"

namespace mtln_fhash_test {

struct mtl_bundle_t {
    std::string name;
    int number_of_conductors = 0;
};

} // namespace mtln_fhash_test

TEST(mtln, fhash) {
    using mtln_fhash_test::mtl_bundle_t;

    fhash_m::fhash_tbl_t tbl;
    mtl_bundle_t bundle;
    bundle.name = "bundle1";
    bundle.number_of_conductors = 5;

    tbl.set(fhash_m::key(bundle.name), std::any(bundle));

    std::any raw;
    int stat = 0;
    tbl.get_raw(fhash_m::key(bundle.name), raw, &stat);
    auto bundle_from_hash_1 = std::any_cast<mtl_bundle_t>(raw);
    bundle_from_hash_1.number_of_conductors = 100;
    tbl.set(fhash_m::key(bundle_from_hash_1.name), std::any(bundle_from_hash_1));

    tbl.get_raw(fhash_m::key(bundle.name), raw, &stat);
    const auto bundle_from_hash_2 = std::any_cast<mtl_bundle_t>(raw);

    EXPECT_EQ(bundle_from_hash_2.number_of_conductors, 100);
}

#endif
