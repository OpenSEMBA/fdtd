#ifndef TEST_BULK_CURRENT_PROBE_H
#define TEST_BULK_CURRENT_PROBE_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <filesystem>
#include <string>

namespace bulk_current_probe_test {

inline std::string casePath(const char* name) {
    return (std::filesystem::path("testData") / "cases" /
            "multipleAssigments" / name).string();
}

} // namespace bulk_current_probe_test

TEST(BulkCurrentProbe, SurfaceImpedanceWritesFortranNamedProbe) {
    const std::string json =
        bulk_current_probe_test::casePath("multipleSurfaceImpedance.fdtd.json");
    const std::string expected =
        "multipleSurfaceImpedance.fdtd_BulkProbeEntry_Jz_49_49_45__50_50_45.dat";
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_bulk_current_probe_output(
        json, expected, 5);
    EXPECT_EQ(err, 0) << "bulkCurrent probe output failed with code " << err;
}

#endif
