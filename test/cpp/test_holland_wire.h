#ifndef TEST_HOLLAND_WIRE_H
#define TEST_HOLLAND_WIRE_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <filesystem>
#include <string>

namespace holland_wire_test {

inline std::string casePath(const char* name) {
    return (std::filesystem::path("testData") / "cases" / "holland" / name).string();
}

} // namespace holland_wire_test

#ifndef CompileWithMTLN
TEST(HollandWire, ShortRunWritesFortranNamedCurrentProbe) {
    const std::string json = holland_wire_test::casePath("holland1981.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));
    const int err = SEMBA_FDTD_m::SEMBA_FDTD_test::test_run_holland_probe_output(json, 10);
    EXPECT_EQ(err, 0) << "Holland short-run probe output failed with code " << err;
}
#endif

#endif
