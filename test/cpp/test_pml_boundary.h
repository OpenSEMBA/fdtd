#ifndef TEST_PML_BOUNDARY_H
#define TEST_PML_BOUNDARY_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <filesystem>
#include <string>

namespace pml_boundary_test {

inline std::string casePath(const std::string& folder,
                            const std::string& filename) {
    return (std::filesystem::path("testData") / "cases" / folder / filename)
        .string();
}

} // namespace pml_boundary_test

TEST(PmlBoundary, PmlAllDoesNotEnableMur) {
    const std::string json = pml_boundary_test::casePath(
        "holland", "holland1981.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));

    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_boundary_mode(json);
    EXPECT_TRUE(info.usePml);
    EXPECT_FALSE(info.useMur);
    EXPECT_TRUE(info.pmlBack);
    EXPECT_TRUE(info.pmlFront);
    EXPECT_TRUE(info.pmlLeft);
    EXPECT_TRUE(info.pmlRight);
    EXPECT_TRUE(info.pmlDown);
    EXPECT_TRUE(info.pmlUp);
    EXPECT_FALSE(info.murBack);
    EXPECT_FALSE(info.murFront);
    EXPECT_FALSE(info.murLeft);
    EXPECT_FALSE(info.murRight);
    EXPECT_FALSE(info.murDown);
    EXPECT_FALSE(info.murUp);
}

TEST(PmlBoundary, MixedPmlFacesDoNotEnableMurFaces) {
    const std::string json = pml_boundary_test::casePath(
        "sgbcShieldingEffectiveness", "shieldingEffectiveness.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));

    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_boundary_mode(json);
    EXPECT_TRUE(info.usePml);
    EXPECT_FALSE(info.useMur);
    EXPECT_FALSE(info.pmlBack);
    EXPECT_FALSE(info.pmlFront);
    EXPECT_FALSE(info.pmlLeft);
    EXPECT_FALSE(info.pmlRight);
    EXPECT_TRUE(info.pmlDown);
    EXPECT_TRUE(info.pmlUp);
    EXPECT_FALSE(info.murDown);
    EXPECT_FALSE(info.murUp);
}

TEST(PmlBoundary, TimestepCallsPmlPhasesForPmlBoundaries) {
    const std::string json = pml_boundary_test::casePath(
        "sgbcShieldingEffectiveness", "shieldingEffectiveness.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));

    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_boundary_mode(
        json, true);
    EXPECT_EQ(info.pmlElectricCalls, 1);
    EXPECT_EQ(info.pmlBodyHCalls, 1);
    EXPECT_EQ(info.pmlMagneticCpmlCalls, 1);
}

TEST(PmlBoundary, TimestepDoesNotCallPmlPhasesForMurBoundaries) {
    const std::string json = pml_boundary_test::casePath(
        "planewave", "pw-in-box.fdtd.json");
    ASSERT_TRUE(std::filesystem::exists(json));

    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_boundary_mode(
        json, true);
    EXPECT_TRUE(info.useMur);
    EXPECT_FALSE(info.usePml);
    EXPECT_EQ(info.pmlElectricCalls, 0);
    EXPECT_EQ(info.pmlBodyHCalls, 0);
    EXPECT_EQ(info.pmlMagneticCpmlCalls, 0);
}

#endif
