#ifndef TEST_PLANEWAVE_INIT_H
#define TEST_PLANEWAVE_INIT_H

#include <gtest/gtest.h>

#include "semba_fdtd.h"

#include <cmath>
#include <filesystem>
#include <string>

TEST(PlanewaveInit, PwInBoxDirectionPolarizationAndBox) {
    const std::string json =
        (std::filesystem::path("testData") / "cases" / "planewave" / "pw-in-box.fdtd.json").string();
    ASSERT_TRUE(std::filesystem::exists(json));
    const auto info = SEMBA_FDTD_m::SEMBA_FDTD_test::test_plane_wave_init(json, 0);

    EXPECT_NEAR(info.px, 0.0, 1e-12);
    EXPECT_NEAR(info.py, 0.0, 1e-12);
    EXPECT_NEAR(info.pz, 1.0, 1e-12);
    EXPECT_NEAR(info.ex, 1.0, 1e-5);
    EXPECT_NEAR(info.ey, 0.0, 1e-5);
    EXPECT_NEAR(info.ez, 0.0, 1e-5);

    const double zvac = std::sqrt(1.2566370614e-6 / 8.854187817e-12);
    EXPECT_NEAR(info.hy, -info.ez * info.px / zvac + info.ex * info.pz / zvac, 1e-6);
    EXPECT_NEAR(info.hz, info.ey * info.px / zvac - info.ex * info.py / zvac, 1e-6);

    EXPECT_EQ(info.esqx1, 2);
    EXPECT_EQ(info.esqx2, 4);
    EXPECT_EQ(info.esqy1, 2);
    EXPECT_EQ(info.esqy2, 4);
    EXPECT_EQ(info.esqz1, 2);
    EXPECT_EQ(info.esqz2, 4);
    EXPECT_TRUE(info.iluminaAb);
    EXPECT_TRUE(info.iluminaAr);

    EXPECT_NEAR(info.distanciaInicial, 0.01, 1e-8);
    EXPECT_NEAR(info.dt, 1.5406665526684904e-11, 1e-17);
    EXPECT_GE(info.numSteps, 1290);
    EXPECT_LE(info.numSteps, 1310);
    EXPECT_GT(info.numSamples, 100);
    EXPECT_NEAR(info.deltaevol, 1.805468626449816074e-13, 1e-20);
}

#endif
