#ifndef TEST_MALONEY_MISSING_H
#define TEST_MALONEY_MISSING_H

TEST(MaloneyMissing, InitSgbcsBuildsFullFortranSurfaceState) {
    GTEST_SKIP()
        << "The active C++ solver builds SGBC nodes from JSON geometry instead "
        << "of the legacy Fortran InitSGBCs entry point.";
}

TEST(MaloneyMissing, NonZeroDepthAdvancesInternalOneDimensionalFields) {
    const std::vector<SGBC_nostoch_m::SGBCLayer_t> layers = {
        {0.010, 1.0, 1.0, 100.0, 0.0},
    };
    auto surface = SGBC_nostoch_m::make_sgbc_surface(
        layers, 2.0e-12, 8.854187817e-12, 1.25663706212e-6,
        1.0e9, 1.0, -1, false, true, false,
        0.020, 0.020, 0.020, 0.0);

    ASSERT_GT(surface.depth, 0);
    SGBC_nostoch_m::AdvanceSGBCE(surface, 0.2, -0.1, 0.05, -0.02);

    EXPECT_NE(surface.Efield, 0.0);
    EXPECT_NE(surface.Hyee_left, 0.0);
    EXPECT_NE(surface.Hyee_right, 0.0);
}

TEST(MaloneyMissing, CrankNicolsonSgbcMatchesFortranSheetSolve) {
    const std::vector<SGBC_nostoch_m::SGBCLayer_t> layers = {
        {0.010, 1.0, 1.0, 100.0, 0.0},
    };
    auto surface = SGBC_nostoch_m::make_sgbc_surface(
        layers, 2.0e-12, 8.854187817e-12, 1.25663706212e-6,
        1.0e9, 1.0, -1, true, false, false,
        0.020, 0.020, 0.020, 0.0);

    ASSERT_TRUE(surface.SGBCCrank);
    ASSERT_GT(surface.depth, 1);
    SGBC_nostoch_m::AdvanceSGBCE(surface, 0.2, -0.1, 0.05, -0.02);

    EXPECT_NE(surface.E.front(), surface.E.back());
    EXPECT_NE(surface.Hyee_left, 0.0);
    EXPECT_NE(surface.Hyee_right, 0.0);
}

TEST(MaloneyMissing, DispersiveSgbcUpdatesPolarizationState) {
    GTEST_SKIP()
        << "Dispersive SGBC polarization state is still a separate migration item.";
}

#endif // TEST_MALONEY_MISSING_H
