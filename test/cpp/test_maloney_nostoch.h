#ifndef TEST_MALONEY_NOSTOCH_H
#define TEST_MALONEY_NOSTOCH_H

#include "maloney_nostoch.h"

TEST(MaloneyNostoch, G1G2MatchesFortranFormula) {
    double g1 = 0.0;
    double g2 = 0.0;

    SGBC_nostoch_m::g1g2(3.0e-12, 8.854187817e-12, 0.0, g1, g2);
    EXPECT_DOUBLE_EQ(g1, 1.0);
    EXPECT_DOUBLE_EQ(g2, 0.3388227200511846);

    SGBC_nostoch_m::g1g2(1.0e-12, 2.0 * 8.854187817e-12, 0.01, g1, g2);
    EXPECT_DOUBLE_EQ(g1, 0.9994354548671793);
    EXPECT_DOUBLE_EQ(g2, 0.056454513282072925);
}

TEST(MaloneyNostoch, G1G2UsesFortranExponentialBranch) {
    double g1 = 0.0;
    double g2 = 0.0;

    SGBC_nostoch_m::g1g2(1.0e-12, 8.854187817e-12, 100.0, g1, g2);
    EXPECT_DOUBLE_EQ(g1, 1.2446256433111056e-05);
    EXPECT_DOUBLE_EQ(g2, 0.009999875537435669);
}

TEST(MaloneyNostoch, GM1GM2MatchesFortranFormula) {
    double gm1 = 0.0;
    double gm2 = 0.0;

    SGBC_nostoch_m::gm1gm2(3.0e-12, 1.25663706212e-6, 0.0, gm1, gm2);
    EXPECT_DOUBLE_EQ(gm1, 1.0);
    EXPECT_DOUBLE_EQ(gm2, 2.387324145078829e-06);

    SGBC_nostoch_m::gm1gm2(1.0e-12, 2.0 * 1.25663706212e-6, 1.0e6, gm1, gm2);
    EXPECT_DOUBLE_EQ(gm1, 0.6681350720946382);
    EXPECT_DOUBLE_EQ(gm2, 3.318649279053619e-07);
}

TEST(MaloneyNostoch, DispersiveG1G2MatchesFortranFormula) {
    double g1 = 0.0;
    double g2 = 0.0;
    std::vector<std::complex<double>> beta;
    std::vector<std::complex<double>> kappa;
    std::vector<std::complex<double>> g3;
    const std::vector<std::complex<double>> a11 = {
        {-1.0e9, 2.0e8},
        {-2.0e9, 0.0},
    };
    const std::vector<std::complex<double>> c11 = {
        {3.0e-3, 1.0e-3},
        {4.0e-3, -2.0e-3},
    };

    SGBC_nostoch_m::g1g2_Dispersive(
        1.0e-12, 2.0 * 8.854187817e-12, 0.01,
        g1, g2, beta, kappa, g3, 2, a11, c11);

    EXPECT_DOUBLE_EQ(g1, 0.9994355663049378);
    EXPECT_DOUBLE_EQ(g2, 0.056443369506215625);
    EXPECT_DOUBLE_EQ(kappa[0].real(), 0.9990004797800952);
    EXPECT_DOUBLE_EQ(kappa[0].imag(), 0.00019980014790405747);
    EXPECT_DOUBLE_EQ(kappa[1].real(), 0.9980019980019981);
    EXPECT_DOUBLE_EQ(kappa[1].imag(), 0.0);
    EXPECT_DOUBLE_EQ(beta[0].real(), 2.9984008195961906e-15);
    EXPECT_DOUBLE_EQ(beta[0].imag(), 9.997999401119037e-16);
    EXPECT_DOUBLE_EQ(beta[1].real(), 3.996003996003997e-15);
    EXPECT_DOUBLE_EQ(beta[1].imag(), -1.9980019980019985e-15);
    EXPECT_DOUBLE_EQ(g3[0].real(), 0.05641516136166512);
    EXPECT_DOUBLE_EQ(g3[0].imag(), 5.638696787772625e-06);
    EXPECT_DOUBLE_EQ(g3[1].real(), 0.05638698252369193);
    EXPECT_DOUBLE_EQ(g3[1].imag(), 0.0);
}

TEST(MaloneyNostoch, TridiagonalSolversMatchFortranThomasAlgorithm) {
    std::vector<double> d = {1.0, 0.0, 1.0};
    std::vector<double> x;

    SGBC_nostoch_m::solve_tridiag_iguales(
        -1.0, 2.0, -1.0,
        0.0, 2.0, -1.0,
        -1.0, 2.0, 0.0,
        d, x, 3);
    ASSERT_EQ(x.size(), 3u);
    EXPECT_DOUBLE_EQ(x[0], 1.0);
    EXPECT_DOUBLE_EQ(x[1], 1.0);
    EXPECT_DOUBLE_EQ(x[2], 1.0);

    const std::vector<double> aa = {0.0, -1.0, 0.0};
    const std::vector<double> bb = {0.0, 2.0, 0.0};
    const std::vector<double> cc = {0.0, -1.0, 0.0};
    SGBC_nostoch_m::solve_tridiag_distintos(
        aa, bb, cc,
        0.0, 2.0, -1.0,
        -1.0, 2.0, 0.0,
        d, x, 3);
    ASSERT_EQ(x.size(), 3u);
    EXPECT_DOUBLE_EQ(x[0], 1.0);
    EXPECT_DOUBLE_EQ(x[1], 1.0);
    EXPECT_DOUBLE_EQ(x[2], 1.0);
}

TEST(MaloneyNostoch, LayerDepthMatchesFortranDepthZeroSingleLayer) {
    const auto result = SGBC_nostoch_m::calculate_sgbc_layer_depth(
        {0.010}, {100.0}, {1.0}, 1.0e9, 1.0, 0);

    EXPECT_EQ(result.depth, 0);
    ASSERT_EQ(result.capa.size(), 1u);
    ASSERT_EQ(result.delta_entreEinterno.size(), 1u);
    EXPECT_EQ(result.capa[0], 1);
    EXPECT_DOUBLE_EQ(result.delta_entreEinterno[0], 0.010);
}

TEST(MaloneyNostoch, LayerDepthMatchesFortranPositiveDepthRounding) {
    auto result = SGBC_nostoch_m::calculate_sgbc_layer_depth(
        {0.010}, {100.0}, {1.0}, 1.0e9, 1.0, 3);

    EXPECT_EQ(result.depth, 2);
    ASSERT_EQ(result.capa.size(), 4u);
    ASSERT_EQ(result.delta_entreEinterno.size(), 4u);
    EXPECT_EQ(result.capa, (std::vector<int32_t>{1, 1, 1, 1}));
    for (double delta : result.delta_entreEinterno) {
        EXPECT_DOUBLE_EQ(delta, 0.0025);
    }

    result = SGBC_nostoch_m::calculate_sgbc_layer_depth(
        {0.010, 0.020}, {100.0, 10.0}, {1.0, 2.0}, 1.0e9, 1.0, 3);

    EXPECT_EQ(result.depth, 3);
    ASSERT_EQ(result.capa.size(), 6u);
    ASSERT_EQ(result.delta_entreEinterno.size(), 6u);
    EXPECT_EQ(result.capa, (std::vector<int32_t>{1, 1, 1, 2, 2, 2}));
    EXPECT_DOUBLE_EQ(result.delta_entreEinterno[0], 0.010 / 3.0);
    EXPECT_DOUBLE_EQ(result.delta_entreEinterno[1], 0.010 / 3.0);
    EXPECT_DOUBLE_EQ(result.delta_entreEinterno[2], 0.010 / 3.0);
    EXPECT_DOUBLE_EQ(result.delta_entreEinterno[3], 0.020 / 3.0);
    EXPECT_DOUBLE_EQ(result.delta_entreEinterno[4], 0.020 / 3.0);
    EXPECT_DOUBLE_EQ(result.delta_entreEinterno[5], 0.020 / 3.0);
}

#endif // TEST_MALONEY_NOSTOCH_H
