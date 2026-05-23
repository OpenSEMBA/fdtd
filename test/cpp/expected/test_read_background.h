#ifndef TEST_READ_BACKGROUND_H
#define TEST_READ_BACKGROUND_H

#include "test_smbjson_helpers.h"

// Scalar checks only (no full Parseador expected builder).

inline constexpr double EXPECTED_BACKGROUND_EPS = 1.7708e-11;
inline constexpr double EXPECTED_BACKGROUND_MU = 2.5133e-6;

inline void expectBackgroundDefaults(const Parseador_t& pr) {
    EXPECT_EQ(pr.Mats->Mats[0].eps, EPSILON_VACUUM);
    EXPECT_EQ(pr.Mats->Mats[0].mu, MU_VACUUM);
}

inline void expectBackgroundSet(const Parseador_t& pr) {
    EXPECT_NEAR(pr.Mats->Mats[0].eps, EXPECTED_BACKGROUND_EPS,
                EXPECTED_BACKGROUND_EPS * 1e-4);
    EXPECT_NEAR(pr.Mats->Mats[0].mu, EXPECTED_BACKGROUND_MU,
                EXPECTED_BACKGROUND_MU * 1e-4);
}

#endif // TEST_READ_BACKGROUND_H
