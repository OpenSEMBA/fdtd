#ifndef TEST_OBSERVATION_H
#define TEST_OBSERVATION_H

#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <sstream>
#include <vector>

#include "nfde_types.h"
#include "observation_types.h"
#include "observation_preprocess.h"
#include "observation_movie_test.h"

using namespace Observa_m;
using NFDETypes_m::mapvtk;

namespace observation_test {

inline std::vector<RKIND_tiempo> create_time_array(int n, RKIND_tiempo dt) {
    std::vector<RKIND_tiempo> arr(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        arr[static_cast<size_t>(i)] = static_cast<RKIND_tiempo>(i) * dt;
    }
    return arr;
}

constexpr RKIND kObsTol = 1e-12;

void expect_vector_size(const std::vector<double>& v, size_t n) {
    EXPECT_EQ(v.size(), n);
}

void expect_vector_size(const std::vector<std::complex<double>>& v, size_t n) {
    EXPECT_EQ(v.size(), n);
}

void expect_vector_size(const std::vector<int>& v, size_t n) {
    EXPECT_EQ(v.size(), n);
}

void check_time_domain_serialized(Serialized_t& serialize, int numberOfSerialized) {
    const size_t n = static_cast<size_t>(numberOfSerialized);
    expect_vector_size(serialize.Valor, n);
    expect_vector_size(serialize.Valor_x, n);
    expect_vector_size(serialize.Valor_y, n);
    expect_vector_size(serialize.Valor_z, n);
    expect_vector_size(serialize.ValorE, n);
    expect_vector_size(serialize.Valor_Ex, n);
    expect_vector_size(serialize.Valor_Ey, n);
    expect_vector_size(serialize.Valor_Ez, n);
    expect_vector_size(serialize.ValorH, n);
    expect_vector_size(serialize.Valor_Hx, n);
    expect_vector_size(serialize.Valor_Hy, n);
    expect_vector_size(serialize.Valor_Hz, n);
}

void check_frequency_domain_serialized(Serialized_t& serialize, int numberOfSerialized) {
    check_time_domain_serialized(serialize, numberOfSerialized);
    const size_t n = static_cast<size_t>(numberOfSerialized);
    expect_vector_size(serialize.ValorComplex_x, n);
    expect_vector_size(serialize.ValorComplex_y, n);
    expect_vector_size(serialize.ValorComplex_z, n);
    expect_vector_size(serialize.ValorComplex_Ex, n);
    expect_vector_size(serialize.ValorComplex_Ey, n);
    expect_vector_size(serialize.ValorComplex_Ez, n);
    expect_vector_size(serialize.ValorComplex_Hx, n);
    expect_vector_size(serialize.ValorComplex_Hy, n);
    expect_vector_size(serialize.ValorComplex_Hz, n);
}

} // namespace observation_test

using observation_test::create_time_array;
using observation_test::check_time_domain_serialized;
using observation_test::check_frequency_domain_serialized;
using observation_test::expect_vector_size;
using observation_test::kObsTol;

TEST(observation, test_allocate_time) {
    Serialized_t serialize;
    const int numberOfSerialized = 4;

    serialize.allocate_for_time_domain(numberOfSerialized);
    check_time_domain_serialized(serialize, numberOfSerialized);
    serialize.deallocate_for_time_domain();

    EXPECT_TRUE(serialize.Valor.empty());
}

TEST(observation, test_allocate_frequency) {
    Serialized_t serialize;
    const int numberOfSerialized = 4;

    serialize.allocate_for_frequency_domain(numberOfSerialized);
    check_frequency_domain_serialized(serialize, numberOfSerialized);
    serialize.deallocate_for_frequency_domain();

    EXPECT_TRUE(serialize.Valor.empty());
    EXPECT_TRUE(serialize.ValorComplex_x.empty());
}

TEST(observation, test_allocate_serialize_current) {
    Serialized_t serialize;
    const int numberOfSerialized = 4;
    const size_t n = static_cast<size_t>(numberOfSerialized);

    serialize.allocate_current_value(numberOfSerialized);
    expect_vector_size(serialize.eI, n);
    expect_vector_size(serialize.eJ, n);
    expect_vector_size(serialize.eK, n);
    expect_vector_size(serialize.currentType, n);
    expect_vector_size(serialize.sggmtag, n);
    serialize.deallocate_current_value();

    EXPECT_TRUE(serialize.eI.empty());
}

TEST(observation, test_preproces_initial_time_less_than_timestep) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 20;
    const RKIND_tiempo dt = 0.1;
    const auto tiempo = create_time_array(100, dt);
    const bool saveall = true;

    obs.TimeStep = 0.5;
    obs.InitialTime = 0.2;
    obs.FinalTime = 2.0;
    obs.nP = 0;
    obs.Volumic = false;
    obs.InitialFreq = 0.0;
    obs.FinalFreq = 1.0;
    obs.FreqStep = 0.1;
    out.SaveAll = false;

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);

    EXPECT_TRUE(approx_equal(obs.InitialTime, 0.0, kObsTol));
}

TEST(observation, test_preproces_timestep_greater_and_mapvtk) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 90;
    const RKIND_tiempo dt = 0.1;
    const auto tiempo = create_time_array(100, dt);
    const bool saveall = false;

    obs.TimeStep = 5.0;
    obs.InitialTime = 0.0;
    obs.FinalTime = 1.0;
    obs.nP = 1;
    obs.P = {{mapvtk}};
    obs.Volumic = false;
    obs.InitialFreq = 0.0;
    obs.FinalFreq = 1.0;
    obs.FreqStep = 0.1;

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);

    EXPECT_TRUE(approx_equal(obs.InitialTime, 0.0, kObsTol));
    EXPECT_TRUE(approx_equal(obs.FinalTime, 0.0, kObsTol));
}

TEST(observation, test_preproces_timestep_greater_not_mapvtk) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 90;
    const RKIND_tiempo dt = 0.1;
    const auto tiempo = create_time_array(100, dt);
    const bool saveall = false;

    obs.TimeStep = 2.0;
    obs.InitialTime = 1.0;
    obs.FinalTime = 1.5;
    obs.nP = 1;
    obs.P = {{999}};
    obs.Volumic = false;
    obs.InitialFreq = 0.0;
    obs.FinalFreq = 1.0;
    obs.FreqStep = 0.1;
    obs.Saveall = false;

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);

    EXPECT_TRUE(approx_equal(obs.FinalTime, obs.InitialTime + obs.TimeStep, kObsTol));
}

TEST(observation, test_preproces_freqstep_zero_or_large) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 90;
    const RKIND_tiempo dt = 0.1;
    const auto tiempo = create_time_array(100, dt);
    const bool saveall = false;

    obs.TimeStep = 0.2;
    obs.InitialTime = 0.0;
    obs.FinalTime = 1.0;
    obs.InitialFreq = 1.0;
    obs.FinalFreq = 3.5;
    obs.FreqStep = 0.0;
    obs.nP = 0;
    obs.Volumic = false;

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);
    EXPECT_TRUE(approx_equal(obs.FreqStep, obs.FinalFreq - obs.InitialFreq, kObsTol));

    obs.FreqStep = 10.0;
    obs.InitialFreq = 0.0;
    obs.FinalFreq = 2.0;
    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);
    EXPECT_TRUE(approx_equal(obs.FreqStep, obs.FinalFreq - obs.InitialFreq, kObsTol));
}

TEST(observation, test_preproces_volumic_false_true_and_saveall) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 90;
    const RKIND_tiempo dt = 0.1;
    auto tiempo = create_time_array(100, dt);
    bool saveall = false;

    obs.Volumic = false;
    obs.Saveall = false;
    obs.TimeStep = 0.2;
    obs.InitialTime = 0.0;
    obs.FinalTime = 1.0;
    obs.InitialFreq = 0.0;
    obs.FinalFreq = 1.0;
    out.SaveAll = false;
    obs.nP = 1;
    obs.P = {{999}};

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);
    const bool ok1 = !out.SaveAll;

    saveall = true;
    obs.Saveall = false;
    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);
    const bool ok2 = obs.Saveall || out.SaveAll;

    EXPECT_TRUE(ok1);
    EXPECT_TRUE(ok2);
}

TEST(observation, test_preproces_saveall_branch) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 90;
    const RKIND_tiempo dt = 0.1;
    const auto tiempo = create_time_array(100, dt);
    const bool saveall = false;

    obs.Volumic = false;
    obs.Saveall = true;
    obs.TimeStep = 0.5;
    obs.InitialTime = 10.0;
    obs.FinalTime = 20.0;
    obs.nP = 1;
    obs.P = {{999}};

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);

    const RKIND expectedFinal =
        tiempo[static_cast<size_t>(finalTimeIndex + 2)];
    EXPECT_EQ(out.Trancos, 1);
    EXPECT_TRUE(approx_equal(obs.InitialTime, 0.0, kObsTol));
    EXPECT_TRUE(approx_equal(obs.FinalTime, expectedFinal, kObsTol));
}

TEST(observation, test_preproces_final_less_than_initial) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 90;
    const RKIND_tiempo dt = 0.1;
    const auto tiempo = create_time_array(100, dt);
    const bool saveall = false;

    obs.Volumic = false;
    obs.Saveall = false;
    obs.TimeStep = 0.2;
    obs.InitialTime = 5.0;
    obs.FinalTime = 1.0;
    obs.nP = 1;
    obs.P = {{999}};

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);

    EXPECT_TRUE(approx_equal(obs.FinalTime, obs.InitialTime + obs.TimeStep, kObsTol));
}

TEST(observation, test_preproces_huge_cap) {
    Obses_t_full obs;
    output_t out;
    const int finalTimeIndex = 90;
    const RKIND_tiempo dt = 0.1;
    const auto tiempo = create_time_array(100, dt);
    const float huge4 = std::numeric_limits<float>::max();
    const bool saveall = false;

    obs.TimeStep = 0.1;
    obs.InitialTime = 0.0;
    obs.FinalTime = obs.InitialTime +
        static_cast<RKIND>(huge4 / 5.0f) * std::min(0.1, obs.TimeStep);
    obs.Volumic = false;
    obs.Saveall = false;
    obs.InitialFreq = 0.0;
    obs.FinalFreq = 1.0;
    obs.FreqStep = 0.1;
    obs.nP = 1;
    obs.P = {{999}};

    preprocess_observation_full(obs, out, tiempo, finalTimeIndex, dt, saveall);

    const RKIND clip = tiempo[static_cast<size_t>(finalTimeIndex + 2)] + dt;
    EXPECT_LE(obs.FinalTime, clip);
}

TEST(observation, test_init_movie_observation) {
    SGGFDTDINFO_t sgg = create_base_sgg();
    sgg.Observation.resize(1);
    sgg.Observation[0] = define_time_movie_observation();

    media_matrices_t media = create_media(sgg.alloc);
    taglist_t tag_numbers = create_tag_list(sgg.alloc);

    bool thereAreObservation = false;
    bool thereAreWires = false;
    bool thereAreFarFields = false;
    int initialtimestep = 1;
    RKIND_tiempo lastexecutedtime = 0.0;

    const std::vector<limit_t> sinpml_fullsize = {
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
    };
    bounds_t bounds;
    const nf2ff_t faces = create_faces_nf2ff(false, false, false, false, false, false);
    sim_control_t control = create_control_flags(
        0, 0, 3, 10, get_temp_dir() + "/entryRoot", "wireflavour", false, false, false, false,
        false, faces);

    InitObservation(sgg, media, tag_numbers, thereAreObservation, thereAreWires, thereAreFarFields,
                    initialtimestep, lastexecutedtime, sinpml_fullsize,
                    NFDETypes_m::EPSILON_VACUUM, NFDETypes_m::MU_VACUUM, bounds, control);

    std::vector<output_movie_t>& output = GetOutput();
    ASSERT_EQ(output.size(), 1u);
    EXPECT_EQ(output[0].timeswritten, 0);
    ASSERT_GE(output[0].item.size(), 1u);
    EXPECT_EQ(output[0].item[0].unit, 1001);

    const observable_t& probe = sgg.Observation[0].P[0];
    const int sweepIdx = NFDETypes_m::iEz - 1;
    const int chark = std::max(sgg.Sweep[static_cast<size_t>(sweepIdx)].ZI, probe.ZI);
    const int chark2 = std::min(sgg.Sweep[static_cast<size_t>(sweepIdx)].ZE, probe.ZE);

    std::ostringstream expected;
    expected << get_temp_dir() << "/entryRoot_timeMovie_ME_" << probe.XI << '_' << probe.YI << '_'
             << chark << "__" << probe.XE << '_' << probe.YE << '_' << chark2 << ".bin";
    EXPECT_EQ(output[0].item[0].path, expected.str());
}

TEST(observation, test_update_movie_observation) {
    SGGFDTDINFO_t sgg = create_base_sgg();
    sgg.Observation.resize(1);
    sgg.Observation[0] = define_time_movie_observation();

    media_matrices_t media = create_media(sgg.alloc);
    taglist_t tag_numbers = create_tag_list(sgg.alloc);

    bool thereAreObservation = false;
    bool thereAreWires = false;
    bool thereAreFarFields = false;
    int initialtimestep = 1;
    RKIND_tiempo lastexecutedtime = 0.0;

    const std::vector<limit_t> sinpml_fullsize = {
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
        create_limit(0, 4, 0, 4, 0, 4, 3, 3, 3),
    };
    bounds_t bounds;
    const nf2ff_t faces = create_faces_nf2ff(false, false, false, false, false, false);
    sim_control_t control = create_control_flags(
        0, 0, 3, 10, get_temp_dir() + "/entryRoot", "wireflavour", false, false, false, false,
        false, faces);

    InitObservation(sgg, media, tag_numbers, thereAreObservation, thereAreWires, thereAreFarFields,
                    initialtimestep, lastexecutedtime, sinpml_fullsize,
                    NFDETypes_m::EPSILON_VACUUM, NFDETypes_m::MU_VACUUM, bounds, control);

    dummyFields_t fields;
    fields.createDummyFields(0, 10, 0.1);

    UpdateObservation(sgg, media, tag_numbers, 5, 0, fields.Ex, fields.Ey, fields.Ez, fields.Hx,
                      fields.Hy, fields.Hz, fields.dxe, fields.dye, fields.dze, fields.dxh,
                      fields.dyh, fields.dzh, control.wiresflavor, sinpml_fullsize, false, false,
                      bounds);

    std::vector<output_movie_t>& output = GetOutput();
    ASSERT_EQ(output.size(), 1u);
    EXPECT_GE(output[0].timeswritten, 1);
}

#endif
