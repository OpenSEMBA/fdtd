#include <gtest/gtest.h>

extern "C" int test_rotate_generate_space_steps();
extern "C" int test_rotate_generate_current_field_sources();
extern "C" int test_rotate_generate_plane_waves();
extern "C" int test_rotate_generate_box_sources();
extern "C" int test_rotate_generate_fronteras();
extern "C" int test_rotate_generate_pecs();
extern "C" int test_rotate_generate_non_metals();
extern "C" int test_rotate_generate_thin_wires();
extern "C" int test_rotate_generate_slanted_wires();
extern "C" int test_rotate_generate_thin_slots();
extern "C" int test_rotate_generate_lossy_thin_surface();
extern "C" int test_rotate_generate_fdms();
extern "C" int test_rotate_generate_sondas();
extern "C" int test_rotate_generate_mas_sondas();
extern "C" int test_rotate_generate_bloque_probes();
extern "C" int test_rotate_generate_volumic_probes();


TEST(rotate, rotate_spaceSteps_test)    { EXPECT_EQ(0, test_rotate_generate_space_steps()); }
TEST(rotate, rotate_currentFieldSources_test) { EXPECT_EQ(0, test_rotate_generate_current_field_sources()); }
TEST(rotate, rotate_planeWaves_test)    { EXPECT_EQ(0, test_rotate_generate_plane_waves()); }
TEST(rotate, rotate_boxSources_test)    { EXPECT_EQ(0, test_rotate_generate_box_sources()); }
TEST(rotate, rotate_fronteras_test)     { EXPECT_EQ(0, test_rotate_generate_fronteras()); }
TEST(rotate, rotate_pecs_test)          { EXPECT_EQ(0, test_rotate_generate_pecs()); }
TEST(rotate, rotate_nonMetals_test)     { EXPECT_EQ(0, test_rotate_generate_non_metals()); }
TEST(rotate, rotate_thinWires_test)     { EXPECT_EQ(0, test_rotate_generate_thin_wires()); }
TEST(rotate, rotate_slantedWires_test)  { EXPECT_EQ(0, test_rotate_generate_slanted_wires()); }
TEST(rotate, rotate_thinSlots_test)     { EXPECT_EQ(0, test_rotate_generate_thin_slots()); }
TEST(rotate, rotate_lossyThinSurface_test) { EXPECT_EQ(0, test_rotate_generate_lossy_thin_surface()); }
TEST(rotate, rotate_freqDependMaterials_test) { EXPECT_EQ(0, test_rotate_generate_fdms()); }
TEST(rotate, rotate_sondas_test)        { EXPECT_EQ(0, test_rotate_generate_sondas()); }
TEST(rotate, rotate_masSondas_test)     { EXPECT_EQ(0, test_rotate_generate_mas_sondas()); }
TEST(rotate, rotate_bloqueProbes_test)  { EXPECT_EQ(0, test_rotate_generate_bloque_probes()); }
TEST(rotate, rotate_volumicProbes_test) { EXPECT_EQ(0, test_rotate_generate_volumic_probes()); }


