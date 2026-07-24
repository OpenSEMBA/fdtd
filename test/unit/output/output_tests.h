#ifdef CompileWithNewOutputModule

#include <gtest/gtest.h>

// Shared utilities and decomposition.
extern "C" int test_count_required_coords();
extern "C" int test_store_required_coords();
extern "C" int test_is_valid_point_current();
extern "C" int test_is_valid_point_field();
extern "C" int test_output_partition_shared_interfaces();
extern "C" int test_output_partition_disjoint_intervals();
extern "C" int test_output_partition_clipping_and_large_shape();
extern "C" int test_output_point_partition_eligibility();
extern "C" int test_output_partition_all_components_cover_volume();

// Distributed publication coordination.
extern "C" int test_output_collective_contract();
extern "C" int test_output_publication_contract();
extern "C" int test_point_publication_plan();
extern "C" int test_planned_metadata_publication();
extern "C" int test_output_transport_serial();

// Artifact metadata and binary layout.
extern "C" int test_output_metadata_publication();
extern "C" int test_output_metadata_fragment_descriptors();
extern "C" int test_portable_binary_output();
extern "C" int test_output_binary_fragment_layout();

// Probe lifecycle behaviour.
extern "C" int test_init_point_probe();
extern "C" int test_update_point_probe();
extern "C" int test_flush_point_probe();
extern "C" int test_multiple_flush_point_probe();
extern "C" int test_init_movie_probe();
extern "C" int test_update_movie_probe();
extern "C" int test_flush_movie_probe();
extern "C" int test_init_frequency_slice_probe();
extern "C" int test_update_frequency_slice_probe();

// Generic output lifecycle and adapters.
extern "C" int test_root_output_manifest();
extern "C" int test_nested_output_path();
extern "C" int test_output_artifact_contract();
extern "C" int test_volumetric_visualisation_output();
extern "C" int test_declared_output_artifacts();
extern "C" int test_output_lifecycle_contract();
extern "C" int test_output_lifecycle_coordination();
extern "C" int test_output_probe_ownership();
extern "C" int test_output_serial_distributed_equivalence();
extern "C" int test_volumetric_output_partition_attachment();

TEST(output, test_volumic_utils_count) { EXPECT_EQ(0, test_count_required_coords()); }
TEST(output, test_volumic_utils_store) { EXPECT_EQ(0, test_store_required_coords()); }
TEST(output, test_volumic_utils_valid_current) { EXPECT_EQ(0, test_is_valid_point_current()); }
TEST(output, test_volumic_utils_valid_field) { EXPECT_EQ(0, test_is_valid_point_field()); }
TEST(output, test_partition_shared_interfaces) { EXPECT_EQ(0, test_output_partition_shared_interfaces()); }
TEST(output, test_partition_disjoint_intervals) { EXPECT_EQ(0, test_output_partition_disjoint_intervals()); }
TEST(output, test_partition_clipping_and_large_shape) { EXPECT_EQ(0, test_output_partition_clipping_and_large_shape()); }
TEST(output, test_point_partition_eligibility) { EXPECT_EQ(0, test_output_point_partition_eligibility()); }
TEST(output, test_partition_all_components_cover_volume) { EXPECT_EQ(0, test_output_partition_all_components_cover_volume()); }
TEST(output, test_collective_contract) { EXPECT_EQ(0, test_output_collective_contract()); }
TEST(output, test_publication_contract) { EXPECT_EQ(0, test_output_publication_contract()); }
TEST(output, test_point_publication_plan) { EXPECT_EQ(0, test_point_publication_plan()); }
TEST(output, test_planned_metadata_publication) { EXPECT_EQ(0, test_planned_metadata_publication()); }
TEST(output, test_transport_serial) { EXPECT_EQ(0, test_output_transport_serial()); }
TEST(output, test_metadata_publication) { EXPECT_EQ(0, test_output_metadata_publication()); }
TEST(output, test_metadata_fragment_descriptors) { EXPECT_EQ(0, test_output_metadata_fragment_descriptors()); }
TEST(output, test_portable_binary_output) { EXPECT_EQ(0, test_portable_binary_output()); }
TEST(output, test_binary_fragment_layout) { EXPECT_EQ(0, test_output_binary_fragment_layout()); }
TEST(output, test_initialize_point_probe) { EXPECT_EQ(0, test_init_point_probe()); }
TEST(output, test_update_point_probe_info) { EXPECT_EQ(0, test_update_point_probe()); }
TEST(output, test_flush_point_probe_info) { EXPECT_EQ(0, test_flush_point_probe()); }
TEST(output, test_flush_multiple_point_probe_info) { EXPECT_EQ(0, test_multiple_flush_point_probe()); }
TEST(output, test_init_movie_probe_for_pec_surface) { EXPECT_EQ(0, test_init_movie_probe()); }
TEST(output, test_update_movie_probe_for_pec_surface) { EXPECT_EQ(0, test_update_movie_probe()); }
TEST(output, test_flush_movie_probe_data) { EXPECT_EQ(0, test_flush_movie_probe()); }
TEST(output, test_init_frequency_slice) { EXPECT_EQ(0, test_init_frequency_slice_probe()); }
TEST(output, test_update_frequency_slice) { EXPECT_EQ(0, test_update_frequency_slice_probe()); }
TEST(output, test_root_output_manifest) { EXPECT_EQ(0, test_root_output_manifest()); }
TEST(output, test_nested_output_path) { EXPECT_EQ(0, test_nested_output_path()); }
TEST(output, test_artifact_contract) { EXPECT_EQ(0, test_output_artifact_contract()); }
TEST(output, test_volumetric_visualisation_output) { EXPECT_EQ(0, test_volumetric_visualisation_output()); }
TEST(output, test_declared_output_artifacts) { EXPECT_EQ(0, test_declared_output_artifacts()); }
TEST(output, test_lifecycle_contract) { EXPECT_EQ(0, test_output_lifecycle_contract()); }
TEST(output, test_lifecycle_coordination) { EXPECT_EQ(0, test_output_lifecycle_coordination()); }
TEST(output, test_probe_ownership) { EXPECT_EQ(0, test_output_probe_ownership()); }
TEST(output, test_serial_distributed_equivalence) { EXPECT_EQ(0, test_output_serial_distributed_equivalence()); }
TEST(output, test_volumetric_partition_attachment) { EXPECT_EQ(0, test_volumetric_output_partition_attachment()); }

#endif
