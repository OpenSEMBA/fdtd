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
extern "C" int test_output_metadata_contract_edges();
extern "C" int test_scalar_metadata_publication();
extern "C" int test_point_publication_plan();
extern "C" int test_planned_metadata_publication();
extern "C" int test_output_transport_serial();

// Artifact metadata and binary layout.
extern "C" int test_output_metadata_publication();
extern "C" int test_output_metadata_fragment_descriptors();
extern "C" int test_atomic_file_replacement();
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
extern "C" int test_output_failure_coordination();
extern "C" int test_output_probe_ownership();
extern "C" int test_output_serial_distributed_equivalence();
extern "C" int test_volumetric_output_partition_attachment();

// Counts valid field coordinates in a volumetric probe region.
TEST(output, test_volumic_utils_count) { EXPECT_EQ(0, test_count_required_coords()); }
// Stores all coordinates selected for a volumetric probe.
TEST(output, test_volumic_utils_store) { EXPECT_EQ(0, test_store_required_coords()); }
// Rejects current probes in empty space without PEC or wires.
TEST(output, test_volumic_utils_valid_current) { EXPECT_EQ(0, test_is_valid_point_current()); }
// Accepts in-bounds field points and rejects out-of-bounds points.
TEST(output, test_volumic_utils_valid_field) { EXPECT_EQ(0, test_is_valid_point_field()); }
// Assigns shared interfaces to exactly one output partition.
TEST(output, test_partition_shared_interfaces) { EXPECT_EQ(0, test_output_partition_shared_interfaces()); }
// Covers disjoint rank intervals without overlap or gaps.
TEST(output, test_partition_disjoint_intervals) { EXPECT_EQ(0, test_output_partition_disjoint_intervals()); }
// Clips requests, preserves large shapes, and rejects invalid components.
TEST(output, test_partition_clipping_and_large_shape) { EXPECT_EQ(0, test_output_partition_clipping_and_large_shape()); }
// Identifies points owned by a non-empty output partition.
TEST(output, test_point_partition_eligibility) { EXPECT_EQ(0, test_output_point_partition_eligibility()); }
// Covers the requested volume exactly once for every field component.
TEST(output, test_partition_all_components_cover_volume) { EXPECT_EQ(0, test_output_partition_all_components_cover_volume()); }
// Selects participants and publication modes for serial and distributed output.
TEST(output, test_collective_contract) { EXPECT_EQ(0, test_output_collective_contract()); }
// Validates artifact identities and complete probe metadata.
TEST(output, test_publication_contract) { EXPECT_EQ(0, test_output_publication_contract()); }
// Rejects incomplete metadata identity and non-relative artifact paths.
TEST(output, test_metadata_contract_edges) { EXPECT_EQ(0, test_output_metadata_contract_edges()); }
// Publishes metadata for a scalar probe from its declared artifacts.
TEST(output, test_scalar_metadata_publication) { EXPECT_EQ(0, test_scalar_metadata_publication()); }
// Selects a canonical point writer and rejects unowned points.
TEST(output, test_point_publication_plan) { EXPECT_EQ(0, test_point_publication_plan()); }
// Publishes metadata only from the planned canonical writer.
TEST(output, test_planned_metadata_publication) { EXPECT_EQ(0, test_planned_metadata_publication()); }
// Preserves serial eligibility, reductions, and flush transfers.
TEST(output, test_transport_serial) { EXPECT_EQ(0, test_output_transport_serial()); }
// Publishes declared and failed metadata states with their artifacts.
TEST(output, test_metadata_publication) { EXPECT_EQ(0, test_output_metadata_publication()); }
// Serializes canonical and fragment metadata descriptors correctly.
TEST(output, test_metadata_fragment_descriptors) { EXPECT_EQ(0, test_output_metadata_fragment_descriptors()); }
// Replaces a completed file without exposing a partial target.
TEST(output, test_atomic_file_replacement) { EXPECT_EQ(0, test_atomic_file_replacement()); }
// Writes portable little-endian real32 binary data.
TEST(output, test_portable_binary_output) { EXPECT_EQ(0, test_portable_binary_output()); }
// Requires identity metadata for binary fragments.
TEST(output, test_binary_fragment_layout) { EXPECT_EQ(0, test_output_binary_fragment_layout()); }
// Registers a point probe and its declared output paths.
TEST(output, test_initialize_point_probe) { EXPECT_EQ(0, test_init_point_probe()); }
// Records point-probe values and their corresponding timesteps.
TEST(output, test_update_point_probe_info) { EXPECT_EQ(0, test_update_point_probe()); }
// Flushes point-probe samples and resets serialized time data.
TEST(output, test_flush_point_probe_info) { EXPECT_EQ(0, test_flush_point_probe()); }
// Preserves point-probe data across consecutive flushes.
TEST(output, test_flush_multiple_point_probe_info) { EXPECT_EQ(0, test_multiple_flush_point_probe()); }
// Allocates a movie probe and creates its output directory.
TEST(output, test_init_movie_probe_for_pec_surface) { EXPECT_EQ(0, test_init_movie_probe()); }
// Samples movie-probe field components for a timestep.
TEST(output, test_update_movie_probe_for_pec_surface) { EXPECT_EQ(0, test_update_movie_probe()); }
// Flushes movie data and creates binary, HDF5, XDMF, and descriptor artifacts.
TEST(output, test_flush_movie_probe_data) { EXPECT_EQ(0, test_flush_movie_probe()); }
// Allocates a frequency-slice probe and its frequency buffers.
TEST(output, test_init_frequency_slice) { EXPECT_EQ(0, test_init_frequency_slice_probe()); }
// Computes frequency-slice values for a field gradient.
TEST(output, test_update_frequency_slice) { EXPECT_EQ(0, test_update_frequency_slice_probe()); }
// Creates a root manifest containing declared output artifacts.
TEST(output, test_root_output_manifest) { EXPECT_EQ(0, test_root_output_manifest()); }
// Creates files and parent directories for nested output paths.
TEST(output, test_nested_output_path) { EXPECT_EQ(0, test_nested_output_path()); }
// Stores binary artifact encoding and lifecycle metadata.
TEST(output, test_artifact_contract) { EXPECT_EQ(0, test_output_artifact_contract()); }
// Writes and verifies volumetric XDMF and HDF5 artifacts.
TEST(output, test_volumetric_visualisation_output) { EXPECT_EQ(0, test_volumetric_visualisation_output()); }
// Declares output artifact kinds and relative paths.
TEST(output, test_declared_output_artifacts) { EXPECT_EQ(0, test_declared_output_artifacts()); }
// Validates output lifecycle terminal states and completeness.
TEST(output, test_lifecycle_contract) { EXPECT_EQ(0, test_output_lifecycle_contract()); }
// Finalizes probes and restricts manifest publication to the root rank.
TEST(output, test_lifecycle_coordination) { EXPECT_EQ(0, test_output_lifecycle_coordination()); }
// Retains failure diagnostics and rejects incomplete finalisation.
TEST(output, test_failure_coordination) { EXPECT_EQ(0, test_output_failure_coordination()); }
// Retains probe participants and validates the scalar writer.
TEST(output, test_probe_ownership) { EXPECT_EQ(0, test_output_probe_ownership()); }
// Keeps artifact identity and partition coverage equivalent across modes.
TEST(output, test_serial_distributed_equivalence) { EXPECT_EQ(0, test_output_serial_distributed_equivalence()); }
// Attaches a volumetric partition and selects serial publication fallback.
TEST(output, test_volumetric_partition_attachment) { EXPECT_EQ(0, test_volumetric_output_partition_attachment()); }

#endif
