#pragma once

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
extern "C" int test_output_transport_serial();

// Artifact metadata and binary layout.
extern "C" int test_atomic_file_replacement();
extern "C" int test_portable_binary_output();
extern "C" int test_output_binary_fragment_layout();
extern "C" int test_output_binary_native_precision();
extern "C" int test_output_binary_mixed_complex_layout();
extern "C" int test_output_binary_append_real64();
extern "C" int test_output_binary_append_empty_real64();

// Probe lifecycle behaviour.
extern "C" int test_init_point_probe();
extern "C" int test_init_point_probe_with_incident();
extern "C" int test_line_probe_integral();
extern "C" int test_line_probe_empty_path();
extern "C" int test_line_probe_dat_output();
extern "C" int test_line_probe_shared_interface_owner();
extern "C" int test_update_point_probe();
extern "C" int test_update_time_probe_ranges();
extern "C" int test_flush_point_probe();
extern "C" int test_flush_wire_probe_dat();
extern "C" int test_flush_bulk_probe_dat();
extern "C" int test_multiple_flush_point_probe();
extern "C" int test_init_movie_probe();
extern "C" int test_update_movie_probe();
extern "C" int test_flush_movie_probe();
extern "C" int test_close_movie_probe();
extern "C" int test_init_frequency_slice_probe();
extern "C" int test_update_frequency_slice_probe();

// Generic output lifecycle and adapters.
extern "C" int test_scalar_probe_has_no_manifest();
extern "C" int test_nested_output_path();
extern "C" int test_output_artifact_contract();
extern "C" int test_declared_output_artifacts();
extern "C" int test_output_lifecycle_contract();
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
// Preserves serial flush transfers.
TEST(output, test_transport_serial) { EXPECT_EQ(0, test_output_transport_serial()); }
// Replaces a completed file without exposing a partial target.
TEST(output, test_atomic_file_replacement) { EXPECT_EQ(0, test_atomic_file_replacement()); }
// Writes portable little-endian real64 binary data through the production append API.
TEST(output, test_portable_binary_output) { EXPECT_EQ(0, test_portable_binary_output()); }
// Requires identity metadata for binary fragments.
TEST(output, test_binary_fragment_layout) { EXPECT_EQ(0, test_output_binary_fragment_layout()); }
// Selects the highest supported native precision for mixed-value records.
TEST(output, test_binary_native_precision) { EXPECT_EQ(0, test_output_binary_native_precision()); }
// Accepts coordinate-plus-complex records at the selected record precision.
TEST(output, test_binary_mixed_complex_layout) { EXPECT_EQ(0, test_output_binary_mixed_complex_layout()); }
// Appends complete double-precision records without replacing prior samples.
TEST(output, test_binary_append_real64) { EXPECT_EQ(0, test_output_binary_append_real64()); }
// Treats an empty double-precision batch as a successful no-op append.
TEST(output, test_binary_append_empty_real64) { EXPECT_EQ(0, test_output_binary_append_empty_real64()); }
// Registers a point probe and its declared output paths.
TEST(output, test_initialize_point_probe) { EXPECT_EQ(0, test_init_point_probe()); }
// Declares the incident field alongside point-probe time samples.
TEST(output, test_initialize_point_probe_with_incident) { EXPECT_EQ(0, test_init_point_probe_with_incident()); }
// Preserves the signed mixed-direction electric-field line integral.
TEST(output, test_line_probe_integral) { EXPECT_EQ(0, test_line_probe_integral()); }
// Keeps an empty line discoverable without fabricating measurements.
TEST(output, test_line_probe_empty_path) { EXPECT_EQ(0, test_line_probe_empty_path()); }
// Publishes every line sample to one flat text artifact.
TEST(output, test_line_probe_dat_output) { EXPECT_EQ(0, test_line_probe_dat_output()); }
TEST(output, test_line_probe_shared_interface_owner) { EXPECT_EQ(0, test_line_probe_shared_interface_owner()); }
// Records point-probe values and their corresponding timesteps.
TEST(output, test_update_point_probe_info) { EXPECT_EQ(0, test_update_point_probe()); }
// Applies explicit time windows and sampling periods to every scalar time probe.
TEST(output, test_update_time_probe_ranges) { EXPECT_EQ(0, test_update_time_probe_ranges()); }
// Flushes point-probe samples and resets serialized time data.
TEST(output, test_flush_point_probe_info) { EXPECT_EQ(0, test_flush_point_probe()); }
// Writes wire-current and wire-charge samples without binary sidecars.
TEST(output, test_flush_wire_probe_dat) { EXPECT_EQ(0, test_flush_wire_probe_dat()); }
// Writes bulk-current samples without a binary sidecar.
TEST(output, test_flush_bulk_probe_dat) { EXPECT_EQ(0, test_flush_bulk_probe_dat()); }
// Preserves point-probe data across consecutive flushes.
TEST(output, test_flush_multiple_point_probe_info) { EXPECT_EQ(0, test_multiple_flush_point_probe()); }
// Allocates a movie probe and creates its output directory.
TEST(output, test_init_movie_probe_for_pec_surface) { EXPECT_EQ(0, test_init_movie_probe()); }
// Samples movie-probe field components for a timestep.
TEST(output, test_update_movie_probe_for_pec_surface) { EXPECT_EQ(0, test_update_movie_probe()); }
// Flushes movie data and creates binary, HDF5, and XDMF artifacts.
TEST(output, test_flush_movie_probe_data) { EXPECT_EQ(0, test_flush_movie_probe()); }
// Closes a movie writer without publishing a JSON sidecar.
TEST(output, test_close_movie_probe_data) { EXPECT_EQ(0, test_close_movie_probe()); }
// Allocates a frequency-slice probe and its frequency buffers.
TEST(output, test_init_frequency_slice) { EXPECT_EQ(0, test_init_frequency_slice_probe()); }
// Computes frequency-slice values for a field gradient.
TEST(output, test_update_frequency_slice) { EXPECT_EQ(0, test_update_frequency_slice_probe()); }
// Does not create metadata or a manifest for a scalar-only run.
TEST(output, test_scalar_probe_has_no_manifest) { EXPECT_EQ(0, test_scalar_probe_has_no_manifest()); }
// Creates files and parent directories for nested output paths.
TEST(output, test_nested_output_path) { EXPECT_EQ(0, test_nested_output_path()); }
// Stores binary artifact encoding and lifecycle metadata.
TEST(output, test_artifact_contract) { EXPECT_EQ(0, test_output_artifact_contract()); }
// Declares output artifact kinds and relative paths.
TEST(output, test_declared_output_artifacts) { EXPECT_EQ(0, test_declared_output_artifacts()); }
// Validates output lifecycle terminal states and completeness.
TEST(output, test_lifecycle_contract) { EXPECT_EQ(0, test_output_lifecycle_contract()); }
// Keeps artifact identity and partition coverage equivalent across modes.
TEST(output, test_serial_distributed_equivalence) { EXPECT_EQ(0, test_output_serial_distributed_equivalence()); }
// Attaches a volumetric partition and selects serial publication fallback.
TEST(output, test_volumetric_partition_attachment) { EXPECT_EQ(0, test_volumetric_output_partition_attachment()); }
