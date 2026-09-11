integer function test_output_publication_contract() bind(c) result(err)
   ! Verifies canonical and fragment identities and metadata completeness.
   use outputTypes_m, only: output_artifact_t, probe_metadata_t, output_artifact_identity_is_valid, &
                            probe_metadata_is_complete, OUTPUT_ARTIFACT_BINARY, &
                            OUTPUT_ARTIFACT_ROLE_FRAGMENT, OUTPUT_LIFECYCLE_COMPLETE
   use assertionTools_m, only: assert_true
   implicit none

   type(output_artifact_t) :: artifact
   type(probe_metadata_t) :: metadata

   err = 0
   artifact%kind = OUTPUT_ARTIFACT_BINARY
   err = err + assert_true(output_artifact_identity_is_valid(artifact), &
                           'Canonical artifact identity is invalid')

   artifact%role = OUTPUT_ARTIFACT_ROLE_FRAGMENT
   artifact%fragment%parent_probe_id = 'movie-001'
   artifact%fragment%contributor_rank = 1
   err = err + assert_true(output_artifact_identity_is_valid(artifact), &
                           'Fragment artifact identity is invalid')

    allocate(metadata%artifacts(1))
    metadata%probe_id = 'movie-001-fragment-1'
    metadata%parent_probe_id = 'movie-001'
    metadata%contributor_rank = 1
    metadata%quantity = 'Ex'
    artifact%relative_path = 'movie-001/rank-1.bin'
    metadata%artifacts(1) = artifact
   metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
   err = err + assert_true(probe_metadata_is_complete(metadata), &
                           'Complete fragment metadata is invalid')

   metadata%artifacts(1)%fragment%contributor_rank = -1
   err = err + assert_true(.not. probe_metadata_is_complete(metadata), &
                           'Fragment metadata accepts an invalid contributor')
end function test_output_publication_contract

integer function test_output_metadata_contract_edges() bind(c) result(err)
   ! Rejects incomplete identity, empty artifacts, and non-relative paths.
   use outputTypes_m, only: probe_metadata_t, OUTPUT_ARTIFACT_TEXT, OUTPUT_LIFECYCLE_COMPLETE, &
                            probe_metadata_is_complete
   use assertionTools_m, only: assert_true
   implicit none

   type(probe_metadata_t) :: metadata

   err = 0
   metadata%probe_id = 'point-001'
   metadata%quantity = 'Ex'
   metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
   allocate(metadata%artifacts(1))
   metadata%artifacts(1)%kind = OUTPUT_ARTIFACT_TEXT
   metadata%artifacts(1)%relative_path = 'point.dat'
   err = err + assert_true(probe_metadata_is_complete(metadata), &
                           'Valid zero-sample metadata is incomplete')

   metadata%artifacts(1)%relative_path = '/absolute/point.dat'
   err = err + assert_true(.not. probe_metadata_is_complete(metadata), &
                           'Absolute artifact path was accepted')

   metadata%artifacts(1)%relative_path = ''
   err = err + assert_true(.not. probe_metadata_is_complete(metadata), &
                           'Empty artifact path was accepted')

   deallocate(metadata%artifacts)
   allocate(metadata%artifacts(0))
   err = err + assert_true(.not. probe_metadata_is_complete(metadata), &
                           'Metadata without artifacts was accepted')
end function test_output_metadata_contract_edges
