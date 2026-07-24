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
   metadata%artifacts(1) = artifact
   metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
   err = err + assert_true(probe_metadata_is_complete(metadata), &
                           'Complete fragment metadata is invalid')

   metadata%artifacts(1)%fragment%contributor_rank = -1
   err = err + assert_true(.not. probe_metadata_is_complete(metadata), &
                           'Fragment metadata accepts an invalid contributor')
end function test_output_publication_contract
