integer function test_output_binary_fragment_layout() bind(c) result(err)
   ! Verifies binary fragment layouts require a valid fragment identity.
   use outputBinary_m, only: validate_binary_layout, BINARY_WRITER_SUCCESS, BINARY_WRITER_INVALID_LAYOUT
   use outputTypes_m, only: output_artifact_t, OUTPUT_ARTIFACT_BINARY, OUTPUT_ARTIFACT_ROLE_FRAGMENT, &
                            BINARY_ENDIAN_LITTLE, BINARY_NUMERIC_REAL32, BINARY_COMPLEX_UNSPECIFIED
   use assertionTools_m, only: assert_integer_equal
   implicit none

   type(output_artifact_t) :: artifact
   integer :: status

   err = 0
   artifact%kind = OUTPUT_ARTIFACT_BINARY
   artifact%role = OUTPUT_ARTIFACT_ROLE_FRAGMENT
   artifact%byte_order = BINARY_ENDIAN_LITTLE
   artifact%numeric_representation = BINARY_NUMERIC_REAL32
    artifact%complex_representation = BINARY_COMPLEX_UNSPECIFIED
    artifact%record_bytes = 4

    call validate_binary_layout(artifact, status)
   err = err + assert_integer_equal(status, BINARY_WRITER_INVALID_LAYOUT, &
                                    'Fragment layout without an identity was accepted')

    artifact%fragment%parent_probe_id = 'movie-001'
    artifact%fragment%contributor_rank = 1
    artifact%component_order = 'Ex'
    call validate_binary_layout(artifact, status)
    err = err + assert_integer_equal(status, BINARY_WRITER_SUCCESS, 'Valid fragment layout was rejected')

    artifact%component_order = ''
    call validate_binary_layout(artifact, status)
    err = err + assert_integer_equal(status, BINARY_WRITER_INVALID_LAYOUT, &
                                     'Binary layout without component order was accepted')
end function test_output_binary_fragment_layout
