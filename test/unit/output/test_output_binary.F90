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

integer function test_output_binary_native_precision() bind(c) result(err)
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use outputTypes_m, only: binary_numeric_representation_for_bytes, formatted_result_tolerance, &
                            BINARY_NUMERIC_REAL32, BINARY_NUMERIC_REAL64, &
                            BINARY_COMPONENTS_SCALAR_TIME, BINARY_COMPONENTS_FAR_FIELD
   use assertionTools_m, only: assert_integer_equal, assert_string_equal, assert_true
   implicit none

   err = 0
   err = err + assert_integer_equal(binary_numeric_representation_for_bytes([4, 8]), &
                                    BINARY_NUMERIC_REAL64, &
                                    'Mixed precision did not select real64')
   err = err + assert_integer_equal(binary_numeric_representation_for_bytes([4, 4]), &
                                    BINARY_NUMERIC_REAL32, &
                                    'Uniform real32 did not retain real32')
   err = err + assert_true(abs(formatted_result_tolerance(6) - 0.5e-6_8) < 1.0e-12_8, &
                            'Six decimal places did not produce a half-unit tolerance')
   err = err + assert_string_equal(BINARY_COMPONENTS_SCALAR_TIME, 'time,value', &
                                   'Scalar record fields changed')
   err = err + assert_string_equal(BINARY_COMPONENTS_FAR_FIELD, &
                                   'frequency,theta,phi,value.real,value.imag', &
                                   'Far-field record fields changed')
end function test_output_binary_native_precision

integer function test_output_binary_mixed_complex_layout() bind(c) result(err)
   use outputBinary_m, only: validate_binary_layout, BINARY_WRITER_SUCCESS
   use outputTypes_m, only: output_artifact_t, OUTPUT_ARTIFACT_BINARY, &
                            BINARY_ENDIAN_LITTLE, BINARY_NUMERIC_REAL64, &
                            BINARY_COMPLEX_REAL_IMAG
   use assertionTools_m, only: assert_integer_equal
   implicit none

   type(output_artifact_t) :: artifact
   integer :: status

   err = 0
   artifact%kind = OUTPUT_ARTIFACT_BINARY
   artifact%byte_order = BINARY_ENDIAN_LITTLE
   artifact%numeric_representation = BINARY_NUMERIC_REAL64
   artifact%complex_representation = BINARY_COMPLEX_REAL_IMAG
   artifact%record_bytes = 24
   artifact%component_order = 'frequency,value.real,value.imag'
   call validate_binary_layout(artifact, status)
   err = err + assert_integer_equal(status, BINARY_WRITER_SUCCESS, &
                                    'Coordinate-plus-complex record was rejected')
end function test_output_binary_mixed_complex_layout

integer function test_output_binary_append_real64() bind(c) result(err)
   use, intrinsic :: iso_fortran_env, only: real64
   use outputBinary_m, only: append_binary_real64, BINARY_WRITER_SUCCESS
   use outputTypes_m, only: output_artifact_t, OUTPUT_ARTIFACT_BINARY, &
                            BINARY_ENDIAN_LITTLE, BINARY_NUMERIC_REAL64, &
                            BINARY_COMPLEX_UNSPECIFIED
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: delete_file, file_exists
   implicit none

   type(output_artifact_t) :: artifact
   integer :: file_size, ios, status
   character(len=*), parameter :: path = 'testing binary/append-real64.bin'

   err = 0
   artifact%kind = OUTPUT_ARTIFACT_BINARY
   artifact%byte_order = BINARY_ENDIAN_LITTLE
   artifact%numeric_representation = BINARY_NUMERIC_REAL64
   artifact%complex_representation = BINARY_COMPLEX_UNSPECIFIED
   artifact%record_bytes = 16
   artifact%component_order = 'time,value'
   call append_binary_real64(path, artifact, [1.0_real64, 2.0_real64], status)
   err = err + assert_integer_equal(status, BINARY_WRITER_SUCCESS, 'First real64 append failed')
   call append_binary_real64(path, artifact, [3.0_real64, 4.0_real64], status)
   err = err + assert_integer_equal(status, BINARY_WRITER_SUCCESS, 'Second real64 append failed')
   inquire(file=path, size=file_size, iostat=ios)
   err = err + assert_true(file_exists(path) .and. ios == 0 .and. file_size == 32, &
                            'Real64 appends did not retain two complete records')
   call delete_file(path, ios)
end function test_output_binary_append_real64

integer function test_output_binary_append_empty_real64() bind(c) result(err)
   use, intrinsic :: iso_fortran_env, only: real64
   use outputBinary_m, only: append_binary_real64, BINARY_WRITER_SUCCESS
   use outputTypes_m, only: output_artifact_t, OUTPUT_ARTIFACT_BINARY, &
                            BINARY_ENDIAN_LITTLE, BINARY_NUMERIC_REAL64, &
                            BINARY_COMPLEX_UNSPECIFIED
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: delete_file, file_exists
   implicit none

   type(output_artifact_t) :: artifact
   real(real64), allocatable :: values(:)
   integer :: file_size, ios, status
   character(len=*), parameter :: path = 'testing binary/append-empty-real64.bin'

   err = 0
   artifact%kind = OUTPUT_ARTIFACT_BINARY
   artifact%byte_order = BINARY_ENDIAN_LITTLE
   artifact%numeric_representation = BINARY_NUMERIC_REAL64
   artifact%complex_representation = BINARY_COMPLEX_UNSPECIFIED
   artifact%record_bytes = 56
   artifact%component_order = 'time,x,y,z,Ex,Ey,Ez'
   allocate(values(0))

   call append_binary_real64(path, artifact, values, status)

   err = err + assert_integer_equal(status, BINARY_WRITER_SUCCESS, 'Empty real64 append failed')
   inquire(file=path, size=file_size, iostat=ios)
   err = err + assert_true(file_exists(path) .and. ios == 0 .and. file_size == 0, &
                            'Empty real64 append did not preserve an empty artifact')
   call delete_file(path, ios)
end function test_output_binary_append_empty_real64
