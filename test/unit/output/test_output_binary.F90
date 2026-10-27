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
