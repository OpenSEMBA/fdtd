integer function test_output_transport_serial() bind(c) result(err)
   use FDETYPES_m, only: RKIND
   use outputTransport_m, only: output_transport_t, init_output_transport, &
                                gather_point_eligibility, reduce_scalar_batch, transfer_flush_batch, &
                                OUTPUT_TRANSPORT_SUCCESS
   use assertionTools_m, only: assert_integer_equal, assert_true
   implicit none

   type(output_transport_t) :: transport
   integer :: status
   integer, allocatable :: counts(:), displacements(:)
   logical, allocatable :: eligibility(:)
   real(kind=RKIND), allocatable :: reduced_values(:), gathered_batch(:)

   err = 0
   call init_output_transport(transport, status=status)
   err = err + assert_integer_equal(status, OUTPUT_TRANSPORT_SUCCESS, 'Serial transport initialisation failed')
   err = err + assert_integer_equal(transport%rank, 0, 'Serial transport rank is not root')
   err = err + assert_integer_equal(transport%rank_count, 1, 'Serial transport has more than one rank')
   err = err + assert_true(.not. transport%distributed, 'Serial transport reports distributed execution')

   call gather_point_eligibility(transport, .true., eligibility, status)
   err = err + assert_integer_equal(status, OUTPUT_TRANSPORT_SUCCESS, 'Serial eligibility gathering failed')
   err = err + assert_true(size(eligibility) == 1 .and. eligibility(1), 'Serial eligibility was not preserved')

   call reduce_scalar_batch(transport, [1.0_RKIND, 2.0_RKIND], reduced_values, status)
   err = err + assert_integer_equal(status, OUTPUT_TRANSPORT_SUCCESS, 'Serial scalar reduction failed')
   err = err + assert_true(all(reduced_values == [1.0_RKIND, 2.0_RKIND]), 'Serial scalar values changed')

   call transfer_flush_batch(transport, [3.0_RKIND, 4.0_RKIND], gathered_batch, counts, displacements, status)
   err = err + assert_integer_equal(status, OUTPUT_TRANSPORT_SUCCESS, 'Serial batch transfer failed')
   err = err + assert_true(all(gathered_batch == [3.0_RKIND, 4.0_RKIND]), 'Serial batch values changed')
   err = err + assert_true(size(counts) == 1 .and. counts(1) == 2, 'Serial batch count is incorrect')
   err = err + assert_true(size(displacements) == 1 .and. displacements(1) == 0, &
                           'Serial batch displacement is incorrect')
end function test_output_transport_serial
