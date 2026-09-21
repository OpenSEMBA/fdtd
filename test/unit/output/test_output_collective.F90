integer function test_output_collective_contract() bind(c) result(err)
   ! Verifies participant selection, ownership, and publication modes.
   use outputCollective_m
   use outputDecomposition_m, only: output_partition_t
   use assertionTools_m, only: assert_integer_equal, assert_true
   implicit none

   type(output_collective_t) :: collective
   type(output_partition_t) :: partition
   integer, allocatable :: participants(:)
   integer :: owner, publication_mode, status
   logical :: local_participates

   err = 0
   call init_output_collective(collective, 0, 1, 0, .false., status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Serial collective initialisation failed')

   call select_output_participants(collective, [.true.], participants, owner, status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Serial participant selection failed')
   err = err + assert_integer_equal(size(participants), 1, 'Serial participant count')
   err = err + assert_integer_equal(owner, 0, 'Serial owner')

   call select_output_publication_mode(collective, publication_mode)
   err = err + assert_integer_equal(publication_mode, OUTPUT_PUBLICATION_ROOT_AGGREGATION, &
                                    'Serial mode is not root aggregation')

   call init_output_collective(collective, 1, 3, 0, .true., status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Distributed collective initialisation failed')
   call select_output_participants(collective, [.false., .true., .true.], participants, owner, status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Participant selection failed')
   err = err + assert_true(all(participants == [1, 2]), 'Participants are not rank ordered')
   err = err + assert_integer_equal(owner, 1, 'Owner is not deterministic')

   call validate_output_ownership(collective, [1, 2], 1, status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Valid ownership was rejected')
   partition%has_data = .true.
   call prepare_output_partition_publication(collective, [1, 2], 1, partition, &
                                             local_participates, publication_mode, status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Valid partition was rejected')
   err = err + assert_true(local_participates, 'Participant was not selected for publication')
   err = err + assert_integer_equal(publication_mode, OUTPUT_PUBLICATION_COLLECTIVE, &
                                    'Distributed collective mode was not selected')
   call validate_output_ownership(collective, [2, 1], 1, status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_INVALID_PARTICIPANTS, &
                                    'Unordered participants were accepted')
   call validate_output_ownership(collective, [1, 2], 2, status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_INVALID_OWNER, &
                                    'Non-deterministic owner was accepted')

end function test_output_collective_contract
