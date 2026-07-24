integer function test_output_metadata_fragment_descriptors() bind(c) result(err)
   ! Verifies canonical and fragment descriptor publication and identity fields.
   use outputMetadata_m, only: publish_final_probe_metadata, OUTPUT_METADATA_SUCCESS
   use outputTypes_m, only: probe_metadata_t, output_artifact_t, OUTPUT_ARTIFACT_METADATA, &
                            OUTPUT_ARTIFACT_ROLE_FRAGMENT, OUTPUT_LIFECYCLE_COMPLETE
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: remove_folder
   implicit none

   type(probe_metadata_t) :: canonical, fragment
   character(len=8192) :: line
   integer :: ios, status, unit
   logical :: canonical_role, fragment_reference, fragment_role, fragment_parent, fragment_contributor

   err = 0
   canonical%probe_id = 'movie-001'
   canonical%quantity = 'Ex'
   canonical%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
   allocate(canonical%artifacts(1), canonical%fragment_descriptors(2))
   canonical%artifacts(1)%kind = OUTPUT_ARTIFACT_METADATA
   canonical%artifacts(1)%relative_path = 'movie-001.json'
   canonical%fragment_descriptors(1)%identity%parent_probe_id = canonical%probe_id
   canonical%fragment_descriptors(1)%identity%contributor_rank = 0
   canonical%fragment_descriptors(1)%relative_path = 'movie-001/rank-0.json'
   canonical%fragment_descriptors(2)%identity%parent_probe_id = canonical%probe_id
   canonical%fragment_descriptors(2)%identity%contributor_rank = 1
   canonical%fragment_descriptors(2)%relative_path = 'movie-001/rank-1.json'

   call publish_final_probe_metadata('testing metadata fragments/canonical.json', canonical, status)
   err = err + assert_integer_equal(status, OUTPUT_METADATA_SUCCESS, 'Canonical descriptor publication failed')
   canonical_role = .false.
   fragment_reference = .false.
   open(newunit=unit, file='testing metadata fragments/canonical.json', status='old', action='read', iostat=ios)
   do while (ios == 0)
      read(unit, '(A)', iostat=ios) line
      if (index(line, '"role":"canonical"') > 0) canonical_role = .true.
      if (index(line, '"contributor_rank":1') > 0 .and. index(line, 'movie-001/rank-1.json') > 0) then
         fragment_reference = .true.
      end if
   end do
   close(unit)
   err = err + assert_true(canonical_role, 'Canonical descriptor has no canonical artifact role')
   err = err + assert_true(fragment_reference, 'Canonical descriptor does not reference its fragment descriptor')

   fragment%probe_id = 'movie-001-rank-1'
   fragment%parent_probe_id = canonical%probe_id
   fragment%contributor_rank = 1
   fragment%quantity = 'Ex'
   fragment%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
   allocate(fragment%artifacts(1))
   fragment%artifacts(1)%kind = OUTPUT_ARTIFACT_METADATA
   fragment%artifacts(1)%role = OUTPUT_ARTIFACT_ROLE_FRAGMENT
   fragment%artifacts(1)%relative_path = 'movie-001/rank-1.bin'
   fragment%artifacts(1)%fragment%parent_probe_id = canonical%probe_id
   fragment%artifacts(1)%fragment%contributor_rank = 1

   call publish_final_probe_metadata('testing metadata fragments/fragment.json', fragment, status)
   err = err + assert_integer_equal(status, OUTPUT_METADATA_SUCCESS, 'Fragment descriptor publication failed')
   fragment_role = .false.
   fragment_parent = .false.
   fragment_contributor = .false.
   open(newunit=unit, file='testing metadata fragments/fragment.json', status='old', action='read', iostat=ios)
   do while (ios == 0)
      read(unit, '(A)', iostat=ios) line
      if (index(line, '"role":"fragment"') > 0) fragment_role = .true.
      if (index(line, '"parent_probe_id":"movie-001"') > 0) fragment_parent = .true.
      if (index(line, '"contributor_rank":1') > 0) fragment_contributor = .true.
   end do
   close(unit)
   err = err + assert_true(fragment_role, 'Fragment descriptor has no fragment artifact role')
   err = err + assert_true(fragment_parent .and. fragment_contributor, 'Fragment descriptor has no parent identity')
   call remove_folder('testing metadata fragments', ios)
end function test_output_metadata_fragment_descriptors
