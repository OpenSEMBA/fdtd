integer function test_output_metadata_fragment_descriptors() bind(c) result(err)
   ! Verifies canonical and fragment descriptor publication and identity fields.
   use outputMetadata_m, only: publish_final_probe_metadata, OUTPUT_METADATA_SUCCESS
   use outputTypes_m, only: probe_metadata_t, output_artifact_t, OUTPUT_ARTIFACT_METADATA, &
                            OUTPUT_ARTIFACT_ROLE_FRAGMENT, OUTPUT_LIFECYCLE_COMPLETE
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: file_exists, remove_folder
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
    err = err + assert_true(.not. file_exists('testing metadata fragments/canonical.json.tmp'), &
                            'Temporary descriptor remains after publication')
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

integer function test_atomic_file_replacement() bind(c) result(err)
   ! Verifies replacement leaves the target complete and removes the temporary file.
   use directoryUtils_m, only: atomic_replace_file, create_file_with_path, file_exists, &
                               delete_file
   use assertionTools_m, only: assert_integer_equal, assert_true, assert_string_equal
   implicit none

   character(len=*), parameter :: target = 'testing atomic replacement/result.json'
   character(len=*), parameter :: temporary = 'testing atomic replacement/result.json.tmp'
   character(len=32) :: line
   integer :: ios, unit

   err = 0
   call create_file_with_path(target, ios)
   call create_file_with_path(temporary, ios)
   open(newunit=unit, file=temporary, status='old', action='write', position='rewind', iostat=ios)
   write(unit, '(a)') 'complete'
   close(unit)

   call atomic_replace_file(temporary, target, ios)
   err = err + assert_integer_equal(ios, 0, 'Atomic replacement failed')
   err = err + assert_true(file_exists(target), 'Replacement target does not exist')
   err = err + assert_true(.not. file_exists(temporary), 'Temporary file remains after replacement')

   open(newunit=unit, file=target, status='old', action='read', iostat=ios)
   read(unit, '(a)', iostat=ios) line
   close(unit)
   err = err + assert_string_equal(line, 'complete', 'Replacement target is incomplete')
   call delete_file(target, ios)
end function test_atomic_file_replacement

integer function test_json_path_escaping() bind(c) result(err)
   ! Preserves native Windows separators through JSON serialization.
   use outputMetadata_m, only: json_escape, json_unescape
   use assertionTools_m, only: assert_string_equal
   implicit none

   character(len=*), parameter :: windows_path = 'C:\temporary\point-probe\sample.dat'

   err = 0
   err = err + assert_string_equal(json_escape(windows_path), &
                                   'C:\\temporary\\point-probe\\sample.dat', &
                                   'Windows path was not JSON escaped')
   err = err + assert_string_equal(json_unescape(json_escape(windows_path)), windows_path, &
                                   'Windows path was not restored after JSON escaping')
end function test_json_path_escaping
