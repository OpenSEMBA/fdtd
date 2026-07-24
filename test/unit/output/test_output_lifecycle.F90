integer function test_point_publication_plan() bind(c) result(err)
   ! Verifies canonical point-writer selection and unowned-point rejection.
   use output_m, only: prepare_point_publication_plan, publication_plan_allows_canonical_write
   use outputTypes_m, only: probe_publication_plan_t
   use outputCollective_m, only: OUTPUT_COLLECTIVE_SUCCESS, OUTPUT_COLLECTIVE_UNOWNED_POINT
   use assertionTools_m, only: assert_integer_equal, assert_true
   implicit none

   type(probe_publication_plan_t) :: plan
   integer :: status

   err = 0
   call prepare_point_publication_plan(plan, 1, 3, [.false., .true., .true.], status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Point plan preparation failed')
   err = err + assert_integer_equal(plan%canonical_writer_rank, 1, 'Unexpected point writer')
   err = err + assert_true(plan%local_eligible, 'Eligible rank was excluded from the point plan')
   err = err + assert_true(publication_plan_allows_canonical_write(plan), &
                            'Point writer is not allowed to publish')

   call prepare_point_publication_plan(plan, 2, 3, [.false., .true., .true.], status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_SUCCESS, 'Non-owner point plan preparation failed')
   err = err + assert_true(.not. publication_plan_allows_canonical_write(plan), &
                            'Non-owner rank is allowed to publish')

   call prepare_point_publication_plan(plan, 0, 3, [.false., .false., .false.], status)
   err = err + assert_integer_equal(status, OUTPUT_COLLECTIVE_UNOWNED_POINT, &
                                    'Unowned point plan is accepted')
end function test_point_publication_plan

integer function test_planned_metadata_publication() bind(c) result(err)
   ! Verifies only the planned canonical writer publishes metadata.
   use output_m, only: publish_planned_probe_metadata, OUTPUT_COORDINATION_SUCCESS
   use outputTypes_m, only: probe_metadata_t, probe_publication_plan_t, OUTPUT_ARTIFACT_TEXT
   use assertionTools_m, only: assert_integer_equal, assert_true
   use directoryUtils_m, only: file_exists, remove_folder
   implicit none

   type(probe_metadata_t) :: metadata
   type(probe_publication_plan_t) :: plan
   integer :: ios, status

   err = 0
   metadata%probe_id = 'point-001'
   metadata%quantity = 'Ex'
   allocate(metadata%artifacts(1))
   metadata%artifacts(1)%kind = OUTPUT_ARTIFACT_TEXT
   metadata%artifacts(1)%relative_path = 'point.dat'
   plan%canonical_writer_rank = 0
   plan%local_is_canonical_writer = .true.

   call publish_planned_probe_metadata('testing planned metadata/point.json', metadata, plan, status)
   err = err + assert_integer_equal(status, OUTPUT_COORDINATION_SUCCESS, 'Canonical metadata publication failed')
   err = err + assert_true(file_exists('testing planned metadata/point.json'), &
                            'Canonical metadata was not created')

   plan%local_is_canonical_writer = .false.
   call publish_planned_probe_metadata('testing planned metadata/non-owner.json', metadata, plan, status)
   err = err + assert_integer_equal(status, OUTPUT_COORDINATION_SUCCESS, 'Non-owner metadata coordination failed')
   err = err + assert_true(.not. file_exists('testing planned metadata/non-owner.json'), &
                            'Non-owner created canonical metadata')
   call remove_folder('testing planned metadata', ios)
end function test_planned_metadata_publication
