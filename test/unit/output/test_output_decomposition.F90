integer function test_output_partition_shared_interfaces() bind(c) result(err)
   ! Verifies unique ownership of shared interfaces across output partitions.
   use iso_fortran_env, only: int64
   use FDETYPES_m, only: iEx, iEy, iHz, limit_t
   use outputTypes_m, only: cell_coordinate_t
   use outputDecomposition_m
   use assertionTools_m, only: assert_integer_equal, assert_true
   implicit none

   type(output_partition_t) :: partition
   type(cell_coordinate_t) :: request_lower, request_upper
   type(limit_t) :: global_bounds, local_sweep
   integer :: components(3), expected_lower(3), expected_upper(3)
   integer :: component_index, coverage(0:10), k, rank, status

   err = 0
   components = [iEx, iEy, iHz]
   expected_lower = [0, 4, 7]
   expected_upper = [3, 6, 10]
   request_lower = cell_coordinate_t(0, 0, 0)
   request_upper = cell_coordinate_t(2, 1, 10)
   global_bounds = limit_t(0, 2, 0, 1, 0, 10, 3, 2, 11)

   do component_index = 1, size(components)
      coverage = 0
      do rank = 0, 2
         select case (rank)
         case (0)
            local_sweep = limit_t(0, 2, 0, 1, 0, 4, 3, 2, 5)
         case (1)
            local_sweep = limit_t(0, 2, 0, 1, 4, 7, 3, 2, 4)
         case (2)
            local_sweep = limit_t(0, 2, 0, 1, 7, 10, 3, 2, 4)
         end select

         call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                                     components(component_index), rank, 3, partition, status)

         err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, "Shared partition status")
         err = err + assert_true(partition%has_data, "Shared partition unexpectedly empty")
         err = err + assert_integer_equal(partition%local_lower%z, expected_lower(rank + 1), &
                                           "Shared partition lower z")
         err = err + assert_integer_equal(partition%local_upper%z, expected_upper(rank + 1), &
                                           "Shared partition upper z")
         err = err + assert_true(partition%file_offset(3) == int(expected_lower(rank + 1), int64), &
                                 "Shared partition file offset")

         do k = partition%local_lower%z, partition%local_upper%z
            coverage(k) = coverage(k) + 1
         end do
      end do

      do k = 0, 10
         err = err + assert_integer_equal(coverage(k), 1, "Shared interface ownership is not unique")
      end do
   end do
end function test_output_partition_shared_interfaces

integer function test_output_partition_disjoint_intervals() bind(c) result(err)
   ! Verifies complete, non-overlapping coverage of disjoint rank intervals.
   use FDETYPES_m, only: iEz, iHx, iHy, limit_t
   use outputTypes_m, only: cell_coordinate_t
   use outputDecomposition_m
   use assertionTools_m, only: assert_integer_equal, assert_true
   implicit none

   type(output_partition_t) :: partition
   type(cell_coordinate_t) :: request_lower, request_upper
   type(limit_t) :: global_bounds, local_sweep
   integer :: components(3), coverage(0:9), component_index, k, rank, status

   err = 0
   components = [iEz, iHx, iHy]
   request_lower = cell_coordinate_t(0, 0, 0)
   request_upper = cell_coordinate_t(2, 1, 9)
   global_bounds = limit_t(0, 2, 0, 1, 0, 9, 3, 2, 10)

   do component_index = 1, size(components)
      coverage = 0
      do rank = 0, 2
         select case (rank)
         case (0)
            local_sweep = limit_t(0, 2, 0, 1, 0, 3, 3, 2, 4)
         case (1)
            local_sweep = limit_t(0, 2, 0, 1, 4, 6, 3, 2, 3)
         case (2)
            local_sweep = limit_t(0, 2, 0, 1, 7, 9, 3, 2, 3)
         end select

         call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                                     components(component_index), rank, 3, partition, status)

         err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, "Disjoint partition status")
         err = err + assert_true(partition%has_data, "Disjoint partition unexpectedly empty")
         do k = partition%local_lower%z, partition%local_upper%z
            coverage(k) = coverage(k) + 1
         end do
      end do

      do k = 0, 9
         err = err + assert_integer_equal(coverage(k), 1, "Disjoint interval ownership is not unique")
      end do
   end do
end function test_output_partition_disjoint_intervals

integer function test_output_partition_clipping_and_large_shape() bind(c) result(err)
   ! Verifies request clipping, empty ranks, large shapes, and invalid components.
   use iso_fortran_env, only: int64
   use FDETYPES_m, only: iEx, limit_t
   use outputTypes_m, only: cell_coordinate_t
   use outputDecomposition_m
   use assertionTools_m, only: assert_integer_equal, assert_true
   implicit none

   type(output_partition_t) :: partition
   type(cell_coordinate_t) :: request_lower, request_upper
   type(limit_t) :: global_bounds, local_sweep
   integer :: status

   err = 0
   request_lower = cell_coordinate_t(-5, 4, 12)
   request_upper = cell_coordinate_t(2, 9, 18)
   global_bounds = limit_t(-2, 5, 3, 7, 10, 20, 8, 5, 11)
   local_sweep = global_bounds

   call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                               iEx, 0, 1, partition, status)

   err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, "Clipped partition status")
   err = err + assert_true(partition%has_data, "Clipped serial partition unexpectedly empty")
   err = err + assert_true(all(partition%global_shape == [5_int64, 4_int64, 7_int64]), &
                           "Clipped global shape")
   err = err + assert_true(all(partition%local_shape == partition%global_shape), &
                           "Serial rank does not own complete partition")
   err = err + assert_true(all(partition%file_offset == 0_int64), "Serial file offset is not zero")

   request_lower = cell_coordinate_t(0, 0, 8)
   request_upper = cell_coordinate_t(2, 1, 10)
   global_bounds = limit_t(0, 2, 0, 1, 0, 10, 3, 2, 11)
   local_sweep = limit_t(0, 2, 0, 1, 0, 4, 3, 2, 5)

   call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                               iEx, 0, 2, partition, status)

   err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, "Empty rank partition status")
   err = err + assert_true(.not. partition%has_data, "Non-intersecting rank owns output data")
   err = err + assert_true(all(partition%global_shape == [3_int64, 2_int64, 3_int64]), &
                           "Empty rank lost the global dataset shape")

   request_lower = cell_coordinate_t(-2000000000, 0, 0)
   request_upper = cell_coordinate_t(2000000000, 0, 0)
   global_bounds = limit_t(-2000000000, 2000000000, 0, 0, 0, 0, 0, 1, 1)
   local_sweep = global_bounds

   call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                               iEx, 0, 1, partition, status)

   err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, "Large partition status")
   err = err + assert_true(partition%global_shape(1) == 4000000001_int64, &
                           "Large partition shape overflowed")

   call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                               -1, 0, 1, partition, status)
   err = err + assert_integer_equal(status, OUTPUT_PARTITION_INVALID_COMPONENT, &
                                    "Invalid component was accepted")
end function test_output_partition_clipping_and_large_shape

integer function test_output_point_partition_eligibility() bind(c) result(err)
   ! Verifies point eligibility at partition boundaries and for empty partitions.
   use outputTypes_m, only: cell_coordinate_t
   use outputDecomposition_m, only: output_partition_t, point_is_in_partition
   use assertionTools_m, only: assert_true
   implicit none

   type(output_partition_t) :: partition

   err = 0
   partition%has_data = .true.
   partition%local_lower = cell_coordinate_t(2, 3, 4)
   partition%local_upper = cell_coordinate_t(5, 6, 7)
   err = err + assert_true(point_is_in_partition(cell_coordinate_t(2, 3, 4), partition), &
                            'Lower partition boundary is not eligible')
   err = err + assert_true(point_is_in_partition(cell_coordinate_t(5, 6, 7), partition), &
                            'Upper partition boundary is not eligible')
   err = err + assert_true(.not. point_is_in_partition(cell_coordinate_t(6, 6, 7), partition), &
                            'Point outside partition is eligible')

   partition%has_data = .false.
   err = err + assert_true(.not. point_is_in_partition(cell_coordinate_t(2, 3, 4), partition), &
                            'Empty partition is eligible for a point')
end function test_output_point_partition_eligibility

integer function test_output_partition_all_components_cover_volume() bind(c) result(err)
   ! Verifies unique volume coverage for every electric and magnetic component.
   use FDETYPES_m, only: iEx, iEy, iEz, iHx, iHy, iHz, limit_t
   use outputTypes_m, only: cell_coordinate_t
   use outputDecomposition_m, only: output_partition_t, build_output_partition, OUTPUT_PARTITION_SUCCESS
   use assertionTools_m, only: assert_integer_equal
   implicit none

   type(output_partition_t) :: partition
   type(cell_coordinate_t) :: request_lower, request_upper
   type(limit_t) :: global_bounds, local_sweep
   integer :: component, components(6), coverage(0:5), rank, status, z

   err = 0
   components = [iEx, iEy, iEz, iHx, iHy, iHz]
   request_lower = cell_coordinate_t(0, 0, 0)
   request_upper = cell_coordinate_t(1, 1, 5)
   global_bounds = limit_t(0, 1, 0, 1, 0, 5, 2, 2, 6)

   do component = 1, size(components)
      coverage = 0
      do rank = 0, 1
         if (rank == 0) then
            if (any(components(component) == [iEx, iEy, iHz])) then
               local_sweep = limit_t(0, 1, 0, 1, 0, 3, 2, 2, 4)
            else
               local_sweep = limit_t(0, 1, 0, 1, 0, 2, 2, 2, 3)
            end if
         else
            local_sweep = limit_t(0, 1, 0, 1, 3, 5, 2, 2, 3)
         end if

         call build_output_partition(request_lower, request_upper, global_bounds, local_sweep, &
                                     components(component), rank, 2, partition, status)
         err = err + assert_integer_equal(status, OUTPUT_PARTITION_SUCCESS, &
                                          'All-component partition status')
         do z = partition%local_lower%z, partition%local_upper%z
            coverage(z) = coverage(z) + 1
         end do
      end do

      do z = 0, 5
         err = err + assert_integer_equal(coverage(z), 1, 'All-component ownership is not unique')
      end do
   end do
end function test_output_partition_all_components_cover_volume
