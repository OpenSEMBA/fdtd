integer function test_output_partition_shared_interfaces() bind(c) result(err)
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
