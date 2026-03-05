!--------------------------------------------------------------------------------
! Test: count_required_coords
!--------------------------------------------------------------------------------
integer function test_count_required_coords() bind(c) result(err)
   use FDETYPES
   use outputTypes
   use mod_volumicProbeUtils
   use mod_assertionTools
   use mod_testOutputUtils
   implicit none

   type(cell_coordinate_t) :: lowerBound, upperBound
   type(problem_info_t)    :: problemInfo
   integer(kind=SINGLE)    :: count
   integer                 :: test_err = 0
   integer, allocatable    :: dummy_coords(:,:)

   ! Setup test case: 3x3x3 domain (1..3)
   lowerBound = cell_coordinate_t(1, 1, 1)
   upperBound = cell_coordinate_t(3, 3, 3)
   
   call setup_dummy_problem_info(problemInfo)

   ! Test Case 1: Field Request (iExC)
   call find_and_store_important_coords(lowerBound, upperBound, iExC, problemInfo, count, dummy_coords)
   
   ! Expected: 3*3*3 = 27 points
   test_err = test_err + assert_integer_equal(count, 27_SINGLE, "Failed count for iExC")

   if (allocated(dummy_coords)) deallocate(dummy_coords)
   call clean_dummy_problem_info(problemInfo)
   err = test_err
end function test_count_required_coords

!--------------------------------------------------------------------------------
! Test: store_required_coords
!--------------------------------------------------------------------------------
integer function test_store_required_coords() bind(c) result(err)
   use FDETYPES
   use outputTypes
   use mod_outputUtils
   use mod_volumicProbeUtils
   use mod_assertionTools
   use mod_testOutputUtils
   implicit none

   type(cell_coordinate_t) :: lowerBound, upperBound
   type(problem_info_t)    :: problemInfo
   integer(kind=SINGLE)    :: nPoints
   integer(kind=SINGLE), allocatable :: stored_coords(:,:)
   integer                 :: test_err = 0

   lowerBound = new_cell_coordinate(1, 1, 1)
   upperBound = new_cell_coordinate(2, 2, 2)
   call setup_dummy_problem_info(problemInfo)

   call find_and_store_important_coords(lowerBound, upperBound, iHyC, problemInfo, nPoints, stored_coords)

   test_err = test_err + assert_integer_equal(nPoints, 8_SINGLE, "Failed nPoints for iHyC")
   
   if (allocated(stored_coords)) then
      test_err = test_err + assert_integer_equal(int(size(stored_coords, 2), SINGLE), 8_SINGLE, "Allocated coords size error")
      ! Verify first coord is (1,1,1)
      test_err = test_err + assert_integer_equal(stored_coords(1,1), 1_SINGLE, "First x coord mismatch")
      test_err = test_err + assert_integer_equal(stored_coords(2,1), 1_SINGLE, "First y coord mismatch")
      test_err = test_err + assert_integer_equal(stored_coords(3,1), 1_SINGLE, "First z coord mismatch")
      deallocate(stored_coords)
   else
      print *, "Coords not allocated."
      test_err = test_err + 1
   end if

   call clean_dummy_problem_info(problemInfo)
   err = test_err
end function test_store_required_coords

!--------------------------------------------------------------------------------
! Test: isValidPointForCurrent
!--------------------------------------------------------------------------------
integer function test_is_valid_point_current() bind(c) result(err)
   use FDETYPES
   use outputTypes
   use mod_volumicProbeUtils
   use mod_testOutputUtils
   implicit none

   type(problem_info_t) :: problemInfo
   integer :: test_err = 0
   logical :: valid

   call setup_dummy_problem_info(problemInfo)

   ! By default, our dummy setup has NO PEC and NO Wires.
   ! So isValidPointForCurrent should be FALSE (as it requires PEC or Wire)
   valid = isValidPointForCurrent(iCur, 1, 1, 1, problemInfo)
   
   if (valid) then
       print *, "Expected False for empty space current probe (no PEC/Wire)"
       test_err = test_err + 1
   end if

   call clean_dummy_problem_info(problemInfo)
   err = test_err
end function test_is_valid_point_current

!--------------------------------------------------------------------------------
! Test: isValidPointForField
!--------------------------------------------------------------------------------
integer function test_is_valid_point_field() bind(c) result(err)
     use FDETYPES
     use outputTypes
     use mod_volumicProbeUtils
     use mod_testOutputUtils
     implicit none

     type(problem_info_t) :: problemInfo
     integer :: test_err = 0
     logical :: valid

     call setup_dummy_problem_info(problemInfo)

     ! Point inside boundary
     valid = isValidPointForField(iEx, 5, 5, 5, problemInfo)
     if (.not. valid) then
         print *, "Expected True for field probe in bounds"
         test_err = test_err + 1
     end if

     ! Point outside boundary (-1)
     valid = isValidPointForField(iEx, -1, 5, 5, problemInfo)
     if (valid) then
          print *, "Expected False for field probe out of bounds"
          test_err = test_err + 1
     end if

     call clean_dummy_problem_info(problemInfo)
     err = test_err
end function test_is_valid_point_field
