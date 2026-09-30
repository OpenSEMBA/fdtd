! Unit tests of the slanted MTLN wire geometry builder.
integer function test_slanted_segments_aligned() bind(C) result(error_cnt)
   use FDETYPES_m, only: rkind
   use NFDETypes_m, only: Desplazamiento_t
   use mtln_types_m, only: segment_t, DIRECTION_X_POS
   use slanted_wire_segments_m
   use smbjson_testingTools

   implicit none

   type(Desplazamiento_t) :: grid
   type(segment_t), dimension(:), allocatable :: segments
   real(kind=rkind), dimension(3, 3) :: points
   real(kind=rkind) :: total_length
   integer :: i

   error_cnt = 0

   call makeUniformGrid(grid, 10, 0.1_rkind)
   ! An axis-aligned wire on a grid line must reduce to the staircase model:
   ! one division per cell, a single unit coupling on the wire component.
   points(:, 1) = [1.0_rkind, 2.0_rkind, 3.0_rkind]
   points(:, 2) = [2.0_rkind, 2.0_rkind, 3.0_rkind]
   points(:, 3) = [4.0_rkind, 2.0_rkind, 3.0_rkind]

   call buildSlantedSegments(points, grid, segments)

   if (size(segments) /= 3) error_cnt = error_cnt + 1
   if (size(segments) == 3) then
      do i = 1, 3
         if (.not. segments(i)%is_slanted) error_cnt = error_cnt + 1
         if (segments(i)%orientation /= DIRECTION_X_POS) error_cnt = error_cnt + 1
         if (.not. expect_near(segments(i)%direction_cosines(1), 1.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
         if (size(segments(i)%couplings) /= 1) then
            error_cnt = error_cnt + 1
            cycle
         end if
         if (segments(i)%couplings(1)%component /= DIRECTION_X_POS) error_cnt = error_cnt + 1
         if (.not. expect_near(segments(i)%couplings(1)%weight, 1.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
         if (segments(i)%couplings(1)%i /= i) error_cnt = error_cnt + 1
         if (segments(i)%couplings(1)%j /= 2) error_cnt = error_cnt + 1
         if (segments(i)%couplings(1)%k /= 3) error_cnt = error_cnt + 1
         if (.not. expect_near(segments(i)%couplings(1)%chord, 0.1_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
      end do
      if (.not. expect_near(norm2(segments(1)%position_begin - [0.1_rkind, 0.2_rkind, 0.3_rkind]), &
                            0.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
      if (.not. expect_near(norm2(segments(3)%position_end - [0.4_rkind, 0.2_rkind, 0.3_rkind]), &
                            0.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
   end if

   ! Reversed direction: same couplings with the opposite orientation.
   points(:, 1) = [4.0_rkind, 2.0_rkind, 3.0_rkind]
   points(:, 2) = [1.0_rkind, 2.0_rkind, 3.0_rkind]
   call buildSlantedSegments(points(:, 1:2), grid, segments)
   if (size(segments) /= 3) error_cnt = error_cnt + 1
   if (size(segments) == 3) then
      do i = 1, 3
         if (segments(i)%orientation /= -DIRECTION_X_POS) error_cnt = error_cnt + 1
         if (size(segments(i)%couplings) /= 1) then
            error_cnt = error_cnt + 1
            cycle
         end if
         if (.not. expect_near(segments(i)%couplings(1)%weight, -1.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
         if (segments(i)%couplings(1)%i /= 4 - i) error_cnt = error_cnt + 1
      end do
   end if

   ! Fractional but axis-aligned coordinates: still one unit coupling per cell.
   points(:, 1) = [1.5_rkind, 2.0_rkind, 3.0_rkind]
   points(:, 2) = [4.5_rkind, 2.0_rkind, 3.0_rkind]
   call buildSlantedSegments(points(:, 1:2), grid, segments)
   if (size(segments) /= 4) error_cnt = error_cnt + 1
   if (size(segments) == 4) then
      total_length = 0.0_rkind
      do i = 1, 4
         if (size(segments(i)%couplings) /= 1) then
            error_cnt = error_cnt + 1
            cycle
         end if
         total_length = total_length + segments(i)%couplings(1)%chord
         if (.not. expect_near(segments(i)%couplings(1)%weight, 1.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
         if (segments(i)%couplings(1)%i /= i) error_cnt = error_cnt + 1
      end do
      if (.not. expect_near(total_length, 0.3_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
      if (.not. expect_near(norm2(segments(1)%position_begin - [0.15_rkind, 0.2_rkind, 0.3_rkind]), &
                            0.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
   end if

   call destroyGrid(grid)
end function

integer function test_slanted_segments_diagonal() bind(C) result(error_cnt)
   use FDETYPES_m, only: rkind
   use NFDETypes_m, only: Desplazamiento_t
   use mtln_types_m, only: segment_t
   use slanted_wire_segments_m
   use smbjson_testingTools

   implicit none

   type(Desplazamiento_t) :: grid
   type(segment_t), dimension(:), allocatable :: segments
   real(kind=rkind), dimension(3, 2) :: points
   real(kind=rkind) :: projection, total_length
   integer :: i, k

   error_cnt = 0

   call makeUniformGrid(grid, 10, 0.1_rkind)
   ! 45 degree wire crossing cell corners: four divisions of 0.7071 and
   ! 1.4142 cell diagonals.
   points(:, 1) = [0.5_rkind, 1.0_rkind, 0.5_rkind]
   points(:, 2) = [3.5_rkind, 1.0_rkind, 3.5_rkind]

   call buildSlantedSegments(points, grid, segments)

   if (size(segments) /= 4) error_cnt = error_cnt + 1
   total_length = 0.0_rkind
   do i = 1, size(segments)
      total_length = total_length + norm2(segments(i)%position_end - segments(i)%position_begin)
      if (.not. segments(i)%is_slanted) error_cnt = error_cnt + 1
      ! The weights project the wire direction: sum(w * t) = |t|^2 = 1.
      projection = 0.0_rkind
      do k = 1, size(segments(i)%couplings)
         projection = projection + segments(i)%couplings(k)%weight* &
                      segments(i)%direction_cosines(segments(i)%couplings(k)%component)
      end do
      if (.not. expect_near(projection, 1.0_rkind, 1e-5_rkind)) error_cnt = error_cnt + 1
   end do
   if (.not. expect_near(total_length, 0.3_rkind*sqrt(2.0_rkind), 1e-6_rkind)) error_cnt = error_cnt + 1

   ! Fully three-dimensional diagonal.
   points(:, 1) = [0.2_rkind, 0.3_rkind, 0.4_rkind]
   points(:, 2) = [3.8_rkind, 3.7_rkind, 3.6_rkind]
   call buildSlantedSegments(points, grid, segments)
   total_length = 0.0_rkind
   do i = 1, size(segments)
      total_length = total_length + norm2(segments(i)%position_end - segments(i)%position_begin)
      projection = 0.0_rkind
      do k = 1, size(segments(i)%couplings)
         if (segments(i)%couplings(k)%component < 1 .or. segments(i)%couplings(k)%component > 3) then
            error_cnt = error_cnt + 1
         end if
         if (segments(i)%couplings(k)%i < 0 .or. segments(i)%couplings(k)%i > 10 .or. &
             segments(i)%couplings(k)%j < 0 .or. segments(i)%couplings(k)%j > 10 .or. &
             segments(i)%couplings(k)%k < 0 .or. segments(i)%couplings(k)%k > 10) then
            error_cnt = error_cnt + 1
         end if
         projection = projection + segments(i)%couplings(k)%weight* &
                      segments(i)%direction_cosines(segments(i)%couplings(k)%component)
      end do
      if (.not. expect_near(projection, 1.0_rkind, 1e-5_rkind)) error_cnt = error_cnt + 1
   end do
   if (.not. expect_near(total_length, 0.1_rkind*sqrt(34.76_rkind), 1e-5_rkind)) error_cnt = error_cnt + 1

   call destroyGrid(grid)
end function

integer function test_slanted_segments_find_node() bind(C) result(error_cnt)
   use FDETYPES_m, only: rkind
   use NFDETypes_m, only: Desplazamiento_t
   use mtln_types_m, only: segment_t
   use slanted_wire_segments_m
   use smbjson_testingTools

   implicit none

   type(Desplazamiento_t) :: grid
   type(segment_t), dimension(:), allocatable :: segments
   real(kind=rkind), dimension(3, 3) :: points
   real(kind=rkind), dimension(3) :: center
   integer :: index, first, last

   error_cnt = 0

   call makeUniformGrid(grid, 20, 0.1_rkind)
   points(:, 1) = [8.5_rkind, 11.0_rkind, 7.669873_rkind]
   points(:, 2) = [11.0_rkind, 11.0_rkind, 12.0_rkind]
   points(:, 3) = [13.5_rkind, 11.0_rkind, 16.330127_rkind]
   center = points(:, 2)

   call buildSlantedSegments(points, grid, segments)

   first = findIndexInSlantedSegments(points(:, 1), grid, segments)
   last = findIndexInSlantedSegments(points(:, 3), grid, segments)
   index = findIndexInSlantedSegments(center, grid, segments)
   if (first /= 1) error_cnt = error_cnt + 1
   if (last /= size(segments) + 1) error_cnt = error_cnt + 1
   if (index < 2 .or. index > size(segments)) then
      error_cnt = error_cnt + 1
   else
      if (norm2(segments(index)%position_begin - center*grid%desX(1)) > 1e-5_rkind) error_cnt = error_cnt + 1
   end if

   call destroyGrid(grid)
end function

integer function test_slanted_segments_graded() bind(C) result(error_cnt)
   use FDETYPES_m, only: rkind
   use NFDETypes_m, only: Desplazamiento_t
   use mtln_types_m, only: segment_t
   use slanted_wire_segments_m
   use smbjson_testingTools

   implicit none

   type(Desplazamiento_t) :: grid
   type(segment_t), dimension(:), allocatable :: segments
   real(kind=rkind), dimension(3, 2) :: points
   real(kind=rkind), dimension(4), parameter :: steps = [0.1_rkind, 0.2_rkind, 0.3_rkind, 0.4_rkind]
   integer :: i

   error_cnt = 0

   grid%nx = 4
   grid%ny = 4
   grid%nz = 4
   grid%mx2 = 4
   grid%my2 = 4
   grid%mz2 = 4
   allocate (grid%desX(0:3), grid%desY(1:1), grid%desZ(1:1))
   grid%desX = steps
   grid%desY = 0.1_rkind
   grid%desZ = 0.1_rkind

   points(:, 1) = [0.0_rkind, 0.0_rkind, 0.0_rkind]
   points(:, 2) = [4.0_rkind, 0.0_rkind, 0.0_rkind]
   call buildSlantedSegments(points, grid, segments)

   if (size(segments) /= 4) error_cnt = error_cnt + 1
   if (size(segments) == 4) then
      do i = 1, 4
         if (.not. expect_near(segments(i)%couplings(1)%chord, steps(i), 1e-6_rkind)) error_cnt = error_cnt + 1
         if (segments(i)%couplings(1)%i /= i - 1) error_cnt = error_cnt + 1
      end do
      if (.not. expect_near(norm2(segments(1)%position_begin), 0.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
      if (.not. expect_near(norm2(segments(4)%position_end - [1.0_rkind, 0.0_rkind, 0.0_rkind]), &
                            0.0_rkind, 1e-6_rkind)) error_cnt = error_cnt + 1
   end if

   deallocate (grid%desX, grid%desY, grid%desZ)
end function

#ifdef CompileWithMTLN
integer function test_read_slanted_wire() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools
   use mtln_types_m, only: cable_t

   implicit none

   character(len=*), parameter :: filename = &
      PATH_TO_TEST_DATA//'cases/holland_slanted/holland1981_slanted.fdtd.json'
   type(Parseador_t) :: problem
   type(parser_t) :: parser
   class(cable_t), pointer :: cable
   integer :: i, slanted_segments

   err = 0

   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   if (.not. allocated(problem%mtln%cables)) then
      err = err + 1
      return
   end if
   if (size(problem%mtln%cables) /= 1) err = err + 1
   cable => problem%mtln%cables(1)%ptr
   if (.not. allocated(cable%segments)) then
      err = err + 1
      return
   end if

   slanted_segments = 0
   do i = 1, size(cable%segments)
      if (cable%segments(i)%is_slanted) slanted_segments = slanted_segments + 1
   end do
   if (slanted_segments /= size(cable%segments)) err = err + 1
   if (slanted_segments == 0) err = err + 1
   ! The rotated wire keeps its 1 m length.
   if (abs(sum(cable%step_size) - 1.0_RKIND) > 1e-4_RKIND) err = err + 1
end function
#endif

subroutine makeUniformGrid(grid, n, step)
   use FDETYPES_m, only: rkind
   use NFDETypes_m, only: Desplazamiento_t
   type(Desplazamiento_t), intent(out) :: grid
   integer, intent(in) :: n
   real(kind=rkind), intent(in) :: step

   grid%nx = n
   grid%ny = n
   grid%nz = n
   grid%mx2 = n
   grid%my2 = n
   grid%mz2 = n
   allocate (grid%desX(1:1), grid%desY(1:1), grid%desZ(1:1))
   grid%desX = step
   grid%desY = step
   grid%desZ = step
end subroutine

subroutine destroyGrid(grid)
   use NFDETypes_m, only: Desplazamiento_t
   type(Desplazamiento_t), intent(inout) :: grid
   deallocate (grid%desX, grid%desY, grid%desZ)
end subroutine
