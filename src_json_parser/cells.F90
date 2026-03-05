module cells_mod

   integer, parameter :: DIR_NULL = 0
   integer, parameter :: DIR_X = 1
   integer, parameter :: DIR_Y = 2
   integer, parameter :: DIR_Z = 3

   integer, parameter :: NO_TAG = -1

   integer, parameter :: CELL_TYPE_PIXEL = 0
   integer, parameter :: CELL_TYPE_LINEL = 1
   integer, parameter :: CELL_TYPE_SURFEL = 2
   integer, parameter :: CELL_TYPE_VOXEL = 3

   ! --- Cells
   type :: cell_t
      integer, dimension(3) :: cell
   end type

   type, extends(cell_t) :: pixel_t
      integer :: tag = NO_TAG
   end type

   type, extends(cell_t) :: linel_t
      integer :: orientation = DIR_NULL                ! DIR_X, DIR_Y, DIR_Z
      integer :: tag = NO_TAG
   end type

   type, extends(cell_t) :: surfel_t
      integer :: orientation = DIR_NULL                ! DIR_X, DIR_Y, DIR_Z
      integer :: tag = NO_TAG
   end type

   type, extends(cell_t) :: voxel_t
      character(len=:), allocatable :: tag
   end type

   type :: cell_interval_t
      type(cell_t) :: ini, end
   contains
      procedure :: getType => cell_interval_getType
      procedure :: getOrientation => cell_interval_getOrientation
      procedure :: getSize => cell_interval_getSize
      procedure, private :: varyingDirections => cell_interval_varyingDirections
   end type

   type :: cell_region_t
      ! Cell regions are defined by semi-open intervals [ini, end).
      ! For linels and surfels, varying directions define orientation.
      type(cell_interval_t), dimension(:), allocatable :: intervals
   contains
      procedure :: toPixels => cell_region_toPixels
      procedure :: getIntervalsOfType => cell_region_getIntervalsOfType
   end type

   interface operator(==)
      module procedure linel_eq
      module procedure pixel_eq
   end interface

contains
   function cell_region_toPixels(this) result(res)
      class(cell_region_t) :: this
      type(pixel_t), dimension(:), allocatable :: res
      type(cell_interval_t), dimension(:), allocatable :: pixelIntervals
      integer :: i
      pixelIntervals = this%getIntervalsOfType(CELL_TYPE_PIXEL)
      allocate(res(size(pixelIntervals)))
      do i = 1, size(res)
         res(i)%cell = pixelIntervals(i)%ini%cell
      end do
   end function

   logical function pixel_eq(a, b)
      type(pixel_t), intent(in) :: a, b
      pixel_eq = &
         all(a%cell == b%cell) .and. &
         a%tag == b%tag
   end function

   logical function linel_eq(a, b)
      type(linel_t), intent(in) :: a, b
      linel_eq = &
         all(a%cell == b%cell) .and. &
         a%tag == b%tag .and. &
         a%orientation == b%orientation
   end function

   elemental function cell_interval_getType(this) result(res)
      class(cell_interval_t), intent(in) :: this
      integer :: res
      res = this%varyingDirections()
   end function

   elemental function cell_interval_getOrientation(this) result(res)
      class(cell_interval_t), intent(in) :: this
      integer :: res
      integer :: i
      select case (this%getType())
       case (CELL_TYPE_LINEL)
         block
            integer :: diff
            do i = DIR_X, DIR_Z
               diff = this%end%cell(i) - this%ini%cell(i)
               if (diff > 0) then
                  res = i
                  return
               else if (diff < 0) then
                  res = -i
                  return
               end if
            end do
         end block
       case (CELL_TYPE_SURFEL)
         block
            integer, dimension(3) :: diff
            diff = this%end%cell - this%ini%cell
            do i = DIR_X, DIR_Z
               if (diff(i) == 0) res = i
            end do
            if ( diff(mod(res,3)+1) < 0 .and. diff(mod(res+1,3)+1) < 0) &
               res = - res
         end block
       case default
         res = DIR_NULL
      end select
   end function

   elemental function cell_interval_getSize(this) result(res)
      class(cell_interval_t), intent(in) :: this
      integer :: res
      integer :: i
      integer, dimension(3) :: diff
      res = 1
      diff = abs(this%end%cell - this%ini%cell)
      do i = DIR_X, DIR_Z
         if (diff(i) /= 0) res = res * diff(i)
      end do
      if (all(diff == 0)) res = 0
   end function

   elemental function cell_interval_varyingDirections(this) result(res)
      class(cell_interval_t), intent(in) :: this
      integer :: res
      integer :: i
      res = 0
      do i = DIR_X, DIR_Z
         if ((this%end%cell(i) - this%ini%cell(i)) /= 0) res = res + 1
      end do
   end function

   function cell_region_getIntervalsOfType(this, cellType) result(res)
      class(cell_region_t), intent(in) :: this
      integer, intent(in) :: cellType
      type(cell_interval_t), dimension(:), allocatable :: res
      integer :: i, j

      allocate(res( count(this%intervals%getType() == cellType) ))
      j = 1
      do i = 1, size(this%intervals)
         if (this%intervals(i)%getType() == cellType) then
            res(j) = this%intervals(i)
            j = j + 1
         end if
      end do
   end function

end module
