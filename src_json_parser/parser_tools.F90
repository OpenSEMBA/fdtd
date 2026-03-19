module parser_tools_m

#ifdef CompileWithSMBJSON
   use mesh_m
   use cells_m
   use json_module
   use json_kinds
   use NFDETypes_m

   use, intrinsic :: iso_fortran_env , only: error_unit

   implicit none

   integer, private, parameter :: J_ERROR_NUMBER = 1

   type :: json_value_ptr_t
      type(json_value), pointer :: p
   end type

#ifdef CompileWithMTLN
   type :: aux_node_t
      type(terminal_node_t) :: node
      integer :: cId
      type(coordinate_t) :: relPos
   end type
#endif

contains

   function getIntervalsInCellRegions(cellRegions, cellType) result (intervals)
      type(cell_region_t), dimension(:), intent(in) :: cellRegions
      integer, intent(in), optional :: cellType
      type(cell_interval_t), dimension(:), allocatable :: intervals, intervalsInRegion
      integer :: numberOfIntervals, copiedIntervals
      integer :: i, j

      numberOfIntervals = 0
      do i = 1, size(cellRegions)
         if (present(cellType)) then
            numberOfIntervals = numberOfIntervals + count(cellRegions(i)%intervals%getType() == cellType)
         else
            numberOfIntervals = numberOfIntervals + size(cellRegions(i)%intervals)
         end if
      end do

      allocate(intervals(numberOfIntervals))
      copiedIntervals = 0
      do i = 1, size(cellRegions)
         if (present(cellType)) then
            intervalsInRegion = cellRegions(i)%getIntervalsOfType(cellType)
         else
            intervalsInRegion = cellRegions(i)%intervals
         end if
         do j = 1, size(intervalsInRegion)
            copiedIntervals = copiedIntervals + 1
            intervals(copiedIntervals) = intervalsInRegion(j)
         end do
      end do

   end function

   function cellRegionToCoords(cellRegion, cellType, tag) result(res)
      type(cell_region_t), intent(in) :: cellRegion
      integer, intent(in), optional :: cellType
      character(len=BUFSIZE), optional, intent(in) :: tag
      type(coords_t), dimension(:), allocatable :: res

      type(cell_interval_t), dimension(:), allocatable :: intervals
      type(coords_t), dimension(:), allocatable :: cs

      intervals = getIntervalsInCellRegions([cellRegion], cellType)
      if (present(tag)) then
         cs = cellIntervalsToCoords(intervals, tag)
      else
         cs = cellIntervalsToCoords(intervals)
      end if
      res = cs
   end

   function coordsToScaledCoords(cs) result(res)
      type(coords_t), intent(in), dimension(:) :: cs
      type(coords_scaled_t), dimension(:), allocatable :: res
      integer :: i

      allocate(res(size(cs)))
      res(:)%Xi = cs(:)%Xi
      res(:)%Xe = cs(:)%Xe
      res(:)%Yi = cs(:)%Yi
      res(:)%Ye = cs(:)%Ye
      res(:)%Zi = cs(:)%Zi
      res(:)%Ze = cs(:)%Ze
      res(:)%Or = cs(:)%Or
      res(:)%tag = cs(:)%tag

      res(:)%xc = 0.0
      res(:)%yc = 0.0
      res(:)%zc = 0.0
      do i = 1, size(cs)
         select case (cs(i)%Or)
          case (iEx)
            res(i)%xc = 1.0
          case (-iEx)
            res(i)%xc = -1.0
          case (iEy)
            res(i)%yc = 1.0
          case (-iEy)
            res(i)%yc = -1.0
          case (iEz)
            res(i)%zc = 1.0
          case (-iEz)
            res(i)%zc = -1.0
         end select
      end do
   end

   subroutine cellRegionsToScaledCoords(res, cellRegions, tag)
      type(coords_scaled_t), dimension(:), pointer :: res
      type(cell_region_t), dimension(:), intent(in) :: cellRegions
      type(cell_interval_t), dimension(:), allocatable :: intervals
      type(coords_t), dimension(:), allocatable :: cs
      type(coords_scaled_t), dimension(:), allocatable :: scaledCoords
      character(len=BUFSIZE), optional, intent(in) :: tag

      intervals = getIntervalsInCellRegions(cellRegions, CELL_TYPE_LINEL)
      if (present(tag)) then
         cs = cellIntervalsToCoords(intervals, tag)
      else
         cs = cellIntervalsToCoords(intervals)
      end if
      scaledCoords = coordsToScaledCoords(cs)
      allocate(res(size(scaledCoords)))
      res = scaledCoords
   end

   function cellIntervalsToCoords(ivls, tag) result(res)
      type(coords_t), dimension(:), pointer :: res
      type(cell_interval_t), dimension(:), intent(in) :: ivls
      integer :: i
      character(len=BUFSIZE), optional, intent(in) :: tag

      allocate(res(size(ivls)))
      do i = 1, size(ivls)
         res(i)%Or = ivls(i)%getOrientation()
         call convertInterval(res(i)%Xi, res(i)%Xe, ivls(i), DIR_X)
         call convertInterval(res(i)%Yi, res(i)%Ye, ivls(i), DIR_Y)
         call convertInterval(res(i)%Zi, res(i)%Ze, ivls(i), DIR_Z)
         if (present(tag)) then
            res(i)%tag = tag
         else
            res(i)%tag = ''
         end if
      end do
   contains
      subroutine convertInterval(xi, xe, interval, dir)
         integer, intent(out) :: xi, xe
         type(cell_interval_t), intent(in) :: interval
         integer, intent(in) :: dir
         integer :: a, b
         a = interval%ini%cell(dir)
         b = interval%end%cell(dir)
         if (a < b) then
            xi = a
            xe = b - 1
         else if (a == b) then
            xi = a
            xe = b
         else
            xi = b
            xe = a - 1
         end if

      end subroutine
   end function

   function getPixelFromElementId(mesh, id) result(res)
      type(pixel_t) :: res
      type(mesh_t), intent(in) :: mesh
      integer, intent(in) :: id

      type(node_t) :: node
      logical :: nodeFound
      logical :: cellRegionFound
      integer :: i

      node = mesh%getNode(id, nodeFound)
      if (nodeFound) then
         res = mesh%nodeToPixel(node)
      else
         stop "Error converting pixel. Node not found."
      end if

   end function


   function vectorToDiagonalMatrix(vector) result(res)
      real(kind=RKIND), dimension(:), intent(in) :: vector
      real, dimension(:, :), allocatable :: res
      integer :: i, n
      n = size(vector, 1)
      allocate(res(n,n), source = 0.0)
      do i = 1, n
         res(i,i) = real(vector(i))
      end do
   end function

   function scalarToMatrix(scalar) result(res)
      real(kind=RKIND), intent(in) :: scalar
      real, dimension(:, :), allocatable :: res
      allocate(res(1,1), source = 0.0)
      res(1,1) = real(scalar)
   end function

   subroutine splitLineIntoWords(line, words)
      character(len=*), intent(in) :: line
      character(len=BUFSIZE), allocatable, dimension(:), intent(out) :: words

      integer :: lenstr, i, start, nwords, wlen

      lenstr = len_trim(line)
      nwords = 0
      i = 1
      ! First pass: count words
      do while (i <= lenstr)
         do while (i <= lenstr .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
            i = i + 1
         end do
         if (i > lenstr) exit
         nwords = nwords + 1
         do while (i <= lenstr .and. (line(i:i) /= ' ' .and. line(i:i) /= achar(9)))
            i = i + 1
         end do
      end do

      if (nwords == 0) then
         allocate(words(0))
         return
      end if

      allocate(words(nwords))
      i = 1
      start = 1
      nwords = 0
      do while (i <= lenstr)
         do while (i <= lenstr .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
            i = i + 1
         end do
         if (i > lenstr) exit
         start = i
         do while (i <= lenstr .and. (line(i:i) /= ' ' .and. line(i:i) /= achar(9)))
            i = i + 1
         end do
         wlen = i - start
         nwords = nwords + 1
         words(nwords) = line(start:start+wlen-1)
      end do
   end subroutine

   pure function to_upper(str) result(res)
      character(len=*), intent(in) :: str
      character(len=len(str)) :: res
      integer :: i
      res = str
      do i = 1, len(str)
         if (res(i:i) >= 'a' .and. res(i:i) <= 'z') then
            res(i:i) = achar(iachar(res(i:i)) - 32)
         end if
      end do
   end function

#endif
end module
