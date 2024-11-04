module parser_tools_mod
    
#ifdef CompileWithSMBJSON
   use labels_mod
   use mesh_mod
   use cells_mod
   use json_module
   use json_kinds
   use NFDETypes

   use, intrinsic :: iso_fortran_env , only: error_unit

   implicit none

   integer, private, parameter :: J_ERROR_NUMBER = 1

   type :: json_value_ptr
      type(json_value), pointer :: p
   end type

   type :: aux_node_t
      type(terminal_node_t) :: node
      integer :: cId
      type(coordinate_t) :: relPos
   end type

   type :: cable_ptr_t
      type(cable_t), pointer :: p
   end type
   
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
   
   subroutine cellRegionsToCoords(res, cellRegions, cellType)
      type(coords), dimension(:), pointer :: res
      type(cell_region_t), dimension(:), intent(in) :: cellRegions
      integer, intent(in), optional :: cellType
      type(cell_interval_t), dimension(:), allocatable :: intervals
      type(coords), dimension(:), allocatable :: cs

      intervals = getIntervalsInCellRegions(cellRegions, cellType)
      cs = cellIntervalsToCoords(intervals)
      allocate(res(size(cs)))
      res = cs
   end

   function coordsToScaledCoords(cs) result(res)
      type(coords), intent(in), dimension(:) :: cs
      type(coords_scaled), dimension(:), allocatable :: res
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

   subroutine cellRegionsToScaledCoords(res, cellRegions)
      type(coords_scaled), dimension(:), pointer :: res
      type(cell_region_t), dimension(:), intent(in) :: cellRegions
      type(cell_interval_t), dimension(:), allocatable :: intervals
      type(coords), dimension(:), allocatable :: cs
      type(coords_scaled), dimension(:), allocatable :: scaledCoords
      
      intervals = getIntervalsInCellRegions(cellRegions, CELL_TYPE_LINEL)
      cs = cellIntervalsToCoords(intervals)
      scaledCoords = coordsToScaledCoords(cs)
      allocate(res(size(scaledCoords)))
      res = scaledCoords
   end

   function cellIntervalsToCoords(ivls) result(res)
      type(coords), dimension(:), pointer :: res
      type(cell_interval_t), dimension(:), intent(in) :: ivls
      integer :: i

      allocate(res(size(ivls)))
      do i = 1, size(ivls)
         res(i)%Or = ivls(i)%getOrientation()
         call convertInterval(res(i)%Xi, res(i)%Xe, ivls(i), DIR_X)
         call convertInterval(res(i)%Yi, res(i)%Ye, ivls(i), DIR_Y)
         call convertInterval(res(i)%Zi, res(i)%Ze, ivls(i), DIR_Z)
         res(i)%tag = ''
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
            xi = b + 1
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
      real, dimension(:), intent(in) :: vector
      real, dimension(:, :), allocatable :: res
      integer :: i, n
      n = size(vector, 1)
      allocate(res(n,n), source = 0.0)
      do i = 1, n
         res(i,i) = vector(i)
      end do
   end function

   function scalarToMatrix(scalar) result(res)
      real, intent(in) :: scalar
      real, dimension(:, :), allocatable :: res
      allocate(res(1,1), source = 0.0)
      res(1,1) = scalar
   end function
#endif
end module
