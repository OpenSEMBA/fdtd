module mesh_mod

#ifdef CompileWithSMBJSON   
   use, intrinsic :: iso_fortran_env , only: error_unit
   
   use fhash, only: fhash_tbl_t, key=>fhash_key
   use cells_mod
   use geometry_mod, only: triangle_t
   integer, private, parameter  ::  MAX_LINE = 256
   integer, parameter :: REGION_TYPE_VOLUME = 3
   integer, parameter :: REGION_TYPE_SURFACE = 2

   type :: element_t
      integer, dimension(:), allocatable :: coordIds
   end type

   type, public, extends(element_t) :: node_t
      ! coordIds must be size 1.
   end type

   type, public, extends(element_t) :: polyline_t
      ! coordIds must be size >1.
   end type

   type, public :: coordinate_t
      real, dimension(3) :: position
   contains
      private
      procedure :: coordinate_diff
      generic, public :: operator(-) => coordinate_diff
      procedure :: coordinate_eq
      generic, public :: operator(==) => coordinate_eq
   end type

   type, public :: mesh_t
      private
      type(fhash_tbl_t) :: coordinates ! Map of CoordinateIds to relative coordinates.
      type(fhash_tbl_t) :: elements    ! Map of ElementIds to elements/cellsRegions.    
   contains
      procedure :: addCoordinate => mesh_addCoordinate
      procedure :: getCoordinate => mesh_getCoordinate
      procedure :: checkCoordinateId => mesh_checkCoordinateId
      procedure :: checkElementId => mesh_checkElementId

      procedure :: addElement => mesh_addElement
      procedure :: addCellRegion  => mesh_addCellRegion
      procedure :: addConformalRegion  => mesh_addConformalRegion
      
      procedure :: getNode => mesh_getNode
      procedure :: getPolyline => mesh_getPolyline
      procedure :: getCellRegion  => mesh_getCellRegion
      procedure :: getCellRegions => mesh_getCellRegions
      procedure :: getConformalRegion => mesh_getConformalRegion
      procedure :: getConformalRegions => mesh_getConformalRegions

      procedure :: countPolylineSegments => mesh_countPolylineSegments
      procedure :: arePolylineSegmentsStructured => mesh_arePolylineSegmentsStructured
      procedure :: polylineToLinels => mesh_polylineToLinels
      procedure :: nodeToPixel => mesh_nodeToPixel

      procedure :: printHashInfo => mesh_printHashInfo
      procedure :: allocateCoordinates => mesh_allocateCoordinates
      procedure :: allocateElements => mesh_allocateElements
   end type

   type, public :: conformal_region_t
      type(triangle_t), dimension(:), allocatable :: triangles
      type(cell_interval_t), dimension(:), allocatable :: intervals
      integer :: type
    end type


contains
   ! __________________________________________________________________
   ! Mesh procedures
   subroutine mesh_allocateCoordinates(this, buck)
      class(mesh_t) :: this
      integer :: buck
      call this%coordinates%allocate(buck)
   end subroutine
   
   subroutine mesh_allocateElements(this, buck)
      class(mesh_t) :: this
      integer :: buck
      call this%elements%allocate(buck)
   end subroutine


   subroutine mesh_printHashInfo(this)
      class(mesh_t) :: this
      integer :: num_buckets, num_items, num_collisions, max_depth
      
      write(*,*) 'Coordinates hash info:'
      call this%coordinates%stats(num_buckets,num_items,num_collisions,max_depth)
      write(*,'(A,T40,I0)') '  Number of buckets allocated: ',num_buckets
      write(*,'(A,T40,I0)') '  Number of key-value pairs stored: ',num_items
      write(*,'(A,T40,I0)') '  Total number of hash-collisions: ',num_collisions
      write(*,'(A,T40,I0)') '  The worst case bucket depth is ',max_depth
      
      write(*,*) 'Elements hash info:'
      call this%elements%stats(num_buckets,num_items,num_collisions,max_depth)
      write(*,'(A,T40,I0)') '  Number of buckets allocated: ',num_buckets
      write(*,'(A,T40,I0)') '  Number of key-value pairs stored: ',num_items
      write(*,'(A,T40,I0)') '  Total number of hash-collisions: ',num_collisions
      write(*,'(A,T40,I0)') '  The worst case bucket depth is ',max_depth

   end subroutine

   function mesh_checkCoordinateId(this, id) result(stat)
      class(mesh_t) :: this
      integer, intent(in) :: id
      integer :: stat
      call this%coordinates%check_key(key(id), stat)
   end function

   
   function mesh_checkElementId(this, id) result(stat)
      class(mesh_t) :: this
      integer, intent(in) :: id
      integer :: stat
      call this%elements%check_key(key(id), stat)
   end function

   subroutine mesh_addCoordinate(this, id, coordinate)
      class(mesh_t) :: this
      integer, intent(in) :: id
      type(coordinate_t), intent(in) :: coordinate
      call this%coordinates%set(key(id), value=coordinate)
   end subroutine

   subroutine mesh_addElement(this, id, e)
      class(mesh_t) :: this
      integer, intent(in) :: id
      class(element_t), intent(in) :: e
      call this%elements%set(key(id), value=e)
   end subroutine

   subroutine mesh_addCellRegion(this, id, e)
      class(mesh_t) :: this
      integer, intent(in) :: id
      class(cell_region_t), intent(in) :: e
      call this%elements%set(key(id), value=e)
   end subroutine

   subroutine mesh_addConformalRegion(this, id, e)
      class(mesh_t) :: this
      integer, intent(in) :: id
      class(conformal_region_t), intent(in) :: e
      call this%elements%set(key(id), value=e)
   end subroutine

   function mesh_getCoordinate(this, id, found) result(res)
      class(mesh_t) :: this
      type(coordinate_t) :: res
      integer, intent(in) :: id
      integer :: stat
      logical, intent(out), optional :: found
      class(*), allocatable :: d

      if (present(found)) found = .false.

      call this%coordinates%get_raw(key(id), d, stat)
      if (stat /= 0) return

      select type(d)
         type is (coordinate_t)
         res = d
         if (present(found)) found = .true.
      end select

   end function

   function mesh_getNode(this, id, found) result(res)
      class(mesh_t) :: this
      type(node_t) :: res
      integer, intent(in) :: id
      logical, optional, intent(out) :: found
      integer :: status
      class(*), allocatable :: d

      if (present(found)) found = .false.
      call this%elements%get_raw(key(id), d, status)
      if (status /= 0) return

      select type(d)
         type is (node_t)
         res = d
         if (present(found)) found = .true.
      end select

   end function

   function mesh_getPolyline(this, id, found) result(res)
      class(mesh_t) :: this
      type(polyline_t) :: res
      integer, intent(in) :: id
      integer :: stat
      logical, intent(out), optional :: found
      class(*), allocatable :: d

      if (present(found)) found = .false.
      call this%elements%get_raw(key(id), d, stat)
      if (stat /= 0) return

      select type(d)
         type is (polyline_t)
         res = d
         if (present(found)) found = .true.
      end select

   end function

   function mesh_getCellRegion(this, id, found) result (res)
      class(mesh_t) :: this
      type(cell_region_t) :: res
      integer, intent(in) :: id
      integer :: stat
      logical, intent(out), optional :: found
      class(*), allocatable :: d

      if (present(found)) found = .false.
      call this%elements%get_raw(key(id), d, stat)
      if (stat /= 0) return

      select type(d)
         type is (cell_region_t)
            res = d
            if (present(found)) found = .true.
         type is (conformal_region_t)
            if (size(d%intervals) == 0) return
            res%intervals = d%intervals
            if (present(found)) found = .true.
      end select

   end function

   function mesh_getCellRegions(this, ids) result (res)
      class(mesh_t) :: this
      type(cell_region_t), dimension(:), allocatable :: res
      integer, dimension(:), intent(in) :: ids
      type(cell_region_t) :: cR
      logical :: found
      integer :: i, j
      integer :: numberOfCellRegions

      ! Precounts
      numberOfCellRegions = 0
      do i = 1, size(ids)
         cR = this%getCellRegion(ids(i), found)
         if (found) then
               numberOfCellRegions = numberOfCellRegions + 1
         end if
      end do     
      
      allocate(res(numberOfCellRegions))
      j = 1
      do i = 1, size(ids)
         cR = this%getCellRegion(ids(i), found)
         if (found) then
               res(j) = cR
               j = j + 1
         end if
      end do

   end function

   function mesh_getConformalRegion(this, id, found) result (res)
      class(mesh_t) :: this
      type(conformal_region_t) :: res
      integer, intent(in) :: id
      integer :: stat
      logical, intent(out), optional :: found
      class(*), allocatable :: d

      if (present(found)) found = .false.
      call this%elements%get_raw(key(id), d, stat)
      if (stat /= 0) return

      select type(d)
         type is (cell_region_t)
            return
         type is (conformal_region_t)
            res = d
            if (present(found)) found = .true.
      end select

   end function

   function mesh_getConformalRegions(this, ids) result (res)
      class(mesh_t) :: this
      type(conformal_region_t), dimension(:), allocatable :: res
      integer, dimension(:), intent(in) :: ids
      type(conformal_region_t) :: cR
      logical :: found
      integer :: i, j
      integer :: numberOfConformalRegions

      ! Precounts
      numberOfConformalRegions = 0
      do i = 1, size(ids)
         cR = this%getConformalRegion(ids(i), found)
         if (found) then
               numberOfConformalRegions = numberOfConformalRegions + 1
         end if
      end do     
      
      allocate(res(numberOfConformalRegions))
      j = 1
      do i = 1, size(ids)
         cR = this%getConformalRegion(ids(i), found)
         if (found) then
               res(j) = cR
               j = j + 1
         end if
      end do

   end function

   function mesh_countPolylineSegments(this, pl) result(res)
      integer :: res
      class(mesh_t) :: this
      type(polyline_t) :: pl
      type(coordinate_t) :: iC, eC
      type(cell_interval_t) :: interval
      
      res = 0
      do i = 1, size(pl%coordIds)-1
         iC = this%getCoordinate(pl%coordIds(i))
         eC = this%getCoordinate(pl%coordIds(i+1))
         interval%ini%cell = int(iC%position)
         interval%end%cell = int(eC%position)
         res = res + interval%getSize()
      end do
      
   end

   function mesh_arePolylineSegmentsStructured(this, pl) result(res)
      logical :: res
      class(mesh_t) :: this
      type(polyline_t) :: pl
      type(coordinate_t) :: iC, eC
      integer :: i, d
      integer :: numberOfVaryingDirections

      do i = 1, size(pl%coordIds)-1
         iC = this%getCoordinate(pl%coordIds(i))
         eC = this%getCoordinate(pl%coordIds(i+1))
         if (any(floor(iC%position) /= iC%position) .or. any(floor(eC%position) /= eC%position)) then
            res = .false.
            return
         end if 

         numberOfVaryingDirections = 0
         do d = DIR_X, DIR_Z
            if (iC%position(d) /= eC%position(d)) then
               numberOfVaryingDirections = numberOfVaryingDirections + 1
            end if
         end do
         if (numberOfVaryingDirections > 1) then
            res = .false.
            return
         end if
      end do

      res = .true.
   end function

   function mesh_polylineToLinels(this, pl) result(res)
      type(linel_t), dimension(:), allocatable :: res
      class(mesh_t), intent(in) :: this
      type(polyline_t), intent(in) :: pl
      type(cell_interval_t) :: interval
      type(coordinate_t) :: iC, eC, mC
      integer :: i, j, lastSegment, nLinelsInSegment
      integer, dimension(3) :: segment

      if (.not. this%arePolylineSegmentsStructured(pl)) then
         allocate(res(0))
         return
      end if

      allocate(res(this%countPolylineSegments(pl)))
      if (size(res) == 0) return

      lastSegment = 1
      do i = 1, size(pl%coordIds)-1
         iC = this%getCoordinate(pl%coordIds(i))
         eC = this%getCoordinate(pl%coordIds(i+1))
         interval%ini%cell = int(iC%position)
         interval%end%cell = int(eC%position)
         if (any(iC%position /= eC%position)) then
            segment = (interval%end%cell - interval%ini%cell) / interval%getSize()
            
            res(lastSegment)%tag = pl%coordIds(i)
            do j = 1, interval%getSize()
               mC%position = iC%position + segment * (real(j-1) + 0.5)
               res(lastSegment)%cell = floor(mc%position)
               res(lastSegment)%orientation = interval%getOrientation()
               lastSegment = lastSegment + 1
            end do
         end if
      end do

      res(1)%tag             = pl%coordIds( 1 )
      res(lastSegment-1)%tag = pl%coordIds( size(pl%coordIds) )
      
   end function

   function mesh_nodeToPixel(this, node) result(res)
      type(pixel_t) :: res
      class(mesh_t), intent(in) :: this
      type(node_t), intent(in) :: node

      type(coordinate_t) :: c
      logical :: coordFound

      c = this%getCoordinate(node%coordIds(1), found=coordFound)
      if (.not. coordFound) then
         write(error_unit, *) "ERROR: converting node to pixel. Coordinate not found."
         return
      end if
      res%cell = c%position
      res%tag = node%coordIds(1)
   end function

   function coordinate_diff(a, b) result(res)
      class(coordinate_t), intent(in) :: a, b
      type(coordinate_t) :: res
      res%position = [(a%position(i) - b%position(i), i = 1, 3)]
   end function

   logical function coordinate_eq(a, b)
      class(coordinate_t), intent(in) :: a, b
      coordinate_eq = all(a%position == b%position)
   end function
#endif
end module
