module vtkAPI_m
   implicit none
   private
   public :: vtk_data_array, vtk_grid, vtk_structured_grid, vtk_unstructured_grid

   !==========================
   ! Data array type
   !==========================
   type :: vtk_data_array
      character(len=:), allocatable :: name
      character(len=:), allocatable :: type
      integer :: num_components
      real, allocatable :: data(:)
   end type vtk_data_array

   !==========================
   ! Base abstract class
   !==========================
   type, abstract :: vtk_grid
      real, allocatable :: points(:, :)
      type(vtk_data_array), allocatable :: scalars(:)
      type(vtk_data_array), allocatable :: vectors(:)
      type(vtk_data_array), allocatable :: cell_scalars(:)
      type(vtk_data_array), allocatable :: cell_vectors(:)
   contains
      procedure(add_points_generic), deferred :: add_points
      procedure(add_scalar_generic), deferred :: add_scalar
      procedure(add_vector_generic), deferred :: add_vector
      procedure(add_cell_scalar_generic), deferred :: add_cell_scalar
      procedure(add_cell_vector_generic), deferred :: add_cell_vector
      procedure(write_file_generic), deferred :: write_file
   end type vtk_grid

   !==========================
   ! Structured Grid (VTS)
   !==========================
   type, extends(vtk_grid) :: vtk_structured_grid
      integer :: nx = 0, ny = 0, nz = 0
   contains
      procedure :: add_points => add_points_structured
      procedure :: add_scalar => add_scalar_structured
      procedure :: add_vector => add_vector_structured
      procedure :: add_cell_scalar => add_cell_scalar_structured
      procedure :: add_cell_vector => add_cell_vector_structured
      procedure :: write_file => write_vts_file
   end type vtk_structured_grid

   !==========================
   ! Unstructured Grid (VTU)
   !==========================
   type, extends(vtk_grid) :: vtk_unstructured_grid
      integer :: num_points = 0
      integer :: num_cells = 0
      integer, allocatable :: connectivity(:)
      integer, allocatable :: offsets(:)
      integer, allocatable :: types(:)
   contains
      procedure :: add_points => add_points_unstructured
      procedure :: add_scalar => add_scalar_unstructured
      procedure :: add_vector => add_vector_unstructured
      procedure :: add_cell_scalar => add_cell_scalar_unstructured
      procedure :: add_cell_vector => add_cell_vector_unstructured
      procedure :: add_cell_connectivity
      procedure :: write_file => write_vtu_file
   end type vtk_unstructured_grid

   !==========================
   ! Generic deferred interfaces
   !==========================
   abstract interface
      subroutine add_points_generic(this, pts)
         import :: vtk_grid
         class(vtk_grid), intent(inout) :: this
         real, intent(in) :: pts(:, :)
      end subroutine add_points_generic

      subroutine add_scalar_generic(this, name, data)
         import :: vtk_grid
         class(vtk_grid), intent(inout) :: this
         character(len=*), intent(in) :: name
         real, intent(in) :: data(:)
      end subroutine add_scalar_generic

      subroutine add_vector_generic(this, name, data)
         import :: vtk_grid
         class(vtk_grid), intent(inout) :: this
         character(len=*), intent(in) :: name
         real, intent(in) :: data(:)
      end subroutine add_vector_generic

      subroutine add_cell_scalar_generic(this, name, data)
         import :: vtk_grid
         class(vtk_grid), intent(inout) :: this
         character(len=*), intent(in) :: name
         real, intent(in) :: data(:)
      end subroutine add_cell_scalar_generic

      subroutine add_cell_vector_generic(this, name, data)
         import :: vtk_grid
         class(vtk_grid), intent(inout) :: this
         character(len=*), intent(in) :: name
         real, intent(in) :: data(:)
      end subroutine add_cell_vector_generic

      subroutine write_file_generic(this, filename)
         import :: vtk_grid
         class(vtk_grid), intent(in) :: this
         character(len=*), intent(in) :: filename
      end subroutine write_file_generic
   end interface

contains
   !==========================
   !==== Structured Grid Methods ====
   !==========================

   subroutine add_points_structured(this, pts)
      class(vtk_structured_grid), intent(inout) :: this
      real, intent(in) :: pts(:, :)
      integer :: npts
      npts = this%nx*this%ny*this%nz
      if (size(pts, 1) /= 3) error stop 'add_points_structured: first dim must be 3'
      if (size(pts, 2) /= npts) error stop 'add_points_structured: wrong number of points'
      if (allocated(this%points)) deallocate (this%points)
      allocate (this%points(3, npts))
      this%points = pts
   end subroutine add_points_structured

   subroutine add_scalar_structured(this, name, data)
      class(vtk_structured_grid), intent(inout) :: this
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      call add_array_generic(this%scalars, name, data, 1)
   end subroutine add_scalar_structured

   subroutine add_vector_structured(this, name, data)
      class(vtk_structured_grid), intent(inout) :: this
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      call add_array_generic(this%vectors, name, data, 3)
   end subroutine add_vector_structured

   subroutine add_cell_scalar_structured(this, name, data)
      class(vtk_structured_grid), intent(inout) :: this
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      call add_array_generic(this%cell_scalars, name, data, 1)
   end subroutine add_cell_scalar_structured

   subroutine add_cell_vector_structured(this, name, data)
      class(vtk_structured_grid), intent(inout) :: this
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      call add_array_generic(this%cell_vectors, name, data, 3)
   end subroutine add_cell_vector_structured

   !==========================
   ! Write VTS
   !==========================
   subroutine write_vts_file(this, filename)
      class(vtk_structured_grid), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer :: iunit
      open (newunit=iunit, file=filename, status='replace', action='write', form='formatted')
      write (iunit, *) '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
      write (iunit, *) '  <StructuredGrid WholeExtent="0 ', this%nx - 1, ' 0 ', this%ny - 1, ' 0 ', this%nz - 1, '">'
      write (iunit, *) '    <Piece Extent="0 ', this%nx - 1, ' 0 ', this%ny - 1, ' 0 ', this%nz - 1, '">'
      call write_pointdata(iunit, this%scalars, this%vectors)
      call write_celldata(iunit, this%cell_scalars, this%cell_vectors)
      call write_points(iunit, this%points)
      write (iunit, *) '    </Piece>'
      write (iunit, *) '  </StructuredGrid>'
      write (iunit, *) '</VTKFile>'
      close (iunit)
   end subroutine write_vts_file

   !==========================
   !==== Unstructured Grid Methods ====
   !==========================
   subroutine add_points_unstructured(this, pts)
      real, intent(in) :: pts(:, :)
      class(vtk_unstructured_grid), intent(inout) :: this
      if (size(pts, 1) /= 3) error stop 'add_points_unstructured: first dim must be 3'
      this%num_points = size(pts, 2)
      if (allocated(this%points)) deallocate (this%points)
      allocate (this%points(3, this%num_points))
      this%points = pts
   end subroutine add_points_unstructured

   subroutine add_scalar_unstructured(this, name, data)
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      class(vtk_unstructured_grid), intent(inout) :: this
      call add_array_generic(this%scalars, name, data, 1)
   end subroutine add_scalar_unstructured

   subroutine add_vector_unstructured(this, name, data)
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      class(vtk_unstructured_grid), intent(inout) :: this
      call add_array_generic(this%vectors, name, data, 3)
   end subroutine add_vector_unstructured

   subroutine add_cell_scalar_unstructured(this, name, data)
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      class(vtk_unstructured_grid), intent(inout) :: this
      call add_array_generic(this%cell_scalars, name, data, 1)
   end subroutine add_cell_scalar_unstructured

   subroutine add_cell_vector_unstructured(this, name, data)
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      class(vtk_unstructured_grid), intent(inout) :: this
      call add_array_generic(this%cell_vectors, name, data, 3)
   end subroutine add_cell_vector_unstructured

   subroutine add_cell_connectivity(this, conn, offsets, types)
      class(vtk_unstructured_grid), intent(inout) :: this
      integer, intent(in) :: conn(:), offsets(:), types(:)
      this%num_cells = size(offsets)
      if (allocated(this%connectivity)) deallocate (this%connectivity)
      if (allocated(this%offsets)) deallocate (this%offsets)
      if (allocated(this%types)) deallocate (this%types)
      allocate (this%connectivity(size(conn)))
      allocate (this%offsets(size(offsets)))
      allocate (this%types(size(types)))
      this%connectivity = conn
      this%offsets = offsets
      this%types = types
   end subroutine add_cell_connectivity

   !==========================
   ! Write VTU
   !==========================
   subroutine write_vtu_file(this, filename)
      character(len=*), intent(in) :: filename
      class(vtk_unstructured_grid), intent(in) :: this
      integer :: iunit
      open (newunit=iunit, file=filename, status='replace', action='write', form='formatted')
      write (iunit, *) '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
      write (iunit, *) '  <UnstructuredGrid>'
      write (iunit, *) '    <Piece NumberOfPoints="', this%num_points, '" NumberOfCells="', this%num_cells, '">'
      call write_pointdata(iunit, this%scalars, this%vectors)
      call write_celldata(iunit, this%cell_scalars, this%cell_vectors)
      call write_cells(iunit, this%connectivity, this%offsets, this%types)
      call write_points(iunit, this%points)
      write (iunit, *) '    </Piece>'
      write (iunit, *) '  </UnstructuredGrid>'
      write (iunit, *) '</VTKFile>'
      close (iunit)
   end subroutine write_vtu_file

!==========================
!==== Shared helper routines ====
!==========================
   subroutine add_array_generic(array_list, name, data, ncomp)
      type(vtk_data_array), allocatable, intent(inout) :: array_list(:)
      character(len=*), intent(in) :: name
      real, intent(in) :: data(:)
      integer, intent(in) :: ncomp
      type(vtk_data_array), allocatable :: tmp(:)
      integer :: n
      if (.not. allocated(array_list)) then
         allocate (array_list(1))
         n = 1
      else
         n = size(array_list) + 1
         allocate (tmp(n))
         tmp(1:n - 1) = array_list
         call move_alloc(tmp, array_list)
      end if
      array_list(n)%name = name
      array_list(n)%type = 'Float32'
      array_list(n)%num_components = ncomp
      array_list(n)%data = data
   end subroutine add_array_generic

   subroutine write_pointdata(iunit, scalars, vectors)
      integer, intent(in) :: iunit
      type(vtk_data_array), allocatable, intent(in) :: scalars(:)
      type(vtk_data_array), allocatable, intent(in) :: vectors(:)
      integer :: i
      character(len=1024) :: tag
      tag = '      <PointData'
      if (allocated(scalars)) then
         if (size(scalars) > 0) tag = trim(tag)//' Scalars="'//trim(scalars(1)%name)//'"'
      end if
      if (allocated(vectors)) then
         if (size(vectors) > 0) tag = trim(tag)//' Vectors="'//trim(vectors(1)%name)//'"'
      end if
      tag = trim(tag)//'>'
      write (iunit, '(A)') trim(tag)
      if (allocated(scalars)) then
         do i = 1, size(scalars)
            write (iunit, '(A)') '        <DataArray type="Float32" Name="'//trim(scalars(i)%name)//'" format="ascii">'
            write (iunit, '(1000(F12.6,1X))') scalars(i)%data
            write (iunit, '(A)') '        </DataArray>'
         end do
      end if
      if (allocated(vectors)) then
         do i = 1, size(vectors)
  write (iunit, '(A)') '        <DataArray type="Float32" Name="'//trim(vectors(i)%name)//'" NumberOfComponents="3" format="ascii">'
            write (iunit, '(1000(F12.6,1X))') vectors(i)%data
            write (iunit, '(A)') '        </DataArray>'
         end do
      end if
      write (iunit, *) '      </PointData>'
   end subroutine write_pointdata

   subroutine write_celldata(iunit, scalars, vectors)
      integer, intent(in) :: iunit
      type(vtk_data_array), allocatable, intent(in) :: scalars(:)
      type(vtk_data_array), allocatable, intent(in) :: vectors(:)
      integer :: i
      if (.not. allocated(scalars) .and. .not. allocated(vectors)) then
         write (iunit, *) '      <CellData>'
         write (iunit, *) '      </CellData>'
         return
      end if
      write (iunit, *) '      <CellData>'
      if (allocated(scalars)) then
         do i = 1, size(scalars)
            write (iunit, '(A)') '        <DataArray type="Float32" Name="'//trim(scalars(i)%name)//'" format="ascii">'
            write (iunit, '(1000(F12.6,1X))') scalars(i)%data
            write (iunit, '(A)') '        </DataArray>'
         end do
      end if
      if (allocated(vectors)) then
         do i = 1, size(vectors)
  write (iunit, '(A)') '        <DataArray type="Float32" Name="'//trim(vectors(i)%name)//'" NumberOfComponents="3" format="ascii">'
            write (iunit, '(1000(F12.6,1X))') vectors(i)%data
            write (iunit, '(A)') '        </DataArray>'
         end do
      end if
      write (iunit, *) '      </CellData>'
   end subroutine write_celldata

   subroutine write_points(iunit, pts)
      integer, intent(in) :: iunit
      real, intent(in) :: pts(:, :)
      write (iunit, *) '      <Points>'
      write (iunit, *) '        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">'
      write (iunit, '(1000(F12.6,1X))') pts
      write (iunit, *) '        </DataArray>'
      write (iunit, *) '      </Points>'
   end subroutine write_points

   subroutine write_cells(iunit, conn, offsets, types)
      integer, intent(in) :: iunit
      integer, intent(in) :: conn(:), offsets(:), types(:)
      write (iunit, *) '      <Cells>'
      write (iunit, *) '        <DataArray type="Int32" Name="connectivity" format="ascii">'
      write (iunit, '(1000(I8,1X))') conn
      write (iunit, *) '        </DataArray>'
      write (iunit, *) '        <DataArray type="Int32" Name="offsets" format="ascii">'
      write (iunit, '(1000(I8,1X))') offsets
      write (iunit, *) '        </DataArray>'
      write (iunit, *) '        <DataArray type="UInt8" Name="types" format="ascii">'
      write (iunit, '(1000(I3,1X))') types
      write (iunit, *) '        </DataArray>'
      write (iunit, *) '      </Cells>'
   end subroutine write_cells

end module vtkAPI_m
