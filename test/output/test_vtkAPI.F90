!==============================
! vtkAPI_m extended testsuite (class-based)
!==============================

!==============================
! Test 1: Structured grid basic allocation
!==============================
integer function test_vtkAPI_points_allocation() bind(C) result(error_cnt)
   use vtkAPI_m
   implicit none
   type(vtk_structured_grid), target :: grid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: points(:, :)
   integer :: nx, ny, nz

   error_cnt = 0
   nx = 2; ny = 2; nz = 2
   grid%nx = nx; grid%ny = ny; grid%nz = nz
   grid_base => grid

   allocate(points(3, nx*ny*nz))
   points = 0.0
   call grid_base%add_points(points)

   if (.not. allocated(grid%points)) error_cnt = error_cnt + 1
   if (size(grid%points,1) /= 3) error_cnt = error_cnt + 1
   if (size(grid%points,2) /= nx*ny*nz) error_cnt = error_cnt + 1
end function

!==============================
! Test 2: Point scalar assignment
!==============================
integer function test_vtkAPI_point_scalar() bind(C) result(error_cnt)
   use vtkAPI_m
   implicit none
   type(vtk_structured_grid), target :: grid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: scalars(:)
   integer :: nx=2, ny=2, nz=2, i

   error_cnt = 0
   grid%nx=nx; grid%ny=ny; grid%nz=nz
   grid_base => grid

   allocate(scalars(nx*ny*nz))
   do i=1,nx*ny*nz
      scalars(i) = real(i)
   end do
   call grid_base%add_scalar('Density', scalars)

   if (.not. allocated(grid%scalars)) error_cnt = error_cnt + 1
   if (grid%scalars(1)%data(1) /= 1.0) error_cnt = error_cnt + 1
end function

!==============================
! Test 3: Point vector assignment
!==============================
integer function test_vtkAPI_point_vector() bind(C) result(error_cnt)
   use vtkAPI_m
   implicit none
   type(vtk_structured_grid), target :: grid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: vec(:)
   integer :: nx=2, ny=2, nz=2, i

   error_cnt = 0
   grid%nx=nx; grid%ny=ny; grid%nz=nz
   grid_base => grid

   allocate(vec(3*nx*ny*nz))
   do i=1,nx*ny*nz
      vec(3*i-2) = real(i)
      vec(3*i-1) = real(i*10)
      vec(3*i)   = real(i*100)
   end do
   call grid_base%add_vector('Momentum', vec)

   if (.not. allocated(grid%vectors)) error_cnt = error_cnt + 1
   if (grid%vectors(1)%data(2) /= 10.0) error_cnt = error_cnt + 1
end function

!==============================
! Test 4: Cell scalar assignment
!==============================
integer function test_vtkAPI_cell_scalar() bind(C) result(error_cnt)
   use vtkAPI_m
   implicit none
   type(vtk_structured_grid), target :: grid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: cell_data(:)
   integer :: nx=2, ny=2, nz=2, n

   error_cnt = 0
   grid%nx=nx; grid%ny=ny; grid%nz=nz
   grid_base => grid

   allocate(cell_data((nx-1)*(ny-1)*(nz-1)))
   do n=1,(nx-1)*(ny-1)*(nz-1)
      cell_data(n) = real(n)
   end do
   call grid_base%add_cell_scalar('Pressure', cell_data)

   if (.not. allocated(grid%cell_scalars)) error_cnt = error_cnt + 1
   if (grid%cell_scalars(1)%data(1) /= 1.0) error_cnt = error_cnt + 1
end function

!==============================
! Test 5: Cell vector assignment
!==============================
integer function test_vtkAPI_cell_vector() bind(C) result(error_cnt)
   use vtkAPI_m
   implicit none
   type(vtk_structured_grid), target :: grid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: vec(:)
   integer :: nx=2, ny=2, nz=2, n

   error_cnt = 0
   grid%nx=nx; grid%ny=ny; grid%nz=nz
   grid_base => grid

   allocate(vec(3*(nx-1)*(ny-1)*(nz-1)))
   do n=1,(nx-1)*(ny-1)*(nz-1)
      vec(3*n-2) = real(n)
      vec(3*n-1) = real(n*10)
      vec(3*n)   = real(n*100)
   end do
   call grid_base%add_cell_vector('Flux', vec)

   if (.not. allocated(grid%cell_vectors)) error_cnt = error_cnt + 1
   if (grid%cell_vectors(1)%data(3) /= 100.0) error_cnt = error_cnt + 1
end function

!==============================
! Test 6: VTS file creation
!==============================
integer function test_vtkAPI_vts_file_creation() bind(C) result(error_cnt)
   use vtkAPI_m
   use directoryUtils_m
   implicit none
   type(vtk_structured_grid), target :: grid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: points(:, :), scalars(:)
   integer :: nx=2, ny=2, nz=2
   integer :: ierr
   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file

   error_cnt = 0
   grid%nx=nx; grid%ny=ny; grid%nz=nz
   grid_base => grid

   allocate(points(3, nx*ny*nz)); points = 0.0
   call grid_base%add_points(points)

   allocate(scalars(nx*ny*nz)); scalars = 1.0
   call grid_base%add_scalar('Density', scalars)

   file = join_path(folder, 'test.vts')
   call create_folder(folder, ierr)
   call grid_base%write_file(file)

   open(unit=10, file=file, status='old', action='read', iostat=ierr)
   if (ierr /= 0) then
      error_cnt = error_cnt + 1
   else
      close(10)
   end if
   call remove_folder(folder, ierr)
end function

!==============================
! Test 7: VTU file creation with cells and point data
!==============================
integer function test_vtkAPI_vtu_file_creation() bind(C) result(error_cnt)
   use vtkAPI_m
   use directoryUtils_m
   implicit none
   type(vtk_unstructured_grid), target :: ugrid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: points(:, :), scalars(:)
   integer, allocatable :: conn(:), offsets(:), types(:)
   integer :: ierr
   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file

   error_cnt = 0
   grid_base => ugrid

   ! Points
   ugrid%num_points = 4
   allocate(points(3,4))
   points(:,1)=(/0.0,0.0,0.0/)
   points(:,2)=(/1.0,0.0,0.0/)
   points(:,3)=(/0.0,1.0,0.0/)
   points(:,4)=(/0.0,0.0,1.0/)
   call grid_base%add_points(points)

   ! Cells
   allocate(conn(4)); conn = (/0,1,2,3/)
   allocate(offsets(1)); offsets = (/4/)
   allocate(types(1)); types = (/10/)
   call ugrid%add_cell_connectivity(conn, offsets, types)

   ! Scalars
   allocate(scalars(4)); scalars = (/1.0,2.0,3.0,4.0/)
   call grid_base%add_scalar('Velocity', scalars)

   file = join_path(folder, 'test.vtu')
   call create_folder(folder, ierr)
   call grid_base%write_file(file)

   open(unit=10, file=file, status='old', action='read', iostat=ierr)
   if (ierr /= 0) then
      error_cnt = error_cnt + 1
   else
      close(10)
   end if
   call remove_folder(folder, ierr)
end function

!==============================
! Test 8: VTU file with cell data
!==============================
integer function test_vtkAPI_vtu_cell_data() bind(C) result(error_cnt)
   use vtkAPI_m
   implicit none
   type(vtk_unstructured_grid), target :: ugrid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: cell_scalars(:), cell_vectors(:)
   integer, allocatable :: conn(:), offsets(:), types(:)

   error_cnt = 0
   grid_base => ugrid

   ! Cells
   ugrid%num_cells = 1
   allocate(conn(4)); conn = (/0,1,2,3/)
   allocate(offsets(1)); offsets = (/4/)
   allocate(types(1)); types = (/10/)
   call ugrid%add_cell_connectivity(conn, offsets, types)

   ! Cell scalar
   allocate(cell_scalars(1)); cell_scalars(1) = 5.0
   call grid_base%add_cell_scalar('Pressure', cell_scalars)
   if (ugrid%cell_scalars(1)%data(1) /= 5.0) error_cnt = error_cnt + 1

   ! Cell vector
   allocate(cell_vectors(3*1)); cell_vectors = (/1.0,2.0,3.0/)
   call grid_base%add_cell_vector('Flux', cell_vectors)
   if (ugrid%cell_vectors(1)%data(3) /= 3.0) error_cnt = error_cnt + 1
end function

!==============================
! Test 9: Verificación de VTS contenido
!==============================
integer function test_vtkAPI_vts_content() bind(C) result(error_cnt)
   use vtkAPI_m
   use directoryUtils_m
   implicit none
   type(vtk_structured_grid), target :: grid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: points(:, :), scalars(:), vectors(:)
   integer :: nx=2, ny=2, nz=2
   integer :: ierr, i
   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file
   character(len=256) :: line
   logical :: found_scalar, found_vector

   error_cnt = 0
   grid%nx=nx; grid%ny=ny; grid%nz=nz
   grid_base => grid

   ! Points
   allocate(points(3, nx*ny*nz)); points=0.0
   call grid_base%add_points(points)

   ! Scalar
   allocate(scalars(nx*ny*nz))
   do i=1,nx*ny*nz
      scalars(i) = real(i)
   end do
   call grid_base%add_scalar('Density', scalars)

   ! Vector
   allocate(vectors(3*nx*ny*nz))
   do i=1,nx*ny*nz
      vectors(3*i-2) = real(i)
      vectors(3*i-1) = real(i*10)
      vectors(3*i)   = real(i*100)
   end do
   call grid_base%add_vector('Momentum', vectors)

   ! Create VTS file
   file = join_path(folder,'test_content.vts')
   call create_folder(folder,ierr)
   call grid_base%write_file(file)

   ! Read file and verify PointData
   found_scalar = .false.; found_vector = .false.
   open(unit=10, file=file, status='old', action='read', iostat=ierr)
   if (ierr /= 0) then
      error_cnt = error_cnt + 1
      return
   end if

   do
      read(10,'(A)', iostat=ierr) line
      if (ierr /= 0) exit
      if (index(line,'Scalars="Density"') /= 0) found_scalar = .true.
      if (index(line,'Vectors="Momentum"') /= 0) found_vector = .true.
      if (found_scalar .and. found_vector) exit
   end do
   close(10)

   if (.not. found_scalar) error_cnt = error_cnt + 1
   if (.not. found_vector) error_cnt = error_cnt + 1

   call remove_folder(folder, ierr)
end function

!==============================
! Test 10: Verificación de VTU contenido
!==============================
integer function test_vtkAPI_vtu_content() bind(C) result(error_cnt)
   use vtkAPI_m
   use directoryUtils_m
   implicit none
   type(vtk_unstructured_grid), target :: ugrid
   class(vtk_grid), pointer :: grid_base
   real, allocatable :: points(:, :), scalars(:), cell_scalars(:)
   integer, allocatable :: conn(:), offsets(:), types(:)
   integer :: ierr
   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file
   character(len=256) :: line
   logical :: found_point_scalar, found_cell_scalar, found_cells, found_points

   error_cnt = 0
   grid_base => ugrid

   ! Points
   ugrid%num_points = 4
   allocate(points(3,4))
   points(:,1)=(/0.0,0.0,0.0/); points(:,2)=(/1.0,0.0,0.0/)
   points(:,3)=(/0.0,1.0,0.0/); points(:,4)=(/0.0,0.0,1.0/)
   call grid_base%add_points(points)

   ! Cells
   ugrid%num_cells = 1
   allocate(conn(4)); conn=(/0,1,2,3/)
   allocate(offsets(1)); offsets=(/4/)
   allocate(types(1)); types=(/10/)
   call ugrid%add_cell_connectivity(conn, offsets, types)

   ! Point scalarll
   allocate(scalars(4)); scalars=(/1.0,2.0,3.0,4.0/)
   call grid_base%add_scalar('Velocity', scalars)

   ! Cell scalar
   allocate(cell_scalars(1)); cell_scalars=(/5.0/)
   call grid_base%add_cell_scalar('Pressure', cell_scalars)

   ! Create VTU file
   file = join_path(folder,'test_content.vtu')
   call create_folder(folder,ierr)
   call grid_base%write_file(file)

   ! Read file
   found_point_scalar = .false.; found_cell_scalar = .false.
   found_cells = .false.; found_points = .false.
   open(unit=10, file=file, status='old', action='read', iostat=ierr)
   if (ierr /= 0) then
      error_cnt = error_cnt + 1
      return
   end if

   do
      read(10,'(A)', iostat=ierr) line
      if (ierr /= 0) exit
      if (index(line,'PointData') /= 0) found_point_scalar = .true.
      if (index(line,'CellData') /= 0) found_cell_scalar = .true.
      if (index(line,'<Cells>') /= 0) found_cells = .true.
      if (index(line,'<Points>') /= 0) found_points = .true.
      if (found_point_scalar .and. found_cell_scalar .and. found_cells .and. found_points) exit
   end do
   close(10)

   if (.not. found_point_scalar) error_cnt = error_cnt + 1
   if (.not. found_cell_scalar) error_cnt = error_cnt + 1
   if (.not. found_cells) error_cnt = error_cnt + 1
   if (.not. found_points) error_cnt = error_cnt + 1

   call remove_folder(folder, ierr)
end function