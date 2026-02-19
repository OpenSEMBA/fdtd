module mod_mapVTKOutput
   use FDETYPES
   use outputTypes
   use mod_outputUtils
   use vtk_fortran
   use mod_directoryUtils
   use mod_allocationUtils
   use Report

   implicit none
contains
   subroutine init_mapvtk_output(this, lowerBound, upperBound, field, outputTypeExtension, mpidir, problemInfo)
      type(mapvtk_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      type(problem_info_t), intent(in) :: problemInfo
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      integer(kind=SINGLE) :: i

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
      this%component = field

      this%path = get_output_path()
      call store_relevant_coordinates(this, problemInfo)

   contains

      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%mainCoords, this%auxCoords, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function

      subroutine store_media_tag()
      end subroutine store_media_tag
   end subroutine init_mapvtk_output

   subroutine store_relevant_coordinates(this, problemInfo)
      type(mapvtk_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo

      integer :: i, j, k, Hfield, counter

      counter = 0
      do k = this%mainCoords%Z, this%auxCoords%Z
      do j = this%mainCoords%Y, this%auxCoords%Y
      do i = this%mainCoords%X, this%auxCoords%X
         do Hfield = iHx, iHz
            if (isWithinBounds(Hfield, i, j, k, problemInfo)) then
               if (isMaterialExceptPML(Hfield, i, j, k, problemInfo)) then
                  counter = counter + 1
               end if
               if (problemInfo%materialTag%getFaceTag(Hfield, i, j, k) < 0 &
                   .and. (btest(iabs(problemInfo%materialTag%getFaceTag(Hfield, i, j, k)), Hfield - 1)) &
                   .and. (.not. isPML(Hfield, i, j, k, problemInfo))) then
                  counter = counter + 1
               end if
            end if
         end do
      end do
      end do
      end do

      this%nPoints = counter
      call alloc_and_init(this%coords, 3, this%nPoints, -99)
      call alloc_and_init(this%materialTag, this%nPoints, -1)

      counter = 0
      do k = this%mainCoords%Z, this%auxCoords%Z
      do j = this%mainCoords%Y, this%auxCoords%Y
      do i = this%mainCoords%X, this%auxCoords%X
         do Hfield = iHx, iHz
            if (isWithinBounds(Hfield, i, j, k, problemInfo)) then
               if (isMaterialExceptPML(Hfield, i, j, k, problemInfo)) then
                  counter = counter + 1
                  call writeFaceTagInfo(this, counter, i, j, k, problemInfo%materialTag%getFaceTag(Hfield, i, j, k))
               end if
               if (problemInfo%materialTag%getFaceTag(Hfield, i, j, k) < 0 &
                   .and. (btest(iabs(problemInfo%materialTag%getFaceTag(Hfield, i, j, k)), Hfield - 1)) &
                   .and. .not. isPML(Hfield, i, j, k, problemInfo)) then
                  counter = counter + 1
                  call writeFaceTagInfo(this, counter, i, j, k, problemInfo%materialTag%getFaceTag(Hfield, i, j, k))
               end if
            end if
         end do
      end do
      end do
      end do
   end subroutine store_relevant_coordinates

   logical function isMaterialExceptPML(field, i, j, k, problemInfo)
      integer, intent(in)             :: field, i, j, k
      type(problem_info_t), intent(in):: problemInfo
      isMaterialExceptPML = .not. isMediaVacuum(field, i, j, k, problemInfo)
      isMaterialExceptPML = isMaterialExceptPML .and. (.not. isPML(field, i, j, k, problemInfo))
   end function isMaterialExceptPML

   subroutine writeFaceTagInfo(this, counter, i, j, k, tag)
      type(mapvtk_output_t), intent(inout) :: this
      integer, intent(in) :: i, j, k, counter, tag
      this%coords(1, counter) = i
      this%coords(2, counter) = j
      this%coords(3, counter) = k
      this%materialTag(counter) = tag
   end subroutine

   subroutine create_geometry_simulation_vtk(this, problemInfo, control)
      implicit none

      type(mapvtk_output_t), intent(in) :: this
      type(problem_info_t), intent(in) :: problemInfo
      type(sim_control_t), intent(in) :: control
      type(vtk_file) :: vtkOutput

      integer :: ierr, i, npts, unit
      real(RKIND), allocatable :: x(:), y(:), z(:), materialTag(:)
      character(len=BUFSIZE) :: info_str
      character(len=BUFSIZE) :: metadata_filename, vtkPath

      !--------------------------------------------- 
      ! Initialize VTK file 
      !--------------------------------------------- 
      call create_folder(this%path, ierr)
      vtkPath = join_path(this%path, get_last_component(this%path))//vtkFileExtension
      ierr = vtkOutput%initialize(format='ASCII', filename=trim(vtkPath), mesh_topology='UnstructuredGrid') 
      if (ierr /= 0) call StopOnError(control%layoutnumber, control%size, 'Error initializing VTK') 

      !---------------------------------------------
      ! Points
      !---------------------------------------------
      npts = this%nPoints
      allocate (x(npts), y(npts), z(npts))

      do i = 1, npts
         x(i) = this%coords(1, i)
         y(i) = this%coords(2, i)
         z(i) = this%coords(3, i)
      end do

      !---------------------------------------------
      ! Metadata: write to .txt file
      !---------------------------------------------
      info_str = 'PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, '// &
                 'WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, '// &
                 'CONF=5/6, OTHER=-1 (ADD +0.5 for borders)'

      ! Create .txt file with same base name as VTK file
      metadata_filename = trim(this%path)//'.txt'
      open (newunit=unit, file=metadata_filename, status='replace', action='write', iostat=ierr)
      if (ierr /= 0) then
         print *, 'Error opening metadata file: ', metadata_filename
      else
         write (unit, '(A)') trim(info_str)
         close (unit)
      end if

      !---------------------------------------------
      ! Write geometry to VTK
      !---------------------------------------------
      ierr = vtkOutput%xml_writer%write_geo(n=npts, x=x, y=y, z=z)

      !---------------------------------------------
      ! Node data: material tags
      !---------------------------------------------
      allocate (materialTag(npts))
      do i = 1, npts
         materialTag(i) = this%materialTag(i)
      end do

      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='open')
      ierr = vtkOutput%xml_writer%write_dataarray(data_name='TagNumber', x=materialTag)
      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='close')

      !---------------------------------------------
      ! Clean up
      !---------------------------------------------
      deallocate (materialTag)
      deallocate (x, y, z)

      !---------------------------------------------
      ! Finalize VTK file
      !---------------------------------------------
      ierr = vtkOutput%xml_writer%finalize()

   end subroutine create_geometry_simulation_vtk

end module mod_mapVTKOutput

