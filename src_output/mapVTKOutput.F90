module mod_mapVTKOutput
   use FDETYPES
   use outputTypes
   use mod_outputUtils
   use vtk_fortran
   use mod_directoryUtils
   use mod_allocationUtils

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
      type(mapvtk_output_t), intent(out) :: this
      type(problem_info_t), intent(in) :: problemInfo

      integer :: i, j, k, Hfield, counter
      counter = 0
      do k = this%mainCoords%Z, this%auxCoords%Z
      do j = this%mainCoords%Y, this%auxCoords%Y
      do i = this%mainCoords%X, this%auxCoords%X
         do Hfield = iHx, iHz
            if (isWithinBounds(Hfield, i, j, k, problemInfo)) then
               if isMaterialExceptPML(Hfield, i, j, k, problemInfo) then
               counter = counter + 1
            end if
            if (problemInfo%tagNumbers%getFaceTag(Hfield, i, j, k) < 0 &
                .and. (btest(iabs(problemInfo%tagNumbers%getFaceTag(Hfield, i, j, k)), Hfield - 1)) &
                .and. .not. isPML(Hfield, i, j, k, problemInfo)) then
               counter = counter + 1
            end if
         end do
      end do
      end do
      end do
      this%nPoints = counter
      call alloc_and_init(this%coords, 3, this%nPoints, -99)
      call alloc_and_init(this%tagNumbers, this%nPoints, -1)

      counter = 0
      do k = this%mainCoords%Z, this%auxCoords%Z
      do j = this%mainCoords%Y, this%auxCoords%Y
      do i = this%mainCoords%X, this%auxCoords%X
         do Hfield = iHx, iHz
            if (isWithinBounds(Hfield, i, j, k, problemInfo)) then
               if isMaterialExceptPML(Hfield, i, j, k, problemInfo) then
               counter = counter + 1
               call writeTagInfo(this, counter, i, j, k, problemInfo%tagNumbers%getFaceTag(Hfield, i, j, k))
            end if
            if (problemInfo%tagNumbers%getFaceTag(Hfield, i, j, k) < 0 &
                .and. (btest(iabs(problemInfo%tagNumbers%getFaceTag(Hfield, i, j, k)), Hfield - 1)) &
                .and. .not. isPML(Hfield, i, j, k, problemInfo)) then
               counter = counter + 1
               call writeTagInfo(this, counter, i, j, k, problemInfo%tagNumbers%getFaceTag(Hfield, i, j, k))
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
         isMaterialExceptPML = isMaterialExceptPML .and. (.not. isPML(Hfield, i, j, k, problemInfo))
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

         type(vtk_file) :: vtk

         integer :: ierr
         integer :: i, nCells, nPoints

         integer, allocatable :: connectivity(:)
         integer, allocatable :: offsets(:)
         integer, allocatable :: celltypes(:)

         real(RKIND), allocatable :: points(:, :)
         real(RKIND), allocatable :: media(:)
         real(RKIND), allocatable :: tags(:)

         character(len=BUFSIZE) :: info_str

         !---------------------------------------------
         ! Initialize
         !---------------------------------------------

         ierr = vtk%initialize( &
                format='ASCII', &
                filename=trim(this%path), &
                mesh_topology='UnstructuredGrid')

         if (ierr /= 0) then
            call StopOnError(control%layoutnumber, control%size, &
                             'Error initializing VTK')
         end if

         !---------------------------------------------
         ! MetaData
         !---------------------------------------------
      info_str = 'PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, CONF=5/6, OTHER=-1 (ADD +0.5 for borders)'
         call vtk%write_field_data('Metadata', info_str)

         !---------------------------------------------
         ! Points
         !---------------------------------------------

         nPoints = this%nPoints + 1

         allocate (points(3, nPoints))

         do i = 0, this%nPoints
            points(1, i + 1) = this%coords(1, i)
            points(2, i + 1) = this%coords(2, i)
            points(3, i + 1) = this%coords(3, i)
         end do

         call vtk%write_points(points)

         deallocate (points)

         !---------------------------------------------
         ! Cells
         !---------------------------------------------

         nCells = numberOfSerialized

         allocate (offsets(nCells))
         allocate (celltypes(nCells))

         ! Maximum possible size
         allocate (connectivity(5*nCells))

         integer :: pos
         pos = 0

         do i = 1, nCells

            if (Elems(i, 3) == -1) then
               !
               ! Edge → VTK_LINE = 3
               !

               connectivity(pos + 1) = Elems(i, 1)
               connectivity(pos + 2) = Elems(i, 2)

               pos = pos + 2

               offsets(i) = pos
               celltypes(i) = 3

            else
               !
               ! Quad → VTK_QUAD = 9
               !

               connectivity(pos + 1) = Elems(i, 1)
               connectivity(pos + 2) = Elems(i, 2)
               connectivity(pos + 3) = Elems(i, 3)
               connectivity(pos + 4) = Elems(i, 4)

               pos = pos + 4

               offsets(i) = pos
               celltypes(i) = 9

            end if

         end do

         call vtk%write_cells(connectivity(1:pos), offsets)
         call vtk%write_cell_types(celltypes)

         deallocate (connectivity, offsets, celltypes)

         !---------------------------------------------
         ! Cell data: mediatype
         !---------------------------------------------

         allocate (media(nCells))

         do i = 1, nCells
            media(i) = Serialized%valor(1, i)
         end do

         call vtk%write_cell_data('mediatype', media)

         deallocate (media)

         !---------------------------------------------
         ! Cell data: tagnumber
         !---------------------------------------------

         allocate (tags(nCells))

         do i = 1, nCells
            tags(i) = real(Serialized%sggMtag(i), RKIND)
         end do

         call vtk%write_cell_data('tagnumber', tags)

         deallocate (tags)

         !---------------------------------------------
         ! Finalize
         !---------------------------------------------

         call vtk%finalize()

      end subroutine create_geometry_simulation_vtk

      end module mod_mapVTKOutput

