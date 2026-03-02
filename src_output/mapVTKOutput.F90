module mod_mapVTKOutput
   use FDETYPES
   use outputTypes
   use mod_outputUtils
   use vtk_fortran
   use mod_directoryUtils
   use mod_allocationUtils
   use mod_vtkAPI
   use mod_volumicProbeUtils
   use Report

   implicit none
contains
   subroutine init_mapvtk_output(this, lowerBound, upperBound, field, outputTypeExtension, mpidir, problemInfo)
      type(mapvtk_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      type(problem_info_t), target ,intent(in) :: problemInfo
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
      type(problem_info_t), pointer, intent(in) :: problemInfo

      integer :: i, j, k, field, counter

      counter = 0
      do k = this%mainCoords%Z, this%auxCoords%Z
      do j = this%mainCoords%Y, this%auxCoords%Y
      do i = this%mainCoords%X, this%auxCoords%X
         do field = iEx, iEz
            if (isEdge(field, i, j, k, problemInfo)) then
               counter = counter + 1
            end if
         end do
         do field = iHx, iHz
            if (isWithinBounds(field, i, j, k, problemInfo)) then
               if (isMaterialExceptPML(field, i, j, k, problemInfo)) then
                  counter = counter + 1
               end if
               if (problemInfo%materialTag%getFaceTag(field, i, j, k) < 0 &
                   .and. (btest(iabs(problemInfo%materialTag%getFaceTag(field, i, j, k)), field - 1)) &
                   .and. (.not. isPML(field, i, j, k, problemInfo))) then
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
      call alloc_and_init(this%currentType, this%nPoints, -1)

      counter = 0
      do k = this%mainCoords%Z, this%auxCoords%Z
      do j = this%mainCoords%Y, this%auxCoords%Y
      do i = this%mainCoords%X, this%auxCoords%X
         do field = iEx, iEz
            if (isEdge(field, i, j, k, problemInfo)) then
               counter = counter + 1
               call writeFaceTagInfo(this, counter, i, j, k, field, problemInfo%materialTag%getFaceTag(field, i, j, k))
            end if
         end do
         do field = iHx, iHz
            if (isWithinBounds(field, i, j, k, problemInfo)) then
               if (isMaterialExceptPML(field, i, j, k, problemInfo)) then
                  counter = counter + 1
                  call writeFaceTagInfo(this, counter, i, j, k, field, problemInfo%materialTag%getFaceTag(field, i, j, k))
               end if
               if (problemInfo%materialTag%getFaceTag(field, i, j, k) < 0 &
                   .and. (btest(iabs(problemInfo%materialTag%getFaceTag(field, i, j, k)), field - 1)) &
                   .and. .not. isPML(field, i, j, k, problemInfo)) then
                  counter = counter + 1
                  call writeFaceTagInfo(this, counter, i, j, k, field, problemInfo%materialTag%getFaceTag(field, i, j, k))
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

   subroutine writeFaceTagInfo(this, counter, i, j, k, field, tag)
      type(mapvtk_output_t), intent(inout) :: this
      integer, intent(in) :: i, j, k, counter, tag, field
      this%coords(1, counter) = i
      this%coords(2, counter) = j
      this%coords(3, counter) = k
      this%currentType(counter) = currentType(field)
      this%materialTag(counter) = tag
   end subroutine

   subroutine create_geometry_simulation_vtu(this, control, realXGrid, realYGrid, realZGrid)
      implicit none

      type(mapvtk_output_t), intent(in) :: this
      type(sim_control_t), intent(in) :: control
      real(KIND=RKIND), pointer, dimension(:), intent(in) :: realXGrid, realYGrid, realZGrid

      type(vtk_file) :: vtkOutput
      type(vtk_unstructured_grid), target :: ugrid

      integer :: ierr, i, npts, unit
      real(RKIND), allocatable :: x(:), y(:), z(:), materialTag(:)
      character(len=BUFSIZE) :: info_str
      character(len=BUFSIZE) :: metadata_filename, vtuPath

      integer, allocatable :: conn(:), offsets(:), types(:)
      integer :: numNodes, numEdges, numQuads
      real(kind=RKIND), allocatable :: nodes(:, :)
      integer(kind=SINGLE), allocatable :: edges(:, :), quads(:, :)

      call create_folder(this%path, ierr)
      vtuPath = join_path(this%path, get_last_component(this%path))//vtuFileExtension

      call createUnstructuredDataForVTU(this%nPoints, this%coords, this%currentType, nodes, edges, quads, numNodes, numEdges, numQuads, control%vtkindex, realXGrid, realYGrid, realZGrid)
      call ugrid%add_points(real(nodes, 4))

      allocate (conn(2*numEdges + 4*numQuads))
      conn(1:2*numEdges) = reshape(edges, [2*numEdges])
      conn(2*numEdges + 1:2*numEdges + 4*numQuads) = reshape(quads, [4*numQuads])

      allocate (offsets(numEdges + numQuads))
      do i = 1, numEdges + numQuads
         if (i <= numEdges) then
            if (i == 1) then
               offsets(i) = 2
            else
               offsets(i) = offsets(i - 1) + 2
            end if
         else
            if (i == 1) then
               offsets(i) = 4
            else
               offsets(i) = offsets(i - 1) + 4
            end if
         end if
      end do

      allocate (types(numEdges + numQuads))
      types(1:numEdges) = 3
      types(numEdges + 1:numEdges + numQuads) = 9

      call ugrid%add_cell_connectivity(conn, offsets, types)
      if (size(offsets) /= numQuads) then
         print *, "Problema con offsets"
      end if
      if (size(types) /= numQuads) then
         print *, "Problema con types"
      end if
      if (offsets(numQuads) /= size(conn)) then
         print *, "Tenemos un problema con conn y offset"
      end if

      call ugrid%write_file(vtuPath)

      !---------------------------------------------
      ! Metadata: write to .txt file
      !---------------------------------------------
      info_str = 'PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, '// &
                 'WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, '// &
                 'CONF=5/6, OTHER=-1 (ADD +0.5 for borders)'

      metadata_filename = trim(this%path)//'.txt'
      open (newunit=unit, file=metadata_filename, status='replace', action='write', iostat=ierr)
      if (ierr /= 0) then
         print *, 'Error opening metadata file: ', metadata_filename
      else
         write (unit, '(A)') trim(info_str)
         close (unit)
      end if

   end subroutine create_geometry_simulation_vtu

   logical function isEdge(campo, iii, jjj, kkk, problemInfo)
      integer(4), intent(in) :: campo, iii, jjj, kkk
      type(problem_info_t), pointer, intent(in) :: problemInfo
      
      type(MediaData_t), pointer, dimension(:) :: mData
      type(limit_t), pointer, dimension(:) :: problemDimension
      
      integer(4) :: imed, imed1, imed2, imed3, imed4, contaborde
      
      mData => problemInfo%materialList
      problemDimension => problemInfo%problemDimension
      isEdge = .false.
      contaborde = 0

      call get_media_from_coord_and_h_neighbours(campo, iii, jjj, kkk,  problemInfo%geometryToMaterialData, imed, imed1, imed2, imed3, imed4)

      if (imed /= 1) then

         if (mData(imed)%is%SGBC) then

            if (mData(imed1)%is%SGBC) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed1)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed1 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed2)%is%SGBC) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed2)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed2 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed3)%is%SGBC) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed3)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed3 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed4)%is%SGBC) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed4)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed4 /= 1) then
               contaborde = contaborde + 1
            end if

         elseif (mData(imed)%is%Multiport) then

            if (mData(imed1)%is%Multiport) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed1)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed1 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed2)%is%Multiport) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed2)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed2 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed3)%is%Multiport) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed3)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed3 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed4)%is%Multiport) then
               if (trim(adjustl(mData(imed)%Multiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed4)%Multiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed4 /= 1) then
               contaborde = contaborde + 1
            end if

         elseif (mData(imed)%is%AnisMultiport) then

            if (mData(imed1)%is%AnisMultiport) then
               if (trim(adjustl(mData(imed)%AnisMultiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed1)%AnisMultiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed1 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed2)%is%AnisMultiport) then
               if (trim(adjustl(mData(imed)%AnisMultiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed2)%AnisMultiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed2 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed3)%is%AnisMultiport) then
               if (trim(adjustl(mData(imed)%AnisMultiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed3)%AnisMultiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed3 /= 1) then
               contaborde = contaborde + 1
            end if

            if (mData(imed4)%is%AnisMultiport) then
               if (trim(adjustl(mData(imed)%AnisMultiport(1)%MultiportFileZ11)) /= &
                   trim(adjustl(mData(imed4)%AnisMultiport(1)%MultiportFileZ11))) contaborde = contaborde + 1
            elseif (imed4 /= 1) then
               contaborde = contaborde + 1
            end if

         else
            if ((imed /= imed1) .and. (imed1 /= 1)) contaborde = contaborde + 1
            if ((imed /= imed2) .and. (imed2 /= 1)) contaborde = contaborde + 1
            if ((imed /= imed3) .and. (imed3 /= 1)) contaborde = contaborde + 1
            if ((imed /= imed4) .and. (imed4 /= 1)) contaborde = contaborde + 1
         end if

         if ((imed1 == 1 .and. imed2 == 1 .and. imed3 == 1 .and. imed4 /= 1) .or. &
             (imed2 == 1 .and. imed3 == 1 .and. imed4 == 1 .and. imed1 /= 1) .or. &
             (imed3 == 1 .and. imed4 == 1 .and. imed1 == 1 .and. imed2 /= 1) .or. &
             (imed4 == 1 .and. imed1 == 1 .and. imed2 == 1 .and. imed3 /= 1) .or. &
             (imed1 == 1 .and. imed2 == 1 .and. imed3 == 1 .and. imed4 == 1) .or. &
             (contaborde > 0)) isEdge = .true.

         if ((iii > problemDimension(campo)%XE) .or. (jjj > problemDimension(campo)%YE) .or. &
             (kkk > problemDimension(campo)%ZE)) isEdge = .false.

         if ((iii < problemDimension(campo)%XI) .or. (jjj < problemDimension(campo)%YI) .or. &
             (kkk < problemDimension(campo)%ZI)) isEdge = .false.

      else
         isEdge = .false.
      end if

   end function

end module mod_mapVTKOutput

