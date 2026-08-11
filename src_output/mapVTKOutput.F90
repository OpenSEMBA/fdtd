module mapVTKOutput_m
   use FDETYPES_m
   use outputTypes_m
   use outputUtils_m
   use directoryUtils_m
   use allocationUtils_m
   use vtkAPI_m
   use volumicProbeUtils_m
    use report_m
#ifdef CompileWithHDF
    use hdf5
#endif

   implicit none
contains
   subroutine init_mapvtk_output(this, lowerBound, upperBound, field, outputTypeExtension, mpidir, problemInfo)
      type(mapvtk_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      type(problem_info_t), target ,intent(in) :: problemInfo
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      character(len=BUFSIZE) :: artifact_paths(2)
      integer :: artifact_kinds(2)

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
       this%component = field

       this%path = get_output_path()
       artifact_paths(1) = trim(join_path(this%path, get_last_component(this%path)))//vtuFileExtension
        artifact_paths(2) = trim(join_path(this%path, get_last_component(this%path)))//'.txt'
       artifact_kinds = [OUTPUT_ARTIFACT_GEOMETRY, OUTPUT_ARTIFACT_TEXT]
       call declare_probe_artifacts(this%artifacts, artifact_paths, artifact_kinds)
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
                call writeFaceTagInfo(this, counter, i, j, k, field, problemInfo%materialTag%getEdgeTag(field, i, j, k))
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

    subroutine create_geometry_simulation_vtu(this, control, realXGrid, realYGrid, realZGrid, problemInfo)
      implicit none

      type(mapvtk_output_t), intent(in) :: this
      type(sim_control_t), intent(in) :: control
       real(KIND=RKIND), pointer, dimension(:), intent(in) :: realXGrid, realYGrid, realZGrid
       type(problem_info_t), intent(in) :: problemInfo

      !type(vtk_file) :: vtkOutput
      type(vtk_unstructured_grid), target :: ugrid

       integer :: ierr, i, npts, unit
       real(RKIND), allocatable :: x(:), y(:), z(:), materialTag(:)
      character(len=BUFSIZE) :: info_str
      character(len=BUFSIZE) :: metadata_filename, vtuPath

      integer, allocatable :: conn(:), offsets(:), types(:)
       integer :: numNodes, numEdges, numQuads
       real(kind=RKIND), allocatable :: nodes(:, :)
       integer(kind=SINGLE), allocatable :: edges(:, :), quads(:, :)
       real, allocatable :: cell_tags(:), cell_media_types(:)

      call create_folder(this%path, ierr)
       vtuPath = this%artifacts(1)%relative_path

      call createUnstructuredDataForVTU(this%nPoints, this%coords, this%currentType, nodes, edges, quads, numNodes, numEdges, numQuads, control%vtkindex, realXGrid, realYGrid, realZGrid)
      call ugrid%add_points(real(nodes, 4))

       if (numEdges + numQuads > 0) then
          allocate (conn(2*numEdges + 4*numQuads))
          if (numEdges > 0) conn(1:2*numEdges) = reshape(edges, [2*numEdges])
          if (numQuads > 0) conn(2*numEdges + 1:2*numEdges + 4*numQuads) = reshape(quads, [4*numQuads])

          allocate (offsets(numEdges + numQuads))
          do i = 1, numEdges + numQuads
             if (i <= numEdges) then
                offsets(i) = 2*i
             else
                offsets(i) = 2*numEdges + 4*(i - numEdges)
             end if
          end do

          allocate (types(numEdges + numQuads))
           if (numEdges > 0) types(1:numEdges) = 3
           if (numQuads > 0) types(numEdges + 1:numEdges + numQuads) = 9

           call ugrid%add_cell_connectivity(conn, offsets, types)
           call build_cell_properties(this, problemInfo, numEdges, numQuads, cell_tags, cell_media_types)
           call ugrid%add_cell_scalar('tagnumber', cell_tags)
           call ugrid%add_cell_scalar('mediatype', cell_media_types)
#ifdef CompileWithHDF
           call write_geometry_xdmf_hdf5(vtuPath, cell_tags, cell_media_types)
#endif
        end if

      call ugrid%write_file(vtuPath)

      !---------------------------------------------
      ! Metadata: write to .txt file
      !---------------------------------------------
      info_str = 'PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, '// &
                 'WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, '// &
                 'CONF=5/6, OTHER=-1 (ADD +0.5 for borders)'

       metadata_filename = this%artifacts(2)%relative_path
      open (newunit=unit, file=metadata_filename, status='replace', action='write', iostat=ierr)
      if (ierr /= 0) then
         print *, 'Error opening metadata file: ', metadata_filename
      else
         write (unit, '(A)') trim(info_str)
         close (unit)
      end if

    end subroutine create_geometry_simulation_vtu

    subroutine build_cell_properties(this, problemInfo, numEdges, numQuads, tags, media_types)
       type(mapvtk_output_t), intent(in) :: this
       type(problem_info_t), intent(in) :: problemInfo
       integer, intent(in) :: numEdges, numQuads
       real, allocatable, intent(out) :: tags(:), media_types(:)
       integer :: edge_index, quad_index, index, field, media

       allocate(tags(numEdges + numQuads), media_types(numEdges + numQuads))
       edge_index = 0
       quad_index = numEdges
       do index = 1, this%nPoints
          select case (this%currentType(index))
          case (iJx, iJy, iJz)
             edge_index = edge_index + 1
             field = electric_field(this%currentType(index))
             media = getMediaIndex(field, this%coords(1, index), this%coords(2, index), this%coords(3, index), &
                                    problemInfo%geometryToMaterialData)
             tags(edge_index) = real(this%materialTag(index))
             media_types(edge_index) = edge_media_type(problemInfo%materialList(media))
          case (iBloqueJx, iBloqueJy, iBloqueJz)
             quad_index = quad_index + 1
             field = magnetic_field(this%currentType(index))
             media = getMediaIndex(field, this%coords(1, index), this%coords(2, index), this%coords(3, index), &
                                    problemInfo%geometryToMaterialData)
             tags(quad_index) = real(this%materialTag(index))
             media_types(quad_index) = surface_media_type(problemInfo%materialList(media))
          end select
       end do
    contains
       real function surface_media_type(material)
          type(MediaData_t), intent(in) :: material

          if (material%is%Pec) then
             surface_media_type = 0.0
          else if (material%is%PMC) then
             surface_media_type = 16.0
          else if (material%is%SGBC .or. material%is%Multiport .or. material%is%AnisMultiport) then
             surface_media_type = 300.0 + material%id
          else if (material%is%EDispersive .or. material%is%MDispersive .or. material%is%EDispersiveAnis .or. &
                   material%is%MDispersiveAnis) then
             surface_media_type = 100.0 + material%id
          else if (material%is%Dielectric .or. material%is%Anisotropic) then
             surface_media_type = 200.0 + material%id
          else if (material%is%ThinSlot) then
             surface_media_type = 400.0 + material%id
          else
             surface_media_type = -1.0
          end if
       end function surface_media_type

       real function edge_media_type(material)
          type(MediaData_t), intent(in) :: material

          if (material%is%Pec) then
             edge_media_type = 0.5
          else if (material%is%PMC) then
             edge_media_type = 16.5
          else if (material%is%SGBC .or. material%is%Multiport .or. material%is%AnisMultiport) then
             edge_media_type = 3.5
          else if (material%is%EDispersive .or. material%is%MDispersive .or. material%is%EDispersiveAnis .or. &
                   material%is%MDispersiveAnis) then
             edge_media_type = 1.5
          else if (material%is%Dielectric .or. material%is%Anisotropic) then
             edge_media_type = 2.5
          else if (material%is%ThinSlot) then
             edge_media_type = 4.5
          else if (material%is%ThinWire) then
             edge_media_type = 7.0
          else if (material%is%Multiwire) then
             edge_media_type = 12.0
          else
             edge_media_type = -0.5
          end if
       end function edge_media_type

       integer function electric_field(current_type)
          integer(kind=SINGLE), intent(in) :: current_type

          select case (current_type)
          case (iJx); electric_field = iEx
          case (iJy); electric_field = iEy
          case (iJz); electric_field = iEz
          end select
       end function electric_field

       integer function magnetic_field(current_type)
          integer(kind=SINGLE), intent(in) :: current_type

          select case (current_type)
          case (iBloqueJx); magnetic_field = iHx
          case (iBloqueJy); magnetic_field = iHy
          case (iBloqueJz); magnetic_field = iHz
          end select
       end function magnetic_field
    end subroutine build_cell_properties

#ifdef CompileWithHDF
    subroutine write_geometry_xdmf_hdf5(vtu_path, tags, media_types)
       character(len=*), intent(in) :: vtu_path
       real, intent(in) :: tags(:), media_types(:)
       character(len=BUFSIZE) :: hdf_path, xdmf_path, hdf_name
       integer(hid_t) :: file_id, space_id, dataset_id = -1
       integer(hsize_t) :: dimensions(1)
       integer :: hdf_error, unit

       hdf_path = trim(vtu_path(:len_trim(vtu_path) - len(vtuFileExtension)))//'.h5'
       xdmf_path = trim(vtu_path(:len_trim(vtu_path) - len(vtuFileExtension)))//'.xdmf'
       hdf_name = get_last_component(hdf_path)
       dimensions = [int(size(tags), hsize_t)]
       call h5open_f(hdf_error)
       call h5fcreate_f(trim(hdf_path), H5F_ACC_TRUNC_F, file_id, hdf_error)
       if (hdf_error /= 0) return
       call h5screate_simple_f(1, dimensions, space_id, hdf_error)
       call h5dcreate_f(file_id, 'tagnumber', H5T_NATIVE_REAL, space_id, dataset_id, hdf_error)
       if (hdf_error == 0) call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, tags, dimensions, hdf_error)
       if (dataset_id >= 0) call h5dclose_f(dataset_id, hdf_error)
       call h5dcreate_f(file_id, 'mediatype', H5T_NATIVE_REAL, space_id, dataset_id, hdf_error)
       if (hdf_error == 0) call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, media_types, dimensions, hdf_error)
       if (dataset_id >= 0) call h5dclose_f(dataset_id, hdf_error)
       call h5sclose_f(space_id, hdf_error)
       call h5fclose_f(file_id, hdf_error)
       call h5close_f(hdf_error)

       open(newunit=unit, file=trim(xdmf_path), status='replace', action='write', iostat=hdf_error)
       if (hdf_error /= 0) return
       write(unit, '(A)') '<Xdmf Version="3.0"><Domain><Grid Name="geometry" GridType="Uniform">'
       write(unit, '(A,I0,A)') '<Topology TopologyType="Polyvertex" NumberOfElements="', size(tags), '"/>'
       write(unit, '(A,I0,A)') '<Attribute Name="tagnumber" Center="Cell"><DataItem Dimensions="', size(tags), &
                                '" Format="HDF">'//trim(hdf_name)//':/tagnumber</DataItem></Attribute>'
       write(unit, '(A,I0,A)') '<Attribute Name="mediatype" Center="Cell"><DataItem Dimensions="', size(tags), &
                                '" Format="HDF">'//trim(hdf_name)//':/mediatype</DataItem></Attribute>'
       write(unit, '(A)') '</Grid></Domain></Xdmf>'
       close(unit)
    end subroutine write_geometry_xdmf_hdf5
#endif

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

end module mapVTKOutput_m
