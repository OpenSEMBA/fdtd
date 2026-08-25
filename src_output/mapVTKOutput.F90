module mapVTKOutput_m
   use FDETYPES_m
   use outputTypes_m
   use outputUtils_m
   use directoryUtils_m
   use allocationUtils_m
    use vtkAPI_m
    use volumicProbeUtils_m
    use report_m
    use HollandWires_m, only: GetHwires
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use xdmf_hdf5_m, only: xdmf_writer_t, xdmf_options_t, xdmf_status_t, xdmf_grid_id_t, &
                            xdmf_attribute_id_t, XDMF_TOPOLOGY_POLYLINE, XDMF_TOPOLOGY_QUADRILATERAL, &
                            XDMF_CENTER_CELL, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64
#ifdef CompileWithMTLN
    use Wire_bundles_mtln_m, only: GetSolverPtr
    use mtln_solver_m, only: mtln_solver_t => mtln_t
#endif
     implicit none (type, external)
    private

    public :: init_mapvtk_output, create_geometry_simulation_vtu, write_geometry_companion
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
#ifdef CompileWithMTLN
        call count_mtln_external_segments(this, counter)
#endif
        call count_holland_wire_segments(this, counter)

       this%nPoints = counter
       call alloc_and_init(this%coords, 3, this%nPoints, -99)
       call alloc_and_init(this%materialTag, this%nPoints, -1)
       call alloc_and_init(this%currentType, this%nPoints, -1)
       call alloc_and_init(this%mediaType, this%nPoints, -1.0_RKIND)

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
#ifdef CompileWithMTLN
        call store_mtln_external_segments(this, counter, problemInfo)
#endif
        call store_holland_wire_segments(this, counter, problemInfo)
    end subroutine store_relevant_coordinates

    subroutine count_holland_wire_segments(this, counter)
       type(mapvtk_output_t), intent(in) :: this
       integer, intent(inout) :: counter

       type(Thinwires_t), pointer :: wires
       integer :: segment_index

       wires => GetHwires()
       do segment_index = 1, wires%NumCurrentSegments
          if (wires%CurrentSegment(segment_index)%i < this%mainCoords%x .or. &
              wires%CurrentSegment(segment_index)%i > this%auxCoords%x .or. &
              wires%CurrentSegment(segment_index)%j < this%mainCoords%y .or. &
              wires%CurrentSegment(segment_index)%j > this%auxCoords%y .or. &
              wires%CurrentSegment(segment_index)%k < this%mainCoords%z .or. &
              wires%CurrentSegment(segment_index)%k > this%auxCoords%z) cycle
          counter = counter + 1
       end do
    end subroutine count_holland_wire_segments

    subroutine store_holland_wire_segments(this, counter, problemInfo)
       type(mapvtk_output_t), intent(inout) :: this
       integer, intent(inout) :: counter
       type(problem_info_t), intent(in) :: problemInfo

       type(Thinwires_t), pointer :: wires
       integer :: segment_index, field

       wires => GetHwires()
       do segment_index = 1, wires%NumCurrentSegments
          if (wires%CurrentSegment(segment_index)%i < this%mainCoords%x .or. &
              wires%CurrentSegment(segment_index)%i > this%auxCoords%x .or. &
              wires%CurrentSegment(segment_index)%j < this%mainCoords%y .or. &
              wires%CurrentSegment(segment_index)%j > this%auxCoords%y .or. &
              wires%CurrentSegment(segment_index)%k < this%mainCoords%z .or. &
              wires%CurrentSegment(segment_index)%k > this%auxCoords%z) cycle

          field = wires%CurrentSegment(segment_index)%tipofield
          counter = counter + 1
          this%coords(:, counter) = [wires%CurrentSegment(segment_index)%i, &
                                     wires%CurrentSegment(segment_index)%j, &
                                     wires%CurrentSegment(segment_index)%k]
          this%currentType(counter) = currentType(field)
          this%materialTag(counter) = problemInfo%materialTag%getEdgeTag(field, this%coords(1, counter), &
                                                                           this%coords(2, counter), &
                                                                           this%coords(3, counter))
          if (wires%CurrentSegment(segment_index)%Is_LeftEnd .or. &
              wires%CurrentSegment(segment_index)%Is_RightEnd) then
             this%mediaType(counter) = 10.0_RKIND
          else if (wires%CurrentSegment(segment_index)%IsEnd_norLeft_norRight) then
             this%mediaType(counter) = 11.0_RKIND
          else
             this%mediaType(counter) = 20.0_RKIND + real(wires%CurrentSegment(segment_index)%NumParallel, RKIND)
          end if
       end do
    end subroutine store_holland_wire_segments

#ifdef CompileWithMTLN
    subroutine count_mtln_external_segments(this, counter)
       type(mapvtk_output_t), intent(in) :: this
       integer, intent(inout) :: counter

       type(mtln_solver_t), pointer :: mtln_local
       integer :: bundle_index, segment_index
       integer :: position(3)

       mtln_local => GetSolverPtr()
       if (.not. associated(mtln_local)) return
       if (.not. allocated(mtln_local%bundles)) return
       do bundle_index = 1, size(mtln_local%bundles)
          if (.not. allocated(mtln_local%bundles(bundle_index)%external_field_segments)) cycle
          do segment_index = 1, size(mtln_local%bundles(bundle_index)%external_field_segments)
             position = mtln_local%bundles(bundle_index)%external_field_segments(segment_index)%position
             if (all(position >= [this%mainCoords%x, this%mainCoords%y, this%mainCoords%z]) .and. &
                 all(position <= [this%auxCoords%x, this%auxCoords%y, this%auxCoords%z])) counter = counter + 1
          end do
       end do
    end subroutine count_mtln_external_segments

    subroutine store_mtln_external_segments(this, counter, problemInfo)
       type(mapvtk_output_t), intent(inout) :: this
       integer, intent(inout) :: counter
       type(problem_info_t), intent(in) :: problemInfo

       type(mtln_solver_t), pointer :: mtln_local
       integer :: bundle_index, segment_index, field
       integer :: position(3)

       mtln_local => GetSolverPtr()
       if (.not. associated(mtln_local)) return
       if (.not. allocated(mtln_local%bundles)) return
       do bundle_index = 1, size(mtln_local%bundles)
          if (.not. allocated(mtln_local%bundles(bundle_index)%external_field_segments)) cycle
          do segment_index = 1, size(mtln_local%bundles(bundle_index)%external_field_segments)
             position = mtln_local%bundles(bundle_index)%external_field_segments(segment_index)%position
             if (.not. all(position >= [this%mainCoords%x, this%mainCoords%y, this%mainCoords%z]) .or. &
                 .not. all(position <= [this%auxCoords%x, this%auxCoords%y, this%auxCoords%z])) cycle
             field = abs(mtln_local%bundles(bundle_index)%external_field_segments(segment_index)%direction)
             counter = counter + 1
             this%coords(:, counter) = position
              this%currentType(counter) = currentType(field)
              this%materialTag(counter) = problemInfo%materialTag%getEdgeTag(field, position(1), position(2), position(3))
              if (segment_index == 1 .or. segment_index == size(mtln_local%bundles(bundle_index)%external_field_segments)) then
                 this%mediaType(counter) = 14.0_RKIND
              else
                 this%mediaType(counter) = 60.0_RKIND + real(sum(mtln_local%bundles(bundle_index)%conductors_in_level), RKIND)
              end if
          end do
       end do
    end subroutine store_mtln_external_segments
#endif

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
      implicit none (type, external)

      type(mapvtk_output_t), intent(in) :: this
      type(sim_control_t), intent(in) :: control
       real(KIND=RKIND), pointer, dimension(:), intent(in) :: realXGrid, realYGrid, realZGrid
       type(problem_info_t), target, intent(in) :: problemInfo

      !type(vtk_file) :: vtkOutput
      type(vtk_unstructured_grid), target :: ugrid

       integer :: ierr, i, unit
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

    subroutine write_geometry_companion(base_path, lower_bound, upper_bound, problemInfo, status, diagnostic)
       character(len=*), intent(in) :: base_path
       type(cell_coordinate_t), intent(in) :: lower_bound, upper_bound
       type(problem_info_t), target, intent(in) :: problemInfo
       integer, intent(out) :: status
       character(len=BUFSIZE), intent(out) :: diagnostic

       type(mapvtk_output_t) :: geometry
       type(xdmf_writer_t) :: writer
       type(xdmf_options_t) :: options
       type(xdmf_status_t) :: writer_status, close_status
       type(xdmf_grid_id_t) :: grid
       type(xdmf_attribute_id_t) :: tags_attribute, media_attribute
       integer :: num_nodes, num_edges, num_quads
       real(RKIND), allocatable :: nodes(:, :)
       integer(kind=SINGLE), allocatable :: edges(:, :), quads(:, :)
       real, allocatable :: tags(:), media_types(:)
       integer(int64), allocatable :: connectivity(:, :)

       geometry%mainCoords = lower_bound
       geometry%auxCoords = upper_bound
       call store_relevant_coordinates(geometry, problemInfo)
       call createUnstructuredDataForVTU(geometry%nPoints, geometry%coords, geometry%currentType, nodes, edges, quads, &
                                         num_nodes, num_edges, num_quads, .true., problemInfo%xSteps, &
                                         problemInfo%ySteps, problemInfo%zSteps)
       call build_cell_properties(geometry, problemInfo, num_edges, num_quads, tags, media_types)

       options%overwrite = .true.
       call writer%create(trim(base_path)//'_geometry', options, writer_status)
       if (writer_status%is_error()) then
          status = 1
          diagnostic = writer_status%message()
          return
       end if

       if (num_edges > 0) then
          allocate(connectivity(2, num_edges))
          connectivity = int(edges, int64) + 1_int64
          call writer%define_unstructured_grid('lines', XDMF_TOPOLOGY_POLYLINE, real(nodes, real64), connectivity, grid, writer_status)
          if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'tagnumber', XDMF_CENTER_CELL, &
               XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, tags_attribute, writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(tags_attribute, real(tags(:num_edges), real64), writer_status)
          if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'mediatype', XDMF_CENTER_CELL, &
               XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, media_attribute, writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(media_attribute, real(media_types(:num_edges), real64), writer_status)
          deallocate(connectivity)
       end if

       if (num_quads > 0 .and. .not. writer_status%is_error()) then
          allocate(connectivity(4, num_quads))
          connectivity = int(quads, int64) + 1_int64
          call writer%define_unstructured_grid('faces', XDMF_TOPOLOGY_QUADRILATERAL, real(nodes, real64), connectivity, grid, writer_status)
          if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'tagnumber', XDMF_CENTER_CELL, &
               XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, tags_attribute, writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(tags_attribute, real(tags(num_edges + 1:), real64), writer_status)
          if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'mediatype', XDMF_CENTER_CELL, &
               XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, media_attribute, writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(media_attribute, real(media_types(num_edges + 1:), real64), writer_status)
          deallocate(connectivity)
       end if
       if (writer_status%is_error()) then
          call writer%close(close_status)
          status = 1
          diagnostic = writer_status%message()
          return
       end if
       call writer%close(writer_status)
       if (writer_status%is_error()) then
          status = 1
          diagnostic = writer_status%message()
       else
          status = 0
          diagnostic = ''
       end if
    end subroutine write_geometry_companion

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
               if (this%mediaType(index) >= 0.0_RKIND) then
                  media_types(edge_index) = this%mediaType(index)
               else
                  media_types(edge_index) = edge_media_type(field, this%coords(:, index), problemInfo, media)
               end if
           case (iBloqueJx, iBloqueJy, iBloqueJz)
             quad_index = quad_index + 1
             field = magnetic_field(this%currentType(index))
              media = getMediaIndex(field, this%coords(1, index), this%coords(2, index), this%coords(3, index), &
                                     problemInfo%geometryToMaterialData)
              tags(quad_index) = real(this%materialTag(index))
              media_types(quad_index) = surface_media_type(problemInfo%materialList(media), media)
          end select
       end do
    contains
       real function surface_media_type(material, media)
           type(MediaData_t), intent(in) :: material
           integer, intent(in) :: media

           if (material%is%Pec) then
              surface_media_type = 0.0
           else if (material%is%PMC) then
              surface_media_type = 16.0
           else if (material%is%ConformalPec) then
              surface_media_type = 1000.0 + media
           else if (material%is%SGBC .or. material%is%Multiport .or. material%is%AnisMultiport) then
              surface_media_type = 300.0 + media
          else if (material%is%EDispersive .or. material%is%MDispersive .or. material%is%EDispersiveAnis .or. &
                   material%is%MDispersiveAnis) then
              surface_media_type = 100.0 + media
          else if (material%is%Dielectric .or. material%is%Anisotropic) then
              surface_media_type = 200.0 + media
           else if (material%is%ThinSlot) then
              surface_media_type = 400.0 + media
           else if (material%is%already_YEEadvanced_byconformal) then
              surface_media_type = 5.0
           else if (material%is%split_and_useless) then
              surface_media_type = 6.0
          else
             surface_media_type = -1.0
          end if
       end function surface_media_type

        real function edge_media_type(field, position, problemInfo, media)
           integer, intent(in) :: field, position(3), media
           type(problem_info_t), intent(in) :: problemInfo
           integer :: candidate_field, candidate_media, candidate_position(3)
           type(MediaData_t), pointer :: material

           material => problemInfo%materialList(media)

           if (material%is%already_YEEadvanced_byconformal) then
              edge_media_type = 5.5
           else if (material%is%split_and_useless) then
              edge_media_type = 6.5
           else if (material%is%Pec) then
             edge_media_type = 0.5
           else if (material%is%PMC) then
              edge_media_type = 16.5
           else if (material%is%ConformalPec) then
              edge_media_type = 2000.0 + media
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
               do candidate_field = iEx, iEz
                  candidate_media = getMediaIndex(candidate_field, position(1), position(2), position(3), &
                                                   problemInfo%geometryToMaterialData)
                  if (candidate_media /= 1 .and. candidate_media /= media) then
                     edge_media_type = 8.0
                     return
                  end if
               end do
               candidate_position = position
               candidate_position(field) = candidate_position(field) + 1
               do candidate_field = iEx, iEz
                  candidate_media = getMediaIndex(candidate_field, candidate_position(1), candidate_position(2), &
                                                   candidate_position(3), problemInfo%geometryToMaterialData)
                  if (candidate_media /= 1 .and. candidate_media /= media) then
                     edge_media_type = 8.0
                     return
                  end if
               end do
            else if (material%is%Multiwire) then
              edge_media_type = 12.0
              do candidate_field = iEx, iEz
                 candidate_media = getMediaIndex(candidate_field, position(1), position(2), position(3), &
                                                  problemInfo%geometryToMaterialData)
                 if (candidate_media /= 1 .and. candidate_media /= media) then
                    edge_media_type = 13.0
                    return
                 end if
              end do
              candidate_position = position
              candidate_position(field) = candidate_position(field) + 1
              do candidate_field = iEx, iEz
                 candidate_media = getMediaIndex(candidate_field, candidate_position(1), candidate_position(2), &
                                                  candidate_position(3), problemInfo%geometryToMaterialData)
                 if (candidate_media /= 1 .and. candidate_media /= media) then
                    edge_media_type = 13.0
                    return
                 end if
              end do
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
