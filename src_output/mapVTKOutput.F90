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
   use outputTransport_m, only: output_transport_t, init_output_transport, transfer_flush_batch, &
                                 OUTPUT_TRANSPORT_SUCCESS
   use xdmf_hdf5_m, only: xdmf_writer_t, xdmf_options_t, xdmf_status_t, xdmf_grid_id_t, &
                          xdmf_attribute_id_t, XDMF_TOPOLOGY_POLYLINE, XDMF_TOPOLOGY_QUADRILATERAL, &
                          XDMF_CENTER_CELL, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64
#ifdef CompileWithMTLN
   use Wire_bundles_mtln_m, only: GetSolverPtr
   use mtln_solver_m, only: mtln_solver_t => mtln_t
#endif
   implicit none(type, external)
   private

   public :: init_mapvtk_output, create_geometry_simulation_vtu, create_parallel_geometry_vtu, &
             write_geometry_companion
contains
   subroutine init_mapvtk_output(this, lowerBound, upperBound, globalLowerBound, globalUpperBound, &
                                 localParticipates, field, outputTypeExtension, mpidir, problemInfo)
      type(mapvtk_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      type(cell_coordinate_t), intent(in) :: globalLowerBound, globalUpperBound
      type(problem_info_t), target, intent(in) :: problemInfo
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      logical, intent(in) :: localParticipates

      character(len=BUFSIZE) :: artifact_paths(1)
      integer :: artifact_kinds(1)

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
      this%component = field
      this%localParticipates = localParticipates
      this%masterPath = trim(get_map_output_path(globalLowerBound, globalUpperBound, field, &
                                                 outputTypeExtension, mpidir))//pvtuFileExtension

      artifact_paths = ''
      if (localParticipates) then
         this%path = get_map_output_path(lowerBound, upperBound, field, outputTypeExtension, mpidir)
         artifact_paths(1) = trim(join_path(this%path, get_last_component(this%path)))//vtuFileExtension
      end if
      artifact_kinds = [OUTPUT_ARTIFACT_GEOMETRY]
      call declare_probe_artifacts(this%artifacts, artifact_paths, artifact_kinds)
      if (localParticipates) call store_relevant_coordinates(this, problemInfo)

   end subroutine init_mapvtk_output

   function get_map_output_path(lowerBound, upperBound, field, outputTypeExtension, mpidir) result(outputPath)
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: field, mpidir
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      character(len=BUFSIZE) :: outputPath
      character(len=BUFSIZE) :: probeBoundsExtension, prefixFieldExtension

      probeBoundsExtension = get_coordinates_extension(lowerBound, upperBound, mpidir)
      prefixFieldExtension = get_prefix_extension(field, mpidir)
      outputPath = trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'// &
                   trim(adjustl(probeBoundsExtension))
   end function get_map_output_path

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
      implicit none(type, external)

      type(mapvtk_output_t), intent(in) :: this
      type(sim_control_t), intent(in) :: control
      real(KIND=RKIND), pointer, dimension(:), intent(in) :: realXGrid, realYGrid, realZGrid
      type(problem_info_t), target, intent(in) :: problemInfo

      !type(vtk_file) :: vtkOutput
      type(vtk_unstructured_grid), target :: ugrid

      integer :: ierr, i
      character(len=BUFSIZE) :: vtuPath

      integer, allocatable :: conn(:), offsets(:), types(:)
      integer :: numNodes, numEdges, numQuads
      real(kind=RKIND), allocatable :: nodes(:, :)
      integer(kind=SINGLE), allocatable :: edges(:, :), quads(:, :)
      real, allocatable :: cell_tags(:), cell_media_types(:)

      call create_folder(this%path, ierr)
      vtuPath = this%artifacts(1)%relative_path

      call createUnstructuredDataForVTU(this%nPoints, this%coords, this%currentType, nodes, edges, quads, numNodes, numEdges, numQuads, control%vtkindex, realXGrid, realYGrid, realZGrid)
      call ugrid%add_points(real(nodes, 4))

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

      call ugrid%write_file(vtuPath)

   end subroutine create_geometry_simulation_vtu

   subroutine create_parallel_geometry_vtu(this, piecePaths)
      type(mapvtk_output_t), intent(in) :: this
      character(len=*), intent(in) :: piecePaths(:)
      character(len=10), parameter :: cellScalarNames(2) = ['tagnumber ', 'mediatype ']

      call write_pvtu_file(this%masterPath, piecePaths, cellScalarNames)
   end subroutine create_parallel_geometry_vtu

    subroutine write_geometry_companion(base_path, lower_bound, upper_bound, problemInfo, communicator, status, diagnostic)
       character(len=*), intent(in) :: base_path
       type(cell_coordinate_t), intent(in) :: lower_bound, upper_bound
       type(problem_info_t), target, intent(in) :: problemInfo
       integer, intent(in) :: communicator
       integer, intent(out) :: status
       character(len=BUFSIZE), intent(out) :: diagnostic

       type(mapvtk_output_t) :: geometry
       type(output_transport_t) :: transport
       type(xdmf_writer_t) :: writer
       type(xdmf_options_t) :: options
       type(xdmf_status_t) :: writer_status, close_status
       integer :: num_nodes, num_edges, num_quads
       real(RKIND), allocatable :: nodes(:, :)
       integer(kind=SINGLE), allocatable :: edges(:, :), quads(:, :)
       real, allocatable :: tags(:), media_types(:)
       real(real64), allocatable :: local_fragment(:), gathered_fragments(:)
       integer, allocatable :: counts(:), displacements(:)
       integer :: transport_status
       logical :: invalid_fragment

       geometry%mainCoords = lower_bound
       geometry%auxCoords = upper_bound
      call store_relevant_coordinates(geometry, problemInfo)
      call createUnstructuredDataForVTU(geometry%nPoints, geometry%coords, geometry%currentType, nodes, edges, quads, &
                                         num_nodes, num_edges, num_quads, .false., problemInfo%xSteps, &
                                         problemInfo%ySteps, problemInfo%zSteps)
       call build_cell_properties(geometry, problemInfo, num_edges, num_quads, tags, media_types)
       call serialise_geometry_fragment(nodes, edges, quads, tags, media_types, local_fragment)

       call init_output_transport(transport, root_rank=0, status=transport_status, communicator=communicator)
       if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
          status = 1
          diagnostic = 'Unable to initialise movie geometry transport'
          return
       end if
       call transfer_flush_batch(transport, local_fragment, gathered_fragments, counts, displacements, transport_status)
       if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
          status = 1
          diagnostic = 'Unable to gather movie geometry fragments'
          return
       end if

       ! Only the transport root opens the companion; every participant supplied its local material fragment above.
       if (transport%rank /= transport%root_rank) then
          status = 0
          diagnostic = ''
          return
       end if

       options%overwrite = .true.
       call writer%create(trim(base_path)//'_geometry', options, writer_status)
      if (writer_status%is_error()) then
         status = 1
          diagnostic = writer_status%message()
          return
       end if
       call write_geometry_fragments(writer, gathered_fragments, counts, displacements, invalid_fragment, writer_status)
       if (invalid_fragment) then
          call writer%close(close_status)
          status = 1
          diagnostic = 'Invalid movie geometry fragment'
          return
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

    subroutine serialise_geometry_fragment(nodes, edges, quads, tags, media_types, fragment)
       real(RKIND), intent(in) :: nodes(:, :)
       integer(kind=SINGLE), intent(in) :: edges(:, :), quads(:, :)
       real, intent(in) :: tags(:), media_types(:)
       real(real64), allocatable, intent(out) :: fragment(:)
       integer :: cursor, num_edges, num_quads

       num_edges = size(edges, 2)
       num_quads = size(quads, 2)
       allocate (fragment(3 + size(nodes) + size(edges) + size(quads) + size(tags) + size(media_types)))
       fragment(1:3) = [real(size(nodes, 2), real64), real(num_edges, real64), real(num_quads, real64)]
       cursor = 4
       if (size(nodes) > 0) then
          fragment(cursor:cursor + size(nodes) - 1) = real(reshape(nodes, [size(nodes)]), real64)
          cursor = cursor + size(nodes)
       end if
       if (size(edges) > 0) then
          fragment(cursor:cursor + size(edges) - 1) = real(reshape(edges, [size(edges)]), real64)
          cursor = cursor + size(edges)
       end if
       if (size(quads) > 0) then
          fragment(cursor:cursor + size(quads) - 1) = real(reshape(quads, [size(quads)]), real64)
          cursor = cursor + size(quads)
       end if
       if (size(tags) > 0) then
          fragment(cursor:cursor + size(tags) - 1) = real(tags, real64)
          cursor = cursor + size(tags)
       end if
       if (size(media_types) > 0) fragment(cursor:) = real(media_types, real64)
    end subroutine serialise_geometry_fragment

    subroutine write_geometry_fragments(writer, fragments, counts, displacements, invalid_fragment, writer_status)
       type(xdmf_writer_t), intent(inout) :: writer
       real(real64), intent(in) :: fragments(:)
       integer, intent(in) :: counts(:), displacements(:)
       logical, intent(out) :: invalid_fragment
       type(xdmf_status_t), intent(inout) :: writer_status
       type(xdmf_grid_id_t) :: grid
       type(xdmf_attribute_id_t) :: tags_attribute, media_attribute
       real(real64), allocatable :: nodes(:, :), tags(:), media_types(:)
       integer(kind=SINGLE), allocatable :: edges(:, :), quads(:, :)
       integer(int64), allocatable :: connectivity(:, :)
       integer :: fragment_index, cursor, num_nodes, num_edges, num_quads, expected_size
       character(len=32) :: grid_name

       invalid_fragment = .false.
       do fragment_index = 1, size(counts)
          if (writer_status%is_error()) return
          if (counts(fragment_index) < 3) then
             invalid_fragment = .true.
             return
          end if
          cursor = displacements(fragment_index) + 1
          num_nodes = nint(fragments(cursor)); num_edges = nint(fragments(cursor + 1)); num_quads = nint(fragments(cursor + 2))
          expected_size = 3 + 3*num_nodes + 2*num_edges + 4*num_quads + 2*(num_edges + num_quads)
          if (num_nodes < 0 .or. num_edges < 0 .or. num_quads < 0 .or. counts(fragment_index) /= expected_size) then
             invalid_fragment = .true.
             return
          end if
          cursor = cursor + 3
          allocate (nodes(3, num_nodes), edges(2, num_edges), quads(4, num_quads), tags(num_edges + num_quads), &
                    media_types(num_edges + num_quads))
          if (num_nodes > 0) then
             nodes = reshape(fragments(cursor:cursor + 3*num_nodes - 1), shape(nodes)); cursor = cursor + 3*num_nodes
          end if
          if (num_edges > 0) then
             edges = reshape(int(nint(fragments(cursor:cursor + 2*num_edges - 1)), kind=SINGLE), shape(edges))
             cursor = cursor + 2*num_edges
          end if
          if (num_quads > 0) then
             quads = reshape(int(nint(fragments(cursor:cursor + 4*num_quads - 1)), kind=SINGLE), shape(quads))
             cursor = cursor + 4*num_quads
          end if
          if (num_edges + num_quads > 0) then
             tags = fragments(cursor:cursor + num_edges + num_quads - 1); cursor = cursor + num_edges + num_quads
             media_types = fragments(cursor:cursor + num_edges + num_quads - 1)
          end if
          if (num_edges > 0) then
             write (grid_name, '("lines_rank_", I0)') fragment_index - 1
             allocate (connectivity(2, num_edges)); connectivity = int(edges, int64) + 1_int64
             call write_geometry_grid(writer, trim(grid_name), XDMF_TOPOLOGY_POLYLINE, nodes, connectivity, tags(:num_edges), &
                                      media_types(:num_edges), grid, tags_attribute, media_attribute, writer_status)
             deallocate (connectivity)
          end if
          if (num_quads > 0 .and. .not. writer_status%is_error()) then
             write (grid_name, '("faces_rank_", I0)') fragment_index - 1
             allocate (connectivity(4, num_quads)); connectivity = int(quads, int64) + 1_int64
             call write_geometry_grid(writer, trim(grid_name), XDMF_TOPOLOGY_QUADRILATERAL, nodes, connectivity, &
                                      tags(num_edges + 1:), media_types(num_edges + 1:), grid, tags_attribute, media_attribute, writer_status)
             deallocate (connectivity)
          end if
          deallocate (nodes, edges, quads, tags, media_types)
       end do
    end subroutine write_geometry_fragments

    subroutine write_geometry_grid(writer, name, topology, nodes, connectivity, tags, media_types, grid, tags_attribute, &
                                   media_attribute, writer_status)
       type(xdmf_writer_t), intent(inout) :: writer
       character(len=*), intent(in) :: name
       integer, intent(in) :: topology
       real(real64), intent(in) :: nodes(:, :), tags(:), media_types(:)
       integer(int64), intent(in) :: connectivity(:, :)
       type(xdmf_grid_id_t), intent(out) :: grid
       type(xdmf_attribute_id_t), intent(out) :: tags_attribute, media_attribute
       type(xdmf_status_t), intent(inout) :: writer_status

       call writer%define_unstructured_grid(name, topology, nodes, connectivity, grid, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'tagnumber', XDMF_CENTER_CELL, &
                                                                         XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, &
                                                                         tags_attribute, writer_status)
       if (.not. writer_status%is_error()) call writer%write_attribute(tags_attribute, tags, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'mediatype', XDMF_CENTER_CELL, &
                                                                         XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, &
                                                                         media_attribute, writer_status)
       if (.not. writer_status%is_error()) call writer%write_attribute(media_attribute, media_types, writer_status)
    end subroutine write_geometry_grid

   subroutine build_cell_properties(this, problemInfo, numEdges, numQuads, tags, media_types)
      type(mapvtk_output_t), intent(in) :: this
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: numEdges, numQuads
      real, allocatable, intent(out) :: tags(:), media_types(:)
      integer :: edge_index, quad_index, index, field

      allocate (tags(numEdges + numQuads), media_types(numEdges + numQuads))
      edge_index = 0
      quad_index = numEdges
      do index = 1, this%nPoints
         select case (this%currentType(index))
         case (iJx, iJy, iJz)
            edge_index = edge_index + 1
            field = electric_field(this%currentType(index))
            tags(edge_index) = real(this%materialTag(index))
            if (this%mediaType(index) >= 0.0_RKIND) then
               media_types(edge_index) = this%mediaType(index)
            else
               media_types(edge_index) = get_output_media_type(field, this%coords(:, index), problemInfo)
            end if
         case (iBloqueJx, iBloqueJy, iBloqueJz)
            quad_index = quad_index + 1
            field = magnetic_field(this%currentType(index))
            tags(quad_index) = real(this%materialTag(index))
            media_types(quad_index) = get_output_media_type(field, this%coords(:, index), problemInfo)
         end select
      end do
   contains
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
