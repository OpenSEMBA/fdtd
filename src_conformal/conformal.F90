module conformal_m

   use geometry_m
   use cell_map_m
   use NFDETypes_m, only: ConformalPECRegions_t, ConformalPECElements_t, ConformalMedia_t, &
                        edge_t, face_t, & 
                        conformal_face_media_t, conformal_edge_media_t, rkind
   
   real (kind=rkind), PARAMETER :: EDGE_RATIO_EQ_TOLERANCE = 1e-5
   real (kind=rkind), PARAMETER :: FACE_RATIO_EQ_TOLERANCE = 1e-3

contains

   function buildSideMaps(elements) result(res)
      type(ConformalPECElements_t), dimension(:), pointer, intent(in) :: elements
      type(side_tris_map_t), dimension(:), allocatable :: res
      integer :: i
      allocate(res(size(elements)))
      do i = 1, size(elements)
         call buildSideMap(res(i), elements(i)%triangles)
      end do
   end function

   subroutine buildConformalMedia(conformalRegs, volumes, surfaces) 
      type(ConformalPECRegions_t), pointer, intent(in) :: conformalRegs
      type(ConformalMedia_t), allocatable, dimension(:), intent(inout) :: volumes, surfaces
      if (associated(conformalRegs%volumes)) then 
         volumes = buildMedia(conformalRegs%volumes)
      else
         allocate(volumes(0))
      end if
      if (associated(conformalRegs%surfaces)) then 
         surfaces = buildSurfaceMedia(conformalRegs%surfaces)
      else
         allocate(surfaces(0))
      end if
   end subroutine

   function buildSurfaceMedia(elements) result(res)
      type(ConformalPECElements_t), dimension(:), pointer, intent(in) :: elements
      type(ConformalPECElements_t) :: canonical_element
      type(ConformalMedia_t), dimension(:), allocatable :: res
      integer :: i

      allocate(res(size(elements)))
      do i = 1, size(elements)
         canonical_element = canonicalClosedSurfaceOrientation(elements(i))
         res(i) = buildMediaFromElement(canonical_element)
      end do
   end function buildSurfaceMedia

   subroutine validateConformalSurface(element, is_valid, message)
      type(ConformalPECElements_t), intent(in) :: element
      logical, intent(out) :: is_valid
      character(len=*), intent(out) :: message

      call validateConformalRegion(element, .false., is_valid, message)
   end subroutine validateConformalSurface

   subroutine validateConformalVolume(element, is_valid, message)
      type(ConformalPECElements_t), intent(in) :: element
      logical, intent(out) :: is_valid
      character(len=*), intent(out) :: message

      call validateConformalRegion(element, .true., is_valid, message)
   end subroutine validateConformalVolume

   subroutine validateConformalRegion(element, require_closed, is_valid, message)
      type(ConformalPECElements_t), intent(in) :: element
      logical, intent(in) :: require_closed
      logical, intent(out) :: is_valid
      character(len=*), intent(out) :: message
      type(cell_map_t) :: cell_map
      type(side_t), dimension(:), allocatable :: sides, sides_on_face
      integer(kind=4), dimension(3) :: cell
      integer :: i, face
      character(len=160) :: detail

      is_valid = .true.
      message = ''
      call validateSurfaceTopology(element%triangles, is_valid, message)
      if (.not. is_valid) return
      if (require_closed) then
         call validateClosedSurfaceTopology(element%triangles, is_valid, message)
         if (.not. is_valid) return
      end if

      call buildCellMap(cell_map, element)
      do i = 1, size(cell_map%keys)
         cell = cell_map%keys(i)%cell
         sides = cell_map%getSidesInCell(cell)
         do face = FACE_X, FACE_Z
            sides_on_face = getSidesOnFace(sides, face)
            if (size(sides_on_face) > 0) then
               call validateSurfaceFaceContour(sides_on_face, face, is_valid, detail)
               if (.not. is_valid) then
                  call setCellFaceValidationError(trim(detail), cell, face, is_valid, message)
                  return
               end if
            end if
         end do
      end do
   end subroutine validateConformalRegion

   function canonicalClosedSurfaceOrientation(element) result(res)
      type(ConformalPECElements_t), intent(in) :: element
      type(ConformalPECElements_t) :: res
      integer, dimension(:), allocatable :: component
      integer :: component_id, triangle_index
      real :: signed_volume, volume_scale, volume_tolerance
      real, dimension(3) :: min_position, max_position, extent

      res = element
      if (.not. allocated(res%triangles)) return
      if (size(res%triangles) == 0) return
      call assignTriangleComponents(res%triangles, component)

      do component_id = 1, maxval(component)
         if (.not. isClosedComponent(res%triangles, component, component_id)) cycle
         signed_volume = 0.0
         min_position = huge(1.0)
         max_position = -huge(1.0)
         do triangle_index = 1, size(res%triangles)
            if (component(triangle_index) /= component_id) cycle
            signed_volume = signed_volume + dot_product(res%triangles(triangle_index)%vertices(1)%position, &
               cross(res%triangles(triangle_index)%vertices(2)%position, res%triangles(triangle_index)%vertices(3)%position))/6.0
            min_position = min(min_position, res%triangles(triangle_index)%vertices(1)%position)
            min_position = min(min_position, res%triangles(triangle_index)%vertices(2)%position)
            min_position = min(min_position, res%triangles(triangle_index)%vertices(3)%position)
            max_position = max(max_position, res%triangles(triangle_index)%vertices(1)%position)
            max_position = max(max_position, res%triangles(triangle_index)%vertices(2)%position)
            max_position = max(max_position, res%triangles(triangle_index)%vertices(3)%position)
         end do
         extent = max_position - min_position
         volume_scale = max(1.0, extent(1)*extent(2)*extent(3))
         volume_tolerance = 100.0*epsilon(1.0)*volume_scale
         if (signed_volume < -volume_tolerance) then
            do triangle_index = 1, size(res%triangles)
               if (component(triangle_index) == component_id) call reverseTriangle(res%triangles(triangle_index))
            end do
         end if
      end do
   end function canonicalClosedSurfaceOrientation

   subroutine validateSurfaceTopology(triangles, is_valid, message)
      type(triangle_t), dimension(:), allocatable, intent(in) :: triangles
      logical, intent(out) :: is_valid
      character(len=*), intent(out) :: message
      integer :: triangle_index, edge, other_triangle, other_edge
      integer :: first_id, second_id, direction, other_first, other_second, other_direction, incidence

      is_valid = .true.
      message = ''
      do triangle_index = 1, size(triangles)
         do edge = 1, 3
            call triangleEdge(triangles(triangle_index), edge, first_id, second_id, direction)
            if (first_id == second_id) then
               write(message, '(A,I0,A,I0)') 'triangle ', triangle_index, ' has a degenerate edge ', edge
               is_valid = .false.
               return
            end if
            incidence = 0
            do other_triangle = 1, size(triangles)
               do other_edge = 1, 3
                  call triangleEdge(triangles(other_triangle), other_edge, other_first, other_second, other_direction)
                  if (first_id == other_first .and. second_id == other_second) then
                     incidence = incidence + 1
                     if (other_triangle /= triangle_index .and. direction == other_direction) then
                        write(message, '(A,I0,A,I0,A,I0,A,I0,A)') 'triangles ', triangle_index, ' and ', other_triangle, &
                           ' traverse shared edge ', first_id, '-', second_id, ' in the same direction'
                        is_valid = .false.
                        return
                     end if
                  end if
               end do
            end do
            if (incidence > 2) then
               write(message, '(A,I0,A,I0,A,I0)') 'non-manifold edge ', first_id, '-', second_id, ' has incidence ', incidence
               is_valid = .false.
               return
            end if
         end do
      end do
   end subroutine validateSurfaceTopology

   subroutine validateClosedSurfaceTopology(triangles, is_valid, message)
      type(triangle_t), dimension(:), allocatable, intent(in) :: triangles
      logical, intent(out) :: is_valid
      character(len=*), intent(out) :: message
      integer :: triangle_index, edge, other_triangle, other_edge, incidence
      integer :: first_id, second_id, direction, other_first, other_second, other_direction

      is_valid = .true.
      message = ''
      if (size(triangles) == 0) then
         message = 'volume has no triangles'
         is_valid = .false.
         return
      end if
      do triangle_index = 1, size(triangles)
         do edge = 1, 3
            call triangleEdge(triangles(triangle_index), edge, first_id, second_id, direction)
            incidence = 0
            do other_triangle = 1, size(triangles)
               do other_edge = 1, 3
                  call triangleEdge(triangles(other_triangle), other_edge, other_first, other_second, other_direction)
                  if (first_id == other_first .and. second_id == other_second) incidence = incidence + 1
               end do
            end do
            if (incidence /= 2) then
               write(message, '(A,I0,A,I0,A,I0)') 'open volume edge ', first_id, '-', second_id, ' has incidence ', incidence
               is_valid = .false.
               return
            end if
         end do
      end do
   end subroutine validateClosedSurfaceTopology

   subroutine validateSurfaceFaceContour(sides, face, is_valid, detail)
      type(side_t), dimension(:), allocatable, intent(in) :: sides
      integer, intent(in) :: face
      logical, intent(out) :: is_valid
      character(len=*), intent(out) :: detail
      type(coord_t), dimension(:), allocatable :: vertices
      integer, dimension(:), allocatable :: incoming, outgoing, queue
      logical, dimension(:), allocatable :: visited
      integer :: side_index, vertex_index, next_vertex, nvertices, head, tail, components

      is_valid = .true.
      detail = ''
      allocate(vertices(2*size(sides)), incoming(2*size(sides)), outgoing(2*size(sides)))
      incoming = 0
      outgoing = 0
      nvertices = 0
      do side_index = 1, size(sides)
         call registerContourVertex(sides(side_index)%init, vertices, incoming, outgoing, nvertices, vertex_index)
         outgoing(vertex_index) = outgoing(vertex_index) + 1
         call registerContourVertex(sides(side_index)%end, vertices, incoming, outgoing, nvertices, vertex_index)
         incoming(vertex_index) = incoming(vertex_index) + 1
      end do
      do vertex_index = 1, nvertices
         if (incoming(vertex_index) > 1 .or. outgoing(vertex_index) > 1) then
            detail = 'branching contour vertex'
            is_valid = .false.
            return
         end if
         if (.not. isOnFaceBoundary(vertices(vertex_index), face) .and. &
             (incoming(vertex_index) /= 1 .or. outgoing(vertex_index) /= 1)) then
            detail = 'dangling contour vertex inside the cell face'
            is_valid = .false.
            return
         end if
      end do

      allocate(visited(size(sides)), queue(size(sides)))
      visited = .false.
      components = 0
      do side_index = 1, size(sides)
         if (visited(side_index)) cycle
         components = components + 1
         head = 1
         tail = 1
         queue(1) = side_index
         visited(side_index) = .true.
         do while (head <= tail)
            next_vertex = queue(head)
            head = head + 1
            do vertex_index = 1, size(sides)
               if (visited(vertex_index)) cycle
               if (sidesTouch(sides(next_vertex), sides(vertex_index))) then
                  tail = tail + 1
                  queue(tail) = vertex_index
                  visited(vertex_index) = .true.
               end if
            end do
         end do
      end do
      if (components /= 1) then
         detail = 'multiple disconnected contours'
         is_valid = .false.
      end if
   end subroutine validateSurfaceFaceContour

   subroutine registerContourVertex(vertex, vertices, incoming, outgoing, nvertices, index)
      type(coord_t), intent(in) :: vertex
      type(coord_t), dimension(:), intent(inout) :: vertices
      integer, dimension(:), intent(inout) :: incoming, outgoing
      integer, intent(inout) :: nvertices
      integer, intent(out) :: index
      integer :: i

      do i = 1, nvertices
         if (sameCoordinate(vertices(i), vertex)) then
            index = i
            return
         end if
      end do
      nvertices = nvertices + 1
      vertices(nvertices) = vertex
      incoming(nvertices) = 0
      outgoing(nvertices) = 0
      index = nvertices
   end subroutine registerContourVertex

   logical function isOnFaceBoundary(vertex, face)
      type(coord_t), intent(in) :: vertex
      integer, intent(in) :: face
      integer :: first_direction, second_direction

      first_direction = mod(face, 3) + 1
      second_direction = mod(face + 1, 3) + 1
      isOnFaceBoundary = isOnGridPlane(vertex%position(first_direction)) .or. &
                         isOnGridPlane(vertex%position(second_direction))
   end function isOnFaceBoundary

   logical function isOnGridPlane(position)
      real, intent(in) :: position
      isOnGridPlane = abs(position - nint(position)) < 1.0e-5
   end function isOnGridPlane

   logical function sidesTouch(first, second)
      type(side_t), intent(in) :: first, second
      sidesTouch = sameCoordinate(first%init, second%init) .or. sameCoordinate(first%init, second%end) .or. &
                   sameCoordinate(first%end, second%init) .or. sameCoordinate(first%end, second%end)
   end function sidesTouch

   logical function sameCoordinate(first, second)
      type(coord_t), intent(in) :: first, second
      sameCoordinate = all(abs(first%position-second%position) < 1.0e-5)
   end function sameCoordinate

   subroutine setCellFaceValidationError(detail, cell, face, is_valid, message)
      character(len=*), intent(in) :: detail
      integer(kind=4), dimension(3), intent(in) :: cell
      integer, intent(in) :: face
      logical, intent(out) :: is_valid
      character(len=*), intent(out) :: message
      write(message, '(A,A,I0,A,I0,A,I0,A,I0)') trim(detail), ' at cell (', cell(1), ',', cell(2), ',', cell(3), '), face ', face
      is_valid = .false.
   end subroutine setCellFaceValidationError

   subroutine assignTriangleComponents(triangles, component)
      type(triangle_t), dimension(:), allocatable, intent(in) :: triangles
      integer, dimension(:), allocatable, intent(out) :: component
      integer, dimension(:), allocatable :: queue
      integer :: component_id, seed, head, tail, current, neighbor

      allocate(component(size(triangles)), queue(size(triangles)))
      component = 0
      component_id = 0
      do seed = 1, size(triangles)
         if (component(seed) /= 0) cycle
         component_id = component_id + 1
         head = 1
         tail = 1
         queue(1) = seed
         component(seed) = component_id
         do while (head <= tail)
            current = queue(head)
            head = head + 1
            do neighbor = 1, size(triangles)
               if (component(neighbor) /= 0) cycle
               if (trianglesShareEdge(triangles(current), triangles(neighbor))) then
                  component(neighbor) = component_id
                  tail = tail + 1
                  queue(tail) = neighbor
               end if
            end do
         end do
      end do
   end subroutine assignTriangleComponents

   logical function isClosedComponent(triangles, component, component_id)
      type(triangle_t), dimension(:), allocatable, intent(in) :: triangles
      integer, dimension(:), intent(in) :: component
      integer, intent(in) :: component_id
      integer :: triangle_index, edge, other_triangle, other_edge, incidence
      integer :: first_id, second_id, direction, other_first, other_second, other_direction

      isClosedComponent = .true.
      do triangle_index = 1, size(triangles)
         if (component(triangle_index) /= component_id) cycle
         do edge = 1, 3
            call triangleEdge(triangles(triangle_index), edge, first_id, second_id, direction)
            incidence = 0
            do other_triangle = 1, size(triangles)
               if (component(other_triangle) /= component_id) cycle
               do other_edge = 1, 3
                  call triangleEdge(triangles(other_triangle), other_edge, other_first, other_second, other_direction)
                  if (first_id == other_first .and. second_id == other_second) incidence = incidence + 1
               end do
            end do
            if (incidence /= 2) then
               isClosedComponent = .false.
               return
            end if
         end do
      end do
   end function isClosedComponent

   logical function trianglesShareEdge(first, second)
      type(triangle_t), intent(in) :: first, second
      integer :: edge, other_edge, first_id, second_id, direction, other_first, other_second, other_direction

      trianglesShareEdge = .false.
      do edge = 1, 3
         call triangleEdge(first, edge, first_id, second_id, direction)
         do other_edge = 1, 3
            call triangleEdge(second, other_edge, other_first, other_second, other_direction)
            if (first_id == other_first .and. second_id == other_second) then
               trianglesShareEdge = .true.
               return
            end if
         end do
      end do
   end function trianglesShareEdge

   subroutine triangleEdge(triangle, edge_index, first_id, second_id, direction)
      type(triangle_t), intent(in) :: triangle
      integer, intent(in) :: edge_index
      integer, intent(out) :: first_id, second_id, direction
      integer :: directed_first, directed_second

      directed_first = triangle%vertices(edge_index)%id
      directed_second = triangle%vertices(mod(edge_index,3)+1)%id
      first_id = min(directed_first, directed_second)
      second_id = max(directed_first, directed_second)
      if (directed_first == first_id) then
         direction = 1
      else
         direction = -1
      end if
   end subroutine triangleEdge

   subroutine reverseTriangle(triangle)
      type(triangle_t), intent(inout) :: triangle
      type(coord_t) :: swapped_vertex

      swapped_vertex = triangle%vertices(2)
      triangle%vertices(2) = triangle%vertices(3)
      triangle%vertices(3) = swapped_vertex
   end subroutine reverseTriangle

   function buildMedia(elements) result(res)
      type(ConformalPECElements_t), dimension(:), pointer :: elements
      type(ConformalMedia_t), dimension(:), allocatable :: res
      integer :: i
      allocate(res(size(elements)))
      do i = 1, size(elements)
         res(i) = buildMediaFromElement(elements(i))
      end do
   end function   

   function buildMediaFromElement(element) result(res)
      type(ConformalPECElements_t), intent(in) :: element
      type(ConformalMedia_t) :: res

      type(cell_map_t) :: cell_map
      real (kind=rkind), dimension(:), allocatable :: edge_ratios, face_ratios
      type(edge_t), dimension(:), allocatable :: edges
      type(face_t), dimension(:), allocatable :: faces
      
      call buildCellMap(cell_map, element)
      call fillElements(cell_map, faces, edges)
      call addNewRatios(edges, faces, edge_ratios, face_ratios)
      res%edge_media => addEdgeMedia(edges, edge_ratios)
      res%face_media => addFaceMedia(faces, face_ratios)

      res%n_edges_media = size(res%edge_media)
      res%n_faces_media = size(res%face_media)

      res%time_step_scale_factor = computeTimeStepScalingFactor(res%edge_media, res%face_media)
      res%tag = element%tag
   end function

   subroutine addNewRatios(edges, faces, edge_ratios, face_ratios)
      type(edge_t), dimension(:), allocatable, intent(in) :: edges
      type(face_t), dimension(:), allocatable, intent(in) :: faces
      real (kind=rkind), dimension(:), allocatable, intent(inout) :: edge_ratios, face_ratios
      integer :: i
      allocate(edge_ratios(0))
      allocate(face_ratios(0))
      do i = 1, size(edges)
         if (isNewRatio(edge_ratios, edges(i)%ratio, EDGE_RATIO_EQ_TOLERANCE)) then 
            call addRatio(edge_ratios, edges(i)%ratio)
         end if
      end do
      do i = 1, size(faces)
         if (isNewRatio(face_ratios, faces(i)%ratio, FACE_RATIO_EQ_TOLERANCE)) then 
            call addRatio(face_ratios, faces(i)%ratio)
         end if
      end do

   end subroutine


   function addEdgeMedia(edges, edge_ratios) result(res)
      real (kind=rkind), dimension(:), allocatable, intent(in) :: edge_ratios
      type(edge_t), dimension(:), allocatable, intent(in) :: edges
      type(edge_t), dimension(:), allocatable :: filtered_edges
      TYPE (conformal_edge_media_t), DIMENSION (:), POINTER :: res
      allocate(res(size(edge_ratios)))
      do i = 1, size(edge_ratios)
         filtered_edges = filterEdgesByMedia(edges, edge_ratios(i))
         allocate(res(i)%edges(size(filtered_edges)))
         res(i)%edges = filtered_edges
         res(i)%ratio = edge_ratios(i)
         res(i)%size = size(filtered_edges)
      end do
   end function

   function addFaceMedia(faces, face_ratios) result(res)
      real (kind=rkind), dimension(:), allocatable, intent(in) :: face_ratios
      type(face_t), dimension(:), allocatable, intent(in) :: faces
      type(face_t), dimension(:), allocatable :: filtered_faces
      TYPE (conformal_face_media_t), DIMENSION (:), POINTER :: res
      allocate(res(size(face_ratios)))
      do i = 1, size(face_ratios)
         filtered_faces = filterFacesByMedia(faces, face_ratios(i))
         allocate(res(i)%faces(size(filtered_faces)))
         res(i)%faces = filtered_faces
         res(i)%ratio = face_ratios(i)
         res(i)%size = size(filtered_faces)
      end do

   end function

   function computeTimeStepScalingFactor(edges_media, faces_media) result(res)
      TYPE (conformal_face_media_t), dimension(:), intent(in), pointer :: faces_media
      TYPE (conformal_edge_media_t), dimension(:), intent(in), pointer :: edges_media
      real (kind=rkind) :: res, l_ratio, area
      type (cell_ratios_map_t) :: cell_ratio_map
      type (cell_ratios_t) :: cell_ratio_info
      integer (kind=4), dimension(3) :: cell, aux_cell
      integer :: idx1, idx2
      integer :: i,j
      res = 1.0
      if (.not. allocated(cell_ratio_map%keys)) allocate(cell_ratio_map%keys(0))
      do i = 1, size(faces_media)
         do j = 1, size(faces_media(i)%faces)
            call cell_ratio_map%addFaceRatio(faces_media(i)%faces(j)%cell, faces_media(i)%faces(j)%direction, faces_media(i)%faces(j)%ratio)
         end do
      end do
      do i = 1, size(edges_media)
         do j = 1, size(edges_media(i)%edges)
            call cell_ratio_map%addEdgeRatio(edges_media(i)%edges(j)%cell, edges_media(i)%edges(j)%direction, edges_media(i)%edges(j)%ratio)
         end do
      end do
      l_ratio = 0.0
      do i = 1, size(cell_ratio_map%keys)
         cell = cell_ratio_map%keys(i)%cell
         aux_cell = cell
         cell_ratio_info = cell_ratio_map%getCellRatiosInCell(cell)
         do j = FACE_X, FACE_Z
            area = cell_ratio_info%area(j)
            idx1 = mod(j,3) + 1
            idx2 = mod(j + 1,3) + 1
            l_ratio = max(cell_ratio_info%length(idx1),cell_ratio_info%length(idx2))
            aux_cell(idx1) = aux_cell(idx1) + 1
            if (cell_ratio_map%hasKey(aux_cell)) then 
               cell_ratio_info = cell_ratio_map%getCellRatiosInCell(aux_cell)
               l_ratio = max(l_ratio, cell_ratio_info%length(idx2))
            else 
               l_ratio = 1.0
            end if
            aux_cell = cell
            aux_cell(idx2) = aux_cell(idx2) + 1
            if (cell_ratio_map%hasKey(aux_cell)) then 
               cell_ratio_info = cell_ratio_map%getCellRatiosInCell(aux_cell)
               l_ratio = max(l_ratio, cell_ratio_info%length(idx1))
            else
               l_ratio = 1.0
            end if
            if (area /= 0.0 .and. l_ratio /= 0.0) then
               res = min(res,sqrt(area/l_ratio))
            end if
         end do
      end do
   end function   


   subroutine fillElements(cell_map, faces, edges)
      type(cell_map_t), intent(in) :: cell_map
      type (face_t), dimension (:), allocatable, intent(inout) :: faces
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      integer (kind=4), dimension(3) :: cell
      integer :: i, edge, face
      type(side_t), dimension(:), allocatable :: sides, sides_on_face, contour, sides_on_edge
      type(triangle_t), dimension(:), allocatable :: tris, tris_on_face
      type(interval_t), dimension(:), allocatable :: intervals
      allocate(faces(0))
      allocate(edges(0))
      do i = 1, size(cell_map%keys)
         cell = cell_map%keys(i)%cell 
         sides = cell_map%getSidesInCell(cell)
         tris =  cell_map%getTrianglesInCell(cell)
         intervals = cell_map%getIntervalsInCell(cell)
         do face = FACE_X, FACE_Z
            sides_on_face = getSidesOnFace(sides, face)
            if (size(sides_on_face) /= 0) then 
               contour = findLargestContour(sides_on_face)
               call fillFaceFromContour(contour, faces)
               call fillEdgesFromContour(contour, edges)
            end if
            tris_on_face = getTrianglesOnFace(tris, face)
            if (size(tris_on_face) /= 0) call fillFullFaces(tris_on_face, faces, edges)
         end do
         call fillIntervals(intervals,edges, faces)
      end do

      do i = 1, size(cell_map%keys)
         cell = cell_map%keys(i)%cell 
         sides = cell_map%getSidesInCell(cell)
         
         do edge = EDGE_X, EDGE_Z
            sides_on_edge = getSidesOnEdge(sides, edge)
            call fillEdges(sides_on_edge, edges)
         end do

         sides = cell_map%getOnSidesInCell(cell)
         do edge = EDGE_X, EDGE_Z
            sides_on_edge = getSidesOnEdge(sides, edge)
            call fillEdges(sides_on_edge, edges)
         end do
      end do
   end subroutine

   function buildSidesFromCellInterval(interval) result(res)
      type(interval_t) :: interval
      integer :: face
      type(side_t), dimension(4) :: res
      type(side_t) :: aux
      type(coord_t), dimension(4) :: cs
      aux%init%position = interval%ini%cell
      aux%end%position = interval%end%cell
      face = aux%getFace()
      cs(1)%position = aux%getCell()
      select case(face)
      case(FACE_X)
         cs(2)%position = cs(1)%position + [0,1,0]
         cs(3)%position = cs(1)%position + [0,1,1]
         cs(4)%position = cs(1)%position + [0,0,1]
      case(FACE_Y)
         cs(2)%position = cs(1)%position + [0,0,1]
         cs(3)%position = cs(1)%position + [1,0,1]
         cs(4)%position = cs(1)%position + [1,0,0]
      case(FACE_Z)
         cs(2)%position = cs(1)%position + [1,0,0]
         cs(3)%position = cs(1)%position + [1,1,0]
         cs(4)%position = cs(1)%position + [0,1,0]
      end select
      res(1)%init = cs(1)
      res(1)%end  = cs(2)
      res(2)%init = cs(2)
      res(2)%end =  cs(3)
      res(3)%init = cs(3)
      res(3)%end  = cs(4)
      res(4)%init = cs(4)
      res(4)%end  = cs(1)
   end function

   subroutine fillIntervals(intervals, edges, faces)
      type(interval_t), dimension(:), allocatable :: intervals
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      type (face_t), dimension (:), allocatable, intent(inout) :: faces
      integer :: i
      type(side_t), dimension(:), allocatable :: contour
      do i = 1, size(intervals) 
         call fillEdgesFromInterval(edges, intervals(i))
         call fillFaceFromInterval(faces, intervals(i))
      end do
   end subroutine

   subroutine fillFullFaces(tris_on_face, faces, edges)
      type(triangle_t), dimension(:), allocatable, intent(in) :: tris_on_face
      type (face_t), dimension (:), allocatable, intent(inout) :: faces
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      type(side_t), dimension(:), allocatable :: tri_sides
      integer :: j, k, s
      real (kind=rkind) :: area, ratio
      integer :: edge, face
      integer, dimension(3) :: cell
      area = 0.0
      ratio = 0.0
      do j = 1, size(tris_on_face)
         area = area + getArea(tris_on_face(j))
      end do
      if (abs(area-1.0) < 1e-4) then 
         cell = tris_on_face(1)%getCell()
         face = tris_on_face(1)%getFace()
         call addFace(faces, cell, face, ratio)
         do k = 1, size(tris_on_face)
            tri_sides = tris_on_face(k)%getSides()
            do s = 1, 3
               if (tri_sides(s)%isOnAnyEdge()) then 
                  cell = tri_sides(s)%getCell()
                  edge = tri_sides(s)%getEdge()
                  if (isNewEdge(edges, cell, edge, ratio)) then 
                     call addEdge(edges, cell, edge, tri_sides(s))
                  end if
               end if
            end do
         end do
      end if

   end subroutine

   subroutine fillEdgesFromContour(contour, edges)
      type(side_t), dimension(:), allocatable, intent(in) :: contour
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      integer :: i, edge
      integer, dimension(3) :: cell
      do i = 1, size(contour)
         edge = contour(i)%getEdge()
         cell = contour(i)%getCell()
         if (edge /= NOT_ON_EDGE) then 
            if (isEdgeFilled(edges, cell, edge)) then 
               call fillSmallerRatio(edges, cell, edge, contour(i))
            else 
               call addEdge(edges, cell, edge, contour(i))
            end if
         end if
      end do
   end subroutine

   subroutine fillEdges(sides, edges)
      type(side_t), dimension(:), allocatable, intent(in) :: sides
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      integer :: i, edge
      integer, dimension(3) :: cell
      do i = 1, size(sides)
         edge = sides(i)%getEdge()
         cell = sides(i)%getCell()
         if (edge /= NOT_ON_EDGE) then 
            if (isEdgeFilled(edges, cell, edge)) then 
               call reduceEdgeRatio(edges, cell, edge, sides(i))
            else 
               call addEdge(edges, cell, edge, sides(i))
            end if
         end if
      end do
   end subroutine

   subroutine fillEdgesFromInterval(edges, interval)
      type (edge_t), dimension (:), allocatable :: edges
      type(interval_t), intent(in) :: interval
      integer :: i, edge
      type(side_t), dimension(4) :: sides
      sides = buildSidesFromCellInterval(interval)
      do i = 1, size(sides)
         edge = sides(i)%getEdge()
         if (edge /= NOT_ON_EDGE) then 
            if (isEdgeFilled(edges, sides(i)%getCell(), edge)) then 
               call fillSmallerRatio(edges, sides(i)%getCell(), edge, sides(i))
            else 
               call addEdge(edges, sides(i)%getCell(), edge, sides(i))
            end if
         end if
      end do
   end subroutine

   subroutine fillFaceFromInterval(faces, interval)
      type (face_t), dimension (:), allocatable :: faces
      type(interval_t), intent(in) :: interval
      type(side_t) :: aux
      real (kind=rkind) :: ratio
      ratio = 0.0
      aux%init%position = interval%ini%cell
      aux%end%position = interval%end%cell
      call addFace(faces, aux%getCell(), aux%getFace(), ratio)
   end subroutine

   subroutine fillFaceFromContour(contour, faces)
      type(side_t), dimension(:), allocatable, intent(in) :: contour
      type (face_t), dimension (:), allocatable :: faces
      real (kind=rkind) :: area
      integer :: face
      integer, dimension(3) :: cell
      cell = findContourCell(contour)
      face = findContourFace(contour)
      if (size(contour) /= 0) then 
         area = 1.0 - contourArea(contour)
         call addFace(faces, cell, face , area)
      end if
   end subroutine

   function findLargestContour(sides) result(res)
      type(side_t), dimension(:), allocatable, intent(in) :: sides
      type(side_t), dimension(:), allocatable :: res
      type(side_t), dimension(:), allocatable :: aux_contour
      type(side_t), dimension(:), allocatable :: aux_side
      integer :: i 
      real :: area, contour_area
      area = 0
      allocate(aux_side(1))
      do i = 1, size(sides)
         aux_side(1) = sides(i)
         aux_contour = buildSidesContour(aux_side)
         contour_area = contourArea(aux_contour)
         if (contour_area > area) then 
            res = aux_contour
            area = contour_area
         end if
      end do

   end function
    
   logical function isNewEdge(edges, cell, edge, ratio)
      type(edge_t), dimension(:), allocatable, intent(inout) :: edges
      integer (kind=4), dimension(3), intent(in) :: cell
      integer (kind=4) :: edge 
      real (kind=rkind) :: ratio
      integer :: i
      isNewEdge = .true.
      do i = 1, size(edges)
         if (all(edges(i)%cell == cell) .and. &
             edges(i)%ratio == ratio .and. &
             edges(i)%direction == edge) then 
               isNewEdge = .false.
               exit
         end if
      end do
   end function

   logical function isEdgeFilled(edges, cell, edge)
      type(edge_t), dimension(:), allocatable, intent(inout) :: edges
      integer (kind=4), dimension(3), intent(in) :: cell
      integer (kind=4) :: edge 
      integer :: i
      isEdgeFilled = .false.
      do i = 1, size(edges)
         if (all(edges(i)%cell == cell) .and. &
             edges(i)%direction == edge) then 
               isEdgeFilled = .true.
               exit
         end if
      end do
   end function

   subroutine reduceEdgeRatio (edges, cell, edge, side)
      type(edge_t), dimension(:), allocatable, intent(inout) :: edges
      integer (kind=4), dimension(3), intent(in) :: cell
      integer (kind=4) :: edge 
      type(side_t), intent(in) :: side
      integer :: i
      do i = 1, size(edges)
         if (all(edges(i)%cell == cell) .and. &
             edges(i)%direction == edge) then 
               if (edges(i)%material_coords(1) /= min(side%init%position(edge), side%end%position(edge)) .and. &
                   edges(i)%material_coords(2) /= max(side%init%position(edge), side%end%position(edge)) .and. &
                   edges(i)%ratio /= 0) then 
                   edges(i)%ratio = edges(i)%ratio - side%length()
               end if
               exit
         end if
      end do
   end subroutine

   subroutine fillSmallerRatio (edges, cell, edge, side)
      type(edge_t), dimension(:), allocatable, intent(inout) :: edges
      integer (kind=4), dimension(3), intent(in) :: cell
      integer (kind=4) :: edge 
      type(side_t), intent(in) :: side
      integer :: i
      real(kind=rkind) :: new_ratio
      do i = 1, size(edges)
         if (all(edges(i)%cell == cell) .and. &
             edges(i)%direction == edge) then 
               new_ratio = 1.0 - side%length()
               if (new_ratio < edges(i)%ratio) then 
                   edges(i)%ratio = new_ratio
               end if
               exit
         end if
      end do
   end subroutine

   subroutine addEdge(edges, cell, edge, side)
      type(edge_t), dimension(:), allocatable, intent(inout) :: edges
      type(edge_t), dimension(:), allocatable :: aux
      integer (kind=4), dimension(3), intent(in) :: cell
      integer (kind=4) :: edge 
      type(side_t), intent(in) :: side
      type(edge_t) :: new_edge
      real (kind=rkind) :: ratio
      real (kind = rkind), dimension(2) :: coords

      ratio = 1.0 - side%length()
      allocate(aux(size(edges) + 1))
      aux(1:size(edges)) = edges
      coords(1) = min(side%init%position(edge), side%end%position(edge))
      coords(2) = max(side%init%position(edge), side%end%position(edge))
      new_edge = edge_t(cell=cell, ratio=ratio, direction=edge, material_coords = coords)
      aux(size(edges) + 1) = new_edge

      deallocate(edges)
      allocate(edges(size(aux)))
      edges = aux
   end subroutine

   subroutine addFace(faces, cell, face, ratio)
      type(face_t), dimension(:), allocatable, intent(inout) :: faces
      type(face_t), dimension(:), allocatable :: aux
      integer (kind=4), dimension(3), intent(in) :: cell
      integer (kind=4) :: face
      type(face_t) :: new_face
      real (kind=rkind) :: ratio
      allocate(aux(size(faces) + 1))
      aux(1:size(faces)) = faces
      new_face = face_t(cell=cell, ratio=ratio, direction=face)
      aux(size(faces) + 1) = new_face

      deallocate(faces)
      allocate(faces(size(aux)))
      faces = aux
   end subroutine

   subroutine addRatio(ratios, ratio)
      real (kind=rkind), dimension(:), allocatable, intent(inout) :: ratios
      real (kind=rkind), dimension(:), allocatable :: aux
      real (kind=rkind) :: ratio
      integer :: i
      logical :: new = .true.
      if (size(ratios) == 0) then 
         deallocate(ratios)
         allocate(ratios(1))
         ratios(1) = ratio
      else
         allocate(aux(size(ratios) + 1))
         aux(1:size(ratios)) = ratios
         aux(size(ratios) + 1) = ratio
         deallocate(ratios)
         allocate(ratios(size(aux)))
         ratios = aux
      end if
   end subroutine

   logical function isNewRatio(ratios, ratio, tol)
      real (kind=rkind), dimension(:), allocatable, intent(in) :: ratios
      real (kind=rkind), intent(in) :: ratio, tol
      integer :: i
      isNewRatio = .true.
      do i = 1, size(ratios)
         if (eq_ratio(ratios(i), ratio, tol)) isNewRatio = .false.
      end do
   end function

   logical function eq_ratio(a,b, tol)
      real (kind=rkind), intent(in) :: a,b, tol
      eq_ratio = abs(a - b) < tol
   end function

   function filterEdgesByMedia(edges, ratio) result(res)
      type(edge_t), dimension(:), allocatable, intent(in) :: edges
      real (kind=rkind) :: ratio
      type(edge_t), dimension(:), allocatable :: res
      integer :: i,n
      n = 0
      do i = 1, size(edges)
         if (eq_ratio(edges(i)%ratio, ratio, EDGE_RATIO_EQ_TOLERANCE)) n = n + 1
      end do
      allocate(res(n))
      n = 0
      do i = 1, size(edges)
         if (eq_ratio(edges(i)%ratio, ratio, EDGE_RATIO_EQ_TOLERANCE)) then 
            n = n + 1
            res(n) = edges(i)
         end if
      end do
   end function

   function filterFacesByMedia(faces, ratio) result(res)
      type(face_t), dimension(:), allocatable, intent(in) :: faces
      real (kind=rkind) :: ratio
      type(face_t), dimension(:), allocatable :: res
      integer :: i,n
      n = 0
      do i = 1, size(faces)
         if (eq_ratio(faces(i)%ratio, ratio, FACE_RATIO_EQ_TOLERANCE)) n = n + 1
      end do
      allocate(res(n))
      n = 0
      do i = 1, size(faces)
         if (eq_ratio(faces(i)%ratio, ratio, FACE_RATIO_EQ_TOLERANCE)) then 
            n = n + 1
            res(n) = faces(i)
         end if
      end do
   end function


end module conformal_m
