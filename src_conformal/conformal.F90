module conformal_mod

   use geometry_mod
   use cell_map_mod
   use NFDETypes, only: ConformalPECRegions, ConformalPECElements, ConformalMedia_t, & 
                        edge_t, face_t, & 
                        conformal_face_media_t, conformal_edge_media_t, rkind
   
   real (kind=rkind), PARAMETER :: EDGE_RATIO_EQ_TOLERANCE = 1e-5
   real (kind=rkind), PARAMETER :: FACE_RATIO_EQ_TOLERANCE = 1e-3

contains

   function buildSideMaps(regions) result(res)
      type(ConformalPECRegions), intent(in) :: regions
      type(side_tris_map_t), dimension(:), allocatable :: res
      integer :: i
      allocate(res(size(regions%volumes)))
      do i = 1, size(regions%volumes)
         call buildSideMap(res(i),regions%volumes(i)%triangles)
      end do
   end function

   function buildConformalMedia(regions) result(res)
      type(ConformalPECRegions), intent(in) :: regions
      type(ConformalMedia_t), dimension(:), allocatable :: res
      integer :: i
      allocate(res(size(regions%volumes)))
      do i = 1, size(regions%volumes)
         res(i) = buildConformalVolume(regions%volumes(i))
      end do
   end function   


   function buildConformalVolume(volume) result(res)
      type(ConformalPECElements), intent(in) :: volume
      type(ConformalMedia_t) :: res

      type(cell_map_t) :: cell_map
      real (kind=rkind), dimension(:), allocatable :: edge_ratios, face_ratios
      type(edge_t), dimension(:), allocatable :: edges, filtered_edges
      type(face_t), dimension(:), allocatable :: faces, filtered_faces
      integer :: i
      
      call buildCellMap(cell_map, volume)
      call fillElements(cell_map, faces, edges)
      call addNewRatios(edges, faces, edge_ratios, face_ratios)
      res%edge_media => addEdgeMedia(edges, edge_ratios)
      res%face_media => addFaceMedia(faces, face_ratios)

      res%n_edges_media = size(res%edge_media)
      res%n_faces_media = size(res%face_media)

      res%time_step_scale_factor = computeTimeStepScalingFactor(res%edge_media, res%face_media)
      res%tag = volume%tag
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
      allocate(faces(0))
      allocate(edges(0))
      block
         type(side_t), dimension(:), allocatable :: sides, sides_on_face, contour, sides_on_edge
         type(triangle_t), dimension(:), allocatable :: tris
         type(cell_interval_t), dimension(:), allocatable :: intervals
         do i = 1, size(cell_map%keys)
            cell = cell_map%keys(i)%cell 
            sides = cell_map%getSidesInCell(cell)
            tris =  cell_map%getTrianglesInCell(cell)
            intervals = cell_map%getIntervalsInCell(cell)
            do face = FACE_X, FACE_Z
               sides_on_face = getSidesOnFace(sides, face)
               contour = findLargestContour(sides_on_face)
               call fillFaceFromContour(contour, faces)
               call fillEdgesFromContour(contour, edges)
               call fillFullFaces(getTrianglesOnFace(tris, face), faces, edges)
            end do
            call fillIntervals(intervals,edges, faces)
         end do
      end block

      block
         type(side_t), dimension(:), allocatable :: sides, sides_on_edge
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
      end block
   end subroutine

   function buildSidesFromCellInterval(interval) result(res)
      type(cell_interval_t) :: interval
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
      type(cell_interval_t), dimension(:), allocatable :: intervals
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
      real :: area
      area = 0.0
      do j = 1, size(tris_on_face)
         area = area + getArea(tris_on_face(j))
      end do
      if (abs(area-1.0) < 1e-4) then 

         call addFace(faces, tris_on_face(1)%getCell(), tris_on_face(1)%getFace(), 0.0)
         do k = 1, size(tris_on_face)
            tri_sides = tris_on_face(k)%getSides()
            do s = 1, 3
               if (tri_sides(s)%isOnAnyEdge()) then 
                  if (isNewEdge(edges, tri_sides(s)%getCell(), tri_sides(s)%getEdge(), 0.0)) then 
                     call addEdge(edges, tri_sides(s)%getCell(), tri_sides(s)%getEdge(), tri_sides(s))
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
      do i = 1, size(contour)
         edge = contour(i)%getEdge()
         if (edge /= NOT_ON_EDGE) then 
            if (isEdgeFilled(edges, contour(i)%getCell(), edge)) then 
               call fillSmallerRatio(edges, contour(i)%getCell(), edge, contour(i))
            else 
               call addEdge(edges, contour(i)%getCell(), edge, contour(i))
            end if
         end if
      end do
   end subroutine

   subroutine fillEdges(sides, edges)
      type(side_t), dimension(:), allocatable, intent(in) :: sides
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      integer :: i, edge
      do i = 1, size(sides)
         edge = sides(i)%getEdge()
         if (edge /= NOT_ON_EDGE) then 
            if (isEdgeFilled(edges, sides(i)%getCell(), edge)) then 
               call reduceEdgeRatio(edges, sides(i)%getCell(), edge, sides(i))
            else 
               call addEdge(edges, sides(i)%getCell(), edge, sides(i))
            end if
         end if
      end do
   end subroutine

   subroutine fillEdgesFromInterval(edges, interval)
      type (edge_t), dimension (:), allocatable :: edges
      type(cell_interval_t), intent(in) :: interval
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
      type(cell_interval_t), intent(in) :: interval
      type(side_t) :: aux
      aux%init%position = interval%ini%cell
      aux%end%position = interval%end%cell
      call addFace(faces, aux%getCell(), aux%getFace(), 0.0)
   end subroutine

   subroutine fillFaceFromContour(contour, faces)
      type(side_t), dimension(:), allocatable, intent(in) :: contour
      type (face_t), dimension (:), allocatable :: faces
      if (size(contour) /= 0) then 
         call addFace(faces, findContourCell(contour), findContourFace(contour), 1.0 - contourArea(contour))
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


end module