module conformal_mod

    use geometry_mod
    use cell_map_mod
    use NFDETypes, only: ConformalPECRegions, ConformalMedia_t, edge_t, face_t, conformal_face_media_t, conformal_edge_media_t, rkind

    real (kind=rkind), PARAMETER :: RATIO_EQ_TOLERANCE = 0.05

contains

   function buildConformalMedia(regions) result(res)
      type(ConformalPECRegions), intent(in) :: regions
      type(ConformalMedia_t) :: res

      type(cell_map_t) :: cell_map
      real (kind=rkind), dimension(:), allocatable :: edge_ratios, face_ratios
      type(edge_t), dimension(:), allocatable :: edges, filtered_edges
      type(face_t), dimension(:), allocatable :: faces, filtered_faces
      integer :: i
      
      do i = 1, size(regions%volumes)
         call buildCellMap(cell_map, regions%volumes(i)%triangles)
      end do
      call fillConformalEdges(cell_map, edges, edge_ratios)
      call fillConformalFaces(cell_map, faces, face_ratios)

      res%edge_media => addEdgeMedia(edges, edge_ratios)
      res%face_media => addFaceMedia(faces, face_ratios)

      res%n_edges_media = size(res%edge_media)
      res%n_faces_media = size(res%face_media)

      res%cfl = computeCFL(res%edge_media, res%face_media)

   end function

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

   function computeCFL(edges_media, faces_media) result(res)
      TYPE (conformal_face_media_t), dimension(:), intent(in), pointer :: faces_media
      TYPE (conformal_edge_media_t), dimension(:), intent(in), pointer :: edges_media
      real (kind=rkind) :: res, a_ratio, l_ratio 
      type(cfl_map_t) :: cfl_map
      type(cfl_info_t) :: cfl_info
      integer, dimension(3) :: cell, aux_cell
      integer :: idx1, idx2
      integer :: i,j
      res = 1.0
      if (.not. allocated(cfl_map%keys)) allocate(cfl_map%keys(0))
      do i = 1, size(faces_media)
         do j = 1, size(faces_media(i)%faces)
            call cfl_map%addFaceRatio(faces_media(i)%faces(j)%cell, faces_media(i)%faces(j)%direction, faces_media(i)%faces(j)%ratio)
         end do
      end do
      do i = 1, size(edges_media)
         do j = 1, size(edges_media(i)%edges)
            call cfl_map%addEdgeRatio(edges_media(i)%edges(j)%cell, edges_media(i)%edges(j)%direction, edges_media(i)%edges(j)%ratio)
         end do
      end do
      l_ratio = 0.0
      do i = 1, size(cfl_map%keys)
         cell = cfl_map%keys(i)%cell
         aux_cell = cell
         cfl_info = cfl_map%getCFLInCell(cell)
         do j = FACE_X, FACE_Z
            area = cfl_info%area(j)
            idx1 = mod(j,3) + 1
            idx2 = mod(j + 1,3) + 1
            l_ratio = max(cfl_info%length(idx1),cfl_info%length(idx2))
            aux_cell(idx1) = aux_cell(idx1) + 1
            if (cfl_map%hasKey(aux_cell)) then 
               cfl_info = cfl_map%getCFLInCell(aux_cell)
               l_ratio = max(l_ratio, cfl_info%length(idx2))
            else 
               l_ratio = 1.0
            end if
            aux_cell = cell
            aux_cell(idx2) = aux_cell(idx2) + 1
            if (cfl_map%hasKey(aux_cell)) then 
               cfl_info = cfl_map%getCFLInCell(aux_cell)
               l_ratio = max(l_ratio, cfl_info%length(idx1))
            else
               l_ratio = 1.0
            end if
            if (area /= 0.0 .and. l_ratio /= 0.0) then
               res = min(res,sqrt(area/l_ratio))
            end if
         end do
      end do
   end function   


   subroutine fillEdgesFromSides(cell, sides, on_sides, edges, edge_ratios)
      integer, dimension(3), intent(in) :: cell
      type(side_t), dimension(:), allocatable, intent(in) :: sides, on_sides
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      real (kind=rkind), dimension(:), allocatable, intent(inout) :: edge_ratios
      type(side_t), dimension(:), allocatable :: sides_in_cell
      real (kind=rkind) :: ratio, delta
      integer :: face, edge
      sides_in_cell = buildCellSideSet(sides, on_sides)
      do j = 1, size(sides_in_cell)
         edge = sides_in_cell(j)%getEdge()
         if (edge /= NOT_ON_EDGE .and. (all(sides_in_cell(j)%getCell() .eq. cell))) then 
            ratio = 1.0 - sides_in_cell(j)%length()
            call addEdge(edges, cell, edge, ratio)
            if (isNewRatio(edge_ratios, ratio)) call addRatio(edge_ratios, ratio)
         end if
      end do
   end subroutine

   subroutine fillFacesFromSides(cell, sides, faces, face_ratios)
      integer, dimension(3), intent(in) :: cell
      type(side_t), dimension(:), allocatable, intent(in) :: sides
      type(side_t), dimension(:), allocatable :: sides_on_face
      type (face_t), dimension (:), allocatable, intent(inout) :: faces
      real (kind=rkind), dimension(:), allocatable, intent(inout) :: face_ratios
      type(side_t), dimension(:), allocatable :: contour
      real (kind=rkind) :: ratio, delta
      integer :: face 
      do face = FACE_X, FACE_Z
         sides_on_face = getSidesOnFace(sides, face)
         contour = buildSidesContour(sidesOnFace)
         if (size(contour) /= 0) then 
            ratio = 1.0 - contourArea(contour)
            call addFace(faces, cell, face, ratio)
            if (isNewRatio(face_ratios, ratio)) call addRatio(face_ratios, ratio)
         end if
      end do
   end subroutine

   subroutine fillConformalFaces(cell_map, faces, face_ratios)
      type(cell_map_t), intent(in) :: cell_map
      type (face_t), dimension (:), allocatable :: faces
      type(side_t), dimension(:), allocatable :: sides
      real (kind=rkind), dimension(:), allocatable :: face_ratios
      integer, dimension(3) :: cell
      integer :: i
      allocate(faces(0))
      allocate(face_ratios(0))
      do i = 1, size(cell_map%keys)
         cell = cell_map%keys(i)%cell 
         sides = cell_map%getSidesInCell(cell)
         call fillFacesFromSides(cell, sides , faces, face_ratios)
      end do
   end subroutine

   subroutine fillConformalEdges(cell_map, edges, edge_ratios)
      type(cell_map_t), intent(in) :: cell_map
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      type (side_t), dimension (:), allocatable :: sides, on_sides
      real (kind=rkind), dimension(:), allocatable, intent(inout) :: edge_ratios
      integer, dimension(3) :: cell
      integer :: i
      allocate(edges(0))
      allocate(edge_ratios(0))
      do i = 1, size(cell_map%keys)
         cell = cell_map%keys(i)%cell 
         sides = cell_map%getSidesInCell(cell) 
         on_sides = cell_map%getOnSidesInCell(cell)
         call fillEdgesFromSides(cell, sides, on_sides, edges, edge_ratios)
      end do
   end subroutine
    
   subroutine addEdge(edges, cell, edge, ratio)
      type(edge_t), dimension(:), allocatable, intent(inout) :: edges
      type(edge_t), dimension(:), allocatable :: aux
      integer, dimension(3), intent(in) :: cell
      integer :: edge 
      type(edge_t) :: new_edge
      real (kind=rkind) :: ratio
      allocate(aux(size(edges) + 1))
      aux(1:size(edges)) = edges
      new_edge = edge_t(cell=cell, ratio=ratio, direction=edge)
      aux(size(edges) + 1) = new_edge

      deallocate(edges)
      allocate(edges(size(aux)))
      edges = aux
   end subroutine

   subroutine addFace(faces, cell, face, ratio)
      type(face_t), dimension(:), allocatable, intent(inout) :: faces
      type(face_t), dimension(:), allocatable :: aux
      integer, dimension(3), intent(in) :: cell
      integer :: face
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

   logical function isNewRatio(ratios, ratio)
      real (kind=rkind), dimension(:), allocatable, intent(in) :: ratios
      real (kind=rkind) :: ratio
      integer :: i
      isNewRatio = .true.
      do i = 1, size(ratios)
         if (eq_ratio(ratios(i), ratio)) isNewRatio = .false.
      end do
   end function

   logical function eq_ratio(a,b)
      real (kind=rkind), intent(in) :: a,b
      eq_ratio = abs(a - b) < RATIO_EQ_TOLERANCE
   end function

   function filterEdgesByMedia(edges, ratio) result(res)
      type(edge_t), dimension(:), allocatable, intent(in) :: edges
      real (kind=rkind) :: ratio
      type(edge_t), dimension(:), allocatable :: res
      integer :: i,n
      n = 0
      do i = 1, size(edges)
         if (eq_ratio(edges(i)%ratio, ratio)) n = n + 1
      end do
      allocate(res(n))
      n = 0
      do i = 1, size(edges)
         if (eq_ratio(edges(i)%ratio, ratio)) then 
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
         if (eq_ratio(faces(i)%ratio, ratio)) n = n + 1
      end do
      allocate(res(n))
      n = 0
      do i = 1, size(faces)
         if (eq_ratio(faces(i)%ratio, ratio)) then 
            n = n + 1
            res(n) = faces(i)
         end if
      end do
   end function


end module