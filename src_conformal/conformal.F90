module conformal_mod

    use geometry_mod
    use cell_map_mod
    use NFDETypes, only: Desplazamiento, ConformalPECRegions, ConformalMedia_t, edge_t, face_t

    REAL, PARAMETER :: RATIO_EQ_TOLERANCE = 0.05

contains

   function buildConformalMedia(regions, grid) result(res)
      type(ConformalPECRegions), intent(in) :: regions
      type(Desplazamiento), intent(in) :: grid
      type(ConformalMedia_t) :: res

      type(cell_map_t) :: cell_map
      real, dimension(:), allocatable :: edge_ratios, face_ratios
      type(edge_t), dimension(:), allocatable :: edges, filtered_edges
      type(face_t), dimension(:), allocatable :: faces, filtered_faces
      integer :: i
      
      do i = 1, size(regions%volumes)
         call buildCellMap(cell_map, regions%volumes(i)%triangles)
      end do
      call fillConformalEdges(cell_map, grid, edges, edge_ratios)
      call fillConformalFaces(cell_map, grid, faces, face_ratios)

      allocate(res%edge_media(size(edge_ratios)))
      do i = 1, size(edge_ratios)
         filtered_edges = filterEdgesByMedia(edges, edge_ratios(i))
         allocate(res%edge_media(i)%edges(size(filtered_edges)))
         res%edge_media(i)%edges = filtered_edges
         res%edge_media(i)%ratio = edge_ratios(i)
      end do

      allocate(res%face_media(size(face_ratios)))
      do i = 1, size(face_ratios)
         filtered_faces = filterFacesByMedia(faces, face_ratios(i))
         allocate(res%face_media(i)%faces(size(filtered_faces)))
         res%face_media(i)%faces = filtered_faces
         res%face_media(i)%ratio = face_ratios(i)
      end do
      res%n_edges_media = 0
      res%n_faces_media = 0
      if (associated(res%edge_media)) res%n_edges_media = size(res%edge_media)
      if (associated(res%face_media)) res%n_faces_media = size(res%face_media)

   end function

   subroutine fillEdgesFromSides(cell, sides, grid, edges, edge_ratios)
      integer, dimension(3), intent(in) :: cell
      type(side_t), dimension(:), allocatable, intent(in) :: sides
      type(Desplazamiento), intent(in) :: grid
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      real, dimension(:), allocatable, intent(inout) :: edge_ratios
      type(side_t), dimension(:), allocatable :: contour
      real :: ratio, delta
      integer :: face, edge
      contour = buildCellSideSet(sides)
      ! contour = buildSidesContour(getSidesOnFace(sides, face))
      do j = 1, size(contour)
         edge = contour(j)%getEdge()
         if (edge /= NOT_ON_EDGE) then 
            if (all(contour(j)%getCell() .eq. cell)) then 
               select case (edge)
               case(EDGE_X) 
                  delta = contour(j)%length()/grid%desX(cell(1))
               case(EDGE_Y)
                  delta = contour(j)%length()/grid%desY(cell(2))
               case(EDGE_Z)
                  delta = contour(j)%length()/grid%desZ(cell(3))
               end select
               ratio = 1.0 - delta
               call addEdge(edges, cell, edge, ratio)
               if (isNewRatio(edge_ratios, ratio)) call addRatio(edge_ratios, ratio)
            end if
         end if
      end do
   end subroutine

   ! subroutine fillEdgesFromTris(cell, tris, grid, edges, edge_ratios)
   !    integer, dimension(3), intent(in) :: cell
   !    type(triangle_t), dimension(:), allocatable, intent(in) :: tris
   !    type(Desplazamiento), intent(in) :: grid
   !    type (edge_t), dimension (:), allocatable, intent(inout) :: edges
   !    real, dimension(:), allocatable, intent(inout) :: edge_ratios
   !    type(side_t), dimension(:), allocatable :: sides
   !    integer :: i, j, edge
   !    do i = 1, size(tris)
   !       sides = tris(i)%getSides()
   !       ! edges on faces w/o anything else
   !       do j = 1, 3
   !          edge = sides(j)%getEdge()
   !          ratio = 1.0 - sides(j)%length()
   !          select case(edge)
   !          case(EDGE_X)
   !             ratio = ratio / grid%desX(cell(1))
   !          case(EDGE_Y)
   !             ratio = ratio / grid%desY(cell(2))
   !          case(EDGE_Z)
   !             ratio = ratio / grid%desZ(cell(3))
   !          end select
   !          if (edge /= NOT_ON_EDGE) then 
   !             call addEdge(edges, cell, edge, ratio)
   !             if (isNewRatio(edge_ratios, ratio)) call addRatio(edge_ratios, ratio)
   !          end if
   !       end do
   !    end do

   ! end subroutine

   subroutine fillFacesFromSides(cell, sides, grid, faces, face_ratios)
      integer, dimension(3), intent(in) :: cell
      type(side_t), dimension(:), allocatable, intent(in) :: sides
      type(Desplazamiento), intent(in) :: grid
      type (face_t), dimension (:), allocatable, intent(inout) :: faces
      real, dimension(:), allocatable, intent(inout) :: face_ratios
      type(side_t), dimension(:), allocatable :: contour
      real :: ratio, delta
      integer :: face 
      do face = FACE_X, FACE_Z
         contour = buildSidesContour(getSidesOnFace(sides, face))
         select case(face)
         case(FACE_X)
            delta = contourArea(contour) / (grid%desY(cell(2)) * grid%desZ(cell(3)))
         case(FACE_Y)
            delta = contourArea(contour) / (grid%desX(cell(1)) * grid%desZ(cell(3)))
         case(FACE_Z)
            delta = contourArea(contour) / (grid%desX(cell(1)) * grid%desY(cell(2)))
         end select
         ratio = 1.0 - delta
         call addFace(faces, cell, face, ratio)
         if (isNewRatio(face_ratios, ratio)) call addRatio(face_ratios, ratio)

      end do

   end subroutine

   ! subroutine fillFacesFromTris(cell, tris, grid, faces, face_ratios)
   !    integer, dimension(3), intent(in) :: cell
   !    type(triangle_t), dimension(:), allocatable, intent(in) :: tris
   !    type(Desplazamiento), intent(in) :: grid
   !    type (face_t), dimension (:), allocatable, intent(inout) :: faces
   !    real, dimension(:), allocatable, intent(inout) :: face_ratios
   ! end subroutine

   subroutine fillConformalFaces(cell_map, grid, faces, face_ratios)
      type(cell_map_t), intent(in) :: cell_map
      type(Desplazamiento), intent(in) :: grid
      type (face_t), dimension (:), allocatable :: faces
      real, dimension(:), allocatable :: face_ratios
      integer, dimension(3) :: cell
      integer :: i
      allocate(faces(0))
      allocate(face_ratios(0))
      do i = 1, size(cell_map%keys)
         cell = cell_map%keys(i)%cell 
         call fillFacesFromSides(cell, cell_map%getSidesInCell(cell), grid, faces, face_ratios)
         ! call fillFacesFromTris(cell, cell_map%getTrianglesInCell(cell), grid, faces, face_ratios)
      end do
   end subroutine

   subroutine fillConformalEdges(cell_map, grid, edges, edge_ratios)
      type(cell_map_t), intent(in) :: cell_map
      ! type(triangle_map_t), intent(in) :: tri_map
      type(Desplazamiento), intent(in) :: grid
      type (edge_t), dimension (:), allocatable, intent(inout) :: edges
      real, dimension(:), allocatable, intent(inout) :: edge_ratios
      integer, dimension(3) :: cell
      integer :: i
      allocate(edges(0))
      allocate(edge_ratios(0))
      do i = 1, size(cell_map%keys)
         cell = cell_map%keys(i)%cell 
         call fillEdgesFromSides(cell, cell_map%getSidesInCell(cell), grid, edges, edge_ratios)
         ! call fillEdgesFromTris(cell, cell_map%getTrianglesInCell(cell), grid, edges, edge_ratios)
      end do

      ! triangles on face are treated later
   end subroutine
  
   ! subroutine buildConformalFaces(cell_map, grid, res, face_ratios)
   !    type(cell_map_t), intent(in) :: cell_map
   !    ! type(triangle_map_t), intent(in) :: tri_map
   !    type(Desplazamiento), intent(in) :: grid

   !    type(side_t), dimension(:), allocatable :: sides, contour
   !    type(triangle_t), dimension(:), allocatable :: tris
   !    type (face_t), dimension (:), allocatable :: res
   !    integer :: i, j, k, face, edge
   !    real :: ratio, area
   !    real, dimension(:), allocatable :: face_ratios
   !    integer, dimension(3) :: cell
   !    allocate(res(0))
   !    allocate(face_ratios(0))
   !    do i = 1, size(cell_map%keys) 
   !       cell = cell_map%keys(i)%cell
   !       sides = cell_map%getSidesInCell(cell)
   !       do face = FACE_X, FACE_Z
   !          ! sides_on_face = getSidesOnFace(sides, face)
   !          contour = buildSidesContour(getSidesOnFace(sides, face))
   !          ! call fillEdges(contour)
   !          ! call fillFace(contour)
   !          do j = 1, size(contour)

   !             edge = contour(j)%getEdge()
   !             if (all(contour(j)%getCell() .eq. cell)) then 
   !                select case (edge)
   !                case(EDGE_X) 
   !                   ratio = contour(j)%length() / grid%desX(cell(1))
   !                   ! call addEdge(res, cell, edge, ratio)
   !                case(EDGE_Y)
   !                   ratio = contour(j)%length() / grid%desY(cell(2))
   !                   ! call addEdge(res, cell, edge, ratio)
   !                case(EDGE_Z)
   !                   ratio = contour(j)%length() / grid%desZ(cell(3))
   !                end select
   !                if (edge /= NOT_ON_EDGE) then
   !                   call addFace(res, cell, edge, ratio)
   !                   if (isNewRatio(edge_ratios, ratio)) call addRatio(edge_ratios, ratio)
   !                end if
   !             end if
   !          end do
   !       end do
   !    ! end do

   !    ! do i = 1, size(tri_map%keys) 
   !       tris = cell_map%getTrianglesInCell(cell)
   !       do j = 1, size(tris)

   !          sides = tris(j)%getSides()
   !          area = contourArea(sides)
   !          face = tris(j)%getFace()
   !          select case(face)
   !          case(FACE_X)
   !             ratio = area / (grid%desY(cell(2)) * grid%desZ(cell(3)))
   !             ! call addFace(faces, cell, face, ratio)
   !          case(FACE_Y)
   !             ratio = area / (grid%desX(cell(1)) * grid%desZ(cell(3)))
   !             ! call addFace(faces, cell, face, ratio)
   !          case(FACE_Z)
   !             ratio = area / (grid%desX(cell(1)) * grid%desY(cell(2)))
   !          end select
   !          if (FACE /= NOT_ON_FACE) then 
   !             call addFace(faces, cell, face, ratio)
   !             if (isNewRatio(face_ratios, ratio)) call addRatio(face_ratios, ratio)
   !          end if

   !          ! edges on faces w/o anything else
   !          do k = 1, 3
   !             edge = sides(edge)%getEdge()
   !             select case(edge)
   !             case(EDGE_X)
   !                ratio = sides(edge)%length() / grid%desX(cell(1))
   !                call addEdge(res, cell, edge, ratio)
   !             case(EDGE_Y)
   !                ratio = sides(edge)%length() / grid%desY(cell(2))
   !                call addEdge(res, cell, edge, ratio)
   !             case(EDGE_Z)
   !                ratio = sides(edge)%length() / grid%desZ(cell(3))
   !                call addEdge(res, cell, edge, ratio)
   !             end select
   !             if (isNewRatio(edge_ratios, ratio)) call addRatio(edge_ratios, ratio)
   !          end do

   !       end do
   !    end do

   !    ! triangles on face are treated later
   ! end subroutine
  
   subroutine addEdge(edges, cell, edge, ratio)
      type(edge_t), dimension(:), allocatable, intent(inout) :: edges
      type(edge_t), dimension(:), allocatable :: aux
      integer, dimension(3), intent(in) :: cell
      integer :: edge 
      type(edge_t) :: new_edge
      real :: ratio
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
      real :: ratio
      allocate(aux(size(faces) + 1))
      aux(1:size(faces)) = faces
      new_face = face_t(cell=cell, ratio=ratio, direction=face)
      aux(size(faces) + 1) = new_face

      deallocate(faces)
      allocate(faces(size(aux)))
      faces = aux
   end subroutine

   subroutine addRatio(ratios, ratio)
      real, dimension(:), allocatable, intent(inout) :: ratios
      real, dimension(:), allocatable :: aux
      real :: ratio
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
      real, dimension(:), allocatable, intent(in) :: ratios
      real :: ratio
      integer :: i
      isNewRatio = .true.
      do i = 1, size(ratios)
         if (eq_ratio(ratios(i), ratio)) isNewRatio = .false.
         ! if ((abs(ratios(i) - ratio) < 0.05)) isNewRatio = .false.
      end do
   end function

   logical function eq_ratio(a,b)
      real, intent(in) :: a,b
      eq_ratio = abs(a - b) < RATIO_EQ_TOLERANCE
   end function

   function filterEdgesByMedia(edges, ratio) result(res)
      type(edge_t), dimension(:), allocatable, intent(in) :: edges
      real :: ratio
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
      real :: ratio
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