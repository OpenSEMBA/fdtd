module conformal

    use geometry_mod
    use cell_map_mod
    use NFDETypes, only: Desplazamiento

contains



    function buildConformalEdgeRegions(side_map, grid) result(res)
        type(side_map_t) :: side_map
        type(Desplazamiento), intent(in) :: grid
  
        type(triangle_t), dimension(:), allocatable :: triangles
        type(side_t) :: side_limit
        type(side_t), dimension(:), allocatable :: sides, sides_on_face, contour
        ! TYPE (conformal_edge_t), dimension (:), allocatable :: res
        integer :: i, j, face, edge
        real :: ratio
        real, dimension(:), allocatable :: media
        ! allocate(res(0))
        allocate(media(0))
        do i = 1, size(side_map%keys) 
           sides = side_map%getSidesInCell(side_map%keys(i)%cell)
           do face = 1, 3
              sides_on_face = getSidesOnFace(sides, face)
              contour = buildSidesContour(sides_on_face)
              ! computeArea(contour)
              ! computeLengths(contour)
              
              do j = 1, size(contour)
                 do edge = 1, 3
                    if (contour(j)%isOnEdge(edge) .and. all(contour(j)%getCell() .eq. side_map%keys(i)%cell) ) then 
                       ! contour(j)%length()
                       if (edge == 1) then 
                          ratio = contour(j)%length() / grid%desX(side_map%keys(i)%cell(1))
                       else if (edge == 2) then 
                          ratio = contour(j)%length() / grid%desY(side_map%keys(i)%cell(2))
                       else if (edge == 3) then 
                          ratio = contour(j)%length() / grid%desZ(side_map%keys(i)%cell(3))
                       end if
  
  
                    end if
                 end do
              end do
  
           end do
        end do
        ! triangles on face are treated later
     end function
  


end module