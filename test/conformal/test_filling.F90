integer function test_conformal_filling_open() bind(C) result(err)
!         /|
!       /  |
!     /    |
!   /______|_______
!   |      |_______|______
!   3     /        |      /
!   |    2         |    /
!   |  /           |  /
!   |/________1____|/

    use conformal_mod
    implicit none

    type(triangle_t) :: t
    type(triangle_t), dimension(:), allocatable :: tris

    type(coord_t) :: c1, c2, c3
    type(ConformalPECRegions) :: cR
    type(ConformalMedia_t) :: cM
    type(cell_map_t) :: cell_map
    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides

    err = 0
    c1 = coord_t(position = [0.6,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.0,0.6,0.0],  id=  2)
    c3 = coord_t(position = [0.0,0.0,0.6], id=  3)

    allocate(tris(1))
    tris(1) = triangle_t(vertices = [c1,c2,c3])

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(1))
    cR%volumes(1)%triangles(1) = tris(1)

    cM = buildConformalMedia(cR)

    if (size(cM%edge_media) /= 1) err = err + 1
    if (abs(cM%edge_media(1)%ratio - 0.4) > 0.01) err = err + 1
    if (size(cM%edge_media(1)%edges) /= 3) err = err + 1
    if (abs(cM%edge_media(1)%edges(1)%ratio-0.4) > 0.01) err = err + 1
    if (abs(cM%edge_media(1)%edges(2)%ratio-0.4) > 0.01) err = err + 1
    if (abs(cM%edge_media(1)%edges(3)%ratio-0.4) > 0.01) err = err + 1

    if (size(cM%face_media) /= 1) err = err + 1
    if (abs(cM%face_media(1)%ratio - 0.82) > 0.01) err = err + 1
    if (size(cM%face_media(1)%faces) /= 3) err = err + 1
    if (abs(cM%face_media(1)%faces(1)%ratio-0.82) > 0.01) err = err + 1
    if (abs(cM%face_media(1)%faces(2)%ratio-0.82) > 0.01) err = err + 1
    if (abs(cM%face_media(1)%faces(3)%ratio-0.82) > 0.01) err = err + 1

    if (abs(cM%cfl - 0.9055) > 0.01) err = err + 1

end function
       

integer function test_conformal_filling_closed() bind(C) result(err)
!         /|
!       /  |
!     /    |
!   /______|_______
!   |      |_______|______
!   3     /        |      /
!   |    2         |    /
!   |  /           |  /
!   4/________1____|/

    use conformal_mod
    implicit none

    type(triangle_t) :: t
    type(triangle_t), dimension(:), allocatable :: tris

    type(coord_t) :: c1, c2, c3, c4
    type(ConformalPECRegions) :: cR
    type(ConformalMedia_t) :: cM
    type(cell_map_t) :: cell_map
    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides

    err = 0
    c1 = coord_t(position = [0.6,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.0,0.6,0.0],  id=  2)
    c3 = coord_t(position = [0.0,0.0,0.6], id=  3)
    c4 = coord_t(position = [0.0,0.0,0.0], id=  4)

    allocate(tris(4))
    tris(1) = triangle_t(vertices = [c1,c2,c3])
    tris(2) = triangle_t(vertices = [c2,c4,c3])
    tris(3) = triangle_t(vertices = [c1,c3,c4])
    tris(4) = triangle_t(vertices = [c1,c4,c2])

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(4))
    cR%volumes(1)%triangles(:) = tris(:)

    cM = buildConformalMedia(cR)

    if (size(cM%edge_media) /= 1) err = err + 1
    if (abs(cM%edge_media(1)%ratio - 0.4) > 0.01) err = err + 1
    if (size(cM%edge_media(1)%edges) /= 3) err = err + 1
    if (abs(cM%edge_media(1)%edges(1)%ratio-0.4) > 0.01) err = err + 1
    if (abs(cM%edge_media(1)%edges(2)%ratio-0.4) > 0.01) err = err + 1
    if (abs(cM%edge_media(1)%edges(3)%ratio-0.4) > 0.01) err = err + 1

    if (size(cM%face_media) /= 1) err = err + 1
    if (abs(cM%face_media(1)%ratio - 0.82) > 0.01) err = err + 1
    if (size(cM%face_media(1)%faces) /= 3) err = err + 1
    if (abs(cM%face_media(1)%faces(1)%ratio-0.82) > 0.01) err = err + 1
    if (abs(cM%face_media(1)%faces(2)%ratio-0.82) > 0.01) err = err + 1
    if (abs(cM%face_media(1)%faces(3)%ratio-0.82) > 0.01) err = err + 1

    if (abs(cM%cfl - 0.9055) > 0.01) err = err + 1

end function

       

integer function  test_conformal_edge_next_cell() bind(C) result(err)
    use conformal_mod
    implicit none

    type(triangle_t) :: t
    type(triangle_t), dimension(:), allocatable :: tris

    type(coord_t) :: c1, c2, c3, c4, c5, c6
    type(ConformalPECRegions) :: cR
    type(ConformalMedia_t) :: cM
    type(cell_map_t) :: cell_map
    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides

    err = 0
    c1 = coord_t(position = [0.0,0.25,0.0], id = 1)
    c2 = coord_t(position = [0.0,1.0,1.0], id=  2)
    c3 = coord_t(position = [0.25,1.0,1.0], id=  3)
    c4 = coord_t(position = [1.0,1.0,0.5], id=  4)
    c5 = coord_t(position = [1.0,1.0,0.0], id=  5)
    c6 = coord_t(position = [0.0,1.0,0.0], id=  6)

    allocate(tris(8))
    tris(1) = triangle_t(vertices = [c1,c3,c2])
    tris(2) = triangle_t(vertices = [c1,c4,c3])
    tris(3) = triangle_t(vertices = [c1,c5,c4])
    tris(4) = triangle_t(vertices = [c1,c2,c6])
    tris(5) = triangle_t(vertices = [c1,c6,c5])
    tris(6) = triangle_t(vertices = [c6,c2,c3])
    tris(7) = triangle_t(vertices = [c6,c3,c4])
    tris(8) = triangle_t(vertices = [c6,c4,c5])

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(8))
    cR%volumes(1)%triangles(:) = tris(:)

    cM = buildConformalMedia(cR)

    if (size(cM%edge_media) /= 4) err = err + 1
    if (abs(cM%edge_media(1)%ratio - 0.25) > 0.01) err = err + 1
    if (abs(cM%edge_media(2)%ratio - 0.0) > 0.01) err = err + 1
    if (abs(cM%edge_media(3)%ratio - 0.75) > 0.01) err = err + 1
    if (abs(cM%edge_media(4)%ratio - 0.5) > 0.01) err = err + 1

    if (size(cM%edge_media(1)%edges) /= 1) err = err + 1
    if (size(cM%edge_media(2)%edges) /= 2) err = err + 1
    if (size(cM%edge_media(3)%edges) /= 1) err = err + 1
    if (size(cM%edge_media(4)%edges) /= 1) err = err + 1


    if (size(cM%face_media) /= 2) err = err + 1
    if (abs(cM%face_media(1)%ratio - 0.625) > 0.01) err = err + 1
    if (abs(cM%face_media(2)%ratio - 0.1875) > 0.01) err = err + 1
    if (size(cM%face_media(1)%faces) /= 2) err = err + 1
    if (size(cM%face_media(2)%faces) /= 1) err = err + 1


end function