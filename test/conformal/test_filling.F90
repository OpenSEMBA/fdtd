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

integer function test_conformal_filling_closed_corner() bind(C) result(err)
!         /6
!       /  |\
!     /    |  \
!   3______|____\__4
!   |      5_______|______
!   |     / \      |      /
!   |    /    \    |    /
!   |  /        \  |  /
!   1/_____________2/

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
    c1 = coord_t(position = [0.0,0.0,0.0], id = 1)
    c2 = coord_t(position = [1.0,0.0,0.0], id=  2)
    c3 = coord_t(position = [0.0,0.0,1.0], id=  3)
    c4 = coord_t(position = [1.0,0.0,1.0], id=  4)
    c5 = coord_t(position = [0.0,1.0,0.0], id=  5)
    c6 = coord_t(position = [0.0,1.0,1.0], id=  6)

    allocate(tris(8))
    tris(1) = triangle_t(vertices = [c1,c2,c3])
    tris(2) = triangle_t(vertices = [c2,c4,c3])
    tris(3) = triangle_t(vertices = [c1,c3,c6])
    tris(4) = triangle_t(vertices = [c1,c6,c5])
    tris(5) = triangle_t(vertices = [c3,c4,c6])
    tris(6) = triangle_t(vertices = [c1,c5,c2])
    tris(7) = triangle_t(vertices = [c2,c5,c4])
    tris(8) = triangle_t(vertices = [c5,c6,c4])

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(8))
    cR%volumes(1)%triangles(:) = tris(:)

    cM = buildConformalMedia(cR)

    if (size(cM%edge_media) /= 1) err = err + 1
    if (abs(cM%edge_media(1)%ratio - 0.0) > 0.01) err = err + 1
    if (size(cM%edge_media(1)%edges) /= 7) err = err + 1
    if (abs(cM%edge_media(1)%edges(1)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(2)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(3)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(4)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(5)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(6)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(7)%ratio-0.0) > 0.00) err = err + 1

    if (size(cM%face_media) /= 2) err = err + 1

    if (abs(cM%face_media(1)%ratio - 0.50) > 0.01) err = err + 1
    if (size(cM%face_media(1)%faces) /= 2) err = err + 1
    if (abs(cM%face_media(1)%faces(1)%ratio-0.0) > 0.50) err = err + 1
    if (abs(cM%face_media(1)%faces(2)%ratio-0.0) > 0.50) err = err + 1

    if (abs(cM%face_media(2)%ratio - 0.00) > 0.01) err = err + 1
    if (size(cM%face_media(2)%faces) /= 2) err = err + 1
    if (abs(cM%face_media(2)%faces(1)%ratio-0.) > 0.00) err = err + 1
    if (abs(cM%face_media(2)%faces(2)%ratio-0.) > 0.00) err = err + 1


    ! if (abs(cM%cfl - 0.9055) > 0.01) err = err + 1

end function

integer function test_conformal_filling_block_and_corner() bind(C) result(err)

!              9 
!            / |\
!          5______6
!         /    |  | 
!       /  |   10 |
!     /    |  /  \|   
!   3______4_/_____
!   |      7______8|______
!   |     /        |      /
!   |    /         |    /
!   |  /           |  /
!   1/______2______ /
    use conformal_mod
    implicit none

    type(triangle_t) :: t
    type(triangle_t), dimension(:), allocatable :: tris

    type(coord_t) :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10
    type(ConformalPECRegions) :: cR
    type(ConformalMedia_t) :: cM
    type(cell_map_t) :: cell_map
    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides

    err = 0
    c1 = coord_t(position = [0.0,0.0,0.0], id = 1)
    c2 = coord_t(position = [0.5,0.0,0.0], id=  2)
    c3 = coord_t(position = [0.0,0.0,1.0], id=  3)
    c4 = coord_t(position = [0.5,0.0,1.0], id=  4)
    c5 = coord_t(position = [0.0,1.0,1.0], id=  5)
    c6 = coord_t(position = [0.5,1.0,1.0], id=  6)
    c7 = coord_t(position = [0.0,1.0,0.0], id=  7)
    c8 = coord_t(position = [0.5,1.0,0.0], id=  8)
    c9 = coord_t(position = [0.0,1.5,1.0], id=  9)
    c10= coord_t(position = [0.0,1.5,0.0], id=  10)

    allocate(tris(16))
    tris(1) = triangle_t(vertices = [c1,c2,c3])
    tris(2) = triangle_t(vertices = [c2,c4,c3])
    tris(3) = triangle_t(vertices = [c1,c3,c5])
    tris(4) = triangle_t(vertices = [c1,c5,c7])
    tris(5) = triangle_t(vertices = [c3,c4,c5])
    tris(6) = triangle_t(vertices = [c4,c6,c5])
    tris(7) = triangle_t(vertices = [c1,c7,c2])
    tris(8) = triangle_t(vertices = [c2,c7,c8])
    tris(9) = triangle_t(vertices = [c2,c8,c6])
    tris(10) = triangle_t(vertices = [c2,c6,c4])
    tris(11) = triangle_t(vertices = [c5,c6,c9])
    tris(12) = triangle_t(vertices = [c7,c10,c8])
    tris(13) = triangle_t(vertices = [c5,c9,c10])
    tris(14) = triangle_t(vertices = [c5,c10,c7])
    tris(15) = triangle_t(vertices = [c9,c6,c10])
    tris(16) = triangle_t(vertices = [c8,c10,c6])

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(16))
    cR%volumes(1)%triangles(:) = tris(:)

    cM = buildConformalMedia(cR)

    if (size(cM%edge_media) /= 2) err = err + 1

    if (abs(cM%edge_media(1)%ratio - 0.5) > 0.01) err = err + 1
    if (size(cM%edge_media(1)%edges) /= 6) err = err + 1
    if (abs(cM%edge_media(1)%edges(1)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(2)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(3)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(4)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(5)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%edge_media(1)%edges(6)%ratio-0.5) > 0.00) err = err + 1

    if (abs(cM%edge_media(2)%ratio - 0.0) > 0.01) err = err + 1
    if (size(cM%edge_media(2)%edges) /= 4) err = err + 1
    if (abs(cM%edge_media(2)%edges(1)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(2)%edges(2)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(2)%edges(3)%ratio-0.0) > 0.00) err = err + 1
    if (abs(cM%edge_media(2)%edges(4)%ratio-0.0) > 0.00) err = err + 1

    if (size(cM%face_media) /= 3) err = err + 1

    if (abs(cM%face_media(1)%ratio - 0.5) > 0.01) err = err + 1
    if (size(cM%face_media(1)%faces) /= 5) err = err + 1
    if (abs(cM%face_media(1)%faces(1)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%face_media(1)%faces(2)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%face_media(1)%faces(3)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%face_media(1)%faces(4)%ratio-0.5) > 0.00) err = err + 1
    if (abs(cM%face_media(1)%faces(5)%ratio-0.5) > 0.00) err = err + 1

    if (abs(cM%face_media(2)%ratio - 0.00) > 0.01) err = err + 1
    if (size(cM%face_media(2)%faces) /= 1) err = err + 1
    if (abs(cM%face_media(2)%faces(1)%ratio-0.0) > 0.50) err = err + 1


    if (abs(cM%face_media(3)%ratio - 0.875) > 0.01) err = err + 1
    if (size(cM%face_media(3)%faces) /= 2) err = err + 1
    if (abs(cM%face_media(3)%faces(1)%ratio-0.875) > 0.00) err = err + 1
    if (abs(cM%face_media(3)%faces(2)%ratio-0.875) > 0.00) err = err + 1

end function

integer function test_conformal_filling_cylinder_base_on_grid_plane() bind(C) result(err)
!         /
!       /  |\
!     /    |  \
!   3______|____\__
!   |      5_______|______
!   |     / \      |      /
!   |    /    \    |    /
!   |  /        \  |  /
!   1/_____________2/

    use conformal_mod
    implicit none

    type(triangle_t) :: t
    type(triangle_t), dimension(:), allocatable :: tris

    type(coord_t), dimension(24) :: c
    type(ConformalPECRegions) :: cR
    type(ConformalMedia_t) :: cM
    type(cell_map_t) :: cell_map

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides

    err = 0
    c(1) = coord_t(position  = [1.0,0.25,0.0], id = 1)
    c(2) = coord_t(position  = [0.0,1.0, 0.0], id = 2)
    c(3) = coord_t(position  = [0.0,2.0, 0.0], id = 3)
    c(4) = coord_t(position  = [1.0,2.25,0.0], id = 4)
    c(5) = coord_t(position  = [2.0,2.25,0.0], id = 5)
    c(6) = coord_t(position  = [2.5,2.0, 0.0], id = 6)
    c(7) = coord_t(position  = [2.5,1.0, 0.0], id = 7)
    c(8) = coord_t(position  = [2.0,0.25,0.0], id = 8)
    c(9) = coord_t(position  = [1.0,1.0, 0.0], id = 9)
    c(10) = coord_t(position = [1.0,2.0, 0.0], id = 10)
    c(11) = coord_t(position = [2.0,2.0, 0.0], id = 11)
    c(12) = coord_t(position = [2.0,1.0, 0.0], id = 12)
    c(13) = coord_t(position = [1.0,0.25,1.0], id = 13)
    c(14) = coord_t(position = [0.0,1.0,   1.0], id = 14)
    c(15) = coord_t(position = [0.0,2.0,   1.0], id = 15)
    c(16) = coord_t(position = [1.0,2.25,1.0], id = 16)
    c(17) = coord_t(position = [2.0,2.25,1.0], id = 17)
    c(18) = coord_t(position = [2.5,2.0, 1.0], id = 18)
    c(19) = coord_t(position = [2.5,1.0, 1.0], id = 19)
    c(20) = coord_t(position = [2.0,0.25,1.0], id = 20)
    c(21) = coord_t(position = [1.0,1.0,   1.0], id = 21)
    c(21) = coord_t(position = [1.0,2.0,   1.0], id = 22)
    c(22) = coord_t(position = [2.0,2.0,   1.0], id = 23)
    c(23) = coord_t(position = [2.0,1.0,   1.0], id = 24)

    allocate(tris(44))

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(44))
    cR%volumes(1)%triangles(:) = tris(:)

    cM = buildConformalMedia(cR)

    ! if (size(cM%edge_media) /= 1) err = err + 1
    ! if (abs(cM%edge_media(1)%ratio - 0.0) > 0.01) err = err + 1
    ! if (size(cM%edge_media(1)%edges) /= 7) err = err + 1
    ! if (abs(cM%edge_media(1)%edges(1)%ratio-0.0) > 0.00) err = err + 1
    ! if (abs(cM%edge_media(1)%edges(2)%ratio-0.0) > 0.00) err = err + 1
    ! if (abs(cM%edge_media(1)%edges(3)%ratio-0.0) > 0.00) err = err + 1
    ! if (abs(cM%edge_media(1)%edges(4)%ratio-0.0) > 0.00) err = err + 1
    ! if (abs(cM%edge_media(1)%edges(5)%ratio-0.0) > 0.00) err = err + 1
    ! if (abs(cM%edge_media(1)%edges(6)%ratio-0.0) > 0.00) err = err + 1
    ! if (abs(cM%edge_media(1)%edges(7)%ratio-0.0) > 0.00) err = err + 1

    ! if (size(cM%face_media) /= 2) err = err + 1

    ! if (abs(cM%face_media(1)%ratio - 0.50) > 0.01) err = err + 1
    ! if (size(cM%face_media(1)%faces) /= 2) err = err + 1
    ! if (abs(cM%face_media(1)%faces(1)%ratio-0.0) > 0.50) err = err + 1
    ! if (abs(cM%face_media(1)%faces(2)%ratio-0.0) > 0.50) err = err + 1

    ! if (abs(cM%face_media(2)%ratio - 0.00) > 0.01) err = err + 1
    ! if (size(cM%face_media(2)%faces) /= 2) err = err + 1
    ! if (abs(cM%face_media(2)%faces(2)%ratio-0.) > 0.00) err = err + 1
    ! if (abs(cM%face_media(2)%faces(1)%ratio-0.) > 0.00) err = err + 1


    ! if (abs(cM%cfl - 0.9055) > 0.01) err = err + 1

end function
