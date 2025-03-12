integer function test_conformal_edges() bind(C) result(err)
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

    type(Desplazamiento) :: grid
    err = 0
    allocate(grid%desX(0:9),grid%desY(0:9),grid%desZ(0:9))
    grid%desX(:) = 1.0
    grid%desY(:) = 1.0
    grid%desZ(:) = 1.0
    c1 = coord_t(position = [0.6,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.0,0.6,0.0],  id=  2)
    c3 = coord_t(position = [0.0,0.0,0.6], id=  3)

    allocate(tris(1))
    tris(1) = triangle_t(vertices = [c1,c2,c3])

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(1))
    cR%volumes(1)%triangles(1) = tris(1)

    cM = buildConformalMedia(cR, grid)

    if (size(cM%edge_media) /= 1) err = err + 1
    if (abs(cM%edge_media(1)%ratio - 0.4) > 0.01) err = err + 1
    if (size(cM%edge_media(1)%edges) /= 3) err = err + 1
    if (abs(cM%edge_media(1)%edges(1)%ratio-0.4) > 0.01) err = err + 1
    if (abs(cM%edge_media(1)%edges(2)%ratio-0.4) > 0.01) err = err + 1
    if (abs(cM%edge_media(1)%edges(3)%ratio-0.4) > 0.01) err = err + 1

end function

integer function test_conformal_faces() bind(C) result(err)
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

    type(Desplazamiento) :: grid
    err = 0
    allocate(grid%desX(0:9),grid%desY(0:9),grid%desZ(0:9))
    grid%desX(:) = 1.0
    grid%desY(:) = 1.0
    grid%desZ(:) = 1.0
    c1 = coord_t(position = [0.6,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.0,0.6,0.0],  id=  2)
    c3 = coord_t(position = [0.0,0.0,0.6], id=  3)

    allocate(tris(1))
    tris(1) = triangle_t(vertices = [c1,c2,c3])

    allocate(cR%volumes(1))
    allocate(cR%volumes(1)%triangles(1))
    cR%volumes(1)%triangles(1) = tris(1)

    cM = buildConformalMedia(cR, grid)

    if (size(cM%face_media) /= 1) err = err + 1
    if (abs(cM%face_media(1)%ratio - 0.82) > 0.01) err = err + 1
    if (size(cM%face_media(1)%faces) /= 3) err = err + 1
    if (abs(cM%face_media(1)%faces(1)%ratio-0.82) > 0.01) err = err + 1
    if (abs(cM%face_media(1)%faces(2)%ratio-0.82) > 0.01) err = err + 1
    if (abs(cM%face_media(1)%faces(3)%ratio-0.82) > 0.01) err = err + 1

end function