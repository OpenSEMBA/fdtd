integer function test_geometry_coord_position() bind(C) result(err)
    use geometry_mod
    implicit none
    type(coord_t) :: c
    err = 0
    c = coord_t(position = [0.0,0.0,0.0], id = 1)
    if (.not. c%isOnVertex()) err = err + 1
    if (c%isOnAnyEdge()) err = err + 1
    if (c%isOnAnyFace()) err = err + 1

    c = coord_t(position = [1.0,0.0,0.0], id = 1)
    if (.not. c%isOnVertex()) err = err + 1
    if (c%isOnAnyEdge()) err = err + 1
    if (c%isOnAnyFace()) err = err + 1

    c = coord_t(position = [0.0,1.0,0.0], id = 1)
    if (.not. c%isOnVertex()) err = err + 1
    if (c%isOnAnyEdge()) err = err + 1
    if (c%isOnAnyFace()) err = err + 1

    c = coord_t(position = [0.0,0.0,1.0], id = 1)
    if (.not. c%isOnVertex()) err = err + 1
    if (c%isOnAnyEdge()) err = err + 1
    if (c%isOnAnyFace()) err = err + 1

    c = coord_t(position = [0.5,0.0,0.0], id = 1)
    if (c%isOnVertex()) err = err + 1
    if (.not. c%isOnEdge(EDGE_X)) err = err + 1
    if (c%isOnAnyFace()) err = err + 1

    c = coord_t(position = [0.0,0.5,0.0], id = 1)
    if (c%isOnVertex()) err = err + 1
    if (.not. c%isOnEdge(EDGE_Y)) err = err + 1
    if (c%isOnANyFace()) err = err + 1

    c = coord_t(position = [0.0,0.0,0.5], id = 1)
    if (c%isOnVertex()) err = err + 1
    if (.not. c%isOnEdge(EDGE_Z)) err = err + 1
    if (c%isOnAnyFace()) err = err + 1

    c = coord_t(position = [0.5,0.5,0.0], id = 1)
    if (c%isOnVertex()) err = err + 1
    if (c%isOnAnyEdge()) err = err + 1
    if (.not. c%isOnFace(FACE_Z)) err = err + 1

    c = coord_t(position = [0.5,0.0,0.5], id = 1)
    if (c%isOnVertex()) err = err + 1
    if (c%isOnAnyEdge()) err = err + 1
    if (.not. c%isOnFace(FACE_Y)) err = err + 1

    c = coord_t(position = [0.0,0.5,0.5], id = 1)
    if (c%isOnVertex()) err = err + 1
    if (c%isOnAnyEdge()) err = err + 1
    if (.not. c%isOnFace(FACE_X)) err = err + 1

end function

integer function test_geometry_side_position() bind(C) result(err)
    use geometry_mod
    implicit none
    type(side_t) :: side
    type(coord_t) :: c1, c2, c3, c4, c5
    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0], id = 1)
    c2 = coord_t(position = [1.0,0.0,0.0], id=  2)
    c3 = coord_t(position = [0.0,0.0,1.0], id=  3)
    c4 = coord_t(position = [0.0,1.0,1.0], id=  4)
    c5 = coord_t(position = [0.0,1.0,0.0], id=  5)

    side = side_t(init=c1, end = c2)
    if (.not. side%isOnEdge(EDGE_X)) err = err + 1
    if (side%isOnAnyFace()) err = err + 1

    side = side_t(init=c1, end = c5)
    if (.not. side%isOnEdge(EDGE_Y)) err = err + 1
    if (side%isOnAnyFace()) err = err + 1

    side = side_t(init=c1, end = c3)
    if (.not. side%isOnEdge(EDGE_Z)) err = err + 1
    if (side%isOnAnyFace()) err = err + 1


    side = side_t(init=c1, end = c4)
    if (side%isOnAnyEdge()) err = err + 1
    if (.not. side%isOnFace(FACE_X)) err = err + 1

end function

integer function test_geometry_triangle_on_face() bind(C) result(err)

    use geometry_mod
    implicit none
    type(triangle_t) :: t
    type(coord_t) :: c1, c2, c3, c4, c5
    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0], id = 1)
    c2 = coord_t(position = [1.0,0.0,0.0], id=  2)
    c3 = coord_t(position = [0.0,0.0,1.0], id=  3)
    c4 = coord_t(position = [0.0,1.0,1.0], id=  4)
    c5 = coord_t(position = [0.0,1.0,0.0], id=  5)
    err = 0

    t = triangle_t(vertices = [c1,c3,c5])
    if (.not. t%isOnAnyFace()) err = err + 1
    if (.not. t%getFace() == FACE_X) err = err + 1

    t = triangle_t(vertices = [c1,c2,c3])
    if (.not. t%isOnAnyFace()) err = err + 1
    if (.not. t%getFace() == FACE_Y) err = err + 1

    t = triangle_t(vertices = [c1,c2,c5])
    if (.not. t%isOnAnyFace()) err = err + 1
    if (.not. t%getFace() == FACE_Z) err = err + 1

    t = triangle_t(vertices = [c1,c2,c4])
    if (t%isOnAnyFace()) err = err + 1
    if (.not. t%getFace() == NOT_ON_FACE) err = err + 1


end function

integer function test_geometry_triangle_normal() bind(C) result(err)

    use geometry_mod
    implicit none
    type(triangle_t) :: t
    err = 0

    
    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,0.0,1.0]; t%vertices(3)%position = [0.0,1.0,0.0]
    if (t%getFace() /= FACE_X) err = err + 1
    
    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,0.0,1.0]; t%vertices(3)%position = [1.0,0.0,0.0]
    if (t%getFace() /= FACE_Y) err = err + 1
    
    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,1.0,0.0]; t%vertices(3)%position = [1.0,0.0,0.0]
    if (t%getFace() /= FACE_Z) err = err + 1

    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,1.0,0.0]; t%vertices(3)%position = [1.0,0.0,1.0]
    if (t%getFace() /= NOT_ON_FACE) err = err + 1


end function

integer function test_geometry_triangle_edges() bind(C) result(err)

    use geometry_mod
    implicit none
    type(triangle_t) :: t
    type(side_t), dimension(3) :: sides
    err = 0
    
    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,0.0,1.0]; t%vertices(3)%position = [0.0,1.0,0.0]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Z) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_Y) err = err + 1
    
    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,0.0,1.0]; t%vertices(3)%position = [1.0,0.0,0.0]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Z) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_X) err = err + 1
    
    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,1.0,0.0]; t%vertices(3)%position = [1.0,0.0,0.0]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Y) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_X) err = err + 1

    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,1.0,0.0]; t%vertices(3)%position = [1.0,0.0,1.0]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Y) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= NOT_ON_EDGE) err = err + 1


end function

integer function test_geometry_triangle_cell() bind(C) result(err)

    use geometry_mod
    implicit none
    type(triangle_t) :: t
    type(side_t), dimension(3) :: sides
    integer, dimension(3) :: cell
    err = 0
    
    t%vertices(1)%position = [0.0,0.0,0.0]; t%vertices(2)%position = [0.0,0.0,1.0]; t%vertices(3)%position = [0.0,1.0,0.0]
    cell = [0.0,0.0,0.0]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1
    
    t%vertices(1)%position = [1.0,0.0,0.0]; t%vertices(2)%position = [1.0,0.0,1.0]; t%vertices(3)%position = [1.0,1.0,0.0]
    cell = [1.0,0.0,0.0]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1
    
    t%vertices(1)%position = [1.0,0.0,1.0]; t%vertices(2)%position = [1.0,0.0,2.0]; t%vertices(3)%position = [1.0,1.0,1.0]
    cell = [1.0,0.0,1.0]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1

    t%vertices(1)%position = [0.0,0.0,1.0]; t%vertices(2)%position = [1.0,0.0,1.0]; t%vertices(3)%position = [1.0,1.0,1.0]
    cell = [0.0,0.0,1.0]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1


end function

integer function test_geometry_elements_in_cell() bind(C) result(err)


    use cell_map_mod
    use geometry_mod
    implicit none 

    type(triangle_map_t) :: tri_map
    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides
    type(triangle_t), dimension(:), allocatable :: triangles, tris
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5
    integer, dimension(3) :: cell

    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0], id = 1)
    c2 = coord_t(position = [1.0,0.0,0.0], id=  2)
    c3 = coord_t(position = [0.0,0.0,1.0], id=  3)
    c4 = coord_t(position = [0.0,1.0,1.0], id=  4)
    c5 = coord_t(position = [0.0,1.0,0.0], id=  5)
    allocate(triangles(3))
    triangles(1) = triangle_t(vertices = [c1,c5,c3])
    triangles(2) = triangle_t(vertices = [c1,c2,c5])
    triangles(3) = triangle_t(vertices = [c1,c2,c4])

    cell = triangles(1)%getCell()
    call buildTriangleMap(tri_map, triangles)
    tris = tri_map%getTrianglesInCell(cell)
    if (size(tris) /= 2) err = err + 1
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    if (size(sides) /= 2) err = err + 1

end function

integer function test_geometry_map_sides() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(triangle_map_t) :: tri_map
    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path
    type(triangle_t), dimension(:), allocatable :: triangles, tris
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6
    integer, dimension(3) :: cell
    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.75,0.0,0.0],  id=  2)
    c3 = coord_t(position = [0.75,0.25,0.0], id=  3)
    c4 = coord_t(position = [0.25,0.75,0.0], id=  4)
    c5 = coord_t(position = [0.0,0.75,0.0],  id=  5)
    c6 = coord_t(position = [0.0,0.0,0.75],  id=  6)
    allocate(triangles(8))
    ! xz
    triangles(1) = triangle_t(vertices = [c1,c2,c6])
    ! yz
    triangles(2) = triangle_t(vertices = [c1,c6,c5])
    ! xy
    triangles(3) = triangle_t(vertices = [c1,c3,c2])
    triangles(4) = triangle_t(vertices = [c1,c4,c3])
    triangles(5) = triangle_t(vertices = [c1,c5,c4])
    ! off faces
    triangles(6) = triangle_t(vertices = [c2,c3,c6])
    triangles(7) = triangle_t(vertices = [c3,c4,c6])
    triangles(8) = triangle_t(vertices = [c4,c5,c6])

    cell = triangles(1)%getCell()
    call buildTriangleMap(tri_map, triangles)
    tris = tri_map%getTrianglesInCell(cell)
    if (size(tris) /= 5) err = err + 1
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    if (size(sides) /= 5) err = err + 1

    if (size(getSidesOnFace(sides, FACE_X)) /= 1) err = err + 1
    if (size(getSidesOnFace(sides, FACE_Y)) /= 1) err = err + 1
    if (size(getSidesOnFace(sides, FACE_Z)) /= 3) err = err + 1

    if (size(getSidesOnEdge(sides, EDGE_X)) /= 0) err = err + 1
    if (size(getSidesOnEdge(sides, EDGE_Y)) /= 0) err = err + 1
    if (size(getSidesOnEdge(sides, EDGE_Z)) /= 0) err = err + 1

end function

integer function test_geometry_path() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(triangle_map_t) :: tri_map
    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path
    type(triangle_t), dimension(:), allocatable :: triangles, tris
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6
    integer, dimension(3) :: cell
    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.75,0.0,0.0],  id=  2)
    c3 = coord_t(position = [0.75,0.25,0.0], id=  3)
    c4 = coord_t(position = [0.25,0.75,0.0], id=  4)
    c5 = coord_t(position = [0.0,0.75,0.0],  id=  5)
    c6 = coord_t(position = [0.0,0.0,0.75],  id=  6)
    allocate(triangles(8))
    ! xz
    triangles(1) = triangle_t(vertices = [c1,c2,c6])
    ! yz
    triangles(2) = triangle_t(vertices = [c1,c6,c5])
    ! xy
    triangles(3) = triangle_t(vertices = [c1,c3,c2])
    triangles(4) = triangle_t(vertices = [c1,c4,c3])
    triangles(5) = triangle_t(vertices = [c1,c5,c4])
    ! off faces
    triangles(6) = triangle_t(vertices = [c2,c3,c6])
    triangles(7) = triangle_t(vertices = [c3,c4,c6])
    triangles(8) = triangle_t(vertices = [c4,c5,c6])

    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    path = getPathOnFace(getSidesOnFace(sides, FACE_Z))
    if (.not. path(1)%init%id == 2) err = err + 1
    if (.not. path(1)%end%id == 3) err = err + 1
    if (.not. path(2)%init%id == 3) err = err + 1
    if (.not. path(2)%end%id == 4) err = err + 1
    if (.not. path(3)%init%id == 4) err = err + 1
    if (.not. path(3)%end%id == 5) err = err + 1

end function

integer function test_geometry_vertex_vertex_contour() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path, contour
    type(triangle_t), dimension(:), allocatable :: triangles
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6
    integer, dimension(3) :: cell
    type(coord_t) :: init, end

    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [1.0,0.0,0.0],  id=  2)
    c3 = coord_t(position = [0.75,0.25,0.0], id=  3)
    c4 = coord_t(position = [0.25,0.75,0.0], id=  4)
    c5 = coord_t(position = [0.0,1.0,0.0],  id=  5)
    c6 = coord_t(position = [0.0,0.0,1.0],  id=  6)
    allocate(triangles(8))
    ! xz
    triangles(1) = triangle_t(vertices = [c1,c2,c6])
    ! yz
    triangles(2) = triangle_t(vertices = [c1,c6,c5])
    ! xy
    triangles(3) = triangle_t(vertices = [c1,c3,c2])
    triangles(4) = triangle_t(vertices = [c1,c4,c3])
    triangles(5) = triangle_t(vertices = [c1,c5,c4])
    ! off faces
    triangles(6) = triangle_t(vertices = [c2,c3,c6])
    triangles(7) = triangle_t(vertices = [c3,c4,c6])
    triangles(8) = triangle_t(vertices = [c4,c5,c6])

    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    
    path = getPathOnFace(getSidesOnFace(sides, FACE_Z))
    init = path(1)%init
    end = path(size(path))%end
    contour = buildSidesContour(path)
    if (size(contour) /= 5) err = err + 1
    if (.not. all(contour(4)%init%position .eq. c5%position)) err = err + 1
    if (.not. all(contour(4)%end%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%init%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%end%position .eq. c2%position)) err = err + 1

end function

integer function test_geometry_vertex_side_contour() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path, contour
    type(triangle_t), dimension(:), allocatable :: triangles
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6
    integer, dimension(3) :: cell
    type(coord_t) :: init, end

    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [1.0,0.0,0.0],  id=  2)
    c3 = coord_t(position = [0.75,0.25,0.0], id=  3)
    c4 = coord_t(position = [0.25,0.75,0.0], id=  4)
    c5 = coord_t(position = [0.0,0.75,0.0],  id=  5)
    c6 = coord_t(position = [0.0,0.0,1.00],  id=  6)
    allocate(triangles(8))
    ! xz
    triangles(1) = triangle_t(vertices = [c1,c2,c6])
    ! yz
    triangles(2) = triangle_t(vertices = [c1,c6,c5])
    ! xy
    triangles(3) = triangle_t(vertices = [c1,c3,c2])
    triangles(4) = triangle_t(vertices = [c1,c4,c3])
    triangles(5) = triangle_t(vertices = [c1,c5,c4])
    ! off faces
    triangles(6) = triangle_t(vertices = [c2,c3,c6])
    triangles(7) = triangle_t(vertices = [c3,c4,c6])
    triangles(8) = triangle_t(vertices = [c4,c5,c6])

    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)

    path = getPathOnFace(getSidesOnFace(sides, FACE_Z))
    init = path(1)%init
    end = path(size(path))%end
    contour = buildSidesContour(path)
    if (size(contour) /= 5) err = err + 1
    if (.not. all(contour(4)%init%position .eq. c5%position)) err = err + 1
    if (.not. all(contour(4)%end%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%init%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%end%position .eq. c2%position)) err = err + 1


end function

integer function test_geometry_side_vertex_contour() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path, contour
    type(triangle_t), dimension(:), allocatable :: triangles
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6
    integer, dimension(3) :: cell
    type(coord_t) :: init, end

    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.75,0.0,0.0],  id=  2)
    c3 = coord_t(position = [0.75,0.25,0.0], id=  3)
    c4 = coord_t(position = [0.25,0.75,0.0], id=  4)
    c5 = coord_t(position = [0.0,1.0,0.0],  id=  5)
    c6 = coord_t(position = [0.0,0.0,1.00],  id=  6)
    allocate(triangles(8))
    ! xz
    triangles(1) = triangle_t(vertices = [c1,c2,c6])
    ! yz
    triangles(2) = triangle_t(vertices = [c1,c6,c5])
    ! xy
    triangles(3) = triangle_t(vertices = [c1,c3,c2])
    triangles(4) = triangle_t(vertices = [c1,c4,c3])
    triangles(5) = triangle_t(vertices = [c1,c5,c4])
    ! off faces
    triangles(6) = triangle_t(vertices = [c2,c3,c6])
    triangles(7) = triangle_t(vertices = [c3,c4,c6])
    triangles(8) = triangle_t(vertices = [c4,c5,c6])

    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)

    path = getPathOnFace(getSidesOnFace(sides, FACE_Z))
    init = path(1)%init
    end = path(size(path))%end
    contour = buildSidesContour(path)
    if (size(contour) /= 5) err = err + 1
    if (.not. all(contour(4)%init%position .eq. c5%position)) err = err + 1
    if (.not. all(contour(4)%end%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%init%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%end%position .eq. c2%position)) err = err + 1


end function

integer function test_geometry_side_side_contour() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path, contour
    type(triangle_t), dimension(:), allocatable :: triangles
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6
    integer, dimension(3) :: cell
    type(coord_t) :: init, end

    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.75,0.0,0.0],  id=  2)
    c3 = coord_t(position = [0.75,0.25,0.0], id=  3)
    c4 = coord_t(position = [0.25,0.75,0.0], id=  4)
    c5 = coord_t(position = [0.0,0.75,0.0],  id=  5)
    c6 = coord_t(position = [0.0,0.0,0.75],  id=  6)
    allocate(triangles(8))
    ! xz
    triangles(1) = triangle_t(vertices = [c1,c2,c6])
    ! yz
    triangles(2) = triangle_t(vertices = [c1,c6,c5])
    ! xy
    triangles(3) = triangle_t(vertices = [c1,c3,c2])
    triangles(4) = triangle_t(vertices = [c1,c4,c3])
    triangles(5) = triangle_t(vertices = [c1,c5,c4])
    ! off faces
    triangles(6) = triangle_t(vertices = [c2,c3,c6])
    triangles(7) = triangle_t(vertices = [c3,c4,c6])
    triangles(8) = triangle_t(vertices = [c4,c5,c6])

    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)

    path = getPathOnFace(getSidesOnFace(sides, FACE_Z))
    init = path(1)%init
    end = path(size(path))%end
    contour = buildSidesContour(path)
    if (size(contour) /= 5) err = err + 1
    if (.not. all(contour(4)%init%position .eq. c5%position)) err = err + 1
    if (.not. all(contour(4)%end%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%init%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(5)%end%position .eq. c2%position)) err = err + 1

end function

integer function test_geometry_side_side_contour_2() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path, contour
    type(triangle_t), dimension(:), allocatable :: triangles
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6, c7
    integer, dimension(3) :: cell
    type(coord_t) :: init, end

    err = 0

    c1 = coord_t(position = [0.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [1.0,0.0,0.0],  id=  2)
    c3 = coord_t(position = [0.0,0.5,0.0], id=  3)
    c4 = coord_t(position = [1.0,0.5,0.0], id=  4)
    c5 = coord_t(position = [0.25,0.25,0.0],  id=  5)
    c6 = coord_t(position = [0.0,0.0,1.0],  id=  6)
    c7 = coord_t(position = [1.0,0.0,1.0],  id=  7)
    allocate(triangles(10))
    ! xz
    triangles(1) = triangle_t(vertices = [c1,c7,c6])
    triangles(2) = triangle_t(vertices = [c1,c2,c7])
    ! yz
    triangles(3) = triangle_t(vertices = [c1,c3,c6])
    triangles(4) = triangle_t(vertices = [c2,c4,c7])
    ! xy
    triangles(5) = triangle_t(vertices = [c1,c3,c5])
    triangles(6) = triangle_t(vertices = [c1,c5,c2])
    triangles(7) = triangle_t(vertices = [c2,c5,c4])
    ! off faces
    triangles(8) = triangle_t(vertices = [c5,c3,c6])
    triangles(9) = triangle_t(vertices = [c4,c5,c7])
    triangles(10) = triangle_t(vertices = [c5,c6,c7])

    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    if (size(sides) /= 3) err = err + 1

    path = getPathOnFace(getSidesOnFace(sides, FACE_Z))
    init = path(1)%init
    end = path(size(path))%end
    contour = buildSidesContour(path)
    if (size(contour) /= 5) err = err + 1
    if (.not. all(contour(3)%init%position .eq. c3%position)) err = err + 1
    if (.not. all(contour(3)%end%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(4)%init%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(4)%end%position .eq. c2%position)) err = err + 1
    if (.not. all(contour(5)%init%position .eq. c2%position)) err = err + 1
    if (.not. all(contour(5)%end%position .eq. c4%position)) err = err + 1

end function

integer function test_geometry_side_side_contour_3() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path, contour
    type(triangle_t), dimension(:), allocatable :: triangles
    type(triangle_t) :: t1
    type(coord_t) :: c1, c2, c3, c4, c5
    integer, dimension(3) :: cell
    type(coord_t) :: init, end

    err = 0

    c1 = coord_t(position = [0.0,0.0,0.25], id = 1)
    c2 = coord_t(position = [1.0,0.0,0.25], id=  2)
    c3 = coord_t(position = [0.0,1.0,0.25], id=  3)
    c4 = coord_t(position = [0.0,0.0,0.0],  id=  4)
    c5 = coord_t(position = [1.0,0.0,0.0],  id=  5)
    allocate(triangles(1))
    triangles(1) = triangle_t(vertices = [c1,c2,c3])

    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    if (size(sides) /= 2) err = err + 1

    path = getPathOnFace(getSidesOnFace(sides, FACE_Y))
    init = path(1)%init
    end = path(size(path))%end
    contour = buildSidesContour(path)
    if (size(contour) /= 4) err = err + 1
    if (.not. all(contour(1)%init%position .eq. c1%position)) err = err + 1
    if (.not. all(contour(1)%end%position .eq. c2%position)) err = err + 1
    if (.not. all(contour(2)%init%position .eq. c2%position)) err = err + 1
    if (.not. all(contour(2)%end%position .eq. c5%position)) err = err + 1
    if (.not. all(contour(3)%init%position .eq. c5%position)) err = err + 1
    if (.not. all(contour(3)%end%position .eq. c4%position)) err = err + 1
    if (.not. all(contour(4)%init%position .eq. c4%position)) err = err + 1
    if (.not. all(contour(4)%end%position .eq. c1%position)) err = err + 1
    
end function


integer function test_geometry_areas() bind(C) result(err)
    use geometry_mod
    use cell_map_mod
    implicit none 

    type(side_map_t) :: side_map
    type(side_t), dimension(:), allocatable :: sides, path, contour
    type(triangle_t), dimension(:), allocatable :: triangles
    type(triangle_t) :: t1,t2,t3
    type(coord_t) :: c1, c2, c3, c4, c5, c6, c7
    integer, dimension(3) :: cell
    type(coord_t) :: init, end
    real :: area

    err = 0

    c1 = coord_t(position = [1.0,0.0,0.0],   id = 1)
    c2 = coord_t(position = [0.0,1.0,0.0],  id=  2)
    c3 = coord_t(position = [0.0,0.0,1.0], id=  3)
    
    allocate(triangles(1))
    triangles(1) = triangle_t(vertices = [c1,c2,c3])
    cell = triangles(1)%getCell()
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    if (size(sides) /= 3) err = err + 1

    path = getPathOnFace(getSidesOnFace(sides, FACE_Z))
    contour = buildSidesContour(path)
    area = contourArea(contour)
    if (area /= 0.5) err = err + 1

    call side_map%unset(key(cell))
    triangles(1)%vertices(1)%position = [0.5,0.0,0.0]
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    if (size(sides) /= 3) err = err + 1

    path = getPathOnFace(getSidesOnFace(sides, FACE_Y))
    contour = buildSidesContour(path)
    area = contourArea(contour)
    if (area /= 0.25) err = err + 1
    
    call side_map%unset(key(cell))
    triangles(1)%vertices(2)%position = [0.0,0.25,0.0]
    call buildSideMap(side_map, triangles)
    sides = side_map%getSidesInCell(cell)
    if (size(sides) /= 3) err = err + 1

    path = getPathOnFace(getSidesOnFace(sides, FACE_X))
    contour = buildSidesContour(path)
    area = contourArea(contour)
    if (area /= 0.125) err = err + 1


end function