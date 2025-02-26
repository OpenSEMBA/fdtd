integer function test_triangle_normal() bind(C) result(err)

    use conformal_types_mod
    implicit none
    type(triangle_t) :: t
    real, dimension(3) :: normal
    integer :: face
    err = 0

    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [0,1,0]
    call t%buildTriangle()
    if (t%getFace() /= FACE_X) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [1,0,0]
    call t%buildTriangle()
    if (t%getFace() /= FACE_Y) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,0]
    call t%buildTriangle()
    if (t%getFace() /= FACE_Z) err = err + 1

    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,1]
    call t%buildTriangle()
    if (t%getFace() /= NOT_ON_FACE) err = err + 1


end function

integer function test_triangle_edges() bind(C) result(err)

    use conformal_types_mod
    implicit none
    type(triangle_t) :: t
    real, dimension(3) :: normal
    type(side_t), dimension(3) :: sides
    err = 0
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [0,1,0]
    call t%buildTriangle()
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Z) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_Y) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [1,0,0]
    call t%buildTriangle()
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Z) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_X) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,0]
    call t%buildTriangle()
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Y) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_X) err = err + 1

    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,1]
    call t%buildTriangle()
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Y) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= NOT_ON_EDGE) err = err + 1


end function