module conformal_types_mod
    implicit none
    type, public :: conformal_t
        type(triangle_t), dimension(:) :: triangles


    end type

    type, public :: triangle_t
        real(kind=RKIND), dimension(3) :: coords
    contains
        procedure :: face
        procedure, private :: centroid
        procedure :: cell
        procedure :: normal
        procedure :: isOnFace
    end type


contains

    function normal(this) result(res)
        class(triangle_t) :: this
        real(kind=rkind) :: v1,v2, cross
        real(kind=rkind) :: res
        v1 = this%coords(2) - this%coords(1)
        v2 = this%coords(3) - this%coords(2)
        cross = [ (v1(2)*v2(3)-v1(3)*v2(2), &
                 -(v1(1)*v2(3)-v(3)*v2(1)), &
                   v1(1)*v2(2)-v1(2)*v2(1))]
        res = cross/norm2(cross)
    end function

    logical function isOnFace(this)
        class(triangle_t) :: this
        real(kind=rkind) :: normal
        normal = this%normal()
        
    end function

    function centroid(this) result(res)
        class(triangle_t) :: this
        real(kind=rkind) :: res(3)
        integer :: i
        res(:) = 0.0
        do i = 1,3
            res(i) = res(i) + this%coords(i)/3
        end do
    end function
    
    function cell(this) result(res)
        class(triangle_t) :: this
        integer(kind=4) :: res(3)
        real(kind=rkind) :: centroid(3)
        centroid = this%getCentrid()
        res = floor(centroid)

