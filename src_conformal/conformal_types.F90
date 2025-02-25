module conformal_types_mod
    use cells_mod
    
    implicit none
    
    type, public :: triangle_t
        real, dimension(3) :: coords
        
    contains
        procedure :: face
        procedure, private :: centroid
        procedure :: cell
        procedure :: normal
        procedure :: isOnFaceDir
        procedure :: isOnAnyFace
    end type
        
    ! type, public :: element_t
    !     integer, dimension(3) :: x, y, z
    ! end type

    ! type, public :: area_t
    !     integer, dimension(3) :: x, y, z
    ! end type

    type, public :: conformal_region_t
        type(triangle_t), dimension(:), allocatable :: triangles
        type(cell_region_t) :: quads
        ! type(element_t), dimension(:), allocatable :: elements
        ! type(area_t), dimension(:), allocatable :: areas
    end type

contains

    integer function face(this)
        class(triangle_t) :: this
        integer :: i
        if (this%isOnAnyFace()) then 
            do i = 1, 3
                if (this%isOnFaceDir(i)) face = i
            end do
        end if
        face = -1
    end function

    function normal(this) result(res)
        class(triangle_t) :: this
        real, dimension(3) :: v1,v2, cross
        real, dimension(3) :: res
        v1 = this%coords(2) - this%coords(1)
        v2 = this%coords(3) - this%coords(2)
        cross = [ v1(2)*v2(3)-v1(3)*v2(2), &
                 -(v1(1)*v2(3)-v1(3)*v2(1)), &
                   v1(1)*v2(2)-v1(2)*v2(1)]
        res = cross/norm2(cross)
    end function

    logical function isOnAnyFace(this)
        class(triangle_t) :: this
        real, dimension(3) :: n
        n = this%normal()
        isOnAnyFace = this%isOnFaceDir(1) .or. &
                      this%isOnFaceDir(2) .or. &
                      this%isOnFaceDir(3)
    end function

    logical function isOnFaceDir(this, dir)
        class(triangle_t) :: this
        integer :: dir
        isOnFaceDir = this%coords(dir) == 1 .and. &
                      this%coords(mod(dir,3) + 1) == 0 .and. &
                      this%coords(mod(dir+1,3) + 1) == 0
    end function

    function centroid(this) result(res)
        class(triangle_t) :: this
        real :: res(3)
        integer :: i
        res(:) = 0.0
        do i = 1,3
            res(i) = res(i) + this%coords(i)/3
        end do
    end function
    
    function cell(this) result(res)
        class(triangle_t) :: this
        integer(kind=4) :: res(3)
        res = floor(this%centroid())
    end function    

end module