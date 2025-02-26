module conformal_types_mod

    use cells_mod

    implicit none
    
    integer, parameter :: FACE_X = 1
    integer, parameter :: FACE_Y = 2
    integer, parameter :: FACE_Z = 3
    integer, parameter :: NOT_ON_FACE = -1
 
    type, public :: coord_t
        real, dimension(3) :: position
        integer :: id
    end type

    type, public :: triangle_t
        ! real, dimension(3) :: vertices(3)
        type(coord_t), dimension(3) :: vertices
        ! private
        real, dimension(3) :: normal
        integer :: face
    contains
        procedure :: getNormal
        procedure :: getFace
        procedure :: centroid
        procedure :: cell
        procedure :: computeNormal
        procedure :: assignFace
        procedure, private :: isOnAnyFace
        procedure, private :: isOnFace
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

    function getNormal(this) result(res)
        class(triangle_t) :: this
        real, dimension(3) :: res
        res = this%normal
    end function

    function getFace(this) result(res)
        class(triangle_t) :: this
        integer :: res
        res = this%face
    end function

    subroutine assignFace(this)
        class(triangle_t) :: this
        integer :: face
        this%face = NOT_ON_FACE
        if (this%isOnAnyFace()) then 
            do face = FACE_X, FACE_Z
                if (this%isOnFace(face)) this%face = face
            end do
        end if
    end subroutine

    function centroid(this) result(res)
        class(triangle_t) :: this
        real, dimension(3) :: res
        integer :: i
        res(:) = 0.0
        do i = 1,3
            res = res + this%vertices(i)%position
        end do
        res = res/3.0
    end function
    
    function cell(this) result(res)
        class(triangle_t) :: this
        integer(kind=4) :: res(3)
        res = floor(this%centroid())
    end function    

    subroutine computeNormal(this)
        class(triangle_t) :: this
        real, dimension(3) :: v1,v2, cross
        v1 = this%vertices(2)%position - this%vertices(1)%position
        v2 = this%vertices(3)%position - this%vertices(2)%position
        cross = [ v1(2)*v2(3)-v1(3)*v2(2), &
                 -(v1(1)*v2(3)-v1(3)*v2(1)), &
                   v1(1)*v2(2)-v1(2)*v2(1)]
        this%normal = cross/norm2(cross)
    end subroutine

    logical function isOnAnyFace(this)
        class(triangle_t) :: this
        isOnAnyFace = this%isOnFace(FACE_X) .or. &
                      this%isOnFace(FACE_Y) .or. &
                      this%isOnFace(FACE_Z)
    end function
    
    logical function isOnFace(this, dir)
        class(triangle_t) :: this
        integer :: dir
        isOnFace = abs(this%normal(dir)) == 1 .and. &
                   abs(this%normal(mod(dir,3) + 1)) == 0 .and. &
                   abs(this%normal(mod(dir+1,3) + 1)) == 0
    end function

end module