module probes_mod

    use mtln_types_mod, only: PROBE_TYPE_CURRENT, PROBE_TYPE_VOLTAGE

    implicit none

    type, public :: probe_t
        integer :: type
        real, allocatable, dimension(:) :: t
        real, allocatable, dimension(:,:) :: val
        real :: dt
        integer :: index, current_frame
        character (len=:), allocatable :: name

    contains
        procedure :: resizeFrames
        procedure :: update
        procedure :: saveFrame


    end type probe_t

    interface probe_t
        module procedure probeCtor
    end interface

contains

    function probeCtor(index, probe_type, dt, name, position) result(res)
        type(probe_t) :: res
        integer, intent(in) :: index
        integer, intent(in) :: probe_type
        real, intent(in) :: dt
        real, dimension(3), optional :: position
        character (len=:), allocatable, optional :: name

        res%type = probe_type
        res%index = index
        res%dt = dt
        res%current_frame = 1

        if (present(name)) then
            res%name = res%name//name//"_"
        end if
        if (probe_type == PROBE_TYPE_VOLTAGE) then
            res%name = res%name//"V"
        else if (probe_type == PROBE_TYPE_CURRENT) then
            res%name = res%name//"I"
        else
            error stop 'Undefined probe'
        end if  
        if (present(position)) then
            block
                character(20) :: a,b,c
                write(a, *) int(position(1))
                write(b, *) int(position(2))
                write(c, *) int(position(3))
                res%name = res%name//"_"//trim(adjustl(a))//"_"//trim(adjustl(b))//"_"//trim(adjustl(c))
                end block
        end if
        write(*,*) 'probe name: ', res%name
        end function

    subroutine resizeFrames(this, num_frames, number_of_conductors)
        class(probe_t) :: this
        integer, intent(in) :: num_frames, number_of_conductors

        allocate(this%t(num_frames + 1))
        allocate(this%val(num_frames + 1, number_of_conductors))
        this%t = 0.0
        this%val = 0.0

    end subroutine

    subroutine update(this, t, v, i)
        class(probe_t) :: this
        real, intent(in) :: t
        real, dimension(:,:), intent(in) :: v
        real, dimension(:,:), intent(in) :: i
        
        if (this%type == PROBE_TYPE_VOLTAGE) then
            call this%saveFrame(t, v(:,this%index))
        else if (this%type == PROBE_TYPE_CURRENT) then
            if (this%index == size(i,2) + 1) then
                call this%saveFrame(t + 0.5*this%dt, i(:,this%index - 1))
            else 
                call this%saveFrame( t+ 0.5*this%dt, i(:,this%index))
            endif
        end if  

    end subroutine

    subroutine saveFrame(this, time, values)
        class(probe_t) :: this
        real, intent(in) :: time
        real, intent(in), dimension(:) :: values

        this%t(this%current_frame) = time
        this%val(this%current_frame,:) = values
        this%current_frame = this%current_frame + 1

    end subroutine


end module probes_mod