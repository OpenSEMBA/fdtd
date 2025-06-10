module probes_mod

    use mtln_types_mod, only: PROBE_TYPE_CURRENT, PROBE_TYPE_VOLTAGE
#ifdef CompileWithMPI
    use FDETYPES, only: SUBCOMM_MPI
#endif

    implicit none

    type, public :: probe_t
        integer :: type
        real, allocatable, dimension(:) :: t
        real, allocatable, dimension(:,:) :: val
        real :: dt
        integer :: index, current_frame
        character (len=:), allocatable :: name
        logical :: in_layer = .true.

    contains
        procedure :: resizeFrames
        procedure :: update
        procedure :: saveFrame


    end type probe_t

    interface probe_t
        module procedure probeCtor
    end interface

contains

    function probeCtor(index, probe_type, dt, name, position, layer_indices) result(res)
        type(probe_t) :: res
        integer, intent(in) :: index
        integer, intent(in) :: probe_type
        real, intent(in) :: dt
        real, dimension(3) :: position
        character (len=:), allocatable :: name
        integer (kind=4), dimension(:,:), intent(in), optional :: layer_indices
        integer :: i, slice
#ifdef CompileWithMPI
        integer :: layer_index, ierr, sizeof
#endif

        res%type = probe_type
        res%index = index
        res%dt = dt
        res%current_frame = 1
        
#ifdef CompileWithMPI
        call MPI_COMM_SIZE(SUBCOMM_MPI, sizeof, ierr)
        if (sizeof > 1) then
            res%in_layer = .false.
            do i = 1, size(layer_indices,1) 
                if (index >= layer_indices(i, 1) .and. index <= layer_indices(i,2)+1) then 
                    res%in_layer = .true.
                    slice = i
                end if
            end do

            layer_index = 0
            if (res%in_layer) then 
                do i = 1, slice - 1
                    layer_index = layer_index + layer_indices(i,2) + 1 - (layer_indices(i,1) - 1)
                end do
                layer_index = layer_index + res%index - layer_indices(i,1) + 1
            end if
            res%index = layer_index
        endif
#endif
        res%name = res%name//name//"_"
        if (probe_type == PROBE_TYPE_VOLTAGE) then
            res%name = res%name//"V"
        else if (probe_type == PROBE_TYPE_CURRENT) then
            res%name = res%name//"I"
        else
            error stop 'Undefined probe'
        end if  
        block
            character(20) :: a,b,c
            write(a, *) int(position(1))
            write(b, *) int(position(2))
            write(c, *) int(position(3))
            res%name = res%name//"_"//trim(adjustl(a))//"_"//trim(adjustl(b))//"_"//trim(adjustl(c))
        end block
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