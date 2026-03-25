module probes_m

    use mtln_types_m, only: PROBE_TYPE_CURRENT, PROBE_TYPE_VOLTAGE
#ifdef CompileWithMPI
    use FDETYPES_m, only: SUBCOMM_MPI, RKIND
#else
    use FDETYPES_m, only: RKIND
#endif
    use FDETYPES_m, only: RKIND, RKIND_TIEMPO

    implicit none

    type, public :: probe_t
        integer :: type
        real(kind=RKIND), allocatable, dimension(:) :: t
        real(kind=RKIND), allocatable, dimension(:,:) :: val
        real(kind=RKIND_TIEMPO) :: dt
        integer :: index, current_frame
        character(len=:), allocatable :: name
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
        real(kind=RKIND_TIEMPO), intent(in) :: dt
        real(kind=RKIND), dimension(3) :: position
        character(len=:), allocatable :: name
        integer(kind=4), dimension(:,:), intent(in), optional :: layer_indices
        integer :: i, slice
#ifdef CompileWithMPI
        integer :: layer_index, ierr, sizeof
#endif

        res%type = probe_type
        res%index = index
        res%dt = dt
        res%current_frame = 1
        
#ifdef CompileWithMPI
        if (present(layer_indices)) then
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
            end if
        end if
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

        allocate(this%t(num_frames))
        allocate(this%val(num_frames, number_of_conductors))
        this%t = 0.0
        this%val = 0.0

    end subroutine

    subroutine update(this, t, v, i)
        class(probe_t) :: this
        real(kind=RKIND_TIEMPO), intent(in) :: t
        real(kind=RKIND), dimension(:,:), intent(in) :: v
        real(kind=RKIND), dimension(:,:), intent(in) :: i
        
        if (this%type == PROBE_TYPE_VOLTAGE) then
            call this%saveFrame(t, v(:,this%index))
        else if (this%type == PROBE_TYPE_CURRENT) then
            if (this%index == size(i,2) + 1) then
                call this%saveFrame(t + 0.5*this%dt, i(:,this%index - 1))
            else 
                call this%saveFrame( t+ 0.5*this%dt, i(:,this%index))
            end if
        end if  

    end subroutine

    subroutine saveFrame(this, time, values)
        class(probe_t) :: this
        real(kind=RKIND_TIEMPO), intent(in) :: time
        real(kind=RKIND), intent(in), dimension(:) :: values
        if (this%current_frame > size(this%t)) return
        this%t(this%current_frame) = time
        this%val(this%current_frame,:) = values
        this%current_frame = this%current_frame + 1

    end subroutine


end module probes_m