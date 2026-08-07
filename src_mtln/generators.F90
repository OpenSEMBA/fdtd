module generators_m

    
    use mtln_types_m, only: SOURCE_TYPE_UNDEFINED, SOURCE_TYPE_VOLTAGE, SOURCE_TYPE_CURRENT   
#ifdef CompileWithMPI
    use FDETYPES_m, only: SUBCOMM_MPI
#endif
    use FDETYPES_m, only: RKIND, RKIND_tiempo

    implicit none

    type generator_t
        integer :: index, conductor
        real(kind=rkind), dimension(:), allocatable :: value
        real(kind=rkind_tiempo), dimension(:), allocatable :: time
        real :: resistance
        integer :: source_type = SOURCE_TYPE_UNDEFINED
        logical :: in_layer = .true.
    contains
        procedure :: initGenerator
        procedure :: interpolate
    end type generator_t

    interface generator_t
        module procedure generatorCtor
    end interface    


contains

    function generatorCtor(index, conductor, gen_type, resistance, path, layer_indices) result(res)
        type(generator_t) :: res
        integer, intent(in) :: index, conductor, gen_type
        real(kind=rkind) :: resistance
        character(*), intent(in) :: path
        integer(kind=4), dimension(:,:), intent(in), optional :: layer_indices
#ifdef CompileWithMPI
        integer :: layer_index, ierr, sizeof, i, slice
#endif

        res%index = index
        res%conductor = conductor
        res%resistance = resistance
        res%source_type = gen_type
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
        if (res%in_layer) call res%initGenerator(path)
    end function

    subroutine initGenerator(this, path)
        class(generator_t) :: this
        character(*), intent(in) :: path
        real(kind=rkind) :: value
        real(kind=rkind_tiempo) :: time
        integer :: io, line_count, i
        
        if (path == "" ) then 
            allocate(this%time(0), this%value(0))
            ! error
        end if
        ! First pass: count the number of lines
        line_count = 0
        open(unit = 1, file = path)
        do
            read(1, *, iostat = io) time, value
            if (io /= 0) exit
            line_count = line_count + 1
        end do
        close(1)
        
        ! Allocate arrays with the exact size needed
        allocate(this%time(line_count))
        allocate(this%value(line_count))
        
        ! Second pass: fill the arrays
        open(unit = 1, file = path)
        do i = 1, line_count
            read(1, *, iostat = io) this%time(i), this%value(i)
            if (io /= 0) then
                ! Handle unexpected read error
                close(1)
                error stop "Error reading excitation file"
            end if
        end do
        close(1)
    end subroutine

    function interpolate(this, t) result(res)
        class(generator_t) :: this
        real(kind=rkind) :: res
        real(kind=rkind_tiempo) :: t, x1, x2
        real(kind=rkind) :: y1, y2
        integer :: index
        real(kind=rkind_tiempo), dimension(:), allocatable :: timediff
        timediff = this%time - t
        index = maxloc(timediff, 1, (timediff) <= 0)
        if (index == 0) index = 1
        x1 = this%time(index)
        y1 = this%value(index)
        if (index+1 > size(this%time)) then
            x2 = x1
            y2 = y1
        else 
            x2 = this%time(index+1)
            y2 = this%value(index+1)
        end if
        res = (t*(y2-y1) + x2*y1 - x1*y2)/(x2-x1)
    end function

end module generators_m