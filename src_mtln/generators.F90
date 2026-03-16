module generators_m

    
    use mtln_types_m, only: SOURCE_TYPE_UNDEFINED, SOURCE_TYPE_VOLTAGE, SOURCE_TYPE_CURRENT   

    implicit none

    type generator_t
        integer :: index, conductor
        real, dimension(:), allocatable :: time, value
        real :: resistance
        integer :: source_type = SOURCE_TYPE_UNDEFINED
    contains
        procedure :: initGenerator
        procedure :: interpolate
    end type generator_t

    interface generator_t
        module procedure generatorCtor
    end interface    


contains

    function generatorCtor(index, conductor, gen_type, resistance, path) result(res)
        type(generator_t) :: res
        integer, intent(in) :: index, conductor, gen_type
        real :: resistance
        character(*), intent(in) :: path

        res%index = index
        res%conductor = conductor
        res%resistance = resistance
        res%source_type = gen_type
        call res%initGenerator(path)
    end function

    subroutine initGenerator(this, path)
        class(generator_t) :: this
        character(*), intent(in) :: path
        real :: time, value
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
        real :: res
        real :: t, x1, x2, y1, y2
        integer :: index
        real, dimension(:), allocatable :: timediff
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