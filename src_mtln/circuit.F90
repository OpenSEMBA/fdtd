module circuit_m

    use ngspice_interface_m
    use mtln_types_m, only: node_source_t, SOURCE_TYPE_CURRENT, SOURCE_TYPE_VOLTAGE
    use Report_m, only: WarnErrReport
    use FDETYPES_m, only: RKIND, RKIND_TIEMPO, SINGLE
    implicit none

    type string_t
        character(len=256) :: name
        integer :: length
    end type string_t

    type source_t
        logical :: has_source = .false.
        real(kind=RKIND_TIEMPO), dimension(:), allocatable :: time
        real(kind=RKIND), dimension(:), allocatable :: value
        integer :: source_type
    contains 
        procedure :: interpolate
    end type

    type VI_t
        real(kind=RKIND) :: voltage
        real(kind=RKIND) :: current
        real(kind=RKIND_TIEMPO) :: time
    end type

    type nodes_t
        type(VI_T), allocatable :: values(:)
        type(source_t), allocatable :: sources(:)
        type(string_t), allocatable :: names(:)
    end type nodes_t

    type, public :: circuit_t
        character(len=:), allocatable :: name
        real(kind=RKIND_TIEMPO) :: time = 0.0, dt = 0.0
        logical :: errorFlag = .false.
        type(nodes_t) :: nodes, saved_nodes   

    contains
        procedure :: init
        procedure :: run
        procedure :: step
        procedure, private :: resume
        procedure :: quit
        procedure, private :: loadNetlist
        procedure :: readInput
        procedure :: setStopTimes
        procedure :: setModStopTimes
        procedure :: getNodeVoltage
        procedure :: getNodeCurrent
        procedure :: updateNodes
        procedure :: getTime
        procedure :: updateNodeCurrent
        procedure :: updateCircuitSources
        procedure :: modifyLineCapacitorValue

        procedure :: printCWD

    end type circuit_t

contains

    real(kind=rkind) function interpolate(this, time, dt) result(res)
        class(source_t) :: this
        real(kind=RKIND_TIEMPO) :: time, dt
        real(kind=RKIND) :: x1,x2, y1, y2
        integer :: index
        real(kind=rkind), dimension(:), allocatable :: timediff
        timediff = this%time - time + dt
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
                
        res = (time*(y2-y1) + x2*y1 - x1*y2)/(x2-x1)
    end function

    subroutine printCWD(this)
        class(circuit_t) :: this
        call command('getcwd' // c_null_char)
    end subroutine

    subroutine init(this, names, sources, netlist)
        class(circuit_t) :: this
        type(string_t), intent(in), dimension(:), optional :: names
        type(node_source_t), intent(in), dimension(:), optional :: sources
        character(len=*), intent(in), optional :: netlist
        integer :: i

        call start()
        if (present(netlist)) then
            write(*,*) 'load netlist'
            call this%loadNetlist(netlist)
        end if

        if (.not.present(names)) then 
            error stop 'Missing node names'
        end if

        allocate(this%nodes%names(size(names)))
        allocate(this%nodes%values(size(names)))

        allocate(this%nodes%sources(size(names)))
        do i = 1, size(names)
            this%nodes%names(i) = names(i)
        end do
        if (present(sources)) then 
            do i = 1, size(sources)
                this%nodes%sources(i) = setSource(sources(i)%path_to_excitation)
                this%nodes%sources(i)%source_type = sources(i)%source_type
            end do
        end if

    end subroutine

    type(source_t) function setSource(source_path) result(res)
        character(*), intent(in) :: source_path
        real(kind=RKIND_TIEMPO) :: time
        real(kind=RKIND) ::value
        integer :: io, line_count, i
        
        if (source_path == "" ) then 
            allocate(res%time(0), res%value(0))
            res%has_source = .false.
            return
        end if
        
        res%has_source = .true.
        
        ! First pass: count the number of lines
        line_count = 0
        open(unit = 1, file = source_path, iostat = io)
        if (io /= 0) then
            call WarnErrReport('Cannot open excitation file: ' // trim(source_path), .true.)
            allocate(res%time(0), res%value(0))
            res%has_source = .false.
            return
        end if
        do
            read(1, *, iostat = io) time, value
            if (io < 0) exit
            if (io > 0) then
                close(1)
                call WarnErrReport('Error reading excitation file: ' // trim(source_path), .true.)
                allocate(res%time(0), res%value(0))
                res%has_source = .false.
                return
            end if
            line_count = line_count + 1
        end do
        close(1)
        
        ! Allocate arrays with the exact size needed
        allocate(res%time(line_count))
        allocate(res%value(line_count))
        
        ! Second pass: fill the arrays (file was verified readable in first pass)
        open(unit = 1, file = source_path)
        do i = 1, line_count
            read(1, *) res%time(i), res%value(i)
        end do
        close(1)
    end function    

    subroutine loadNetlist(this, netlist)
        class(circuit_t) :: this
        character(len=*, kind=c_char), intent(in) :: netlist
        call command(c_char_'source ' // trim(netlist) // c_null_char)
    end subroutine

    subroutine step(this)
        class(circuit_t) :: this
        call this%updateCircuitSources(this%time)
        if (this%time == 0) then
            call this%run()
        else
            call this%resume()
        end if
        call this%updateNodes()

    end subroutine

    subroutine run(this)
        class(circuit_t) :: this
        call command('run ' // c_null_char)
    end subroutine



    subroutine setStopTimes(this, finalTime, dt)
        class(circuit_t) :: this
        real(kind=RKIND_TIEMPO), intent(in) :: finalTime, dt
        character(50) :: charTime
        real(kind=rkind) :: time

        time = 0.0_rkind
        do while (time < finalTime)
            time = time + dt
            write(charTime, *) time
            call command('stop when time = '//charTime // c_null_char)
        end do
    end subroutine

    subroutine setModStopTimes(this, dt)
        class(circuit_t) :: this
        real(kind=RKIND_TIEMPO), intent(in) :: dt
        character(50) :: charTime
        write(charTime, *) real(dt, SINGLE)
        call command('stop when time mod '//charTime // c_null_char)
    end subroutine

    subroutine resume(this)
        class(circuit_t) :: this
        call command('resume ' // c_null_char)
    end subroutine

    subroutine quit(this)
        class(circuit_t) :: this
        call command('quit 0' // c_null_char)
    end subroutine

    subroutine readInput(this, input, printInput) 
        class(circuit_t) :: this
        character(*), intent(in) :: input(:)
        logical, optional :: printInput
        type(c_ptr) :: argv_c(size(input))
        integer :: i   

        type c_string_t
            character(len=:,kind=c_char), allocatable :: item
        end type c_string_t
        type(c_string_t), target :: tmp(size(input))

        if (present(printInput)) then
            if (printInput .eqv. .true.) then 
                do i = 1 , size(input)
                    write(*,*) input(i)
                end do
            end if
        end if

        do i = 1, size(input)
            tmp(i)%item = trim(input(i)) // c_null_char
            argv_c(i) = c_loc(tmp(i)%item)
        end do
        call circ(argv_c)
    end subroutine

    function getName(cName) result(res)
        type(c_ptr) :: cName
        type(string_t) :: res
        character, pointer :: f_output(:) => null()
        integer :: i
        res%name = ""
        res%length = 0
        call c_f_pointer(cName, f_output,[100])
        do i = 1,100
            if (f_output(i) == c_null_char) exit
            res%name(i:i) = f_output(i)
        end do
        res%length = i-1

    end function

    subroutine updateCircuitSources(this, time)
        class(circuit_t) :: this
        real(kind=RKIND_TIEMPO), intent(in) :: time
        real(kind=RKIND) :: interp
        character(50) :: source_value
        integer :: i, index
        do i = 1, size(this%nodes%sources)
            if (this%nodes%sources(i)%has_source) then
                if (this%nodes%sources(i)%source_type == SOURCE_TYPE_VOLTAGE) then 
                    interp = this%nodes%sources(i)%interpolate(time, 0.0_RKIND_TIEMPO) 
                    write(source_value, *) interp
                    call command("alter @V"//trim(this%nodes%names(i)%name)//"_s[dc] = "//trim(source_value) // c_null_char)
                else if (this%nodes%sources(i)%source_type == SOURCE_TYPE_CURRENT) then 
                    interp = this%nodes%sources(i)%interpolate(time, 0.0_RKIND_TIEMPO) 
                    write(source_value, *) interp
                    call command("alter @I"//trim(this%nodes%names(i)%name)//"_s[dc] = "//trim(source_value) // c_null_char)
                end if
            end if
        end do
    end subroutine

    subroutine modifyLineCapacitorValue(this, name, c)
        class(circuit_t) :: this
        character(*), intent(in) :: name
        real(kind=rkind), intent(in) :: c
        character(50) :: sC

        write(sC, *) c
        call command("alter @CL"//trim(name)//" = "//trim(sC) // c_null_char)

    end subroutine

    subroutine updateNodeCurrent(this, node_name, current)
        class(circuit_t) :: this
        real(kind=rkind) :: current
        character(50) :: sCurrent
        character(*) :: node_name
        if (index(node_name, "initial") /= 0) then
            write(sCurrent, *) current
        else if (index(node_name, "end") /= 0) then
            write(sCurrent, *) -current
        end if
        call command("alter @I"//trim(node_name)//"[dc] = "//trim(sCurrent) // c_null_char)
    end subroutine

    subroutine updateNodes(this) 
        class(circuit_t) :: this
        integer :: i
        type(vectorInfo_t), pointer :: info
        real(kind=c_double), pointer :: values(:)
        do i = 1, size(this%nodes%names)
            call c_f_pointer(get_vector_info(trim(this%nodes%names(i)%name)//c_null_char), info)
            call c_f_pointer(info%vRealData, values,shape=[info%vLength])
            if (this%nodes%names(i)%name /= "time") then 
                this%nodes%values(i)%voltage = values(ubound(values,1))
            else 
                this%nodes%values(i)%time = values(ubound(values,1))
            end if
        end do
    end subroutine

    function getNodeVoltage(this, name) result(res)
        class(circuit_t) :: this
        character(len=*), intent(in) :: name
        real(kind=rkind) :: res
        res = this%nodes%values(findVoltageIndexByName(this%nodes%names, name))%voltage
    end function

    function getNodeCurrent(this, name) result(res)
        class(circuit_t) :: this
        character(len=*), intent(in) :: name
        real(kind=rkind) :: res
        res = this%nodes%values(findVoltageIndexByName(this%nodes%names, name))%current
    end function

    function getTime(this) result(res)
        class(circuit_t) :: this
        real(kind=rkind_tiempo) :: res
        res = this%nodes%values(findIndexByName(this%nodes%names, "time"))%time
    end function

    function findIndexByName(names, name) result(res)
        type(string_t) :: names(:)
        character(len=*), intent(in) :: name
        integer :: res, i
        res = 0
        do i = 1, size(names)
            if ( names(i)%name(1:names(i)%length) == trim(name)) then 
                res = i
                exit
            end if
        end do
    end function    

    function findVoltageIndexByName(names, name) result(res)
        type(string_t) :: names(:)
        character(len=*), intent(in) :: name
        integer :: res, i
        res = 0
        do i = 1, size(names)
            if ( names(i)%name(1:names(i)%length) == 'V('//trim(name)//')') then 
                res = i
                exit
            else if ( names(i)%name(1:names(i)%length) == trim(name)) then 
                res = i
                exit
            end if
        end do
    end function    

end module 