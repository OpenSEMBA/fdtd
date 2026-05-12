module benchmark_m

    use FDETYPES_m, only: RKIND
    implicit none

    private
    public :: benchmark_t

    type benchmark_t
        real(kind=rkind) :: total_time
        real(kind=rkind) :: min_time
        real(kind=rkind) :: max_time
        real(kind=rkind) :: last_time
        integer :: count
        logical :: running
        real(kind=rkind) :: start_time
        character(len=128) :: name
    contains
        procedure :: benchmark_init
        procedure :: benchmark_start
        procedure :: benchmark_stop
        procedure :: benchmark_get_elapsed
        procedure :: benchmark_reset
        procedure :: benchmark_report
    end type benchmark_t

contains

    subroutine benchmark_init(this, name)
        class(benchmark_t) :: this
        character(*), intent(in), optional :: name

        this%total_time = 0.0_rkind
        this%min_time = huge(1.0_rkind)
        this%max_time = 0.0_rkind
        this%last_time = 0.0_rkind
        this%count = 0
        this%running = .false.
        this%start_time = 0.0_rkind
        if (present(name)) then
            this%name = trim(name)
        else
            this%name = 'benchmark'
        end if
    end subroutine benchmark_init

    subroutine benchmark_start(this)
        class(benchmark_t) :: this
        integer :: i, n

        call system_clock(count=i, count_rate=n)
        this%start_time = real(i, kind=rkind) / real(n, kind=rkind)
        this%running = .true.
    end subroutine benchmark_start

    subroutine benchmark_stop(this)
        class(benchmark_t) :: this
        integer :: i, n
        real(kind=rkind) :: elapsed

        call system_clock(count=i, count_rate=n)
        elapsed = real(i, kind=rkind) / real(n, kind=rkind) - this%start_time

        if (elapsed < 0.0_rkind) elapsed = 0.0_rkind

        this%last_time = elapsed
        this%total_time = this%total_time + elapsed
        this%count = this%count + 1

        if (elapsed < this%min_time) this%min_time = elapsed
        if (elapsed > this%max_time) this%max_time = elapsed

        this%running = .false.
    end subroutine benchmark_stop

    function benchmark_get_elapsed(this) result(elapsed)
        class(benchmark_t) :: this
        real(kind=rkind) :: elapsed
        integer :: i, n

        if (this%running) then
            call system_clock(count=i, count_rate=n)
            elapsed = real(i, kind=rkind) / real(n, kind=rkind) - this%start_time
            if (elapsed < 0.0_rkind) elapsed = 0.0_rkind
        else
            elapsed = this%last_time
        end if
    end function benchmark_get_elapsed

    subroutine benchmark_reset(this)
        class(benchmark_t) :: this

        this%total_time = 0.0_rkind
        this%min_time = huge(1.0_rkind)
        this%max_time = 0.0_rkind
        this%last_time = 0.0_rkind
        this%count = 0
        this%running = .false.
        this%start_time = 0.0_rkind
    end subroutine benchmark_reset

    subroutine benchmark_report(this, unit)
        class(benchmark_t) :: this
        integer, intent(in), optional :: unit
        integer :: u

        u = 6
        if (present(unit)) u = unit

        if (this%count > 0) then
            write(u, '(A, A, I5)') 'Benchmark: ', trim(this%name), this%count
            write(u, '(A, ES12.4)') '  Total time (s): ', this%total_time
            write(u, '(A, ES12.4)') '  Average time (s): ', this%total_time / this%count
            write(u, '(A, ES12.4)') '  Min time (s): ', this%min_time
            write(u, '(A, ES12.4)') '  Max time (s): ', this%max_time
            write(u, '(A, F10.2, A)') '  Calls per second: ', &
                 real(this%count, kind=rkind) / this%total_time, ' Hz'
        else
            write(u, '(A, A)') 'No measurements for benchmark: ', trim(this%name)
        end if
    end subroutine benchmark_report

end module benchmark_m
