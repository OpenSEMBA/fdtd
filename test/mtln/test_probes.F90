integer function test_probes_number_of_frames() bind(C) result(error_cnt)
    use probes_m
    use mtln_types_m, only: PROBE_TYPE_VOLTAGE, PROBE_TYPE_CURRENT
    implicit none

    type(probe_t) :: probe
    integer, parameter :: num_frames = 20
    integer, parameter :: num_conductors = 2

    real, dimension(3) :: position = [0.0, 0.0, 0.0]
    character(len=:), allocatable :: name

    error_cnt = 0
    name = "test_probe"

    probe = probe_t(index = 1, probe_type = PROBE_TYPE_VOLTAGE, &
                    dt = 1.0e-9, name = name, position = position)
    call probe%resizeFrames(num_frames, num_conductors)

    ! After resizeFrames, probe arrays must have exactly num_frames entries.
    ! Previously they were allocated with num_frames+1, leaving the last slot
    ! zeroed and causing a spurious (0,0) data point in output (issue #268).
    if (size(probe%t) /= num_frames) then
        error_cnt = error_cnt + 1
    end if
    if (size(probe%val, 1) /= num_frames) then
        error_cnt = error_cnt + 1
    end if

end function
