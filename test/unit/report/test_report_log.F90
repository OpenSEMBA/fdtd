! Unit tests for the leveled log API in Report_m

integer function test_log_line_format() bind(C, name="test_log_line_format") result(status)
    use Report_m
    implicit none

    status = 0
    call set_log_step(-1_4)
    call expect_line(log_line(LOG_WARN, '', 'plain message'), 'WARN: plain message')
    call expect_line(log_line(LOG_INFO, 'MTLN', 'with origin'), '[MTLN] INFO: with origin')

    call set_log_step(42_4)
    call expect_line(log_line(LOG_ERROR, 'WIRES', 'in loop'), '[n=42] [WIRES] ERROR: in loop')
    call expect_line(log_line(LOG_WARN, '', 'no origin'), '[n=42] WARN: no origin')
    call set_log_step(-1_4)

contains

    subroutine expect_line(actual, expected)
        character(len=*), intent(in) :: actual, expected
        if (trim(actual) /= expected) then
            print *, 'test_log_line_format FAILED: got "', trim(actual), '" expected "', expected, '"'
            status = 1
        end if
    end subroutine expect_line

end function test_log_line_format

integer function test_log_msg_sets_fatal_error() bind(C, name="test_log_msg_sets_fatal_error") result(status)
    use Report_m
    implicit none

    status = 0
    call resetFatalError()
    call log_msg(LOG_INFO, 'TEST', 'info does not flag')
    call log_msg(LOG_WARN, 'TEST', 'warning does not flag')
    if (isFatalError()) then
        print *, 'test_log_msg_sets_fatal_error FAILED: INFO/WARN flagged a fatal error'
        status = 1
    end if

    call log_msg(LOG_ERROR, 'TEST', 'error flags')
    if (.not. isFatalError()) then
        print *, 'test_log_msg_sets_fatal_error FAILED: ERROR did not flag a fatal error'
        status = 1
    end if
    call resetFatalError()

    call resetFatalError()
end function test_log_msg_sets_fatal_error
