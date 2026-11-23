module ngspice_interface_m

    use iso_c_binding
    implicit none

    type, bind(c) :: vectorInfo_t
        type(c_ptr) :: vName
        integer(c_int) :: vType
        integer(c_short) :: vFlags
        type(c_ptr) :: vRealData !real
        type(c_ptr) :: vCompData !ngcomplex
        integer(c_int) :: vLength
    end type

    interface
        subroutine command(input) bind (C, name = "command")
            use iso_c_binding, only: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: input
        end subroutine

        integer(c_int) function circ(input) bind(C, name = "circ")
            use iso_c_binding, only: c_ptr, c_int
            type(c_ptr), intent(in) :: input(*)
        end function

        subroutine start() bind (C, name = "start")
            use iso_c_binding, only: c_int, c_ptr
        end subroutine

        type(c_ptr) function get_vector_info(name) bind (C, name="get_vector_info")
            use iso_c_binding, only: c_char, c_ptr
            character(kind=c_char), dimension(*), intent(in) :: name
        end function
            
        type(c_ptr) function get_all_plots() bind (C, name="get_all_plots")
            use iso_c_binding, only: c_ptr
        end function

        integer(c_int) function has_error() bind (C, name="has_error")
            use iso_c_binding, only: c_int
        end function

        integer(c_int) function codemodels_report_text(buffer, size) bind (C, name="codemodels_report_text")
            use iso_c_binding, only: c_char, c_int
            character(kind=c_char), dimension(*), intent(out) :: buffer
            integer(c_int), value, intent(in) :: size
        end function
       
    end interface
end module