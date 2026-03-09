module ngspice_interface_mod

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

        subroutine circ(input) bind(C, name = "circ")
            use iso_c_binding, only: c_ptr
            type(c_ptr), intent(in) :: input(*)
        end subroutine

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

    end interface


end module