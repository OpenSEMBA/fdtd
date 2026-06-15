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

        integer(c_int) function has_error() bind (C, name="has_error")
            use iso_c_binding, only: c_int
        end function

        function curplot() bind(C, name="ngspice_cur_plot") result(res)
            use iso_c_binding
            type(c_ptr) :: res
        end function curplot

        function allvecs_current() bind(C, name="ngspice_all_vecs_current") result(res)
            use iso_c_binding
            type(c_ptr) :: res
        end function allvecs_current

        function allplots() bind(c, name="ngspice_all_plots") result(res)
            use iso_c_binding
            type(c_ptr) :: res
        end function allplots


    end interface
end module