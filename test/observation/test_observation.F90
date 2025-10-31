module test_observation_m
    use Observa

    implicit none

    contains
        subroutine check_shape_real(arr, n_expected, test_err, name)
            real(kind=8), intent(in), dimension(:,:) :: arr
            integer, intent(in) :: n_expected
            integer, intent(inout) :: test_err
            character(len=*), intent(in), optional :: name

            integer :: rank_arr
            integer, allocatable :: shp(:)
            character(len=:), allocatable :: nm

            nm = merge(name, "array", present(name))

            rank_arr = rank(arr)
            shp = shape(arr)

            ! Expect exactly two dimensions [1, n_expected]
            if (rank_arr /= 2 .or. any(shp /= [1, n_expected])) then
              test_err = test_err + 1
            end if
        end subroutine check_shape_real

        subroutine check_shape_complex(arr, n_expected, test_err, name)
            complex(kind=8), intent(in), dimension(:,:) :: arr
            integer, intent(in) :: n_expected
            integer, intent(inout) :: test_err
            character(len=*), intent(in), optional :: name

            integer :: rank_arr
            integer, allocatable :: shp(:)
            character(len=:), allocatable :: nm

            nm = merge(name, "array", present(name))

            rank_arr = rank(arr)
            shp = shape(arr)

            ! Expect exactly two dimensions [1, n_expected]
            if (rank_arr /= 2 .or. any(shp /= [1, n_expected])) then
              test_err = test_err + 1
            end if
        end subroutine check_shape_complex

        subroutine check_size(arr, n_expected, test_err, name)
            integer, intent(in), dimension(:) :: arr
            integer, intent(in) :: n_expected
            integer, intent(inout) :: test_err
            character(len=*), intent(in), optional :: name

            integer :: rank_arr
            integer :: siz
            character(len=:), allocatable :: nm

            nm = merge(name, "array", present(name))

            rank_arr = rank(arr)
            siz = size(arr)

            if (rank_arr /= 1 .or. siz /= n_expected) then
                test_err = test_err + 1
            end if
        end subroutine check_size

        integer function test_allocate_serialize_for_time_domain() bind(C) result(err)
            type(Serialized_t) :: serialize
            integer(kind=4) :: numberOfSerialized = 4
            integer :: test_err = 0

            call serialize%allocate_for_time_domain(numberOfSerialized)

            call check_shape_real(serialize%Valor, numberOfSerialized, test_err, "valor")
            call check_shape_real(serialize%Valor_x, numberOfSerialized, test_err, "Valor_x")
            call check_shape_real(serialize%Valor_y, numberOfSerialized, test_err, "Valor_y")
            call check_shape_real(serialize%Valor_z, numberOfSerialized, test_err, "Valor_z")
            call check_shape_real(serialize%ValorE, numberOfSerialized, test_err, "ValorE")
            call check_shape_real(serialize%Valor_Ex, numberOfSerialized, test_err, "Valor_Ex")
            call check_shape_real(serialize%Valor_Ey, numberOfSerialized, test_err, "Valor_Ey")
            call check_shape_real(serialize%Valor_Ez, numberOfSerialized, test_err, "Valor_Ez")
            call check_shape_real(serialize%ValorH, numberOfSerialized, test_err, "ValorH")
            call check_shape_real(serialize%Valor_Hx, numberOfSerialized, test_err, "Valor_Hx")
            call check_shape_real(serialize%Valor_Hy, numberOfSerialized, test_err, "Valor_Hy")
            call check_shape_real(serialize%Valor_Hz, numberOfSerialized, test_err, "Valor_Hz")

            call serialize%deallocate_for_time_domain()

            err = test_err
        end function test_allocate_serialize_for_time_domain

        integer function test_allocate_serialize_for_frequency_domain() bind(C) result(err)
            type(Serialized_t) :: serialize
            integer(kind=4) :: numberOfSerialized = 4
            integer :: test_err = 0

            call serialize%allocate_for_frequency_domain(numberOfSerialized)

            call check_shape_real(serialize%Valor, numberOfSerialized, test_err, "valor")
            call check_shape_real(serialize%Valor_x, numberOfSerialized, test_err, "Valor_x")
            call check_shape_real(serialize%Valor_y, numberOfSerialized, test_err, "Valor_y")
            call check_shape_real(serialize%Valor_z, numberOfSerialized, test_err, "Valor_z")
            call check_shape_real(serialize%ValorE, numberOfSerialized, test_err, "ValorE")
            call check_shape_real(serialize%Valor_Ex, numberOfSerialized, test_err, "Valor_Ex")
            call check_shape_real(serialize%Valor_Ey, numberOfSerialized, test_err, "Valor_Ey")
            call check_shape_real(serialize%Valor_Ez, numberOfSerialized, test_err, "Valor_Ez")
            call check_shape_real(serialize%ValorH, numberOfSerialized, test_err, "ValorH")
            call check_shape_real(serialize%Valor_Hx, numberOfSerialized, test_err, "Valor_Hx")
            call check_shape_real(serialize%Valor_Hy, numberOfSerialized, test_err, "Valor_Hy")
            call check_shape_real(serialize%Valor_Hz, numberOfSerialized, test_err, "Valor_Hz")

            call check_shape_complex(serialize%ValorComplex_x, numberOfSerialized, test_err, "ValorComplex_x")
            call check_shape_complex(serialize%ValorComplex_y, numberOfSerialized, test_err, "ValorComplex_y")
            call check_shape_complex(serialize%ValorComplex_z, numberOfSerialized, test_err, "ValorComplex_z")
            call check_shape_complex(serialize%ValorComplex_Ex, numberOfSerialized, test_err, "ValorComplex_Ex")
            call check_shape_complex(serialize%ValorComplex_Ey, numberOfSerialized, test_err, "ValorComplex_Ey")
            call check_shape_complex(serialize%ValorComplex_Ez, numberOfSerialized, test_err, "ValorComplex_Ez")
            call check_shape_complex(serialize%ValorComplex_Hx, numberOfSerialized, test_err, "ValorComplex_Hx")
            call check_shape_complex(serialize%ValorComplex_Hy, numberOfSerialized, test_err, "ValorComplex_Hy")
            call check_shape_complex(serialize%ValorComplex_Hz, numberOfSerialized, test_err, "ValorComplex_Hz")

            call serialize%deallocate_for_frequency_domain()

            err = test_err

        end function test_allocate_serialize_for_frequency_domain

        integer function test_allocate_current() bind(C) result(err)
            type(Serialized_t) :: serialize
            integer(kind=4) :: numberOfSerialized = 4
            integer :: test_err = 0

            call serialize%allocate_current_value(numberOfSerialized)

            call check_size(serialize%eI, numberOfSerialized, test_err, "eI")
            call check_size(serialize%eJ, numberOfSerialized, test_err, "eJ")
            call check_size(serialize%eK, numberOfSerialized, test_err, "eK")
            call check_size(serialize%currentType, numberOfSerialized, test_err, "currentType")
            call check_size(serialize%sggMtag, numberOfSerialized, test_err, "sggMtag")

            call serialize%deallocate_current_value()

            err = test_err
        end function test_allocate_current 



end module test_observation_m