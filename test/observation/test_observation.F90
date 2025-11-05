integer function test_allocate_serialize_for_time_domain() bind(C) result(err)
  use observation_testingTools
  use OBSERVA
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
  use observation_testingTools
  use OBSERVA
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
  use observation_testingTools
  use OBSERVA
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
