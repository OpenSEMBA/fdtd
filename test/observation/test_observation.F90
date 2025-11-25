integer function test_allocate_serialize_for_time_domain() bind(C) result(err)
  use observation_testingTools
  use OBSERVA
  use FDETYPES
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
  use FDETYPES
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
  use FDETYPES
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

   function test_create_vtk_file() bind(C) result(error)
      use vtk_fortran
      use FDETYPES
      
      type(vtk_file)     :: a_vtk_file                             ! A VTK file.
      integer, parameter :: nx1=0_4                              ! X lower bound extent.
      integer, parameter :: nx2=3_4                              ! X upper bound extent.
      integer, parameter :: ny1=0_4                              ! Y lower bound extent.
      integer, parameter :: ny2=2_4                              ! Y upper bound extent.
      integer, parameter :: nz1=0_4                              ! Z lower bound extent.
      integer, parameter :: nz2=1_4                              ! Z upper bound extent.
      integer, parameter :: nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1) ! Number of elements.
      real               :: x(nx1:nx2,ny1:ny2,nz1:nz2) = 0.1_RKIND             ! X coordinates.
      real               :: y(nx1:nx2,ny1:ny2,nz1:nz2) = 0.2_RKIND             ! Y coordinates.
      real               :: z(nx1:nx2,ny1:ny2,nz1:nz2) = 0.3_RKIND             ! Z coordinates.
      real               :: v(nx1:nx2,ny1:ny2,nz1:nz2) = 0.4_RKIND             ! Variable at coordinates.
      integer            :: error                                  ! Error status.

      ! initialize the data...

      error = a_vtk_file%initialize(format='ASCII', filename='XML_STRG-binary.vts', &
                                    mesh_topology='StructuredGrid',                  &
                                    nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)
      error = a_vtk_file%xml_writer%write_piece(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)
      error = a_vtk_file%xml_writer%write_geo(n=nn, x=x, y=y, z=z)
      error = a_vtk_file%xml_writer%write_dataarray(location='node', action='open')
      error = a_vtk_file%xml_writer%write_dataarray(data_name='float64_scalar', x=v, one_component=.true.)
      error = a_vtk_file%xml_writer%write_dataarray(location='node', action='close')
      error = a_vtk_file%xml_writer%write_piece()
      error = a_vtk_file%finalize()
   end function test_create_vtk_file
