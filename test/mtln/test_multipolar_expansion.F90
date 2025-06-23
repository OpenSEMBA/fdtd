integer function test_multipolar_expansion_of_dipole() bind(C) result(error_cnt)    
    use mtln_types_mod
    use multipolar_expansion_mod

    real, dimension(2) :: expansionCenter = [0.0, 0.0]
    real :: d=0.1, r=1.0
    real :: vComputed, vExpected
    real, dimension(2) :: pos
    type(multipolar_coefficient_t), dimension(2) :: ab
    
    ab(1) = [0.0, d]
    ab(2) = [0.0, 0.0]
    
    ! First test
    block
        pos = [r, 0.0]
        vComputed = multipolarExpansion(pos, ab, expansionCenter)
        vExpected = 1.0 / (2.0 * pi()) * log((r + d / 2.0) / (r - d / 2.0))
        if (.not. expect_near(vExpected, vComputed, 1e-4)) then
            error_cnt = error_cnt + 1
        end if
    end block

    ! Second test
    block
        pos = [0.0, r]
        vComputed = multipolarExpansion(pos, ab, expansionCenter)
        vExpected = 0.0
        if (.not. expect_near(vExpected, vComputed, 1e-4)) then
            error_cnt = error_cnt + 1
        end if
    end block

end function

integer function test_multipolar_expansion_for_single_wire_with_dielectric() bind(C) result(error_cnt)
! type :: external_dielectric_t
!     real :: radius = 0.0
!     real :: relative_permittivity = 1.0
!     real :: effective_relative_permittivity = 1.0
! end type

! type :: external_field_segment_t
!     integer, dimension(3) :: position
!     integer :: direction = 0
!     real :: radius = 0.0
!     logical :: has_dielectric = .false.
!     type(external_dielectric_t) :: dielectric
!     real (kind=rkind) , pointer  ::  field => null()
! contains
!     private
!     procedure :: external_field_segments_eq
!     generic, public :: operator(==) => external_field_segments_eq
! end type
    error_cnt = 1
end function