integer function test_multipolar_expansion_of_dipole() bind(C) result(error_cnt)    
    use mtln_types_mod
    use multipolar_expansion_mod
    use mtln_testingTools_mod

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
        if (.not. checkNear(vExpected, vComputed, 1e-4)) then
            error_cnt = error_cnt + 1
        end if
    end block

    ! Second test
    block
        pos = [0.0, r]
        vComputed = multipolarExpansion(pos, ab, expansionCenter)
        vExpected = 0.0
        if (.not. checkNear(vExpected, vComputed, 1e-4)) then
            error_cnt = error_cnt + 1
        end if
    end block

end function


integer function test_multipolar_expansion_for_lansink2024_wire_with_dielectric() bind(C) result(error_cnt)
    use mtln_types_mod
    use multipolar_expansion_mod

    type(multipolar_expansion_t) :: multipolarExpansion


end function

integer function test_multipolar_expansion_for_lansink2024_two_wires() bind(C) result(error_cnt)
    use mtln_types_mod
    use multipolar_expansion_mod

    type(multipolar_expansion_t) :: mE
    type(box_2d_t) :: fdtdCell

    real, dimension(:,:), allocatable :: computedL, computedC

    mE%inner_region%min = [-0.004, -0.004]
    mE%inner_region%max = [ 0.004,  0.004]

    allocate(mE%electric(1))
    mE%electric(1)%inner_region_average_potential = 0.90407844239490087
    mE%electric(1)%expansion_center = [0.0, 0.0]
    mE%electric(1)%conductor_potentials = [1.0]
    allocate(mE%electric(1)%ab(1))
    mE%electric(1)%ab(1)%a = 0.97667340898489752
    mE%electric(1)%ab(1)%b = 0.0

    allocate(mE%magnetic(1))
    mE%magnetic(1)%inner_region_average_potential = 0.84903792056711014
    mE%magnetic(1)%expansion_center = [0.0, 0.0]
    mE%magnetic(1)%conductor_potentials = [1.0]
    allocate(mE%magnetic(1)%ab(1))
    mE%magnetic(1)%ab(1)%a = 0.90929569352397666
    mE%magnetic(1)%ab(1)%b = 0.0


    fdtdCell%min = [-0.0075, -0.0075]
    fdtdCell%max = [ 0.0075,  0.0075]

    computedC = getCellCapacitanceOnBox(mE, fdtdCell)
    if (.not. checkNear(49e-12, computedC, 1e-4)) then
        error_cnt = error_cnt + 1
    end if

    computedL = getCellInductanceOnBox(mE, fdtdCell)
    if (.not. checkNear(320e-9, computedL, 1e-4)) then
        error_cnt = error_cnt + 1
    end if

end function