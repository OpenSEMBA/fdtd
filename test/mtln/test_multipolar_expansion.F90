integer function test_multipolar_expansion_for_dipole() bind(C) result(error_cnt)    
    use mtln_types_mod
    use multipolar_expansion_mod
    use mtln_testingTools_mod
    
    real, dimension(2) :: expansionCenter = [0.0, 0.0]
    real :: d=0.1, r=1.0 
    real :: vComputed, vExpected
    real, dimension(2) :: pos
    type(multipolar_coefficient_t), dimension(2) :: ab

    error_cnt = 0
    
    ab(1)%a = 0.0; ab(1)%b = 0.0
    ab(2)%a = d; ab(2)%b = 0.0
    
    ! First test
    block
        pos = [r, 0.0]
        vComputed = multipolarExpansion2D(pos, ab, expansionCenter)
        vExpected = 1.0 / (2.0 * pi) * log((r + d / 2.0) / (r - d / 2.0))
        if (abs(vExpected - vComputed) > 1e-4) then
            error_cnt = error_cnt + 1
        end if
    end block

    ! Second test
    block
        pos = [0.0, r]
        vComputed = multipolarExpansion2D(pos, ab, expansionCenter)
        vExpected = 0.0
        if (abs(vExpected - vComputed) > 1e-4) then
            error_cnt = error_cnt + 1
        end if
    end block

end function


integer function test_multipolar_expansion_for_lansink_two_wires() bind(C) result(error_cnt)
	! From:
	! Rotgerink, J.L. et al. (2024, September).
	! Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	! In 2024 International Symposium on Electromagnetic Compatibility
	! EMC Europe(pp. 334 - 339). IEEE.
    use mtln_types_mod
    use multipolar_expansion_mod
    use mtln_testingTools_mod

    type(multipolar_expansion_t) :: mE
    type(box_2d_t) :: fdtdCell

    real, dimension(:,:), allocatable :: computedL, computedC

    error_cnt = 0

    mE%inner_region%min = [-0.0265, -0.031]
    mE%inner_region%max = [ 0.0355,  0.031]

    allocate(mE%electric(2))
    mE%electric(1)%inner_region_average_potential = 0.56086362615993235
    mE%electric(1)%expansion_center = [-0.00497, 0.0]
    mE%electric(1)%conductor_potentials = [1.0, 0.59092722039872780]
    allocate(mE%electric(1)%ab(3))
    mE%electric(1)%ab(1)%a = 0.94888369862564237
    mE%electric(1)%ab(1)%b = 0.0
    mE%electric(1)%ab(2)%a = -4.5026111057062017e-07
    mE%electric(1)%ab(2)%b = 0.0
    mE%electric(1)%ab(3)%a = 7.1226466480672610e-05
    mE%electric(1)%ab(3)%b = 7.4826684830590268e-09

    mE%electric(2)%inner_region_average_potential = 0.80708482435726114
    mE%electric(2)%expansion_center = [0.0099205134400286565, 0.0]
    mE%electric(2)%conductor_potentials = [0.84971105674469871, 1.0]
    allocate(mE%electric(2)%ab(3))
    mE%electric(2)%ab(1)%a = 1.3644011168458479
    mE%electric(2)%ab(1)%b = 0.0
    mE%electric(2)%ab(2)%a = 1.7915912102364171e-06
    mE%electric(2)%ab(2)%b = 0.0
    mE%electric(2)%ab(3)%a = 1.4620553866347293e-06
    mE%electric(2)%ab(3)%b = -1.4363492844460606e-09

    mE%magnetic = mE%electric

    fdtdCell%min = [-0.100, -0.100]
    fdtdCell%max = [ 0.100,  0.100]

    computedC = getCellCapacitanceOnBox(mE, fdtdCell)
    if (.not. checkNear(14.08e-12, computedC(1,1), 6e-2)) then
        error_cnt = error_cnt + 1
    end if
    if (.not. checkNear(43.99e-12, computedC(1,2), 6e-2)) then
        error_cnt = error_cnt + 1
    end if
    if (.not. checkNear(44.31e-12, computedC(2,1), 6e-2)) then
        error_cnt = error_cnt + 1
    end if
    if (.not. checkNear(28.79e-12, computedC(2,2), 6e-2)) then
        error_cnt = error_cnt + 1
    end if

    computedL = getCellInductanceOnBox(mE, fdtdCell)
    if (.not. checkNear(791e-9, computedL(1,1), 6e-2)) then
        error_cnt = error_cnt + 1
    end if
    if (.not. checkNear(253e-9, computedL(1,2), 6e-2)) then
        error_cnt = error_cnt + 1
    end if
    if (.not. checkNear(251e-9, computedL(2,1), 6e-2)) then
        error_cnt = error_cnt + 1
    end if
    if (.not. checkNear(387e-9, computedL(2,2), 6e-2)) then
        error_cnt = error_cnt + 1
    end if


end function

integer function test_multipolar_expansion_for_lansink_wire_with_dielectric() bind(C) result(error_cnt)
	! From:
	! Rotgerink, J.L. et al. (2024, September).
	! Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	! In 2024 International Symposium on Electromagnetic Compatibility
	! EMC Europe(pp. 334 - 339). IEEE.
    use mtln_types_mod
    use multipolar_expansion_mod
    use mtln_testingTools_mod

    type(multipolar_expansion_t) :: mE
    type(box_2d_t) :: fdtdCell

    real, dimension(:,:), allocatable :: computedL, computedC

    error_cnt = 0

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
    if (.not. checkNear(49e-12, computedC(1,1), 1e-4)) then
        error_cnt = error_cnt + 1
    end if

    computedL = getCellInductanceOnBox(mE, fdtdCell)
    if (.not. checkNear(320e-9, computedL(1,1), 1e-4)) then
        error_cnt = error_cnt + 1
    end if

end function