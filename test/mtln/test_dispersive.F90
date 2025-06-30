integer function test_dispersive_init_1_pole() bind(C) result(error_cnt)    
    use dispersive_mod
    use mtln_testingTools_mod
    use mtl_bundle_mod

    implicit none
    type(transfer_impedance_per_meter_t) :: zt
    complex, dimension(1) :: residues, poles
    type(pol_res_t):: connector
    type(mtl_t) :: line_out, line_in

    type(transmission_line_level_t), dimension(2) :: levels
    type(transmission_line_level_t) :: level1
    type(transmission_line_level_t) :: level2
    type(mtl_bundle_t) :: bundle
    
    integer :: ndiv, nc, npoles
    error_cnt = 0

    residues(1)%re = 1e4
    residues(1)%im = 1e3
    poles(1)%re = -1e6
    poles(1)%im = 1e1
    !complex conjugate
    ! residues(2)%re = 1e4
    ! residues(2)%im = -1e3
    ! poles(2)%re = -1e6
    ! poles(2)%im = -1e1
    
    zt%direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH
    zt%resistive_term = 1e-2
    zt%inductive_term = 1e-6
    zt%poles = poles
    zt%residues = residues
    
    
    line_out = buildLineWithNConductors(1, 'line_out')
    line_in = buildLineWithNConductors(1, 'line_in', parent_name = "line_out", conductor_in_parent = 1)
    level1%lines = [line_out]
    level2%lines = [line_in]
    levels = [level1, level2]
    
    bundle = mtl_bundle_t(levels, name='bundle')
    call bundle%addTransferImpedance(1, [2], zt)

    ndiv = 5
    nc = 2
    npoles = 1

    if (size(bundle%transfer_impedance%q1) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%phi) /= ndiv*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%d) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%e) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q1_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3_phi) /= ndiv*nc) then
        error_cnt = error_cnt + 1
    end if
    
    if(.not.all(bundle%transfer_impedance%q1(:,1,2,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,2,1,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,2,2,:) == 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,1,1,:) == 0)) then
        error_cnt = error_cnt + 1
    end if


end function

integer function test_dispersive_init_1_pole_lines_with_lumped() bind(C) result(error_cnt)    
    use dispersive_mod
    use mtln_testingTools_mod
    use mtl_bundle_mod

    implicit none

    type(transfer_impedance_per_meter_t) :: zt
    complex, dimension(1) :: residues, poles
    type(pol_res_t):: connector
    type(mtl_t) :: line_out, line_in

    type(transmission_line_level_t), dimension(2) :: levels
    type(transmission_line_level_t) :: level1
    type(transmission_line_level_t) :: level2
    type(mtl_bundle_t) :: bundle
    
    integer :: ndiv, nc, npoles
    
    error_cnt = 0

    residues(1)%re = 1e4
    residues(1)%im = 1e3
    poles(1)%re = -1e6
    poles(1)%im = 1e1
    !complex conjugate
    ! residues(2)%re = 1e4
    ! residues(2)%im = -1e3
    ! poles(2)%re = -1e6
    ! poles(2)%im = -1e1
    
    zt%direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH
    zt%resistive_term = 1e-2
    zt%inductive_term = 1e-6
    zt%poles = poles
    zt%residues = residues
    
    
    line_out = buildLineWithNConductors(1, 'line_out')
    line_in = buildLineWithNConductors(1, 'line_in', parent_name = "line_out", conductor_in_parent = 1)

    call line_out%lumped_elements%addDispersiveLumped(1, 1, zt)
    call line_in%lumped_elements%addDispersiveLumped(5, 1, zt)

    level1%lines = [line_out]
    level2%lines = [line_in]
    levels = [level1, level2]
    
    bundle = mtl_bundle_t(levels, name='bundle')
    call bundle%addTransferImpedance(1, [2], zt)

    ndiv = 5
    nc = 2
    npoles = 1

    if (size(bundle%transfer_impedance%q1) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%phi) /= ndiv*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%d) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%e) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q1_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3_phi) /= ndiv*nc) then
        error_cnt = error_cnt + 1
    end if
    
    if(.not.all(bundle%transfer_impedance%q1(:,1,2,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,2,1,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    block 
        integer :: i
        integer, dimension(2) :: pos = [1,5]
        do i=2,4
            if(.not.all(bundle%transfer_impedance%q1(i,2,2,:) == 0)) then
                error_cnt = error_cnt + 1
            end if
            if(.not.all(bundle%transfer_impedance%q1(i,1,1,:) == 0)) then
                error_cnt = error_cnt + 1
            end if
        end do
        if(.not.all(bundle%transfer_impedance%q1(5,2,2,:) /= 0)) then
            error_cnt = error_cnt + 1
        end if
        if(.not.all(bundle%transfer_impedance%q1(1,1,1,:) /= 0)) then
            error_cnt = error_cnt + 1
        end if
        if(.not.all(bundle%transfer_impedance%q1(5,1,1,:) == 0)) then
            error_cnt = error_cnt + 1
        end if
        if(.not.all(bundle%transfer_impedance%q1(1,2,2,:) == 0)) then
            error_cnt = error_cnt + 1
        end if
    end block

end function

integer function test_dispersive_init_2_poles() bind(C) result(error_cnt)    
    use dispersive_mod
    use mtln_testingTools_mod
    use mtl_bundle_mod

    implicit none
    type(transfer_impedance_per_meter_t) :: zt
    complex, dimension(2) :: residues, poles
    type(pol_res_t):: connector
    type(mtl_t) :: line_out, line_in

    type(transmission_line_level_t), dimension(2) :: levels
    type(transmission_line_level_t) :: level1
    type(transmission_line_level_t) :: level2
    type(mtl_bundle_t) :: bundle
    
    integer :: ndiv, nc, npoles
    error_cnt = 0

    residues(1)%re = 1e4
    residues(1)%im = 1e3
    poles(1)%re = -1e6
    poles(1)%im = 1e1
    !complex conjugate
    residues(2)%re = 1e4
    residues(2)%im = -1e3
    poles(2)%re = -1e6
    poles(2)%im = -1e1
    
    zt%direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH
    zt%resistive_term = 1e-2
    zt%inductive_term = 1e-6
    zt%poles = poles
    zt%residues = residues
    
    
    line_out = buildLineWithNConductors(1, 'line_out')
    line_in = buildLineWithNConductors(3, 'line_in', parent_name = "line_out", conductor_in_parent = 1)
    level1%lines = [line_out]
    level2%lines = [line_in]
    levels = [level1, level2]
    
    bundle = mtl_bundle_t(levels, name='bundle')
    call bundle%addTransferImpedance(1, [2,3,4], zt)

    ndiv = 5
    nc = 4
    npoles = 2

    if (size(bundle%transfer_impedance%q1) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%phi) /= ndiv*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%d) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%e) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q1_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3_phi) /= ndiv*nc) then
        error_cnt = error_cnt + 1
    end if
    
    if(.not.all(bundle%transfer_impedance%q1(:,1,2:4,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,2:4,1,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,2:4,2:4,:) == 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,1,1,:) == 0)) then
        error_cnt = error_cnt + 1
    end if


end function

integer function test_dispersive_init_1_pole_3_levels() bind(C) result(error_cnt)    
    use dispersive_mod
    use mtln_testingTools_mod
    use mtl_bundle_mod

    implicit none
    type(transfer_impedance_per_meter_t) :: zt
    complex, dimension(1) :: residues, poles
    type(pol_res_t):: connector
    type(mtl_t) :: line_out, line_in, line_in_2

    type(transmission_line_level_t), dimension(3) :: levels
    type(transmission_line_level_t) :: level1, level2, level3
    type(mtl_bundle_t) :: bundle
    
    integer :: ndiv, nc, npoles
    error_cnt = 0

    residues(1)%re = 1e4
    residues(1)%im = 1e3
    poles(1)%re = -1e6
    poles(1)%im = 1e1
    
    zt%direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
    zt%resistive_term = 1e-2
    zt%inductive_term = 1e-6
    zt%poles = poles
    zt%residues = residues
    
    
    line_out = buildLineWithNConductors(1, 'line_out')
    line_in = buildLineWithNConductors(2, 'line_in', parent_name = "line_out", conductor_in_parent = 1)
    line_in_2 = buildLineWithNConductors(2, 'line_in_2', parent_name = "line_in", conductor_in_parent = 1)
    level1%lines = [line_out]
    level2%lines = [line_in]
    level3%lines = [line_in_2]
    levels = [level1, level2, level3]
    
    bundle = mtl_bundle_t(levels, name='bundle')
    call bundle%addTransferImpedance(1, [2,3], zt)
    call bundle%addTransferImpedance(2, [4,5], zt)

    ndiv = 5
    nc = 5
    npoles = 1

    if (size(bundle%transfer_impedance%q1) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3) /= ndiv*nc*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%phi) /= ndiv*nc*npoles) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%d) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%e) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q1_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q2_sum) /= ndiv*nc*nc) then
        error_cnt = error_cnt + 1
    end if
    if (size(bundle%transfer_impedance%q3_phi) /= ndiv*nc) then
        error_cnt = error_cnt + 1
    end if
    
    if(.not.all(bundle%transfer_impedance%q1(:,2:3,1,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,4:5,2,:) /= 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,3:5,3:5,:) == 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,1,1:3,:) == 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,2:3,2,:) == 0)) then
        error_cnt = error_cnt + 1
    end if
    if(.not.all(bundle%transfer_impedance%q1(:,4:5,1,:) == 0)) then
        error_cnt = error_cnt + 1
    end if


end function