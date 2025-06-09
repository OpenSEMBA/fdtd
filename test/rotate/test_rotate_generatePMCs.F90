integer function test_rotate_generate_pmcs() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0  ! Error counter for expect_eq calls
    
    ! Test case 1: X->Y->Z->X rotation (mpidir=2)
    mpidir = 2
    allocate(this%pmcRegs)
    
    ! Initialize test data for volumes
    this%pmcRegs%nvols = 1
    allocate(this%pmcRegs%Vols(1))
    this%pmcRegs%Vols(1)%XI = 1
    this%pmcRegs%Vols(1)%YI = 2
    this%pmcRegs%Vols(1)%ZI = 3
    this%pmcRegs%Vols(1)%XE = 4
    this%pmcRegs%Vols(1)%YE = 5
    this%pmcRegs%Vols(1)%ZE = 6
    this%pmcRegs%Vols(1)%OR = iEx
    
    ! Initialize test data for surfaces
    this%pmcRegs%nsurfs = 1
    allocate(this%pmcRegs%Surfs(1))
    this%pmcRegs%Surfs(1)%XI = 7
    this%pmcRegs%Surfs(1)%YI = 8
    this%pmcRegs%Surfs(1)%ZI = 9
    this%pmcRegs%Surfs(1)%XE = 10
    this%pmcRegs%Surfs(1)%YE = 11
    this%pmcRegs%Surfs(1)%ZE = 12
    this%pmcRegs%Surfs(1)%OR = iEy
    
    ! Initialize test data for lines
    this%pmcRegs%nlins = 1
    allocate(this%pmcRegs%Lins(1))
    this%pmcRegs%Lins(1)%XI = 13
    this%pmcRegs%Lins(1)%YI = 14
    this%pmcRegs%Lins(1)%ZI = 15
    this%pmcRegs%Lins(1)%XE = 16
    this%pmcRegs%Lins(1)%YE = 17
    this%pmcRegs%Lins(1)%ZE = 18
    this%pmcRegs%Lins(1)%OR = iEz
    
    ! Call the routine
    call rotate_generatePMCs(this, mpidir)
    
    ! Verify results for volumes
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%XI, 3, "rotate_generatePMCs: Vols XI should be 3")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%YI, 1, "rotate_generatePMCs: Vols YI should be 1")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%ZI, 2, "rotate_generatePMCs: Vols ZI should be 2")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%XE, 6, "rotate_generatePMCs: Vols XE should be 6")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%YE, 4, "rotate_generatePMCs: Vols YE should be 4")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%ZE, 5, "rotate_generatePMCs: Vols ZE should be 5")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%OR, iEy, "rotate_generatePMCs: Vols OR should be iEy")
    
    ! Verify results for surfaces
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%XI, 9, "rotate_generatePMCs: Surfs XI should be 9")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%YI, 7, "rotate_generatePMCs: Surfs YI should be 7")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%ZI, 8, "rotate_generatePMCs: Surfs ZI should be 8")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%XE, 12, "rotate_generatePMCs: Surfs XE should be 12")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%YE, 10, "rotate_generatePMCs: Surfs YE should be 10")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%ZE, 11, "rotate_generatePMCs: Surfs ZE should be 11")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%OR, iEz, "rotate_generatePMCs: Surfs OR should be iEz")
    
    ! Verify results for lines
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%XI, 15, "rotate_generatePMCs: Lins XI should be 15")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%YI, 13, "rotate_generatePMCs: Lins YI should be 13")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%ZI, 14, "rotate_generatePMCs: Lins ZI should be 14")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%XE, 18, "rotate_generatePMCs: Lins XE should be 18")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%YE, 16, "rotate_generatePMCs: Lins YE should be 16")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%ZE, 17, "rotate_generatePMCs: Lins ZE should be 17")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%OR, iEx, "rotate_generatePMCs: Lins OR should be iEx")
    
    deallocate(this%pmcRegs%Vols)
    deallocate(this%pmcRegs%Surfs)
    deallocate(this%pmcRegs%Lins)
    deallocate(this%pmcRegs)
    
    ! Test case 2: X->Z->Y->X rotation (mpidir=1)
    mpidir = 1
    allocate(this%pmcRegs)
    
    ! Initialize test data (same as before)
    this%pmcRegs%nvols = 1
    allocate(this%pmcRegs%Vols(1))
    this%pmcRegs%Vols(1)%XI = 1
    this%pmcRegs%Vols(1)%YI = 2
    this%pmcRegs%Vols(1)%ZI = 3
    this%pmcRegs%Vols(1)%XE = 4
    this%pmcRegs%Vols(1)%YE = 5
    this%pmcRegs%Vols(1)%ZE = 6
    this%pmcRegs%Vols(1)%OR = iEx
    
    this%pmcRegs%nsurfs = 1
    allocate(this%pmcRegs%Surfs(1))
    this%pmcRegs%Surfs(1)%XI = 7
    this%pmcRegs%Surfs(1)%YI = 8
    this%pmcRegs%Surfs(1)%ZI = 9
    this%pmcRegs%Surfs(1)%XE = 10
    this%pmcRegs%Surfs(1)%YE = 11
    this%pmcRegs%Surfs(1)%ZE = 12
    this%pmcRegs%Surfs(1)%OR = iEy
    
    this%pmcRegs%nlins = 1
    allocate(this%pmcRegs%Lins(1))
    this%pmcRegs%Lins(1)%XI = 13
    this%pmcRegs%Lins(1)%YI = 14
    this%pmcRegs%Lins(1)%ZI = 15
    this%pmcRegs%Lins(1)%XE = 16
    this%pmcRegs%Lins(1)%YE = 17
    this%pmcRegs%Lins(1)%ZE = 18
    this%pmcRegs%Lins(1)%OR = iEz
    
    ! Call the routine
    call rotate_generatePMCs(this, mpidir)
    
    ! Verify results for volumes
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%XI, 2, "rotate_generatePMCs: Vols XI should be 2")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%YI, 3, "rotate_generatePMCs: Vols YI should be 3")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%ZI, 1, "rotate_generatePMCs: Vols ZI should be 1")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%XE, 5, "rotate_generatePMCs: Vols XE should be 5")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%YE, 6, "rotate_generatePMCs: Vols YE should be 6")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%ZE, 4, "rotate_generatePMCs: Vols ZE should be 4")
    call expect_eq_int(test_err, this%pmcRegs%Vols(1)%OR, iEz, "rotate_generatePMCs: Vols OR should be iEz")
    
    ! Verify results for surfaces
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%XI, 8, "rotate_generatePMCs: Surfs XI should be 8")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%YI, 9, "rotate_generatePMCs: Surfs YI should be 9")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%ZI, 7, "rotate_generatePMCs: Surfs ZI should be 7")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%XE, 11, "rotate_generatePMCs: Surfs XE should be 11")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%YE, 12, "rotate_generatePMCs: Surfs YE should be 12")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%ZE, 10, "rotate_generatePMCs: Surfs ZE should be 10")
    call expect_eq_int(test_err, this%pmcRegs%Surfs(1)%OR, iEx, "rotate_generatePMCs: Surfs OR should be iEx")
    
    ! Verify results for lines
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%XI, 14, "rotate_generatePMCs: Lins XI should be 14")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%YI, 15, "rotate_generatePMCs: Lins YI should be 15")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%ZI, 13, "rotate_generatePMCs: Lins ZI should be 13")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%XE, 17, "rotate_generatePMCs: Lins XE should be 17")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%YE, 18, "rotate_generatePMCs: Lins YE should be 18")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%ZE, 16, "rotate_generatePMCs: Lins ZE should be 16")
    call expect_eq_int(test_err, this%pmcRegs%Lins(1)%OR, iEy, "rotate_generatePMCs: Lins OR should be iEy")
    
    deallocate(this%pmcRegs%Vols)
    deallocate(this%pmcRegs%Surfs)
    deallocate(this%pmcRegs%Lins)
    deallocate(this%pmcRegs)
    
    err = test_err  ! Set the function result to the accumulated error count
end function