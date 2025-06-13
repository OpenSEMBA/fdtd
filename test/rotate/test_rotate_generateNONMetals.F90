integer function test_rotate_generate_non_metals() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0  
    
    
    mpidir = 2
    allocate(this%DielRegs)
    
    
    this%DielRegs%nvols = 1
    allocate(this%DielRegs%vols(1))
    this%DielRegs%vols(1)%n_c1P = 2
    this%DielRegs%vols(1)%n_c2P = 2
    allocate(this%DielRegs%vols(1)%C1P(2))
    allocate(this%DielRegs%vols(1)%C2P(2))
    
    
    this%DielRegs%vols(1)%C1P(1)%XI = 1
    this%DielRegs%vols(1)%C1P(1)%YI = 2
    this%DielRegs%vols(1)%C1P(1)%ZI = 3
    this%DielRegs%vols(1)%C1P(1)%OR = iEx
    this%DielRegs%vols(1)%C1P(2)%XI = 4
    this%DielRegs%vols(1)%C1P(2)%YI = 5
    this%DielRegs%vols(1)%C1P(2)%ZI = 6
    this%DielRegs%vols(1)%C1P(2)%OR = iEy
    
    
    this%DielRegs%vols(1)%C2P(1)%XI = 7
    this%DielRegs%vols(1)%C2P(1)%YI = 8
    this%DielRegs%vols(1)%C2P(1)%ZI = 9
    this%DielRegs%vols(1)%C2P(1)%OR = iEz
    this%DielRegs%vols(1)%C2P(2)%XI = 10
    this%DielRegs%vols(1)%C2P(2)%YI = 11
    this%DielRegs%vols(1)%C2P(2)%ZI = 12
    this%DielRegs%vols(1)%C2P(2)%OR = iEx
    
    
    this%DielRegs%nsurfs = 1
    allocate(this%DielRegs%surfs(1))
    this%DielRegs%surfs(1)%n_c1P = 2
    this%DielRegs%surfs(1)%n_c2P = 2
    allocate(this%DielRegs%surfs(1)%C1P(2))
    allocate(this%DielRegs%surfs(1)%C2P(2))
    
    
    this%DielRegs%surfs(1)%C1P(1)%XI = 13
    this%DielRegs%surfs(1)%C1P(1)%YI = 14
    this%DielRegs%surfs(1)%C1P(1)%ZI = 15
    this%DielRegs%surfs(1)%C1P(1)%OR = iEy
    this%DielRegs%surfs(1)%C1P(2)%XI = 16
    this%DielRegs%surfs(1)%C1P(2)%YI = 17
    this%DielRegs%surfs(1)%C1P(2)%ZI = 18
    this%DielRegs%surfs(1)%C1P(2)%OR = iEz
    
    
    this%DielRegs%surfs(1)%C2P(1)%XI = 19
    this%DielRegs%surfs(1)%C2P(1)%YI = 20
    this%DielRegs%surfs(1)%C2P(1)%ZI = 21
    this%DielRegs%surfs(1)%C2P(1)%OR = iEx
    this%DielRegs%surfs(1)%C2P(2)%XI = 22
    this%DielRegs%surfs(1)%C2P(2)%YI = 23
    this%DielRegs%surfs(1)%C2P(2)%ZI = 24
    this%DielRegs%surfs(1)%C2P(2)%OR = iEy
    
    
    this%DielRegs%nlins = 1
    allocate(this%DielRegs%lins(1))
    this%DielRegs%lins(1)%n_c1P = 2
    this%DielRegs%lins(1)%n_c2P = 2
    allocate(this%DielRegs%lins(1)%C1P(2))
    allocate(this%DielRegs%lins(1)%C2P(2))
    
    
    this%DielRegs%lins(1)%C1P(1)%XI = 25
    this%DielRegs%lins(1)%C1P(1)%YI = 26
    this%DielRegs%lins(1)%C1P(1)%ZI = 27
    this%DielRegs%lins(1)%C1P(1)%OR = iEz
    this%DielRegs%lins(1)%C1P(2)%XI = 28
    this%DielRegs%lins(1)%C1P(2)%YI = 29
    this%DielRegs%lins(1)%C1P(2)%ZI = 30
    this%DielRegs%lins(1)%C1P(2)%OR = iEx
    
    
    this%DielRegs%lins(1)%C2P(1)%XI = 31
    this%DielRegs%lins(1)%C2P(1)%YI = 32
    this%DielRegs%lins(1)%C2P(1)%ZI = 33
    this%DielRegs%lins(1)%C2P(1)%OR = iEy
    this%DielRegs%lins(1)%C2P(2)%XI = 34
    this%DielRegs%lins(1)%C2P(2)%YI = 35
    this%DielRegs%lins(1)%C2P(2)%ZI = 36
    this%DielRegs%lins(1)%C2P(2)%OR = iEz
    
    
    call rotate_generateNONMetals(this, mpidir)
    
    
    call expect_eq_int(test_err, 3, this%DielRegs%vols(1)%C1P(1)%XI, "rotate_generateNONMetals: Vols C1P(1) XI should be 3")
    call expect_eq_int(test_err, 1, this%DielRegs%vols(1)%C1P(1)%YI, "rotate_generateNONMetals: Vols C1P(1) YI should be 1")
    call expect_eq_int(test_err, 2, this%DielRegs%vols(1)%C1P(1)%ZI, "rotate_generateNONMetals: Vols C1P(1) ZI should be 2")
    call expect_eq_int(test_err, iEy, this%DielRegs%vols(1)%C1P(1)%OR, "rotate_generateNONMetals: Vols C1P(1) OR should be iEy")
    
    call expect_eq_int(test_err, 6, this%DielRegs%vols(1)%C1P(2)%XI, "rotate_generateNONMetals: Vols C1P(2) XI should be 6")
    call expect_eq_int(test_err, 4, this%DielRegs%vols(1)%C1P(2)%YI, "rotate_generateNONMetals: Vols C1P(2) YI should be 4")
    call expect_eq_int(test_err, 5, this%DielRegs%vols(1)%C1P(2)%ZI, "rotate_generateNONMetals: Vols C1P(2) ZI should be 5")
    call expect_eq_int(test_err, iEz, this%DielRegs%vols(1)%C1P(2)%OR, "rotate_generateNONMetals: Vols C1P(2) OR should be iEz")
    
    
    call expect_eq_int(test_err, 9, this%DielRegs%vols(1)%C2P(1)%XI, "rotate_generateNONMetals: Vols C2P(1) XI should be 9")
    call expect_eq_int(test_err, 7, this%DielRegs%vols(1)%C2P(1)%YI, "rotate_generateNONMetals: Vols C2P(1) YI should be 7")
    call expect_eq_int(test_err, 8, this%DielRegs%vols(1)%C2P(1)%ZI, "rotate_generateNONMetals: Vols C2P(1) ZI should be 8")
    call expect_eq_int(test_err, iEx, this%DielRegs%vols(1)%C2P(1)%OR, "rotate_generateNONMetals: Vols C2P(1) OR should be iEx")
    
    call expect_eq_int(test_err, 12, this%DielRegs%vols(1)%C2P(2)%XI, "rotate_generateNONMetals: Vols C2P(2) XI should be 12")
    call expect_eq_int(test_err, 10, this%DielRegs%vols(1)%C2P(2)%YI, "rotate_generateNONMetals: Vols C2P(2) YI should be 10")
    call expect_eq_int(test_err, 11, this%DielRegs%vols(1)%C2P(2)%ZI, "rotate_generateNONMetals: Vols C2P(2) ZI should be 11")
    call expect_eq_int(test_err, iEy, this%DielRegs%vols(1)%C2P(2)%OR, "rotate_generateNONMetals: Vols C2P(2) OR should be iEy")
    
    
    call expect_eq_int(test_err, 15, this%DielRegs%surfs(1)%C1P(1)%XI, "rotate_generateNONMetals: Surfs C1P(1) XI should be 15")
    call expect_eq_int(test_err, 13, this%DielRegs%surfs(1)%C1P(1)%YI, "rotate_generateNONMetals: Surfs C1P(1) YI should be 13")
    call expect_eq_int(test_err, 14, this%DielRegs%surfs(1)%C1P(1)%ZI, "rotate_generateNONMetals: Surfs C1P(1) ZI should be 14")
    call expect_eq_int(test_err, iEz, this%DielRegs%surfs(1)%C1P(1)%OR, "rotate_generateNONMetals: Surfs C1P(1) OR should be iEz")
    
    call expect_eq_int(test_err, 18, this%DielRegs%surfs(1)%C1P(2)%XI, "rotate_generateNONMetals: Surfs C1P(2) XI should be 18")
    call expect_eq_int(test_err, 16, this%DielRegs%surfs(1)%C1P(2)%YI, "rotate_generateNONMetals: Surfs C1P(2) YI should be 16")
    call expect_eq_int(test_err, 17, this%DielRegs%surfs(1)%C1P(2)%ZI, "rotate_generateNONMetals: Surfs C1P(2) ZI should be 17")
    call expect_eq_int(test_err, iEx, this%DielRegs%surfs(1)%C1P(2)%OR, "rotate_generateNONMetals: Surfs C1P(2) OR should be iEx")
    
    
    call expect_eq_int(test_err, 21, this%DielRegs%surfs(1)%C2P(1)%XI, "rotate_generateNONMetals: Surfs C2P(1) XI should be 21")
    call expect_eq_int(test_err, 19, this%DielRegs%surfs(1)%C2P(1)%YI, "rotate_generateNONMetals: Surfs C2P(1) YI should be 19")
    call expect_eq_int(test_err, 20, this%DielRegs%surfs(1)%C2P(1)%ZI, "rotate_generateNONMetals: Surfs C2P(1) ZI should be 20")
    call expect_eq_int(test_err, iEy, this%DielRegs%surfs(1)%C2P(1)%OR, "rotate_generateNONMetals: Surfs C2P(1) OR should be iEy")
    
    call expect_eq_int(test_err, 24, this%DielRegs%surfs(1)%C2P(2)%XI, "rotate_generateNONMetals: Surfs C2P(2) XI should be 24")
    call expect_eq_int(test_err, 22, this%DielRegs%surfs(1)%C2P(2)%YI, "rotate_generateNONMetals: Surfs C2P(2) YI should be 22")
    call expect_eq_int(test_err, 23, this%DielRegs%surfs(1)%C2P(2)%ZI, "rotate_generateNONMetals: Surfs C2P(2) ZI should be 23")
    call expect_eq_int(test_err, iEz, this%DielRegs%surfs(1)%C2P(2)%OR, "rotate_generateNONMetals: Surfs C2P(2) OR should be iEz")
    
    
    call expect_eq_int(test_err, 27, this%DielRegs%lins(1)%C1P(1)%XI, "rotate_generateNONMetals: Lins C1P(1) XI should be 27")
    call expect_eq_int(test_err, 25, this%DielRegs%lins(1)%C1P(1)%YI, "rotate_generateNONMetals: Lins C1P(1) YI should be 25")
    call expect_eq_int(test_err, 26, this%DielRegs%lins(1)%C1P(1)%ZI, "rotate_generateNONMetals: Lins C1P(1) ZI should be 26")
    call expect_eq_int(test_err, iEx, this%DielRegs%lins(1)%C1P(1)%OR, "rotate_generateNONMetals: Lins C1P(1) OR should be iEx")
    
    call expect_eq_int(test_err, 30, this%DielRegs%lins(1)%C1P(2)%XI, "rotate_generateNONMetals: Lins C1P(2) XI should be 30")
    call expect_eq_int(test_err, 28, this%DielRegs%lins(1)%C1P(2)%YI, "rotate_generateNONMetals: Lins C1P(2) YI should be 28")
    call expect_eq_int(test_err, 29, this%DielRegs%lins(1)%C1P(2)%ZI, "rotate_generateNONMetals: Lins C1P(2) ZI should be 29")
    call expect_eq_int(test_err, iEy, this%DielRegs%lins(1)%C1P(2)%OR, "rotate_generateNONMetals: Lins C1P(2) OR should be iEy")
    
    
    call expect_eq_int(test_err, 33, this%DielRegs%lins(1)%C2P(1)%XI, "rotate_generateNONMetals: Lins C2P(1) XI should be 33")
    call expect_eq_int(test_err, 31, this%DielRegs%lins(1)%C2P(1)%YI, "rotate_generateNONMetals: Lins C2P(1) YI should be 31")
    call expect_eq_int(test_err, 32, this%DielRegs%lins(1)%C2P(1)%ZI, "rotate_generateNONMetals: Lins C2P(1) ZI should be 32")
    call expect_eq_int(test_err, iEz, this%DielRegs%lins(1)%C2P(1)%OR, "rotate_generateNONMetals: Lins C2P(1) OR should be iEz")
    
    call expect_eq_int(test_err, 36, this%DielRegs%lins(1)%C2P(2)%XI, "rotate_generateNONMetals: Lins C2P(2) XI should be 36")
    call expect_eq_int(test_err, 34, this%DielRegs%lins(1)%C2P(2)%YI, "rotate_generateNONMetals: Lins C2P(2) YI should be 34")
    call expect_eq_int(test_err, 35, this%DielRegs%lins(1)%C2P(2)%ZI, "rotate_generateNONMetals: Lins C2P(2) ZI should be 35")
    call expect_eq_int(test_err, iEx, this%DielRegs%lins(1)%C2P(2)%OR, "rotate_generateNONMetals: Lins C2P(2) OR should be iEx")
    
    
    deallocate(this%DielRegs%vols(1)%C1P)
    deallocate(this%DielRegs%vols(1)%C2P)
    deallocate(this%DielRegs%vols)
    deallocate(this%DielRegs%surfs(1)%C1P)
    deallocate(this%DielRegs%surfs(1)%C2P)
    deallocate(this%DielRegs%surfs)
    deallocate(this%DielRegs%lins(1)%C1P)
    deallocate(this%DielRegs%lins(1)%C2P)
    deallocate(this%DielRegs%lins)
    deallocate(this%DielRegs)
    
    
    mpidir = 1
    allocate(this%DielRegs)
    
    
    this%DielRegs%nvols = 1
    allocate(this%DielRegs%vols(1))
    this%DielRegs%vols(1)%n_c1P = 2
    this%DielRegs%vols(1)%n_c2P = 2
    allocate(this%DielRegs%vols(1)%C1P(2))
    allocate(this%DielRegs%vols(1)%C2P(2))
    
    
    this%DielRegs%vols(1)%C1P(1)%XI = 1
    this%DielRegs%vols(1)%C1P(1)%YI = 2
    this%DielRegs%vols(1)%C1P(1)%ZI = 3
    this%DielRegs%vols(1)%C1P(1)%OR = iEx
    this%DielRegs%vols(1)%C1P(2)%XI = 4
    this%DielRegs%vols(1)%C1P(2)%YI = 5
    this%DielRegs%vols(1)%C1P(2)%ZI = 6
    this%DielRegs%vols(1)%C1P(2)%OR = iEy
    
    
    this%DielRegs%vols(1)%C2P(1)%XI = 7
    this%DielRegs%vols(1)%C2P(1)%YI = 8
    this%DielRegs%vols(1)%C2P(1)%ZI = 9
    this%DielRegs%vols(1)%C2P(1)%OR = iEz
    this%DielRegs%vols(1)%C2P(2)%XI = 10
    this%DielRegs%vols(1)%C2P(2)%YI = 11
    this%DielRegs%vols(1)%C2P(2)%ZI = 12
    this%DielRegs%vols(1)%C2P(2)%OR = iEx
    
    
    this%DielRegs%nsurfs = 1
    allocate(this%DielRegs%surfs(1))
    this%DielRegs%surfs(1)%n_c1P = 2
    this%DielRegs%surfs(1)%n_c2P = 2
    allocate(this%DielRegs%surfs(1)%C1P(2))
    allocate(this%DielRegs%surfs(1)%C2P(2))
    
    
    this%DielRegs%surfs(1)%C1P(1)%XI = 13
    this%DielRegs%surfs(1)%C1P(1)%YI = 14
    this%DielRegs%surfs(1)%C1P(1)%ZI = 15
    this%DielRegs%surfs(1)%C1P(1)%OR = iEy
    this%DielRegs%surfs(1)%C1P(2)%XI = 16
    this%DielRegs%surfs(1)%C1P(2)%YI = 17
    this%DielRegs%surfs(1)%C1P(2)%ZI = 18
    this%DielRegs%surfs(1)%C1P(2)%OR = iEz
    
    
    this%DielRegs%surfs(1)%C2P(1)%XI = 19
    this%DielRegs%surfs(1)%C2P(1)%YI = 20
    this%DielRegs%surfs(1)%C2P(1)%ZI = 21
    this%DielRegs%surfs(1)%C2P(1)%OR = iEx
    this%DielRegs%surfs(1)%C2P(2)%XI = 22
    this%DielRegs%surfs(1)%C2P(2)%YI = 23
    this%DielRegs%surfs(1)%C2P(2)%ZI = 24
    this%DielRegs%surfs(1)%C2P(2)%OR = iEy
    
    
    this%DielRegs%nlins = 1
    allocate(this%DielRegs%lins(1))
    this%DielRegs%lins(1)%n_c1P = 2
    this%DielRegs%lins(1)%n_c2P = 2
    allocate(this%DielRegs%lins(1)%C1P(2))
    allocate(this%DielRegs%lins(1)%C2P(2))
    
    
    this%DielRegs%lins(1)%C1P(1)%XI = 25
    this%DielRegs%lins(1)%C1P(1)%YI = 26
    this%DielRegs%lins(1)%C1P(1)%ZI = 27
    this%DielRegs%lins(1)%C1P(1)%OR = iEz
    this%DielRegs%lins(1)%C1P(2)%XI = 28
    this%DielRegs%lins(1)%C1P(2)%YI = 29
    this%DielRegs%lins(1)%C1P(2)%ZI = 30
    this%DielRegs%lins(1)%C1P(2)%OR = iEx
    
    
    this%DielRegs%lins(1)%C2P(1)%XI = 31
    this%DielRegs%lins(1)%C2P(1)%YI = 32
    this%DielRegs%lins(1)%C2P(1)%ZI = 33
    this%DielRegs%lins(1)%C2P(1)%OR = iEy
    this%DielRegs%lins(1)%C2P(2)%XI = 34
    this%DielRegs%lins(1)%C2P(2)%YI = 35
    this%DielRegs%lins(1)%C2P(2)%ZI = 36
    this%DielRegs%lins(1)%C2P(2)%OR = iEz
    
    
    call rotate_generateNONMetals(this, mpidir)
    
    
    call expect_eq_int(test_err, 2, this%DielRegs%vols(1)%C1P(1)%XI, "rotate_generateNONMetals: Vols C1P(1) XI should be 2")
    call expect_eq_int(test_err, 3, this%DielRegs%vols(1)%C1P(1)%YI, "rotate_generateNONMetals: Vols C1P(1) YI should be 3")
    call expect_eq_int(test_err, 1, this%DielRegs%vols(1)%C1P(1)%ZI, "rotate_generateNONMetals: Vols C1P(1) ZI should be 1")
    call expect_eq_int(test_err, iEz, this%DielRegs%vols(1)%C1P(1)%OR, "rotate_generateNONMetals: Vols C1P(1) OR should be iEz")
    
    call expect_eq_int(test_err, 5, this%DielRegs%vols(1)%C1P(2)%XI, "rotate_generateNONMetals: Vols C1P(2) XI should be 5")
    call expect_eq_int(test_err, 6, this%DielRegs%vols(1)%C1P(2)%YI, "rotate_generateNONMetals: Vols C1P(2) YI should be 6")
    call expect_eq_int(test_err, 4, this%DielRegs%vols(1)%C1P(2)%ZI, "rotate_generateNONMetals: Vols C1P(2) ZI should be 4")
    call expect_eq_int(test_err, iEx, this%DielRegs%vols(1)%C1P(2)%OR, "rotate_generateNONMetals: Vols C1P(2) OR should be iEx")
    
    
    call expect_eq_int(test_err, 8, this%DielRegs%vols(1)%C2P(1)%XI, "rotate_generateNONMetals: Vols C2P(1) XI should be 8")
    call expect_eq_int(test_err, 9, this%DielRegs%vols(1)%C2P(1)%YI, "rotate_generateNONMetals: Vols C2P(1) YI should be 9")
    call expect_eq_int(test_err, 7, this%DielRegs%vols(1)%C2P(1)%ZI, "rotate_generateNONMetals: Vols C2P(1) ZI should be 7")
    call expect_eq_int(test_err, iEy, this%DielRegs%vols(1)%C2P(1)%OR, "rotate_generateNONMetals: Vols C2P(1) OR should be iEy")
    
    call expect_eq_int(test_err, 11, this%DielRegs%vols(1)%C2P(2)%XI, "rotate_generateNONMetals: Vols C2P(2) XI should be 11")
    call expect_eq_int(test_err, 12, this%DielRegs%vols(1)%C2P(2)%YI, "rotate_generateNONMetals: Vols C2P(2) YI should be 12")
    call expect_eq_int(test_err, 10, this%DielRegs%vols(1)%C2P(2)%ZI, "rotate_generateNONMetals: Vols C2P(2) ZI should be 10")
    call expect_eq_int(test_err, iEz, this%DielRegs%vols(1)%C2P(2)%OR, "rotate_generateNONMetals: Vols C2P(2) OR should be iEz")
    
    
    call expect_eq_int(test_err, 14, this%DielRegs%surfs(1)%C1P(1)%XI, "rotate_generateNONMetals: Surfs C1P(1) XI should be 14")
    call expect_eq_int(test_err, 15, this%DielRegs%surfs(1)%C1P(1)%YI, "rotate_generateNONMetals: Surfs C1P(1) YI should be 15")
    call expect_eq_int(test_err, 13, this%DielRegs%surfs(1)%C1P(1)%ZI, "rotate_generateNONMetals: Surfs C1P(1) ZI should be 13")
    call expect_eq_int(test_err, iEx, this%DielRegs%surfs(1)%C1P(1)%OR, "rotate_generateNONMetals: Surfs C1P(1) OR should be iEx")
    
    call expect_eq_int(test_err, 17, this%DielRegs%surfs(1)%C1P(2)%XI, "rotate_generateNONMetals: Surfs C1P(2) XI should be 17")
    call expect_eq_int(test_err, 18, this%DielRegs%surfs(1)%C1P(2)%YI, "rotate_generateNONMetals: Surfs C1P(2) YI should be 18")
    call expect_eq_int(test_err, 16, this%DielRegs%surfs(1)%C1P(2)%ZI, "rotate_generateNONMetals: Surfs C1P(2) ZI should be 16")
    call expect_eq_int(test_err, iEy, this%DielRegs%surfs(1)%C1P(2)%OR, "rotate_generateNONMetals: Surfs C1P(2) OR should be iEy")
    
    
    call expect_eq_int(test_err, 20, this%DielRegs%surfs(1)%C2P(1)%XI, "rotate_generateNONMetals: Surfs C2P(1) XI should be 20")
    call expect_eq_int(test_err, 21, this%DielRegs%surfs(1)%C2P(1)%YI, "rotate_generateNONMetals: Surfs C2P(1) YI should be 21")
    call expect_eq_int(test_err, 19, this%DielRegs%surfs(1)%C2P(1)%ZI, "rotate_generateNONMetals: Surfs C2P(1) ZI should be 19")
    call expect_eq_int(test_err, iEz, this%DielRegs%surfs(1)%C2P(1)%OR, "rotate_generateNONMetals: Surfs C2P(1) OR should be iEz")
    
    call expect_eq_int(test_err, 23, this%DielRegs%surfs(1)%C2P(2)%XI, "rotate_generateNONMetals: Surfs C2P(2) XI should be 23")
    call expect_eq_int(test_err, 24, this%DielRegs%surfs(1)%C2P(2)%YI, "rotate_generateNONMetals: Surfs C2P(2) YI should be 24")
    call expect_eq_int(test_err, 22, this%DielRegs%surfs(1)%C2P(2)%ZI, "rotate_generateNONMetals: Surfs C2P(2) ZI should be 22")
    call expect_eq_int(test_err, iEx, this%DielRegs%surfs(1)%C2P(2)%OR, "rotate_generateNONMetals: Surfs C2P(2) OR should be iEx")
    
    
    call expect_eq_int(test_err, 26, this%DielRegs%lins(1)%C1P(1)%XI, "rotate_generateNONMetals: Lins C1P(1) XI should be 26")
    call expect_eq_int(test_err, 27, this%DielRegs%lins(1)%C1P(1)%YI, "rotate_generateNONMetals: Lins C1P(1) YI should be 27")
    call expect_eq_int(test_err, 25, this%DielRegs%lins(1)%C1P(1)%ZI, "rotate_generateNONMetals: Lins C1P(1) ZI should be 25")
    call expect_eq_int(test_err, iEy, this%DielRegs%lins(1)%C1P(1)%OR, "rotate_generateNONMetals: Lins C1P(1) OR should be iEy")
    
    call expect_eq_int(test_err, 29, this%DielRegs%lins(1)%C1P(2)%XI, "rotate_generateNONMetals: Lins C1P(2) XI should be 29")
    call expect_eq_int(test_err, 30, this%DielRegs%lins(1)%C1P(2)%YI, "rotate_generateNONMetals: Lins C1P(2) YI should be 30")
    call expect_eq_int(test_err, 28, this%DielRegs%lins(1)%C1P(2)%ZI, "rotate_generateNONMetals: Lins C1P(2) ZI should be 28")
    call expect_eq_int(test_err, iEz, this%DielRegs%lins(1)%C1P(2)%OR, "rotate_generateNONMetals: Lins C1P(2) OR should be iEz")
    
    
    call expect_eq_int(test_err, 32, this%DielRegs%lins(1)%C2P(1)%XI, "rotate_generateNONMetals: Lins C2P(1) XI should be 32")
    call expect_eq_int(test_err, 33, this%DielRegs%lins(1)%C2P(1)%YI, "rotate_generateNONMetals: Lins C2P(1) YI should be 33")
    call expect_eq_int(test_err, 31, this%DielRegs%lins(1)%C2P(1)%ZI, "rotate_generateNONMetals: Lins C2P(1) ZI should be 31")
    call expect_eq_int(test_err, iEx, this%DielRegs%lins(1)%C2P(1)%OR, "rotate_generateNONMetals: Lins C2P(1) OR should be iEx")
    
    call expect_eq_int(test_err, 35, this%DielRegs%lins(1)%C2P(2)%XI, "rotate_generateNONMetals: Lins C2P(2) XI should be 35")
    call expect_eq_int(test_err, 36, this%DielRegs%lins(1)%C2P(2)%YI, "rotate_generateNONMetals: Lins C2P(2) YI should be 36")
    call expect_eq_int(test_err, 34, this%DielRegs%lins(1)%C2P(2)%ZI, "rotate_generateNONMetals: Lins C2P(2) ZI should be 34")
    call expect_eq_int(test_err, iEy, this%DielRegs%lins(1)%C2P(2)%OR, "rotate_generateNONMetals: Lins C2P(2) OR should be iEy")
    
    
    deallocate(this%DielRegs%vols(1)%C1P)
    deallocate(this%DielRegs%vols(1)%C2P)
    deallocate(this%DielRegs%vols)
    deallocate(this%DielRegs%surfs(1)%C1P)
    deallocate(this%DielRegs%surfs(1)%C2P)
    deallocate(this%DielRegs%surfs)
    deallocate(this%DielRegs%lins(1)%C1P)
    deallocate(this%DielRegs%lins(1)%C2P)
    deallocate(this%DielRegs%lins)
    deallocate(this%DielRegs)
    
    err = test_err  
end function test_rotate_generate_non_metals