module test_rotate_mtln_m
    use mtln_types_m
    use smbjson_m
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_mtln() bind(C) result(err)
    type(Parseador_t) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    type(cable_t), target :: cable
    allocate(this%mtln)
    allocate(this%mtln%cables(1))
    allocate(cable%segments(1))

    cable%segments(1)%x = 1
    cable%segments(1)%y = 2 
    cable%segments(1)%z = 3
    cable%segments(1)%orientation = 1
    this%mtln%cables(1)%ptr => cable

    mpidir = 2
    call rotate_mtln(this, mpidir)
    call expect_eq_int(test_err, 3, this%mtln%cables(1)%ptr%segments(1)%x, "rotate_mtln: x should be rotated with mpidir 2")
    call expect_eq_int(test_err, 1, this%mtln%cables(1)%ptr%segments(1)%y, "rotate_mtln: y should be rotated with mpidir 2")
    call expect_eq_int(test_err, 2, this%mtln%cables(1)%ptr%segments(1)%z, "rotate_mtln: z should be rotated with mpidir 2")
    call expect_eq_int(test_err, 2, this%mtln%cables(1)%ptr%segments(1)%orientation, "rotate_mtln: orientation should be rotated with mpidir 2")

    cable%segments(1)%x = 1
    cable%segments(1)%y = 2 
    cable%segments(1)%z = 3
    cable%segments(1)%orientation = 1
    this%mtln%cables(1)%ptr => cable

    mpidir = 1
    call rotate_mtln(this, mpidir)
    call expect_eq_int(test_err, 2, this%mtln%cables(1)%ptr%segments(1)%x, "rotate_mtln: x should be rotated with mpidir 1")
    call expect_eq_int(test_err, 3, this%mtln%cables(1)%ptr%segments(1)%y, "rotate_mtln: y should be rotated with mpidir 1")
    call expect_eq_int(test_err, 1, this%mtln%cables(1)%ptr%segments(1)%z, "rotate_mtln: z should be rotated with mpidir 1")
    call expect_eq_int(test_err, 3, this%mtln%cables(1)%ptr%segments(1)%orientation, "rotate_mtln: orientation should be rotated with mpidir 1")

    deallocate(cable%segments)
    deallocate(this%mtln%cables)
    deallocate(this%mtln)
    err = test_err

end function

end module