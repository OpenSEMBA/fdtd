integer function test_mtl_wrong_dt() bind(C) result(error_cnt)

    use mtl_m
    use mtln_testingTools_mod
    implicit none


    type(mtl_t) :: line
    real(kind=RKIND_tiempo) :: dt = 1.0 
    line = buildLineWithNConductors(2,'line0', dt = dt, type = "shielded")
    error_cnt = 0
    if (line%dt == dt) then 
        error_cnt = error_cnt + 1
    end if

end function

integer function test_mtl_init_homogeneous() bind(C) result(error_cnt) 
    use mtl_m
    use mtln_testingTools_mod
    implicit none

    character(len=*), parameter :: name = 'line0'
    integer :: i,j

    
    real(kind=rkind),dimension(2,2) :: lpul = reshape( &
        source = [ 4.4712610E-07_rkind, 1.4863653E-07_rkind, 1.4863653E-07_rkind, 4.4712610E-07_rkind ], shape = [ 2,2 ] )
    real(kind=rkind),dimension(2,2) :: cpul = reshape( &
        source = [ 2.242e-10_rkind, -7.453e-11_rkind,-7.453e-11_rkind, 2.242e-10_rkind ], shape = [ 2,2 ] )
    real(kind=rkind),dimension(2,2) :: rpul = reshape( &
        source = [ 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind ], shape = [ 2,2 ] )
    real(kind=rkind),dimension(2,2) :: gpul = reshape( &
        source = [ 0.0_rkind, 0.0_rkind, 0.0_rkind, 0.0_rkind ], shape = [ 2,2 ] )
    real(kind=rkind), dimension(5) :: step_size = [20.0_rkind, 20.0_rkind, 20.0_rkind, 20.0_rkind, 20.0_rkind]
    type(segment_t), dimension(:), allocatable :: segments

    type(mtl_t) :: line 
    type(transfer_impedance_per_meter_t):: Zt
    type(multipolar_expansion_t), dimension(:), allocatable:: mE

    Zt%inductive_term = 0.0
    Zt%resistive_term = 0.0
    allocate(Zt%poles(0), Zt%residues(0))
    allocate(mE(0))

    allocate(segments(5))
    do i = 1, 5
        segments(i)%x = i
        segments(i)%y = 1
        segments(i)%z = 1
        segments(i)%orientation = DIRECTION_X_POS
    end do

    error_cnt = 0
    line = mtl_shielded(lpul, cpul, rpul, gpul, step_size, name, segments=segments, dt = 1e-12_RKIND_TIEMPO, parent_name ="p", conductor_in_parent = 1, transfer_impedance = Zt)
    call comparePULMatrices(error_cnt, line%lpul, lpul)
    call comparePULMatrices(error_cnt, line%cpul, cpul)
    call comparePULMatrices(error_cnt, line%rpul, rpul)
    call comparePULMatrices(error_cnt, line%gpul, gpul)
    line = mtl_unshielded(lpul, cpul, rpul, gpul, step_size, name, segments=segments, dt = 1e-12_RKIND_TIEMPO, multipolar_expansion = mE, radius = 0.0_rkind)
    call comparePULMatrices(error_cnt, line%lpul, lpul)
    call comparePULMatrices(error_cnt, line%cpul, cpul)
    call comparePULMatrices(error_cnt, line%rpul, rpul)
    call comparePULMatrices(error_cnt, line%gpul, gpul)

end function

integer function test_mtl_time_step() bind(C) result(error_cnt)    

    use mtl_m
    use mtln_testingTools_mod

    implicit none

    real(kind=rkind), dimension(5,2) :: phase_velocities
    real(kind=rkind) :: time_step, max_vel


    type(mtl_t) :: line 
    line = buildLineWithNConductors(2, "line0", dt = 1e-6_rkind_tiempo, type = "unshielded")

    error_cnt = 0

    phase_velocities = line%getPhaseVelocities()
    max_vel = maxval(phase_velocities)
    time_step = line%getMaxTimeStep()
    !expected
    if (.not.(checkNear(phase_velocities(1,1),1.05900008e+08_rkind, 0.01_rkind))) then
        error_cnt = error_cnt +1
    end if
    if (.not.(checkNear(phase_velocities(1,2), 1.05900010e+08_rkind, 0.01_rkind))) then
        error_cnt = error_cnt +1
    end if
    if (.not.(checkNear(time_step, 1.888573951383424e-07_rkind, 0.01_rkind))) then
        error_cnt = error_cnt +1
    end if

end function

integer function test_mtl_inactive_mpi_slice() bind(C) result(error_cnt)
    use mtl_m
    use mtln_types_m, only: DIRECTION_X_POS
    implicit none

    type(mtl_t) :: shielded_line, unshielded_line
    type(transfer_impedance_per_meter_t) :: zt
    type(multipolar_expansion_t), dimension(:), allocatable :: multipolar_expansion
    real(kind=rkind), dimension(1,1) :: lpul, cpul, rpul, gpul
    real(kind=rkind), dimension(2) :: step_size
    type(segment_t), dimension(:), allocatable :: segments
#ifdef CompileWithMPI
    integer(kind=4), dimension(:,:), allocatable :: layer_indices
#endif
    integer :: i
    real(kind=rkind_tiempo), parameter :: dt = 1e-11_rkind_tiempo

    error_cnt = 0
#ifdef CompileWithMPI
    lpul = reshape([4.4712610e-7_rkind], [1,1])
    cpul = reshape([2.242e-10_rkind], [1,1])
    rpul = 0.0_rkind
    gpul = 0.0_rkind
    step_size = [20.0_rkind, 20.0_rkind]
    zt%inductive_term = 0.0_rkind
    zt%resistive_term = 0.0_rkind
    allocate(zt%poles(0), zt%residues(0), multipolar_expansion(0), segments(2), layer_indices(0,2))
    do i = 1, size(segments)
        segments(i)%x = i
        segments(i)%y = 1
        segments(i)%z = 1
        segments(i)%orientation = DIRECTION_X_POS
    end do

    shielded_line = mtl_shielded(lpul, cpul, rpul, gpul, step_size, 'inactive', segments, dt, &
                                 'parent', 1, zt, layer_indices, .false.)
    unshielded_line = mtl_unshielded(lpul, cpul, rpul, gpul, step_size, 'inactive', segments, dt, &
                                     multipolar_expansion, 0.0_rkind, layer_indices, .false.)

    if (shielded_line%bundle_in_layer) error_cnt = error_cnt + 1
    if (size(shielded_line%step_size) /= 0) error_cnt = error_cnt + 1
    if (shielded_line%dt /= dt) error_cnt = error_cnt + 1
    if (unshielded_line%bundle_in_layer) error_cnt = error_cnt + 1
    if (size(unshielded_line%step_size) /= 0) error_cnt = error_cnt + 1
    if (unshielded_line%dt /= dt) error_cnt = error_cnt + 1
#endif
end function

integer function test_mtl_replicated_mpi_probe() bind(C) result(error_cnt)
    use probes_m, only: probe_t, probeCtor
    use mtln_types_m, only: PROBE_TYPE_CURRENT
    use FDETYPES_m, only: RKIND, RKIND_TIEMPO
    implicit none

    type(probe_t) :: probe
    character(len=:), allocatable :: name
#ifdef CompileWithMPI
    integer(kind=4), dimension(:,:), allocatable :: layer_indices
#endif

    error_cnt = 0
#ifdef CompileWithMPI
    allocate(layer_indices(0,0))
    name = "replicated"
    probe = probeCtor(2, PROBE_TYPE_CURRENT, 1.0e-12_RKIND_TIEMPO, name, &
                      [0.0_RKIND, 0.0_RKIND, 0.0_RKIND], layer_indices)
    if (.not. probe%in_layer) error_cnt = error_cnt + 1
    if (probe%index /= 2) error_cnt = error_cnt + 1
#endif
end function
