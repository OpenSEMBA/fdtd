integer function test_mtl_wrong_dt() bind(C) result(error_cnt)

    use mtl_mod
    use mtln_testingTools_mod
    implicit none


    type(mtl_t) :: line
    real :: dt = 1.0 
    line = buildLineWithNConductors(2,'line0', dt = dt, type = "shielded")
    error_cnt = 0
    if (line%dt == dt) then 
        error_cnt = error_cnt + 1
    end if

end function

integer function test_mtl_init_homogeneous() bind(C) result(error_cnt) 
    use mtl_mod
    use mtln_testingTools_mod
    implicit none

    character(len=*), parameter :: name = 'line0'
    integer :: i,j

    
    real,dimension(2,2) :: lpul = reshape( &
        source = [ 4.4712610E-07, 1.4863653E-07, 1.4863653E-07, 4.4712610E-07 ], shape = [ 2,2 ] )
    real,dimension(2,2) :: cpul = reshape( &
        source = [ 2.242e-10, -7.453e-11,-7.453e-11, 2.242e-10 ], shape = [ 2,2 ] )
    real,dimension(2,2) :: rpul = reshape( &
        source = [ 0.0, 0.0, 0.0, 0.0 ], shape = [ 2,2 ] )
    real,dimension(2,2) :: gpul = reshape( &
        source = [ 0.0, 0.0, 0.0, 0.0 ], shape = [ 2,2 ] )
    real, dimension(5) :: step_size = [20.0, 20.0, 20.0, 20.0, 20.0]
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
    line = mtl_shielded(lpul, cpul, rpul, gpul, step_size, name, segments=segments, dt = 1e-12, parent_name ="p", conductor_in_parent = 1, transfer_impedance = Zt)
    call comparePULMatrices(error_cnt, line%lpul, lpul)
    call comparePULMatrices(error_cnt, line%cpul, cpul)
    call comparePULMatrices(error_cnt, line%rpul, rpul)
    call comparePULMatrices(error_cnt, line%gpul, gpul)
    line = mtl_unshielded(lpul, cpul, rpul, gpul, step_size, name, segments=segments, dt = 1e-12, multipolar_expansion = mE)
    call comparePULMatrices(error_cnt, line%lpul, lpul)
    call comparePULMatrices(error_cnt, line%cpul, cpul)
    call comparePULMatrices(error_cnt, line%rpul, rpul)
    call comparePULMatrices(error_cnt, line%gpul, gpul)

end function

integer function test_mtl_time_step() bind(C) result(error_cnt)    

    use mtl_mod
    use mtln_testingTools_mod

    implicit none

    real, dimension(5,2) :: phase_velocities
    real :: time_step, max_vel


    type(mtl_t) :: line 
    line = buildLineWithNConductors(2, "line0", dt = 1e-6, type = "unshielded")

    error_cnt = 0

    phase_velocities = line%getPhaseVelocities()
    max_vel = maxval(phase_velocities)
    time_step = line%getMaxTimeStep()
    !expected
    if (.not.(checkNear(phase_velocities(1,1),1.05900008e+08, 0.01))) then
        error_cnt = error_cnt +1
    end if
    if (.not.(checkNear(phase_velocities(1,2), 1.05900010e+08, 0.01))) then
        error_cnt = error_cnt +1
    end if
    if (.not.(checkNear(time_step, 1.888573951383424e-07, 0.01))) then
        error_cnt = error_cnt +1
    end if

end function