integer function test_mtl_wrong_dt() bind(C) result(error_cnt)

    use mtl_mod
    use mtln_testingTools_mod
    implicit none


    type(mtl_t) :: line
    real :: dt = 1.0 
    line = buildLineWithNConductors(2,'line0', dt = dt)
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
    type(external_field_segment_t), dimension(5) :: external_field_segments

    type(mtl_t) :: line 
    do i = 1, 5
        external_field_segments(i)%position =(/i,1,1/)            
        external_field_segments(i)%direction = DIRECTION_X_POS  
        external_field_segments(i)%field => null()
    end do

    error_cnt = 0
    line = mtl_t(lpul, cpul, rpul, gpul, step_size, name, external_field_segments=external_field_segments)
    call comparePULMatrices(error_cnt, line%lpul, lpul)
    call comparePULMatrices(error_cnt, line%cpul, cpul)
    call comparePULMatrices(error_cnt, line%rpul, rpul)
    call comparePULMatrices(error_cnt, line%gpul, gpul)
    


end function

integer function test_mtl_init_inhomogeneous() bind(C) result(error_cnt)    
    use mtl_mod
    use mtln_testingTools_mod

    implicit none

    character(len=*), parameter :: name = 'line0'
    integer :: i,j

    
    real,dimension(2,2) :: lpul1 = reshape( &
        source = [ 4.4712610E-07, 1.4863653E-07, 1.4863653E-07, 4.4712610E-07 ], shape = [ 2,2 ] )
    real,dimension(2,2) :: lpul2 = reshape( &
        source = [ 2.4712610E-07, 0.4863653E-07, 0.4863653E-07, 2.4712610E-07 ], shape = [ 2,2 ] )
    real,dimension(2,2,2) :: lpul
    
    real,dimension(2,2) :: cpul1 = reshape( &
        source = [ 1.242e-10, -6.453e-11,-6.453e-11, 1.242e-10 ], shape = [ 2,2 ] )
    real,dimension(2,2) :: cpul2 = reshape( &
        source = [ 2.242e-10, -7.453e-11,-7.453e-11, 2.242e-10 ], shape = [ 2,2 ] )
    real,dimension(2,2) :: cpul3 = reshape( &
        source = [ 3.242e-10, -8.453e-11,-8.453e-11, 3.242e-10 ], shape = [ 2,2 ] )
    real,dimension(3,2,2) :: cpul
    
    real,dimension(2,2) :: rpul0 = reshape( source = [ 0.0, 0.0, 0.0, 0.0 ], shape = [ 2,2 ] )
    real,dimension(2,2,2) :: rpul
    
    real,dimension(2,2) :: gpul0 = reshape( source = [ 0.0, 0.0, 0.0, 0.0 ], shape = [ 2,2 ] )
    real,dimension(3,2,2) :: gpul
    
    real, dimension(2) :: step_size = [50.0, 50.0]
    type(mtl_t) :: line 
    
    lpul(1,:,:) = lpul1
    lpul(2,:,:) = lpul2
    cpul(1,:,:) = cpul1
    cpul(2,:,:) = cpul2
    cpul(3,:,:) = cpul3
    rpul(1,:,:) = rpul0
    rpul(2,:,:) = rpul0
    gpul(1,:,:) = gpul0
    gpul(2,:,:) = gpul0
    gpul(3,:,:) = gpul0

    error_cnt = 0
    line = mtl_t(lpul, cpul, rpul, gpul, step_size, name)
    call comparePULMatricesIH(error_cnt, line%lpul, lpul)
    call comparePULMatricesIH(error_cnt, line%cpul, cpul)
    call comparePULMatricesIH(error_cnt, line%rpul, rpul)
    call comparePULMatricesIH(error_cnt, line%gpul, gpul)
    


end function

integer function test_mtl_time_step() bind(C) result(error_cnt)    

    use mtl_mod
    use mtln_testingTools_mod

    implicit none

    real, dimension(5,2) :: phase_velocities
    real :: time_step, max_vel


    type(mtl_t) :: line 
    line = buildLineWithNConductors(2, "line0")

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