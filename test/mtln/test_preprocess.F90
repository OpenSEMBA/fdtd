integer function test_preprocess_conductors_before_cable() bind(C) result(error_cnt)

    use mtl_mod 
    use preprocess_mod, only: conductorsInLevel, findConductorsBeforeCable
    use mtln_testingTools_mod

    type(mtl_t) :: line1, line2, line3_1, line3_2, line4
    type(transmission_line_bundle_t) :: line_bundle
    error_cnt = 0
  
    line1   = buildLineWithNConductors(1, "line1")
    line2   = buildLineWithNConductors(2, "line2",   parent_name = "line1",   conductor_in_parent = 1)
    line3_1 = buildLineWithNConductors(2, "line3_1", parent_name = "line2",   conductor_in_parent = 1)
    line3_2 = buildLineWithNConductors(2, "line3_2", parent_name = "line2",   conductor_in_parent = 2)
    line4   = buildLineWithNConductors(2, "line4",   parent_name = "line3_2", conductor_in_parent = 2)

    allocate(line_bundle%levels(4))
    line_bundle%levels(1)%lines = [line1]
    line_bundle%levels(2)%lines = [line2]
    line_bundle%levels(3)%lines = [line3_1, line3_2]
    line_bundle%levels(4)%lines = [line4]

    if (findConductorsBeforeCable("line2", line_bundle%levels(2)) /= 0) then 
        error_cnt = error_cnt + 1
    end if
    if (findConductorsBeforeCable("line3_1", line_bundle%levels(3)) /= 0) then 
        error_cnt = error_cnt + 1
    end if
    if (findConductorsBeforeCable("line3_2", line_bundle%levels(3)) /= 2) then 
        error_cnt = error_cnt + 1
    end if
    if (findConductorsBeforeCable("line4", line_bundle%levels(4)) /= 0) then 
        error_cnt = error_cnt + 1
    end if

end function

integer function test_preprocess_conductors_in_level() bind(C) result(error_cnt)

    use mtl_mod 
    use preprocess_mod, only: conductorsInLevel, findConductorsBeforeCable
    use mtln_testingTools_mod

    type(mtl_t) :: line1, line2, line3_1, line3_2, line4
    type(transmission_line_bundle_t) :: line_bundle

    error_cnt = 0
    
    line1   = buildLineWithNConductors(1, name= "line1")
    line2   = buildLineWithNConductors(2, name= "line2"  , parent_name = "line1",   conductor_in_parent = 1)
    line3_1 = buildLineWithNConductors(2, name= "line3_1", parent_name = "line2",   conductor_in_parent = 1)
    line3_2 = buildLineWithNConductors(2, name= "line3_2", parent_name = "line2",   conductor_in_parent = 2)
    line4   = buildLineWithNConductors(2, name= "line4"  , parent_name = "line3_2", conductor_in_parent = 2)

    allocate(line_bundle%levels(4))
    line_bundle%levels(1)%lines = [line1]
    line_bundle%levels(2)%lines = [line2]
    line_bundle%levels(3)%lines = [line3_1, line3_2]
    line_bundle%levels(4)%lines = [line4]

    
    if (all(conductorsInLevel(line_bundle) /= [1,2,4,2])) then 
        error_cnt = error_cnt + 1
    end if

end function

integer function test_preprocess_zt_conductor_ranges() bind(C) result(error_cnt)

    use mtl_mod 
    use mtl_bundle_mod 
    use preprocess_mod, only: conductorsInLevel, findOuterConductorNumber, findInnerConductorRange
    use mtln_testingTools_mod

    type :: range
        integer, dimension(:), allocatable :: idx
    end type

    type(mtl_t) :: line1, line2, line3_1, line3_2, line4
    type(transmission_line_bundle_t) :: line_bundle
    type(mtl_bundle_t) :: mtl_bundle

    integer, dimension(4) :: conductors_in_level
    integer :: i, j, cnt

    integer :: conductor_out
    integer, dimension(:), allocatable :: range_in

    integer, dimension(4) :: expected_out
    type(range), dimension(4) :: expected_in

    expected_out = [1,2,3,7]
    expected_in(1)%idx = [2,3]
    expected_in(2)%idx = [4,5]
    expected_in(3)%idx = [6,7]
    expected_in(4)%idx = [8,9]

    error_cnt = 0
    line1 =   buildLineWithNConductors(1, "line1")
    line2 =   buildLineWithNConductors(2, "line2",  parent_name= "line1", conductor_in_parent= 1)
    line3_1 = buildLineWithNConductors(2, "line3_1",parent_name= "line2", conductor_in_parent= 1)
    line3_2 = buildLineWithNConductors(2, "line3_2",parent_name= "line2", conductor_in_parent= 2)
    line4 =   buildLineWithNConductors(2, "line4",  parent_name= "line3_2", conductor_in_parent= 2)

    allocate(line_bundle%levels(4))
    line_bundle%levels(1)%lines = [line1]
    line_bundle%levels(2)%lines = [line2]
    line_bundle%levels(3)%lines = [line3_1, line3_2]
    line_bundle%levels(4)%lines = [line4]

    conductors_in_level = conductorsInLevel(line_bundle)

    mtl_bundle = mtl_bundle_t(line_bundle%levels, name="bundle1")

    cnt = 1
    do i = 2, size(line_bundle%levels)
        do j = 1, size(line_bundle%levels(i)%lines)

            conductor_out = findOuterConductorNumber(line_bundle%levels(i)%lines(j), &
                                                     line_bundle%levels(i-1), &
                                                     sum(conductors_in_level(1:i-2)))

            range_in = findInnerConductorRange(line_bundle%levels(i)%lines(j), &
                                               line_bundle%levels(i), &
                                               sum(conductors_in_level(1:i-1)))
            
            if (expected_out(cnt) /= conductor_out) then 
                error_cnt = error_cnt + 1
            end if
            if (.not.all(expected_in(cnt)%idx.eq.range_in)) then 
                error_cnt = error_cnt + 1
            end if
            cnt = cnt + 1
            ! write(*,*) 'out: ', conductor_out, ' in: [', range_in,']'
        end do
    end do  

end function

integer function test_preprocess_zt_conductor_ranges_2() bind(C) result(error_cnt)

    use mtl_mod 
    use mtl_bundle_mod 
    use preprocess_mod, only: conductorsInLevel, findOuterConductorNumber, findInnerConductorRange
    use mtln_testingTools_mod

    type :: range
        integer, dimension(:), allocatable :: idx
    end type

    type(mtl_t) :: line1, line2, line3_1, line3_2, line3_3, line4_1, line4_2, line4_3
    type(transmission_line_bundle_t) :: line_bundle
    type(mtl_bundle_t) :: mtl_bundle

    integer, dimension(4) :: conductors_in_level
    integer :: i, j, cnt

    integer :: conductor_out
    integer, dimension(:), allocatable :: range_in

    integer, dimension(7) :: expected_out
    type(range), dimension(7) :: expected_in

    expected_out = [1,2,3,4,7,8,9]
    expected_in(1)%idx = [2,3,4]
    expected_in(2)%idx = [5]
    expected_in(3)%idx = [6,7]
    expected_in(4)%idx = [8,9]
    expected_in(5)%idx = [10,11]
    expected_in(6)%idx = [12]
    expected_in(7)%idx = [13]

    error_cnt = 0
    line1   = buildLineWithNConductors(1, name = "line1")
    line2   = buildLineWithNConductors(3, name = "line2",   parent_name = "line1",   conductor_in_parent = 1)
    line3_1 = buildLineWithNConductors(1, name = "line3_1", parent_name = "line2", conductor_in_parent = 1)
    line3_2 = buildLineWithNConductors(2, name = "line3_2", parent_name = "line2", conductor_in_parent = 2)
    line3_3 = buildLineWithNConductors(2, name = "line3_3", parent_name = "line2", conductor_in_parent = 3)
    line4_1 = buildLineWithNConductors(2, name = "line4_1", parent_name = "line3_2", conductor_in_parent = 2)
    line4_2 = buildLineWithNConductors(1, name = "line4_2", parent_name = "line3_3", conductor_in_parent = 1)
    line4_3 = buildLineWithNConductors(1, name = "line4_3", parent_name = "line3_3", conductor_in_parent = 2)

    allocate(line_bundle%levels(4))
    line_bundle%levels(1)%lines = [line1]
    line_bundle%levels(2)%lines = [line2]
    line_bundle%levels(3)%lines = [line3_1, line3_2, line3_3]
    line_bundle%levels(4)%lines = [line4_1, line4_2, line4_3]

    conductors_in_level = conductorsInLevel(line_bundle)

    mtl_bundle = mtl_bundle_t(line_bundle%levels, name="bundle1")

    cnt = 1
    do i = 2, size(line_bundle%levels)
        do j = 1, size(line_bundle%levels(i)%lines)

            conductor_out = findOuterConductorNumber(line_bundle%levels(i)%lines(j), &
                                                     line_bundle%levels(i-1), &
                                                     sum(conductors_in_level(1:i-2)))

            range_in = findInnerConductorRange(line_bundle%levels(i)%lines(j), &
                                               line_bundle%levels(i), &
                                               sum(conductors_in_level(1:i-1)))
            
            if (expected_out(cnt) /= conductor_out) then 
                error_cnt = error_cnt + 1
            end if
            if (.not.all(expected_in(cnt)%idx.eq.range_in)) then 
                error_cnt = error_cnt + 1
            end if
            cnt = cnt + 1
            ! write(*,*) 'out: ', conductor_out, ' in: [', range_in,']'
        end do
    end do  


end function

