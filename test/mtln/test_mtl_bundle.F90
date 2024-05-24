integer function test_mtl_bundle_init() bind(C) result(error_cnt)
    use mtl_mod
    use mtl_bundle_mod
    implicit none

    type(mtl_t) :: mtl_out, mtl_in
    type(mtl_bundle_t) :: bundle
    type(mtl_array_t), dimension(2) :: levels

    real,dimension(1,1) :: l1 = reshape( source = [ 4.4712610E-07 ], shape = [ 1,1 ] )
    real,dimension(1,1) :: c1 = reshape( source = [ 2.242e-10 ], shape = [ 1,1 ] )
    real,dimension(1,1) :: r1 = reshape( source = [ 0.0 ], shape = [ 1,1 ] )
    real,dimension(1,1) :: g1 = reshape( source = [ 0.0 ], shape = [ 1,1 ] )


    real, dimension(5) :: step_size = [20.0, 20.0, 20.0, 20.0, 20.0]
    type(external_field_segment_t), dimension(5) :: external_field_segments

    error_cnt = 0
    block
        integer :: i
        do i = 1, 5
            external_field_segments(i)%position = (/i, 1, 1/)
            external_field_segments(i)%direction = 1
            external_field_segments(i)%field => null()
        end do
    end block

    mtl_out   = mtl_t(l1, c1, r1, g1, step_size, "line_out", external_field_segments = external_field_segments)
    mtl_in   =  mtl_t(l1, c1, r1, g1, step_size, "line_in",  &
                      parent_name = "line_out", &
                      conductor_in_parent = 1, &
                      external_field_segments = external_field_segments)

    allocate(levels(1)%lines(1))
    allocate(levels(2)%lines(1))
    levels(1)%lines = [mtl_out]
    levels(2)%lines = [mtl_in]
    bundle = mtl_bundle_t(levels, name="bundle")

    if ((size(bundle%lpul,1) /= 5) .or. &
        (size(bundle%lpul,2) /= 2) .or. &
        (size(bundle%lpul,3) /= 2)) then 
        error_cnt = error_cnt + 1
    end if

    if ((bundle%lpul(1,1,1) /= mtl_out%lpul(1,1,1)) .or. &
        bundle%lpul(1,2,2) /= mtl_in%lpul(1,1,1)) then
        error_cnt = error_cnt + 1
    end if
    !check size of pul matrices and V I vectors


    

end function