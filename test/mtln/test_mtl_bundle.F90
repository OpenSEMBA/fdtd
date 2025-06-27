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

    integer :: i
    real, dimension(5) :: step_size = [20.0, 20.0, 20.0, 20.0, 20.0]
    type(direction_t), allocatable, dimension(:) :: segments

    error_cnt = 0
    allocate(segments(5))
    do i = 1, 5
        segments(i)%x = i
        segments(i)%y = 1 
        segments(i)%z = 1
        segments(i)%orientation = 1
    end do

    mtl_in   =  mtl_t(l1, c1, r1, g1, step_size, "line_in", conductor_in_parent = 1, segments = segments)
    mtl_out   = mtl_t(l1, c1, r1, g1, step_size, "line_out", segments = segments)

    allocate(levels(1)%lines(1))
    allocate(levels(2)%lines(1))
    levels(1)%lines = mtl_out
    levels(2)%lines = mtl_in
    bundle = mtldCtor(levels, name="bundle")

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