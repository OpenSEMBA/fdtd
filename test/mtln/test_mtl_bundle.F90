integer function test_mtl_bundle_init() bind(C) result(error_cnt)
    use mtl_m
    use mtl_bundle_m
    implicit none

    type(mtl_t) :: mtl_out, mtl_in
    type(mtl_bundle_t) :: bundle
    type(transmission_line_level_t), dimension(2) :: levels

    real(kind=rkind),dimension(1,1) :: l1 = reshape( source = [ 4.4712610E-07_rkind ], shape = [ 1,1 ] )
    real(kind=rkind),dimension(1,1) :: c1 = reshape( source = [ 2.242e-10_rkind ], shape = [ 1,1 ] )
    real(kind=rkind),dimension(1,1) :: r1 = reshape( source = [ 0.0_rkind ], shape = [ 1,1 ] )
    real(kind=rkind),dimension(1,1) :: g1 = reshape( source = [ 0.0_rkind ], shape = [ 1,1 ] )

    integer :: i
    real(kind=rkind), dimension(5) :: step_size = [20.0_rkind, 20.0_rkind, 20.0_rkind, 20.0_rkind, 20.0_rkind]
    type(segment_t), allocatable, dimension(:) :: segments

    type(transfer_impedance_per_meter_t):: Zt
    type(multipolar_expansion_t), dimension(:), allocatable:: mE
    Zt%inductive_term = 0.0
    Zt%resistive_term = 0.0
    allocate(Zt%poles(0), Zt%residues(0))
    allocate(mE(0))

    error_cnt = 0
    allocate(segments(5))
    do i = 1, 5
        segments(i)%x = i
        segments(i)%y = 1 
        segments(i)%z = 1
        segments(i)%orientation = 1
    end do

    mtl_in   =  mtl_shielded(&
                    l1, c1, r1, g1, &
                    step_size, &
                    name = "line_in", &
                    segments = segments, &
                    dt = 1e-11_RKIND_TIEMPO, & 
                    parent_name = "line_out", &
                    conductor_in_parent = 1, &
                    transfer_impedance = Zt)
    mtl_out   = mtl_unshielded(&
                    l1, c1, r1, g1, &
                    step_size, &
                    name = "line_out", &
                    segments = segments, &
                    dt = 1e-11_RKIND_TIEMPO, &
                    multipolar_expansion = mE, &
                    radius = 0.0_rkind )

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

integer function test_mtl_bundle_generator() bind(C) result(error_cnt)
    use FDETYPES_m, only: RKIND
    use mtl_bundle_m, only: mtl_bundle_t
    use mtln_types_m, only: SOURCE_TYPE_CURRENT
    implicit none

    type(mtl_bundle_t) :: bundle

    error_cnt = 0
    allocate(bundle%generators(0))
    allocate(bundle%rpul(3, 1, 1), source=0.0_rkind)
    allocate(bundle%du(3, 1, 1), source=2.0_rkind)

    call bundle%addGenerator(2, 1, SOURCE_TYPE_CURRENT, 4.0_rkind, &
                             './testData/cases/planewave/gauss_1GHz.exc')

    if (size(bundle%generators) /= 1) error_cnt = error_cnt + 1
    if (bundle%rpul(2, 1, 1) /= 2.0_rkind) error_cnt = error_cnt + 1
end function
