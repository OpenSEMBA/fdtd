! Unit tests for src_main_pub/healer.F90
! Tests SortInitEndWithIncreasingOrder and Readjust.

integer function test_sort_all_swapped() bind(C, name="test_sort_all_swapped") result(status)
    use CreateMatrices_m, only: SortInitEndWithIncreasingOrder
    use FDETYPES_m, only: XYZlimit_t
    implicit none
    type(XYZlimit_t) :: p

    status = 0
    p%XI = 10; p%XE = 5
    p%YI = 20; p%YE = 2
    p%ZI = 15; p%ZE = 3

    call SortInitEndWithIncreasingOrder(p)

    if (p%XI /= 5 .or. p%XE /= 10) then
        print *, "test_sort_all_swapped FAILED: X not sorted"
        status = 1
    end if
    if (p%YI /= 2 .or. p%YE /= 20) then
        print *, "test_sort_all_swapped FAILED: Y not sorted"
        status = 1
    end if
    if (p%ZI /= 3 .or. p%ZE /= 15) then
        print *, "test_sort_all_swapped FAILED: Z not sorted"
        status = 1
    end if
end function test_sort_all_swapped

integer function test_sort_already_ordered() bind(C, name="test_sort_already_ordered") result(status)
    use CreateMatrices_m, only: SortInitEndWithIncreasingOrder
    use FDETYPES_m, only: XYZlimit_t
    implicit none
    type(XYZlimit_t) :: p

    status = 0
    p%XI = 1; p%XE = 5
    p%YI = 2; p%YE = 20
    p%ZI = 3; p%ZE = 15

    call SortInitEndWithIncreasingOrder(p)

    if (p%XI /= 1 .or. p%XE /= 5) then
        print *, "test_sort_already_ordered FAILED: X changed unexpectedly"
        status = 1
    end if
    if (p%YI /= 2 .or. p%YE /= 20) then
        print *, "test_sort_already_ordered FAILED: Y changed unexpectedly"
        status = 1
    end if
    if (p%ZI /= 3 .or. p%ZE /= 15) then
        print *, "test_sort_already_ordered FAILED: Z changed unexpectedly"
        status = 1
    end if
end function test_sort_already_ordered

integer function test_sort_x_reversed_yz_ok() bind(C, name="test_sort_x_reversed_yz_ok") result(status)
    use CreateMatrices_m, only: SortInitEndWithIncreasingOrder
    use FDETYPES_m, only: XYZlimit_t
    implicit none
    type(XYZlimit_t) :: p

    status = 0
    p%XI = 9; p%XE = 1   ! reversed
    p%YI = 3; p%YE = 7   ! already ordered
    p%ZI = 4; p%ZE = 4   ! equal

    call SortInitEndWithIncreasingOrder(p)

    if (p%XI /= 1 .or. p%XE /= 9) then
        print *, "test_sort_x_reversed_yz_ok FAILED: X not sorted"
        status = 1
    end if
    if (p%YI /= 3 .or. p%YE /= 7) then
        print *, "test_sort_x_reversed_yz_ok FAILED: Y changed unexpectedly"
        status = 1
    end if
    if (p%ZI /= 4 .or. p%ZE /= 4) then
        print *, "test_sort_x_reversed_yz_ok FAILED: Z changed unexpectedly"
        status = 1
    end if
end function test_sort_x_reversed_yz_ok

integer function test_readjust_grow() bind(C, name="test_readjust_grow") result(status)
    use CreateMatrices_m, only: Readjust
    use FDETYPES_m, only: MediaData_t, set_priorities
    implicit none
    type(MediaData_t), pointer, dimension(:) :: med
    integer(kind=4) :: numMedia, newNumMedia

    status = 0
    call set_priorities(.false., .false., .false.)

    numMedia = 2
    newNumMedia = 5
    allocate(med(0:numMedia))
    med(0)%Epr = 1.0
    med(1)%Epr = 2.0
    med(2)%Epr = 3.0

    call Readjust(numMedia, med, newNumMedia)

    if (numMedia /= 5) then
        print *, "test_readjust_grow FAILED: NumMedia not updated, got", numMedia
        status = 1
    end if
    if (size(med) /= 6) then  ! 0:5 = 6 elements
        print *, "test_readjust_grow FAILED: med size wrong, got", size(med)
        status = 1
    end if
    ! Original entries preserved
    if (med(0)%Epr /= 1.0 .or. med(1)%Epr /= 2.0 .or. med(2)%Epr /= 3.0) then
        print *, "test_readjust_grow FAILED: original entries corrupted"
        status = 1
    end if
    ! New entries initialised to sentinel -1
    if (med(3)%Epr /= -1.0 .or. med(4)%Epr /= -1.0 .or. med(5)%Epr /= -1.0) then
        print *, "test_readjust_grow FAILED: new entries not initialised to -1"
        status = 1
    end if

    deallocate(med)
end function test_readjust_grow

integer function test_readjust_shrink() bind(C, name="test_readjust_shrink") result(status)
    use CreateMatrices_m, only: Readjust
    use FDETYPES_m, only: MediaData_t, set_priorities
    implicit none
    type(MediaData_t), pointer, dimension(:) :: med
    integer(kind=4) :: numMedia, newNumMedia

    status = 0
    call set_priorities(.false., .false., .false.)

    numMedia = 5
    newNumMedia = 2
    allocate(med(0:numMedia))
    med(0)%Epr = 10.0
    med(1)%Epr = 11.0
    med(2)%Epr = 12.0
    med(3)%Epr = 13.0
    med(4)%Epr = 14.0
    med(5)%Epr = 15.0

    call Readjust(numMedia, med, newNumMedia)

    if (numMedia /= 2) then
        print *, "test_readjust_shrink FAILED: NumMedia not updated, got", numMedia
        status = 1
    end if
    if (size(med) /= 3) then  ! 0:2 = 3 elements
        print *, "test_readjust_shrink FAILED: med size wrong, got", size(med)
        status = 1
    end if
    ! First entries preserved
    if (med(0)%Epr /= 10.0 .or. med(1)%Epr /= 11.0 .or. med(2)%Epr /= 12.0) then
        print *, "test_readjust_shrink FAILED: remaining entries corrupted"
        status = 1
    end if

    deallocate(med)
end function test_readjust_shrink

integer function test_readjust_same_size() bind(C, name="test_readjust_same_size") result(status)
    use CreateMatrices_m, only: Readjust
    use FDETYPES_m, only: MediaData_t, set_priorities
    implicit none
    type(MediaData_t), pointer, dimension(:) :: med
    integer(kind=4) :: numMedia, newNumMedia

    status = 0
    call set_priorities(.false., .false., .false.)

    numMedia = 3
    newNumMedia = 3
    allocate(med(0:numMedia))
    med(0)%Epr = 5.0
    med(1)%Epr = 6.0
    med(2)%Epr = 7.0
    med(3)%Epr = 8.0

    call Readjust(numMedia, med, newNumMedia)

    if (numMedia /= 3) then
        print *, "test_readjust_same_size FAILED: NumMedia changed unexpectedly, got", numMedia
        status = 1
    end if
    if (size(med) /= 4) then  ! 0:3 = 4 elements
        print *, "test_readjust_same_size FAILED: med size wrong, got", size(med)
        status = 1
    end if
    if (med(0)%Epr /= 5.0 .or. med(1)%Epr /= 6.0 .or. med(2)%Epr /= 7.0 .or. med(3)%Epr /= 8.0) then
        print *, "test_readjust_same_size FAILED: entries corrupted"
        status = 1
    end if

    deallocate(med)
end function test_readjust_same_size
