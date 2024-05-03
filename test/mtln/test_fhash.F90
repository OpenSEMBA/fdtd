integer function test_fhash() bind(C) result(error_cnt)

    use mtl_bundle_mod
    use fhash, only: fhash_tbl_t, key=>fhash_key
    

    implicit none

    type(fhash_tbl_t) :: tbl
    type(mtl_bundle_t) :: bundle
    class(*), allocatable :: bundle_from_hash_1, bundle_from_hash_2
    class(*), pointer :: ptr 
    integer :: stat, n_cond
    error_cnt = 0

    bundle%name = "bundle1"
    bundle%number_of_conductors = 5

    call tbl%set(key(bundle%name), value = bundle)

    call tbl%get_raw(key(bundle%name), bundle_from_hash_1,stat)
    
    select type(bundle_from_hash_1)
    type is (mtl_bundle_t)
        bundle_from_hash_1%number_of_conductors = 100
        call tbl%set(key(bundle_from_hash_1%name), value = bundle_from_hash_1)
    end select

    call tbl%get_raw(key(bundle%name), bundle_from_hash_2,stat)
    select type(bundle_from_hash_2)
    type is (mtl_bundle_t)
        n_cond = bundle_from_hash_2%number_of_conductors
    end select
    
    write(*,*) 'number of conductors ', n_cond
    if (n_cond /= 100) then
        error_cnt = error_cnt +1
    end if

end function


