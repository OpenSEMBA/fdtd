! Unit tests for preprocess_geom.F90 functions
! Following the pattern of test_geometry.F90

integer function test_searchtag() bind(C, name="test_searchtag") result(status)
    use, intrinsic :: iso_c_binding
    use FDETYPES_m
    use Preprocess_m
    implicit none
    type(tagtype_t) :: test_tags
    integer(c_int) :: result_idx
    
    status = 0
    
    ! Initialize tag array with test data
    test_tags%numertags = 5
    allocate(test_tags%tag(5))
    test_tags%tag(1) = 'tag_alpha'
    test_tags%tag(2) = 'tag_beta'
    test_tags%tag(3) = 'tag_gamma'
    test_tags%tag(4) = 'tag_delta'
    test_tags%tag(5) = 'tag_epsilon'
    
    ! Test 1: Search for non-existent tag (should return -1)
    result_idx = searchtag(test_tags, 'nonexistent')
    if (result_idx /= -1) then
        print *, 'test_searchtag FAILED: Expected -1 for non-existent tag, got', result_idx
        status = 1
    end if
    
    ! Test 2: Search for existing tag at index 1
    result_idx = searchtag(test_tags, 'tag_alpha')
    if (result_idx /= 1) then
        print *, 'test_searchtag FAILED: Expected 1 for tag_alpha, got', result_idx
        status = 1
    end if
    
    ! Test 3: Search for existing tag at middle index
    result_idx = searchtag(test_tags, 'tag_gamma')
    if (result_idx /= 3) then
        print *, 'test_searchtag FAILED: Expected 3 for tag_gamma, got', result_idx
        status = 1
    end if
    
    ! Test 4: Search for existing tag at last index
    result_idx = searchtag(test_tags, 'tag_epsilon')
    if (result_idx /= 5) then
        print *, 'test_searchtag FAILED: Expected 5 for tag_epsilon, got', result_idx
        status = 1
    end if
    
    ! Test 5: Search with leading/trailing whitespace (should be trimmed)
    result_idx = searchtag(test_tags, '  tag_beta  ')
    if (result_idx /= 2) then
        print *, 'test_searchtag FAILED: Expected 2 for tag_beta with whitespace, got', result_idx
        status = 1
    end if
    
    ! Test 6: Empty tag search
    result_idx = searchtag(test_tags, '')
    if (result_idx /= -1) then
        print *, 'test_searchtag FAILED: Expected -1 for empty tag, got', result_idx
        status = 1
    end if
    
    ! Test 7: Case sensitivity test
    result_idx = searchtag(test_tags, 'TAG_ALPHA')
    if (result_idx /= -1) then
        print *, 'test_searchtag FAILED: Expected -1 for TAG_ALPHA (case sensitive), got', result_idx
        status = 1
    end if
    
    ! Cleanup
    deallocate(test_tags%tag)
    
    if (status == 0) then
        print *, 'test_searchtag PASSED'
    end if
end function test_searchtag

! Test with empty tagtype
integer function test_searchtag_empty() bind(C, name="test_searchtag_empty") result(status)
    use, intrinsic :: iso_c_binding
    use FDETYPES_m
    use Preprocess_m
    implicit none
    type(tagtype_t) :: empty_tags
    integer(c_int) :: result_idx
    
    status = 0
    
    ! Initialize empty tagtype
    empty_tags%numertags = 0
    
    ! Search should return -1 for empty tagtype
    result_idx = searchtag(empty_tags, 'any_tag')
    if (result_idx /= -1) then
        print *, 'test_searchtag_empty FAILED: Expected -1 for empty tagtype, got', result_idx
        status = 1
    end if
    
    if (status == 0) then
        print *, 'test_searchtag_empty PASSED'
    end if
end function test_searchtag_empty

! Test with single tag
integer function test_searchtag_single() bind(C, name="test_searchtag_single") result(status)
    use, intrinsic :: iso_c_binding
    use FDETYPES_m
    use Preprocess_m
    implicit none
    type(tagtype_t) :: single_tags
    integer(c_int) :: result_idx
    
    status = 0
    
    ! Initialize with single tag
    single_tags%numertags = 1
    allocate(single_tags%tag(1))
    single_tags%tag(1) = 'only_tag'
    
    ! Search for the tag
    result_idx = searchtag(single_tags, 'only_tag')
    if (result_idx /= 1) then
        print *, 'test_searchtag_single FAILED: Expected 1, got', result_idx
        status = 1
    end if
    
    ! Search for different tag
    result_idx = searchtag(single_tags, 'other_tag')
    if (result_idx /= -1) then
        print *, 'test_searchtag_single FAILED: Expected -1 for other_tag, got', result_idx
        status = 1
    end if
    
    ! Cleanup
    deallocate(single_tags%tag)
    
    if (status == 0) then
        print *, 'test_searchtag_single PASSED'
    end if
end function test_searchtag_single

! Test with tags containing special characters
integer function test_searchtag_special_chars() bind(C, name="test_searchtag_special_chars") result(status)
    use, intrinsic :: iso_c_binding
    use FDETYPES_m
    use Preprocess_m
    implicit none
    type(tagtype_t) :: special_tags
    integer(c_int) :: result_idx
    
    status = 0
    
    ! Initialize with special character tags
    special_tags%numertags = 4
    allocate(special_tags%tag(4))
    special_tags%tag(1) = 'tag_with_underscores'
    special_tags%tag(2) = 'tag-with-dashes'
    special_tags%tag(3) = 'tag.with.dots'
    special_tags%tag(4) = 'Tag123'
    
    ! Test each tag
    result_idx = searchtag(special_tags, 'tag_with_underscores')
    if (result_idx /= 1) then
        print *, 'test_searchtag_special_chars FAILED: Expected 1, got', result_idx
        status = 1
    end if
    
    result_idx = searchtag(special_tags, 'tag-with-dashes')
    if (result_idx /= 2) then
        print *, 'test_searchtag_special_chars FAILED: Expected 2, got', result_idx
        status = 1
    end if
    
    result_idx = searchtag(special_tags, 'tag.with.dots')
    if (result_idx /= 3) then
        print *, 'test_searchtag_special_chars FAILED: Expected 3, got', result_idx
        status = 1
    end if
    
    result_idx = searchtag(special_tags, 'Tag123')
    if (result_idx /= 4) then
        print *, 'test_searchtag_special_chars FAILED: Expected 4, got', result_idx
        status = 1
    end if
    
    ! Cleanup
    deallocate(special_tags%tag)
    
    if (status == 0) then
        print *, 'test_searchtag_special_chars PASSED'
    end if
end function test_searchtag_special_chars
