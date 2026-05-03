integer function test_check_dielectric_tags_c1p_no_duplicates() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(Dielectric_t) :: comp
   type(Dielectric_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c1P(2))
   comp%n_C1P = 2
   comp%c1P(1)%tag = 'tagA'
   comp%c1P(2)%tag = 'tagB'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkDielectricComponentTags(comp, prev, 0, 'c1P', numertag, tagtype, 1, 'error')

   if (numertag /= 2) err = err + 1
   if (trim(tagtype%tag(1)) /= 'tagA') err = err + 1
   if (trim(tagtype%tag(2)) /= 'tagB') err = err + 1

   deallocate(comp%c1P, prev, tagtype%tag)
end function

integer function test_check_dielectric_tags_c1p_self_duplicate() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(Dielectric_t) :: comp
   type(Dielectric_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c1P(3))
   comp%n_C1P = 3
   comp%c1P(1)%tag = 'tagA'
   comp%c1P(2)%tag = 'tagA'   ! duplicate of first
   comp%c1P(3)%tag = 'tagB'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkDielectricComponentTags(comp, prev, 0, 'c1P', numertag, tagtype, 1, 'error')

   ! tagA counted once (duplicate skipped), tagB counted once → numertag = 2
   if (numertag /= 2) err = err + 1

   deallocate(comp%c1P, prev, tagtype%tag)
end function

integer function test_check_dielectric_tags_precounting_zero() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(Dielectric_t) :: comp
   type(Dielectric_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c1P(2))
   comp%n_C1P = 2
   comp%c1P(1)%tag = 'tagA'
   comp%c1P(2)%tag = 'tagB'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0
   tagtype%tag(:) = ''

   ! precounting=0 means tags are not stored, only counted
   call checkDielectricComponentTags(comp, prev, 0, 'c1P', numertag, tagtype, 0, 'error')

   if (numertag /= 2) err = err + 1
   ! tags should NOT be written since precounting=0
   if (trim(tagtype%tag(1)) /= '') err = err + 1

   deallocate(comp%c1P, prev, tagtype%tag)
end function

integer function test_check_dielectric_tags_c2p() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(Dielectric_t) :: comp
   type(Dielectric_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c2P(2))
   comp%n_C2P = 2
   comp%c2P(1)%tag = 'c2tag1'
   comp%c2P(2)%tag = 'c2tag2'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkDielectricComponentTags(comp, prev, 0, 'c2P', numertag, tagtype, 1, 'error')

   if (numertag /= 2) err = err + 1
   if (trim(tagtype%tag(1)) /= 'c2tag1') err = err + 1
   if (trim(tagtype%tag(2)) /= 'c2tag2') err = err + 1

   deallocate(comp%c2P, prev, tagtype%tag)
end function
