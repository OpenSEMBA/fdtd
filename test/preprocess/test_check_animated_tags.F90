integer function test_check_animated_tags_c1p_no_duplicates() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(ANISOTROPICbody_t) :: comp
   type(ANISOTROPICbody_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c1P(2))
   comp%n_C1P = 2
   comp%c1P(1)%tag = 'animA'
   comp%c1P(2)%tag = 'animB'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkAnimatedComponentTags(comp, prev, 0, 'c1P', numertag, tagtype, 1, 'error')

   if (numertag /= 2) err = err + 1
   if (trim(tagtype%tag(1)) /= 'animA') err = err + 1
   if (trim(tagtype%tag(2)) /= 'animB') err = err + 1

   deallocate(comp%c1P, prev, tagtype%tag)
end function

integer function test_check_animated_tags_c1p_self_duplicate() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(ANISOTROPICbody_t) :: comp
   type(ANISOTROPICbody_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c1P(3))
   comp%n_C1P = 3
   comp%c1P(1)%tag = 'animA'
   comp%c1P(2)%tag = 'animA'   ! duplicate
   comp%c1P(3)%tag = 'animC'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkAnimatedComponentTags(comp, prev, 0, 'c1P', numertag, tagtype, 1, 'error')

   ! animA once + animC = numertag 2
   if (numertag /= 2) err = err + 1

   deallocate(comp%c1P, prev, tagtype%tag)
end function

integer function test_check_animated_tags_c2p() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(ANISOTROPICbody_t) :: comp
   type(ANISOTROPICbody_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c2P(1))
   comp%n_C2P = 1
   comp%c2P(1)%tag = 'c2animX'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkAnimatedComponentTags(comp, prev, 0, 'c2P', numertag, tagtype, 1, 'error')

   if (numertag /= 1) err = err + 1
   if (trim(tagtype%tag(1)) /= 'c2animX') err = err + 1

   deallocate(comp%c2P, prev, tagtype%tag)
end function
