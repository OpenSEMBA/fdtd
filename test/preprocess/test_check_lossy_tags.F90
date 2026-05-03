integer function test_check_lossy_tags_no_duplicates() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(LossyThinSurface_t) :: comp
   type(LossyThinSurface_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c(2))
   comp%nc = 2
   comp%c(1)%tag = 'lossyA'
   comp%c(2)%tag = 'lossyB'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkLossyTags(comp, prev, 0, numertag, tagtype, 1)

   if (numertag /= 2) err = err + 1
   if (trim(tagtype%tag(1)) /= 'lossyA') err = err + 1
   if (trim(tagtype%tag(2)) /= 'lossyB') err = err + 1

   deallocate(comp%c, prev, tagtype%tag)
end function

integer function test_check_lossy_tags_self_duplicate() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(LossyThinSurface_t) :: comp
   type(LossyThinSurface_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c(3))
   comp%nc = 3
   comp%c(1)%tag = 'lossyA'
   comp%c(2)%tag = 'lossyA'   ! duplicate of first
   comp%c(3)%tag = 'lossyC'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkLossyTags(comp, prev, 0, numertag, tagtype, 1)

   ! lossyA (once) + lossyC = 2
   if (numertag /= 2) err = err + 1

   deallocate(comp%c, prev, tagtype%tag)
end function

integer function test_check_lossy_tags_precounting_zero() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(LossyThinSurface_t) :: comp
   type(LossyThinSurface_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c(2))
   comp%nc = 2
   comp%c(1)%tag = 'lossyX'
   comp%c(2)%tag = 'lossyY'
   allocate(prev(0))

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0
   tagtype%tag(:) = ''

   call checkLossyTags(comp, prev, 0, numertag, tagtype, 0)

   if (numertag /= 2) err = err + 1
   ! tags should NOT be written since precounting=0
   if (trim(tagtype%tag(1)) /= '') err = err + 1

   deallocate(comp%c, prev, tagtype%tag)
end function

integer function test_check_lossy_tags_with_prev_duplicate() bind(C) result(err)
   use preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   type(LossyThinSurface_t) :: comp
   type(LossyThinSurface_t), allocatable :: prev(:)
   type(tagtype_t) :: tagtype
   integer :: numertag

   err = 0
   allocate(comp%c(2))
   comp%nc = 2
   comp%c(1)%tag = 'lossyA'
   comp%c(2)%tag = 'newTag'
   allocate(prev(1))
   allocate(prev(1)%c(1))
   prev(1)%nc = 1
   prev(1)%c(1)%tag = 'lossyA'   ! duplicate with prev

   numertag = 0
   allocate(tagtype%tag(10))
   tagtype%numertags = 0

   call checkLossyTags(comp, prev, 1, numertag, tagtype, 1)

   ! lossyA is duplicate with prev (foundDuplicate → numertag stays 1 effectively)
   ! newTag is new → numertag = 1
   if (numertag /= 1) err = err + 1
   if (trim(tagtype%tag(1)) /= 'newTag') err = err + 1

   deallocate(comp%c, tagtype%tag)
   deallocate(prev(1)%c)
   deallocate(prev)
end function
