integer function test_searchtag_found() bind(C) result(err)
   use preprocess_tags_m
   use FDETYPES_m
   implicit none
   type(tagtype_t) :: tagtype
   integer :: idx

   err = 0
   tagtype%numertags = 3
   allocate(tagtype%tag(3))
   tagtype%tag(1) = 'alpha'
   tagtype%tag(2) = 'beta'
   tagtype%tag(3) = 'gamma'

   idx = searchtag(tagtype, 'beta')
   if (idx /= 2) err = err + 1

   idx = searchtag(tagtype, 'alpha')
   if (idx /= 1) err = err + 1

   idx = searchtag(tagtype, 'gamma')
   if (idx /= 3) err = err + 1

   deallocate(tagtype%tag)
end function

integer function test_searchtag_notfound() bind(C) result(err)
   use preprocess_tags_m
   use FDETYPES_m
   implicit none
   type(tagtype_t) :: tagtype
   integer :: idx

   err = 0
   tagtype%numertags = 2
   allocate(tagtype%tag(2))
   tagtype%tag(1) = 'alpha'
   tagtype%tag(2) = 'beta'

   idx = searchtag(tagtype, 'delta')
   if (idx /= -1) err = err + 1

   deallocate(tagtype%tag)
end function

integer function test_searchtag_empty() bind(C) result(err)
   use preprocess_tags_m
   use FDETYPES_m
   implicit none
   type(tagtype_t) :: tagtype
   integer :: idx

   err = 0
   tagtype%numertags = 0
   allocate(tagtype%tag(0))

   idx = searchtag(tagtype, 'alpha')
   if (idx /= -1) err = err + 1

   deallocate(tagtype%tag)
end function
