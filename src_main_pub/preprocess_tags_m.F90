module preprocess_tags_m
   use NFDETypes_m
   use FDETYPES_m
   implicit none
   private
   public :: checkDielectricComponentTags, checkAnimatedComponentTags, checkLossyTags, searchtag

contains

   subroutine checkDielectricComponentTags (component, prev_components, n_prev, coord_type, numertag, tagtype, precounting, error_msg)
      type(Dielectric_t), intent(in) :: component        ! Current component
      type(Dielectric_t), intent(in) :: prev_components(:) ! Array of previous components
      integer, intent(in) :: n_prev                            ! Number of previous components
      character(len=*), intent(in) :: coord_type               ! 'c1P' or 'c2P'
      integer, intent(inout) :: numertag
      type(tagtype_t), intent(inout) :: tagtype
      integer, intent(in) :: precounting
      character(len=*), intent(in) :: error_msg

      logical :: foundDuplicate
      character(len=BUFSIZE) :: tagToCheck
      integer :: i, j, k, m, tama2, prev_size

      if (coord_type == 'c1P') then
         tama2 = component%n_c1P
      else
         tama2 = component%n_c2P
      end if

      check_tags: do j = 1, tama2
         numertag = numertag + 1
         foundDuplicate = .false.

         ! Get tag to check based on coord_type
         if (coord_type == 'c1P') then
            tagToCheck = trim(adjustl(component%c1P(j)%tag))
         else
            tagToCheck = trim(adjustl(component%c2P(j)%tag))
         end if

         ! Check current component up to j-1
         if (j > 1) then
            check_current: do k = 1, j-1
               if (coord_type == 'c1P') then
                  if (tagToCheck == trim(adjustl(component%c1P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  end if
               else
                  if (tagToCheck == trim(adjustl(component%c2P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  end if
               end if
            end do check_current
         end if

         ! If not found, check all previous components
         if (.not. foundDuplicate) then
            check_previous: do m = 1, n_prev
               ! Check c1P of previous component
               if (prev_components(m)%n_c1P > 0) then
                  do k = 1, prev_components(m)%n_c1P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c1P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     end if
                  end do
               end if

               ! Check c2P of previous component
               if (prev_components(m)%n_c2P > 0) then
                  do k = 1, prev_components(m)%n_c2P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c2P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     end if
                  end do
               end if
            end do check_previous
         end if

         if (foundDuplicate) then
            numertag = numertag - 1
         else if (precounting == 1) then
            tagtype%tag(numertag) = tagToCheck
         end if
      end do check_tags
   end subroutine checkDielectricComponentTags


   subroutine checkAnimatedComponentTags (component, prev_components, n_prev, coord_type, numertag, tagtype, precounting, error_msg)
      type(ANISOTROPICbody_t), intent(in) :: component        ! Current component
      type(ANISOTROPICbody_t), intent(in) :: prev_components(:) ! Array of previous components
      integer, intent(in) :: n_prev                            ! Number of previous components
      character(len=*), intent(in) :: coord_type               ! 'c1P' or 'c2P'
      integer, intent(inout) :: numertag
      type(tagtype_t), intent(inout) :: tagtype
      integer, intent(in) :: precounting
      character(len=*), intent(in) :: error_msg

      logical :: foundDuplicate
      character(len=BUFSIZE) :: tagToCheck
      integer :: i, j, k, m, tama2, prev_size

      if (coord_type == 'c1P') then
         tama2 = component%n_c1P
      else
         tama2 = component%n_c2P
      end if

      check_tags: do j = 1, tama2
         numertag = numertag + 1
         foundDuplicate = .false.

         ! Get tag to check based on coord_type
         if (coord_type == 'c1P') then
            tagToCheck = trim(adjustl(component%c1P(j)%tag))
         else
            tagToCheck = trim(adjustl(component%c2P(j)%tag))
         end if

         ! Check current component up to j-1
         if (j > 1) then
            check_current: do k = 1, j-1
               if (coord_type == 'c1P') then
                  if (tagToCheck == trim(adjustl(component%c1P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  end if
               else
                  if (tagToCheck == trim(adjustl(component%c2P(k)%tag))) then
                     foundDuplicate = .true.
                     exit check_current
                  end if
               end if
            end do check_current
         end if

         ! If not found, check all previous components
         if (.not. foundDuplicate) then
            check_previous: do m = 1, n_prev
               ! Check c1P of previous component
               if (prev_components(m)%n_c1P > 0) then
                  do k = 1, prev_components(m)%n_c1P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c1P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     end if
                  end do
               end if

               ! Check c2P of previous component
               if (prev_components(m)%n_c2P > 0) then
                  do k = 1, prev_components(m)%n_c2P
                     if (tagToCheck == trim(adjustl(prev_components(m)%c2P(k)%tag))) then
                        print *, error_msg
                        print *, 'Duplicate tag found:', tagToCheck
                        stop
                     end if
                  end do
               end if
            end do check_previous
         end if

         if (foundDuplicate) then
            numertag = numertag - 1
         else if (precounting == 1) then
            tagtype%tag(numertag) = tagToCheck
         end if
      end do check_tags
   end subroutine checkAnimatedComponentTags


   subroutine checkLossyTags(component, prev_components, n_prev, numertag, tagtype, precounting)
      type(LossyThinSurface_t), intent(in) :: component        ! Current component
      type(LossyThinSurface_t), intent(in) :: prev_components(:) ! Array of previous components
      integer, intent(in) :: n_prev                         ! Number of previous components
      integer, intent(inout) :: numertag
      type(tagtype_t), intent(inout) :: tagtype
      integer, intent(in) :: precounting

      logical :: foundDuplicate
      character(len=BUFSIZE) :: tagToCheck
      integer :: i, j, k, m, tama2

      tama2 = component%nc
      if (tama2 == 0) then
         print *, 'Bug in LossyThinSurf Tags. Missing coordinates'
         stop
      end if

      check_tags: do j = 1, tama2
         numertag = numertag + 1
         foundDuplicate = .false.

         tagToCheck = trim(adjustl(component%C(j)%tag))

         ! Check current component up to j-1
         if (j > 1) then
            check_current: do k = 1, j-1
               if (tagToCheck == trim(adjustl(component%C(k)%tag))) then
                  foundDuplicate = .true.
                  exit check_current
               end if
            end do check_current
         end if

         ! If not found, check all previous components
         
         if ((.not. foundDuplicate) .and. (n_prev>0)) then
            check_previous: do m = 1, n_prev
               if (prev_components(m)%nc > 0) then
                  do k = 1, prev_components(m)%nc
                     if (tagToCheck == trim(adjustl(prev_components(m)%C(k)%tag))) then
                        foundDuplicate = .true.
                     end if
                  end do
               end if
            end do check_previous
         end if

         if (foundDuplicate) then
            numertag = numertag - 1
         else if (precounting == 1) then
            tagtype%tag(numertag) = tagToCheck
         end if
      end do check_tags
   end subroutine checkLossyTags


   function searchtag(tagtype,tag) result(numertag)


      character(len=BUFSIZE) :: tag
      integer(Kind=4) :: i,numertag
      type(tagtype_t) :: tagtype

      numertag=-1
      busca: do i=1,tagtype%numertags
         if (trim(adjustl(tagtype%tag(i)))==trim(adjustl(tag))) then
            numertag=i
            exit busca
         end if
      end do busca
      return
   end function searchtag


end module preprocess_tags_m
