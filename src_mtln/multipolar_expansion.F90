module multipolar_expansion_mod

   use mtln_types_mod

contains
   function multipolarExpansion(position, ab, expansionCenter) result(res)
    
    ! 2D multipolar expansion from:
      ! TSOGTGEREL GANTUMUR, MULTIPOLE EXPANSIONS IN THE PLANE.
      ! 2016-04-16 lecture notes.
      
    implicit none

      real, intent(in) :: position(2)
      real, intent(in) :: expansionCenter(2)
      type_t(multipolar_coefficient_t), dimension(:), intent(in) :: ab
      real :: res

      real :: rVec(2), r, phi
      integer :: n

      rVec = position - expansionCenter
      r = sqrt(rVec(1)**2 + rVec(2)**2)
      phi = atan2(rVec(2), rVec(1))

      res = 0.0
      do n = 1, size(ab)
         if (n == 1) then
            res = res - ab(n)%a * log(r)
            ! b0 should always be zero
         else
            res = res + (ab(n)%a * cos((n-1)*phi) + ab(n)%b * sin((n-1)*phi)) / r**(n-1)
         end if
      end do

      res = res / (2.0 * pi())
   end function multipolarExpansion

end module multipolar_expansion_mod
