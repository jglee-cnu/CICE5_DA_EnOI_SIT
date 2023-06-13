module m_modstate_point
contains

      real function modstate_point(cvar,mem,ix,iy)
      use mod_dimensions
      use mod_states
      implicit none

      character(len=*), intent(in) :: cvar
      type(gstates)   , intent(in) :: mem
      integer         , intent(in) :: ix,iy

      select case (trim(cvar))

      case ('var')
         modstate_point=mem%var(ix,iy)

      case default
         print *,'No match in modstate_point'
         print *,cvar
         stop
      end select
   end function modstate_point

end module m_modstate_point
