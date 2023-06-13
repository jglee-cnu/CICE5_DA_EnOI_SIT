module m_Generate_element_Sij 
contains
! DESCRIPTION
! Function which generates the S matix element by element. The function operates 
! for different observation types: temp, sal, u, v and SSH/SLA data. 
! The four first data types are found as progonostic parameters in the model
! While the SSH/SLA data has to be generated from the prognositc data. This
! Is done by a separte function in the case of this type of measurments.
! For the given location of the observation we combine the four surrounding 
! vertical model data using the biilinear coeffisients given for each obs(j).
! This is done for all types of observations. 
! In the  case of horizontal data we do not need to do more with the data
! While in case of vertical measurment profiles (data for depth .ne. 0) 
! the data has to be interporlated to a denser gridd in the vertical and for
! This we use a sline interpolation.

real function Generate_element_Sij(obs,mem)
   use mod_states
   use mod_measurement
   use mod_dimensions
   use ice_domain, only: nblocks

   implicit none
   type(measurement),    intent(in) :: obs
   type(gstates),        intent(in) :: mem
   integer                          :: i,j,ip1,jp1,k,m,km1
   integer                          :: ii,jj, imin,imax,jmin,jmax, cnt
   integer                          :: iblk
   real                             :: x0,y0,model_upper,model_lower
   real, dimension(nz)              :: zt


   real, parameter :: undef=-9.99     ! land points have value huge()

! Get model gridcell
   i   = obs%ipiv
   j   = obs%jpiv
   ip1 = min(i+1,nx)
   jp1 = min(j+1,ny)

! interpolate model depths
10 continue

      Generate_element_Sij = mem%var(i  ,j  )*obs%a1 &
                           + mem%var(ip1,j  )*obs%a2 &
                           + mem%var(ip1,jp1)*obs%a3 &
                           + mem%var(i  ,jp1)*obs%a4

end function Generate_element_Sij
end module m_Generate_element_Sij

