module mod_dimensions

   use ice_domain, only: nblocks

   implicit none
   integer, parameter :: nx=320                 ! x-dimension of analysis region
   integer, parameter :: cnz=1                  ! z-dimension of analysis region
   integer, parameter :: nz=1                   ! z-dimension of all model grid
   integer, parameter :: any=384                ! y-dimension of all model grid
   integer, save      :: ny
end module mod_dimensions

