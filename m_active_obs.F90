module m_active_obs

  use ice_kinds_mod
  character (char_len_long), public :: comm_dir

contains
subroutine active_obs(nrobs,gA,nrens,work1,work2,work_g1,work_g2,TLON,TLAT,TLON_g,TLAT_g,cvar)

  use ice_blocks, only: nx_block, ny_block
  use ice_domain_size, only: nx_global, ny_global
  use ice_constants, only: c0, c12, c1000, secday, rad_to_deg, c360, c180, &
                           field_loc_center, field_type_scalar, c2
  use ice_communicate, only: my_task, master_task
  use ice_domain, only: nblocks, distrb_info
  use ice_domain_size, only: max_blocks
  use ice_fileunits, only: nu_diag, nu_forcing
  use ice_state, only: aice, vice
  use ice_gather_scatter, only: scatter_global, gather_global

  use mod_dimensions
  use mod_measurement
  use mod_states
  use m_obsgrid
  use m_obs

  implicit none
  character(len=*), intent(in) :: cvar
  type(gstates)                 :: gA(nrens)
  real,dimension(any) :: yt
  real,dimension(nx) :: xt
  real,dimension(nz) :: zt

  real tmpdist
  real olat,olon, lat, lon
  integer       ::  m,joff,l,iobs,nrobs,nst
  integer       ::  i,j,nrens,kk,k,jd,ntemp
  integer       ::  istr,iend,ims,ime,lengthx,lengthy
  integer       ::  iblk, nstation, ios
  integer, dimension(2) :: a1, a2, a3, a4, a11, a22, a33, a44, a111, a222, a333, a444
  real :: b, b1, b2, b3, b4
  real,dimension(ininrobs) :: TTLON, TTLAT
  character(100) ::  cdummy
  character (char_len_long) :: station_file
  real, dimension(nx_block,ny_block), intent(in) :: &          
          work1, work2, TLON, TLAT
  real, dimension(nx_global,ny_global), intent(in) :: &
          work_g1, work_g2, TLON_g, TLAT_g
  integer :: op

  ntemp = 0

if (allocated(obs)) deallocate(obs); allocate(obs(ininrobs))

  do m = 1,ininrobs

   if (iniobs(m)%lon .lt. 0) then
   iniobs(m)%lon = iniobs(m)%lon + 360.
   endif

    if(iniobs(m)%id .eq. trim(cvar)) then
    ntemp = ntemp + 1

       obs(ntemp)%d = iniobs(m)%d
       obs(ntemp)%var = iniobs(m)%var
       obs(ntemp)%id = iniobs(m)%id

       obs(ntemp)%lon = iniobs(m)%lon
       obs(ntemp)%lat = iniobs(m)%lat
       TTLON(ntemp) = iniobs(m)%lon
       TTLAT(ntemp) = iniobs(m)%lat

         a1 = minloc(sqrt((work_g1-TTLON(ntemp))**2+(work_g2-TTLAT(ntemp))**2) &
              ,TTLON(ntemp)-work_g1.ge.0.and.TTLAT(ntemp)-work_g2.ge.0)
         a2 = minloc(sqrt((work_g1-TTLON(ntemp))**2+(work_g2-TTLAT(ntemp))**2) &
              ,TTLON(ntemp)-work_g1.le.0.and.TTLAT(ntemp)-work_g2.ge.0)
         a3 = minloc(sqrt((work_g1-TTLON(ntemp))**2+(work_g2-TTLAT(ntemp))**2) &
              ,TTLON(ntemp)-work_g1.lt.0.and.TTLAT(ntemp)-work_g2.lt.0)
         a4 = minloc(sqrt((work_g1-TTLON(ntemp))**2+(work_g2-TTLAT(ntemp))**2) &
              ,TTLON(ntemp)-work_g1.ge.0.and.TTLAT(ntemp)-work_g2.lt.0)

         b1 = sqrt((work_g1(a1(1),a1(2))-TTLON(ntemp))**2+(work_g2(a1(1),a1(2))-TTLAT(ntemp))**2)
         b2 = sqrt((work_g1(a2(1),a2(2))-TTLON(ntemp))**2+(work_g2(a2(1),a2(2))-TTLAT(ntemp))**2)
         b3 = sqrt((work_g1(a3(1),a3(2))-TTLON(ntemp))**2+(work_g2(a3(1),a3(2))-TTLAT(ntemp))**2)
         b4 = sqrt((work_g1(a4(1),a4(2))-TTLON(ntemp))**2+(work_g2(a4(1),a4(2))-TTLAT(ntemp))**2)

         b = min(b1, b2, b3, b4)

         if (b.eq.b1) then
            obs(ntemp)%ipiv = a1(1) !--> TLON_g grid location
            obs(ntemp)%jpiv = a1(2) !--> TLAT_g grid location
         else if (b.eq.b2) then
            obs(ntemp)%ipiv = a2(1)
            obs(ntemp)%jpiv = a2(2)
         else if (b.eq.b3) then
            obs(ntemp)%ipiv = a3(1)
            obs(ntemp)%jpiv = a3(2)
         else if (b.eq.b4) then
            obs(ntemp)%ipiv = a4(1)
            obs(ntemp)%jpiv = a4(2)
         endif

       else
       cycle
       endif ! if iniobs(m)%id

       obs(ntemp)%a1 = 1
       obs(ntemp)%a2 = 0
       obs(ntemp)%a3 = 0
       obs(ntemp)%a4 = 0

  enddo ! do m = 1,ininrobs

  nrobs = ntemp

end subroutine active_obs
end module m_active_obs
