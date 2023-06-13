      module m_read_obs

      use ice_kinds_mod

      implicit none

      character (char_len_long), public :: &
              obs_data_dir

!-------------------------------------------------------------

      contains

!-------------------------------------------------------------

      subroutine read_obs(fyear, month, mday)

      use mod_measurement
      use m_obs
      use m_random

      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag, nu_forcing
      use ice_read_write, only: ice_open, ice_read
      use ice_forcing, only: hi_data_type

      implicit none

      integer :: m, l, iobs, nrobs, k, jj, reclf
      integer :: fyear, month, mday, tmpyear
      integer :: sit_ininrobs, sit_ininrobs_cs2, sit_ininrobs_smos, tmpininrobs
      integer :: ntemp_sit, ntemp_sit_cs2, ntemp_sit_smos, sit_nrobs
      character (char_len_long) :: obs_file_sit, obs_file_sit2
      character*4 :: ayr
      character*3 :: sit
      character*2 :: amo, ady
      type(inimeasurement4) :: sit_obs
      logical :: lexist

      ininrobs = 0; sit_ininrobs = 0
      sit_ininrobs_cs2 = 0; sit_ininrobs_smos = 0
      ntemp_sit = 0; ntemp_sit_cs2 = 0; ntemp_sit_smos = 0
      sit_nrobs = 0

      write(ayr, '(i4.4)') fyear
      write(amo, '(i2.2)') month
      write(ady, '(i2.2)') mday

      if (trim(hi_data_type) == 'enoi') then
      obs_file_sit = trim(obs_data_dir)//'sit.'//ayr//amo//ady//'.dat'
      obs_file_sit2 = trim(obs_data_dir)//'sit.smos.'&
                                        //ayr//amo//ady//'.dat'

         open (nu_forcing, file=obs_file_sit, form='unformatted', &
                           access='sequential')
         inquire(nu_forcing, exist=lexist) 
         inquire(nu_forcing, size=reclf)
         if (lexist .eq. .true. .and. reclf .ne. 0) then
         do jj = 1, reclf/28
         read(nu_forcing) sit_obs
         if (sit_obs%d .ne. -9999 ) then
         sit_ininrobs_cs2 = sit_ininrobs_cs2 + 1
         endif;enddo
         endif
         close(nu_forcing)

         open (nu_forcing, file=obs_file_sit2, form='unformatted', &
                           access='sequential')
         inquire(nu_forcing, exist=lexist) 
         inquire(nu_forcing, size=reclf)
         if (lexist .eq. .true. .and. reclf .ne. 0) then
         do jj = 1, reclf/28
         read(nu_forcing) sit_obs
         if (sit_obs%d .ne. -999 .and. sit_obs%d .le. 1.0 .and. sit_obs%d .ne. 0.0) then
         sit_ininrobs_smos = sit_ininrobs_smos + 1
         endif;enddo
         endif
         close(nu_forcing)

         open (nu_forcing, file=obs_file_sit, form='unformatted', &
                           access='sequential')
         inquire(nu_forcing, exist=lexist)
         inquire(nu_forcing, size=reclf)
         if (lexist .eq. .true. .and. reclf .ne. 0) then
         allocate(sit_var_cs2(sit_ininrobs_cs2))
         do jj = 1, reclf/28
         read(nu_forcing) sit_obs
         if (sit_obs%d .ne. -9999) then
         ntemp_sit_cs2 = ntemp_sit_cs2 + 1
         sit_var_cs2(ntemp_sit_cs2) = sit_obs
         endif;enddo
         endif
         close(nu_forcing)

         open (nu_forcing, file=obs_file_sit2, form='unformatted', &
                           access='sequential')
         inquire(nu_forcing, exist=lexist)
         inquire(nu_forcing, size=reclf)
         if (lexist .eq. .true. .and. reclf .ne. 0.0) then
         allocate(sit_var_smos(sit_ininrobs_smos))
         do jj = 1, reclf/28
         read(nu_forcing) sit_obs
         if (sit_obs%d .ne. -999 .and. sit_obs%d .le. 1.0 .and. sit_obs%d .ne. 0.0) then
         ntemp_sit_smos = ntemp_sit_smos + 1
         sit_var_smos(ntemp_sit_smos) = sit_obs
         endif;enddo
         endif
         close(nu_forcing)

         sit_ininrobs = ntemp_sit_cs2 + ntemp_sit_smos
         ntemp_sit = ntemp_sit_cs2 + ntemp_sit_smos

         allocate (sit_var(sit_ininrobs))
         do jj = 1, sit_ininrobs
         if (jj .le. sit_ininrobs_cs2) then
         sit_var(jj) = sit_var_cs2(jj)
         else if (jj .gt. sit_ininrobs_cs2 .and. &
                  jj .le. sit_ininrobs_cs2+sit_ininrobs_smos) then
         sit_var(jj) = sit_var_smos(jj-sit_ininrobs_cs2)
         endif;enddo

       endif

      sit_ininrobs = ntemp_sit
      tmpininrobs = sit_ininrobs
      ininrobs = tmpininrobs

      if (allocated (iniobs)) deallocate(iniobs); allocate(iniobs(ininrobs))
      if (allocated (gE))     deallocate(gE);     allocate(gE(ininrobs))

      gE = 0.
      do m=1,ininrobs
        if (m .le. sit_ininrobs) then
        sit_nrobs = sit_nrobs + 1
        iniobs(m) = sit_var(sit_nrobs)
        endif
        call random(gE(m),1)
        gE(m) = gE(m)*sqrt(iniobs(m)%var)
      enddo

      if (trim(hi_data_type) == 'enoi') then
      if (allocated(sit_var)) deallocate(sit_var)
      if (allocated(sit_var_cs2)) deallocate(sit_var_cs2)
      if (allocated(sit_var_smos)) deallocate(sit_var_smos)
      endif

      end subroutine read_obs

      end module m_read_obs
