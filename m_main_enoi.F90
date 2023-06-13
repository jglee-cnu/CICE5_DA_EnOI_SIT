      module m_main_enoi

      use ice_kinds_mod

      implicit none

      character (char_len_long), public :: &
         err_data_dir

      integer (kind=int_kind), public :: nrens
!=======================================================================

          contains
           
!=======================================================================

        subroutine main_enoi(iblk)

        use ice_constants, only: c0, c10, c180, field_loc_center, field_type_scalar, &
                                 rad_to_deg, puny, puny2, c1, c360, c2, p01, p001, &
                                 p2, Tocnfrz
        use ice_domain_size, only: ncat, max_blocks, nx_global, ny_global, nilyr, nslyr

        use ice_domain, only: nblocks, distrb_info
        use ice_domain_size, only: nx_global, ny_global
        use ice_blocks, only: nx_block, ny_block
        use ice_state, only: aice, aicen, vice_up, vicen_up, vice, &
                             vice_org, vicen_org, vicen, &
                             aice_org, aicen_org, aicen_up, &
                             hice_org, hicen_org, &
                             hsno_org, hsnon_org, &
                             aice_up, vsnon, vsno, uvel, vvel, &!, gaice3, ogaice3
                             trcrn, nt_Tsfc,nt_qice,nt_qsno,ainc_sit
        use ice_flux, only: strocnx, strocny, strocnxT, strocnyT, sst
        use ice_calendar, only: istep, istep1, time, time_forc, year_init, &
                                sec, mday, month, nyr, yday, daycal, dayyr, &
                                daymo,  days_per_year
        use ice_fileunits, only: nu_diag, nu_forcing
        use ice_blocks, only: nx_block, ny_block, nblocks_tot

        use ice_grid, only: TLON, TLAT
        use ice_forcing, only: hi_data_type, fyear
        use ice_gather_scatter, only: scatter_global, gather_global
        use ice_calendar, only: yday, dt,new_day,new_month,new_year
        use ice_calendar, only: year_init, month, mday, daymo 
        use ice_communicate, only: my_task, master_task
        use ice_read_write, only: ice_open, ice_read
        use ice_itd, only: hin_max
        use EnOI
        use m_read_obs
        use m_obs

       integer :: i,j,k, year, day, d, nbits, iens, cate, task
       integer (kind=int_kind), intent(in) :: iblk
       integer (kind=int_kind) :: iyear, imonth, iday
       real (kind=dbl_kind) :: yday1
       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,nrens) :: &
               ogaice, gaice, ogvice, gvice, err33_sit
       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
               err33_sitm
       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
               err
       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,nrens) :: &
               gathv
       real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,nrens) :: &
                ogathv
       real (kind=dbl_kind), dimension(nx_global,ny_global) :: &
               err1_sit, err11_sit
       real (kind=dbl_kind), dimension(nx_global,ny_global,nrens) :: &
               err3_sit, err2_sit, err22_sit, ggathv
       real (kind=dbl_kind), dimension(nx_global,ny_global) :: &
               err3_sitm
       real :: work1(nx_block,ny_block,max_blocks), &
               work2(nx_block,ny_block,max_blocks)
       real, dimension(:,:), allocatable :: &
             work_g1, work_g2, TLON_g, TLAT_g
       logical owrite
       character(len=3) :: num1, num2
       character (char_len_long) :: ens_sit_file1, ens_sit_file2

         
       character*4 :: ayr
       character*2 :: amo, ady
       logical     :: exist
       real :: dmiss, dmiss2, we
       data dmiss/-9.99e+08/
       data dmiss2/1e+30/

      nbits = 64

      if (trim(hi_data_type) == 'enoi') then
      if (nx_global == 320) then ! gx1

      iyear = nyr + year_init - 1 ! set year_init=1 in ice_in to get iyear=nyr
      imonth = month
      iday = mday
      yday1 = yday

      if (new_year) then
       iyear = iyear - 1
       imonth = 12
       iday = daymo(imonth)
      elseif (new_month) then
       imonth = month - 1
       iday = daymo(imonth)
      elseif (new_day) then
       iday = iday - 1
       yday1 = yday1 - 2
      endif
      if (yday1 .eq. 0.0) yday1=1.0

      write(ayr, '(i4.4)') iyear
      write(amo, '(i2.2)') imonth
      write(ady, '(i2.2)') iday

      d = int(yday1)
      write(num1,'(i3.3)') d
      ens_sit_file1 = trim(err_data_dir)//'ensemble.sit.'//num1//'.dat'
      endif

      if (my_task == master_task) then
      write (nu_diag,*) ' '
      write (nu_diag,*) 'err sit file1: ', trim(ens_sit_file1)
      endif

      if (my_task == master_task) then
      open (nu_forcing,file=ens_sit_file1,recl=nx_global*ny_global*nbits/8, &
                     form='unformatted',access='direct')
      endif
      if (my_task == master_task) then

      do iens = 1, nrens
      read (nu_forcing,rec=iens) err1_sit

      do j=1,ny_global; do i=1,nx_global
      err2_sit(i,j,iens) = 0.0
      err2_sit(i,j,iens) = err1_sit(i,j)
      enddo;enddo
      enddo ! iens
      endif ! my_task
      if (my_task == master_task) close(nu_forcing)

      if (my_task == master_task) then
        err3_sit=0.0; err3_sitm=0.0
        err33_sit=0.0; err33_sitm=0.0
        do iens = 1, nrens
        err3_sitm(:,:) = err3_sitm(:,:) + err2_sit(:,:,iens)/nrens
        enddo
        do iens = 1, nrens;do j = 1, ny_global;do i = 1, nx_global
        err3_sit(i,j,iens) = err2_sit(i,j,iens) - err3_sitm(i,j)
        enddo;enddo;enddo
      endif ! my_task

      do iens = 1, nrens
      call scatter_global(err33_sit(:,:,:,iens), err3_sit(:,:,iens), &
                          master_task, distrb_info, &
                          field_loc_center, field_type_scalar)
      enddo

      do iens = 1, nrens
      err33_sitm(:,:,:) = err33_sitm(:,:,:) + err33_sit(:,:,:,iens)/nrens
      enddo

      do iens = 1, nrens
      err33_sit(:,:,:,iens) = err33_sit(:,:,:,iens) - err33_sitm(:,:,:)
      enddo

      do iens = 1, nrens; do j = 1, ny_block;do i = 1, nx_block
      gathv(i,j,iblk,iens) = vice(i,j,iblk)+err33_sit(i,j,iblk,iens)
      enddo;enddo;enddo
       
      do task = 0, nblocks_tot-1; do iens = 1, nrens   
      call gather_global(ggathv(:,:,iens), gathv(:,:,:,iens), &
                         task, distrb_info, spc_val=c0)
      enddo;enddo

      do j=1,ny_block;do i=1,nx_block
      work1(i,j,iblk) = TLON(i,j,iblk)*rad_to_deg + c360
      if (work1(i,j,iblk) .gt. c360) work1(i,j,iblk) = work1(i,j,iblk)-c360
      if (work1(i,j,iblk) .lt. c0)   work1(i,j,iblk) = work1(i,j,iblk)+c360
      work2(i,j,iblk) = TLAT(i,j,iblk)*rad_to_deg
      enddo;enddo

      allocate(work_g1(nx_global,ny_global))
      allocate(work_g2(nx_global,ny_global))
      allocate(TLON_g(nx_global,ny_global))
      allocate(TLAT_g(nx_global,ny_global))

      do task = 0, nblocks_tot-1
      call gather_global(work_g1(:,:), work1(:,:,:), &
                         task, distrb_info, spc_val=c0)
      call gather_global(work_g2(:,:), work2(:,:,:), &
                         task, distrb_info, spc_val=c0)
      call gather_global(TLON_g(:,:), TLON(:,:,:), &
                         task, distrb_info, spc_val=c0)
      call gather_global(TLAT_g(:,:), TLAT(:,:,:), &
                         task, distrb_info, spc_val=c0)
      enddo
   
      do j=1,ny_block;do i=1,nx_block
      vice_org(i,j,iblk) = vice(i,j,iblk)
      aice_org(i,j,iblk) = aice(i,j,iblk)
      hice_org(i,j,iblk) = vice(i,j,iblk) / aice(i,j,iblk)
      hsno_org(i,j,iblk) = vsno(i,j,iblk) / aice(i,j,iblk)
      do cate=1,5
      vicen_org(i,j,cate,iblk) = vicen(i,j,cate,iblk)
      aicen_org(i,j,cate,iblk) = aicen(i,j,cate,iblk)
      hicen_org(i,j,cate,iblk) = vicen(i,j,cate,iblk) / aicen(i,j,cate,iblk)
      hsnon_org(i,j,cate,iblk) = vsnon(i,j,cate,iblk) / aicen(i,j,cate,iblk)
      enddo;enddo;enddo
 
      call read_obs(iyear,imonth,iday)
      
      call core_enoi(ogathv(:,:,iblk,:),gathv(:,:,iblk,:)  ,&
                ggathv(:,:,:), nrens                 ,&
                work1(:,:,iblk),work2(:,:,iblk)       ,&
                work_g1(:,:), work_g2(:,:)            ,&
                TLON(:,:,iblk), TLAT(:,:,iblk)        ,&
                TLON_g(:,:),TLAT_g(:,:),'sit')

      ogathv(:,:,iblk,:) = ogathv(:,:,iblk,:)

      deallocate(iniobs)
      deallocate(gE)
      deallocate(work_g1)
      deallocate(work_g2)
      deallocate(TLON_g)
      deallocate(TLAT_g)
  
      aicen_up = aicen_org
      aice_up = aice_org
      vice_up = 0.
      do j=1,ny_block
      do i=1,nx_block
      do iens = 1,nrens
      vice_up(i,j,iblk) = vice_up(i,j,iblk) + ogathv(i,j,iblk,iens)/nrens
      enddo

      vice_up(i,j,iblk) = max(vice_up(i,j,iblk),c0) !6th

! for category

      do cate = 1, 5
      if (aice_org(i,j,iblk) .ne. c0 .and. vice_org(i,j,iblk) .ne. c0) then
         aicen_up(i,j,cate,iblk) = aicen_org(i,j,cate,iblk)
         vicen_up(i,j,cate,iblk) = vice_up(i,j,iblk)* &
                                   (vicen_org(i,j,cate,iblk)/vice_org(i,j,iblk))
      else if (aice_org(i,j,iblk) .eq. c0 .and. vice_org(i,j,iblk) .eq. c0) then
        if (aice_up(i,j,iblk) .eq. c0 .and. vice_up(i,j,iblk) .ge. puny2) then
          if (cate.eq.1) then
          vicen_up(i,j,cate,iblk) = vice_up(i,j,iblk) * 0.6
          aicen_up(i,j,cate,iblk) = vicen_up(i,j,cate,iblk) &
                                  / (hin_max(cate) - hin_max(cate-1))/2
          else if (cate.eq.2) then
          vicen_up(i,j,cate,iblk) = vice_up(i,j,iblk) * 0.2
          aicen_up(i,j,cate,iblk) = vicen_up(i,j,cate,iblk) &
                                  / (hin_max(cate) - hin_max(cate-1))/2
          else if (cate.eq.3) then
          vicen_up(i,j,cate,iblk) = vice_up(i,j,iblk) * 0.1
          aicen_up(i,j,cate,iblk) = vicen_up(i,j,cate,iblk) &
                                  / (hin_max(cate) - hin_max(cate-1))/2
          else if (cate.ge.4) then
          vicen_up(i,j,cate,iblk) = vice_up(i,j,iblk) * 0.05
          aicen_up(i,j,cate,iblk) = vicen_up(i,j,cate,iblk) &
                                  / (hin_max(cate) - hin_max(cate-1))/2
          endif
        endif
      endif
      enddo

      aice_up(i,j,iblk) = sum(aicen_up(i,j,:,iblk), dim=1)

      do cate = 1, 5

      if (aice_up(i,j,iblk) .gt. c1) then
      aicen_up(i,j,cate,iblk) = aicen_up(i,j,cate,iblk) &
                               -((aice_up(i,j,iblk)-c1) &
                               *(aicen_up(i,j,cate,iblk)/aice_up(i,j,iblk)))
      endif

        aicen_up(i,j,cate,iblk) = min(max(aicen_up(i,j,cate,iblk),c0),c1)

      if (trim(hi_data_type) ==  'enoi') then
      if (vicen_org(i,j,cate,iblk) .ge. puny2) then
          aicen_up(i,j,cate,iblk) = aicen_up(i,j,cate,iblk) &
                                   * (vicen_up(i,j,cate,iblk)/vicen_org(i,j,cate,iblk))
      else if (vicen_org(i,j,cate,iblk) .eq. c0 .and. vicen_up(i,j,cate,iblk) .eq. c0) then
        aicen_up(i,j,cate,iblk) = c0 ! 2nd
        vicen_up(i,j,cate,iblk) = c0
      endif
      endif
      enddo

      aice_up(i,j,iblk) = sum(aicen_up(i,j,:,iblk), dim=1)

      do cate = 1, 5

      if (aice_up(i,j,iblk) .gt. c1) then
      aicen_up(i,j,cate,iblk) = aicen_up(i,j,cate,iblk) &
                               -((aice_up(i,j,iblk)-c1) &
                               *(aicen_up(i,j,cate,iblk)/aice_up(i,j,iblk)))
      endif

      vicen_up(i,j,cate,iblk) = max(vicen_up(i,j,cate,iblk),c0) !5th

      vicen_up(i,j,cate,iblk) = &
      max(vicen_up(i,j,cate,iblk), hin_max(cate-1)*aicen_up(i,j,cate,iblk))
      vicen_up(i,j,cate,iblk) = &
      min(vicen_up(i,j,cate,iblk), 9.33*aicen_up(i,j,cate,iblk))

      if (aicen_org(i,j,cate,iblk) .ge. puny2) then
      vsnon(i,j,cate,iblk) = vsnon(i,j,cate,iblk) &
                            * min(aicen_up(i,j,cate,iblk)/aicen_org(i,j,cate,iblk), c10)
      endif

      trcrn(i,j,nt_Tsfc,cate,iblk) = min(trcrn(i,j,nt_Tsfc,cate,iblk),c0)
      enddo

      do cate = 1, 5
        if (aicen_org(i,j,cate,iblk) .lt. puny2) then
          vsnon(i,j,cate,iblk) = p2 * vicen_up(i,j,cate,iblk)
        endif
        if (aicen_up(i,j,cate,iblk) .eq. c0) then
          trcrn(i,j,nt_Tsfc,cate,iblk) = Tocnfrz
          vicen_up(i,j,cate,iblk) = c0 !3rd
          vsnon(i,j,cate,iblk) = c0
        endif
      enddo

      aice_up(i,j,iblk) = sum(aicen_up(i,j,:,iblk), dim=1)

      if (aice_up(i,j,iblk) .eq. c0) then
          uvel(i,j,iblk) = c0  !4th
          vvel(i,j,iblk) = c0
          strocnx(i,j,iblk) = c0
          strocny(i,j,iblk) = c0
          strocnxT(i,j,iblk) = c0      
          strocnyT(i,j,iblk) = c0
      endif

      vice_up(i,j,iblk) = sum(vicen_up(i,j,:,iblk),dim=1)

      enddo;enddo

      ainc_sit = 0.0
      ainc_sit=vice_up-vice_org
      aice = aice_up
      aicen = aicen_up
      vice = vice_up
      vicen = vicen_up

      endif ! if (trim(hi_data_type) == 'enoi') 

      return
      end subroutine main_enoi

      end module m_main_enoi
