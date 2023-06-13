      module EnOI

!=======================================================================

          contains
           
!=======================================================================

          subroutine core_enoi(ogath,gath,ggath,nrens,work1,work2,work_g1,work_g2,TLON,TLAT,TLON_g,TLAT_g, cvar)

          use ice_kinds_mod
          use ice_blocks, only: nx_block, ny_block
          use ice_domain_size, only: nx_global, ny_global
          use ice_constants, only: c0, c12, c1000, secday, depressT, &
              field_loc_center, field_type_scalar
          use ice_domain, only: nblocks, distrb_info
          use ice_state, only: vice
          use ice_calendar, only: year_init, month, mday
          use ice_forcing, only: fyear
          use ice_gather_scatter, only: scatter_global, gather_global

          use mod_dimensions
          use mod_states
          use mod_measurement
          use m_active_obs
          use m_ensmean
          use m_ensvar
          use m_prep_4_EnOI
          use m_obs
          use m_local_analysis
          use common_mpi

          character(len=*), intent(in) :: cvar
          integer :: nrens, nrens2, midy, nrobs, task
          integer :: iglobal_ndim
          integer :: m, i, j, ii, jj, o
          type(gstates) :: gA(nrens)
          type(states)  :: A((nx_block)*ny_block,nrens)
          real,dimension(nx_block,ny_block,nrens) :: gath, ogath
          real,dimension(nx_block,ny_block),intent(in) :: work1, work2, TLON, TLAT
          real,dimension(nx_global,ny_global),intent(in) :: work_g1,work_g2,TLON_g,TLAT_g
          real,dimension(nx_global,ny_global,nrens),intent(in) :: ggath
          real,dimension(nx) :: xt
          real,dimension(any) :: yt
          real,dimension(nz) :: zt
          real, allocatable, dimension(:,:) :: S,D
          real, allocatable, dimension(:)   :: meanD
          real :: lon,lat

          real, allocatable, dimension(:) ::  ave_sit
          real, allocatable, dimension(:,:) :: A4_sit, ST, HBHT
          real :: radius
          real, dimension(nx_block, ny_block, max_blocks, nrens, nrens) :: X4_blk
          real, dimension(nx_global, ny_global, nrens, nrens) :: X4
          integer,dimension(nx_block,ny_block,max_blocks) :: nobs_array_blk
          integer,dimension(nx_global,ny_global) :: nobs_array

          ! model variables part
          do m=1,nrens;do j=1,ny_global;do i=1,nx_global
          gA(m)%var(i,j) = ggath(i,j,m)
          enddo;enddo;enddo

          do m=1,nrens
          do j=1,ny_block
          do i=1,nx_block
          A(i+(j-1)*(nx_block),m)%var = gath(i,j,m)
          enddo;enddo;enddo

          ! activate observations..
          call active_obs(nrobs,gA,nrens,work1,work2,work_g1,work_g2,TLON,TLAT,TLON_g,TLAT_g, cvar)

          if(nrobs.eq.0) then
          go to 100
          endif

          ! D,S

          allocate(D(nrobs,nrens))
          allocate(meanD(nrobs))
          allocate(S(nrobs,nrens))
          call prep_4_EnOI(gA,nrens,nrobs,D,S)

! Each Proccessors

          radius = 800
          call pre_local_analysis(nrens,nrobs,obs,work1,work2, &
               D, S, radius, nobs_array_blk, nx_block, ny_block, X4_blk)

          iglobal_ndim = (nx_block)*ny_block

          ! mean innovation 
          meanD=0.
          do jj=1,nrens
          do ii=1,nrobs
                meanD(ii)=meanD(ii)+D(ii,jj)
          enddo
          enddo
          meanD(:)=meanD(:)/float(nrens)

          call local_analysis(nobs_array_blk(:,:,1), A%var, nrens, nx_block, ny_block, X4_blk(:,:,1,:,:), cvar)

          deallocate(D, meanD, S)

          100 continue
 
          ! OUTPUT = A(xa)

          ogath = 0
          do m=1,nrens;do j=1,ny_block;do i=1,nx_block
          ogath(i,j,m) = ogath(i,j,m) + A(i+(j-1)*(nx_block),m)%var
          enddo;enddo;enddo

          if(allocated(obs)) deallocate(obs)
           
        call barrier()

        return
        end subroutine core_enoi

      end module EnOI
