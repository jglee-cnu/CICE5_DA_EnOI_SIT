module m_local_analysis

  use ice_kinds_mod
  use ice_blocks, only: nblocks_tot

  implicit none


  public local_analysis
  public pre_local_analysis

  integer, parameter, private :: LOCFUN_NONE = 1
  integer, parameter, private :: LOCFUN_STEP = 2
  integer, parameter, private :: LOCFUN_GASPARI_COHN = 3

  integer, parameter, public :: LOCFUN_USED = LOCFUN_GASPARI_COHN

contains 

  subroutine local_analysis(nobs_array, fld, nrens, ni, nj, X4, cvar)
    use ice_communicate, only: my_task, master_task
    use mod_measurement
    use mod_states
    implicit none
    character(len=*), intent(in) :: cvar
    real, dimension(ni,nj,nrens,nrens), intent(in) :: X4
    integer, intent(in) :: ni, nj ! size of grid
    integer, intent(in) :: nrens ! size of ensemble
    real, dimension(ni * nj, nrens), intent(inout) :: fld ! field - to be updated
    integer, dimension(ni, nj), intent(in) :: nobs_array! Local no obs

    real, dimension(nrens) :: sub_fld
    real, dimension(nrens, nrens) :: X4tmp
    integer :: m, i, j, v, op
    integer :: irecl

    do j = 1, nj
       do i = 1, ni
          if (nobs_array(i, j) == 0) then
             cycle
          end if

          X4tmp = X4(i,j,:,:)

          sub_fld = fld((j - 1) * ni + i, :)
          ! Final matrix multiplication - this is the alternative to simplified
          ! analysis2()
          fld((j - 1) * ni + i, :) = 0.0
          do m = 1, nrens
             fld((j - 1) * ni + i, :) = fld((j - 1) * ni + i, :) + sub_fld(m) * X4tmp(m, :)
          end do

       enddo
    enddo
  end subroutine local_analysis

  subroutine pre_local_analysis(nrens, nrobs, obs, modlon, modlat, &
                                d, S, radius, nobs_array, ni, nj, X4)

    use common_mpi
    use mod_measurement
    use mod_states
    use m_spherdist
    use m_random
    use ice_communicate, only: my_task, master_task
    use m_obs, only: tmp_DD, tmp_S, tmp_obs

    integer, intent(in) :: nrens
    integer, intent(in) :: nrobs
    type(measurement), intent(in) :: obs(nrobs)
    real, dimension(ni, nj), intent(in) :: modlon, modlat
    real, dimension(nrobs), intent(in) :: d ! innovations
    real, dimension(nrobs, nrens), intent(inout) :: S ! HA
    real, intent(in) :: radius ! localisation radius in km
    integer, dimension(ni, nj), intent(out) :: nobs_array ! # of local obs for&
    real, dimension(ni,nj,nrens,nrens), intent(out) :: X4
    integer, intent(in) :: ni, nj ! horizontal grid size

    real, dimension(nrens, nrens) :: X4tmp
    real, dimension(nrobs, nrens) :: DD1, DD     ! ensemble innovations
    real, dimension(nrobs, nrens) :: DD_send, DD_recv, S_send, S_recv
    real, dimension(nrobs*nrens)  :: DD_send2, S_send2, DD_recv2, S_recv2
    real, allocatable, dimension(:,:) :: X1 ! nobs x nobs
    real, allocatable, dimension(:,:) :: subD, subS ! nobs x nrens
    real, allocatable, dimension(:,:) :: G
    real, allocatable, dimension(:) :: x
    real :: sqrtm

    integer :: iostatus

    integer :: lapack_info

    integer :: nobs ! # of local obs
    integer :: m, i, j, o, jj
    integer :: irecl
    integer :: nobs_max ! maximal number of local obs
    real :: dist, lfactor
    type(measurement) :: obs0
    integer :: task

    if (my_task == master_task) then
       do o = 1, nrobs
          call randn(nrens, DD1(o, :))
          DD1(o, :) = DD1(o, :) * sqrt(obs(o) % var)
       end do
       do m = 1, nrens
          DD1(:, m) = d + DD1(:, m) - S(:, m)
       end do
   
    sqrtm = sqrt(real(nrens) - 1.0d0)
    do o = 1, nrobs
       DD_send(o, :) = DD1(o, :) / (sqrt(obs(o) % var) * sqrtm)
        S_send(o, :) = S(o, :) / (sqrt(obs(o) % var) * sqrtm)
    end do

    do o = 1, nrobs; do m = 1, nrens
       DD_send2(nrens*(o-1)+m) = DD_send(o, m)
        S_send2(nrens*(o-1)+m) =  S_send(o, m)
    enddo;enddo

       do task = 1, nblocks_tot-1
       call common_mpi_send_1D(DD_send2, task, 1)
       call common_mpi_send_1D( S_send2, task, 2)
       enddo

    else
    
       call common_mpi_recv_1D(DD_recv2, master_task, 1)
       call common_mpi_recv_1D( S_recv2, master_task, 2)
    endif
 
    call barrier()

    do o = 1, nrobs; do m = 1, nrens
    DD_recv(o,m) = DD_recv2(nrens*(o-1)+m)
     S_recv(o,m) =  S_recv2(nrens*(o-1)+m)
    enddo;enddo

    DD=0; S=0
    do o = 1, nrobs; do m = 1, nrens
    if (my_task == master_task ) then
       DD(o,m) = DD_send(o,m)
        S(o,m) =  S_send(o,m)
    else
       DD(o,m) = DD_recv(o,m)
        S(o,m) =  S_recv(o,m)
    endif
    enddo;enddo

    nobs_array = 0
    do j = 1, nj
       do i = 1, ni
       nobs_array(i,j) = 0
       X4(i,j,:,:) = 0.0
         do m = 1, nrens
         X4(i,j,m,m) = 1.0
         enddo

             nobs = 0 ! no upper limit on the number of local observations
             call get_local_obs(i, j, nrobs, nrens, obs, radius * 1000.0,&
                  modlon, modlat, ni, nj, nobs, S, DD)!, tmp_S, tmp_DD)
             nobs_array(i, j) = nobs

          if (nobs == 0) then
             ! just in case
             X4(i,j,:,:) = 0.0
             do m = 1, nrens
                X4(i,j,m,m) = 1.0
             enddo
             deallocate(tmp_S, tmp_DD, tmp_obs)
             cycle
          end if

          ! Allocate local arrays
          allocate(subS(nobs, nrens))
          allocate(subD(nobs, nrens))
          allocate(G(nrens, nobs))
          if (nobs < nrens) then
             allocate(X1(nobs, nobs))
          else
             allocate(X1(nrens, nrens))
          end if

          subD = tmp_DD(:,:)
          subS = tmp_S(:,:)

          ! taper ensemble observation anomalies and innovations
          !
          if (LOCFUN_USED /= LOCFUN_NONE) then
             do o = 1, nobs
                obs0 = tmp_obs(o)
                dist = spherdist(modlon(i, j), modlat(i, j),&
                     obs0 % lon, obs0 % lat)
                lfactor = max(locfun(dist / radius / 1000.0d0), 0.0d0)

                subS(o, :) = subS(o, :) * lfactor
                subD(o, :) = subD(o, :) * lfactor
             end do
          end if

          if (nobs < nrens) then
             ! Construct matrix (S * S' + I) - to be inverted
             !
             X1 = matmul(subS, transpose(subS))
             do o = 1, nobs
                X1(o, o) = X1(o, o) + 1.0d0
             end do
             
             ! Inversion via Cholesky decomposition, done in two stages
             call dpotrf('U', nobs, X1, nobs, lapack_info)
             call dpotri('U', nobs, X1, nobs, lapack_info)
             
             ! fill the lower triangular part of (symmetric) X1
             do o = 2, nobs
                X1(o, 1 :  o - 1) = X1(1 : o - 1, o)
             end do

             ! X4tmp is subS^T * X3              
             G = matmul(transpose(subS), X1)
          else ! nobs >= nrens
             X1 = matmul(transpose(subS), subS)
             do m = 1, nrens
                X1(m, m) = X1(m, m) + 1.0d0
             end do

             call dpotrf('U', nrens, X1, nrens, lapack_info)
             call dpotri('U', nrens, X1, nrens, lapack_info)
             
             do m = 2, nrens
                X1(m, 1 :  m - 1) = X1(1 : m - 1, m)
             end do
             G = matmul(X1, transpose(subS))
          end if
          X4tmp = matmul(G, subD)

          deallocate(tmp_DD)
          deallocate(tmp_S)
          deallocate(tmp_obs)
          deallocate(X1, subD, subS, G)

          ! Add I: A^a = A^f + A^f * X4 = A^f * (I + X4)          
          do m = 1, nrens
             X4tmp(m, m) = X4tmp(m, m) + 1.0d0
          enddo

          X4(i,j,:,:) = real(X4tmp, 4)
       end do ! i = 1, ni
    end do ! j = my_first_iteration, my_last_iteration

  end subroutine pre_local_analysis

  real function locfun(x)
    real, intent(in) :: x

    real :: xx, xx2, xx3

    select case(LOCFUN_USED)

    case (LOCFUN_NONE)
       locfun = 1.0
    case (LOCFUN_STEP)
       if (x > 1.0) then
          locfun = 0.0
       else
          locfun = 1.0
       end if
    case (LOCFUN_GASPARI_COHN)
       if (x > 1.0) then
          locfun = 0.0
       else
          xx = x * 2.0
          xx2 = xx * xx
          xx3 = xx2 * xx
          if (xx < 1.0) then
             locfun = 1.0 + xx2 * (- xx3 / 4.0 + xx2 / 2.0)&
                  + xx3 * (5.0 / 8.) - xx2 * (5.0 / 3.0)
          else
             locfun = xx2 * (xx3 / 12.0 - xx2 / 2.0)&
                  + xx3 * (5.0 / 8.0) + xx2 * (5.0 / 3.0)&
                  - xx * 5.0 + 4.0 - (2.0 / 3.0) / xx
          end if
       end if
    case default
       print *, 'ERROR: m_local_analysis.F90: locfun(): LOCFUN_USED =', LOCFUN_USED, 'is unknown'
       stop
    end select
  end function locfun

  !
  subroutine get_local_obs(i, j, nrobs, nrens, obs, rmax, modlon, modlat, &
       ni, nj, nobs, S, DD)!, tmp_S, tmp_DD)
    use mod_measurement
    use m_spherdist
    use m_obs, only: tmp_S, tmp_DD, tmp_obs

    implicit none
    real, dimension(nrobs, nrens), intent(in) :: S, DD
!    real, allocatable, dimension(:,:), intent(out) :: tmp_S, tmp_DD
    integer, intent(in) :: i, j
    integer, intent(in) :: nrobs
    integer, intent(in) :: nrens
    type(measurement), intent(in) :: obs(nrobs)
    real, intent(in) :: rmax ! maximal allowed distance
    real, intent(in) :: modlon(ni, nj)
    real, intent(in) :: modlat(ni, nj)
    integer, intent(in) :: ni, nj
    integer, intent(inout) :: nobs ! input : max allowed # of local obs
                                   ! output: actual # of local obs for this
                                   !         point

    integer :: ngood, ngood2
    integer :: sorted(nrobs)
    real :: dist(nrobs)
    integer :: o

    ngood = 0
    do o = 1, nrobs
       dist(o) = spherdist(obs(o) % lon, obs(o) % lat, modlon(i, j), modlat(i, j))
       if (dist(o) <= rmax) then
          ngood = ngood + 1
       end if
    end do

    allocate (tmp_DD(ngood,nrens))
    allocate (tmp_S(ngood,nrens))
    allocate (tmp_obs(ngood))

    ngood = 0
    do o = 1, nrobs
       dist(o) = spherdist(obs(o) % lon, obs(o) % lat, modlon(i, j), modlat(i, j))
       if (dist(o) <= rmax) then
          ngood = ngood + 1
          tmp_DD(ngood,:) = DD(o,:)
          tmp_S(ngood,:) = S(o,:)
          tmp_obs(ngood) = obs(o)
       end if
    enddo

    if (nobs <= 0 .or. nobs >= ngood) then
       nobs = ngood
    end if
  end subroutine get_local_obs

end module m_local_analysis
