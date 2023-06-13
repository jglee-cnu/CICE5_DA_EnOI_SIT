module m_prep_4_EnOI

contains
!subroutine prep_4_EnOI(obs,gA,nrens,nrobs,E,D,S,zt)
subroutine prep_4_EnOI(gA,nrens,nrobs,D,S)
   use ice_kinds_mod
   use ice_communicate, only: my_task, master_task
   use ice_fileunits, only: nu_diag, nu_forcing
   use ice_read_write, only: ice_open
   use mod_dimensions
   use mod_states
   use mod_measurement
   use m_obs
! Functions and subroutines
   use m_obs_pert
   use m_Generate_element_Sij

   implicit none

! In put variables
   integer,           intent(in) :: nrobs         ! Number of measurements
   integer,           intent(in) :: nrens         ! Size of ensemble
   type(gstates),     intent(in) :: gA(nrens)     ! Ensemble of model states 
   real, dimension(nz) :: zt
   real :: tt

! Output variables
   real, intent(inout) :: D(nrobs,nrens)
   real, intent(inout) :: S(nrobs,nrens)

! Local variables 
   real meanS(nrobs)    ! Automatic array
   integer i,j,m,iens,o

   character (char_len_long) :: &
           err_sit_file

   character (char_len_long) :: &
           err_data_dir

   logical :: isnan

! Observe ensemble to construct the matrix S=HA (Hxb)
   do iens =1, nrens
    do m =1, nrobs
     S(m,iens)      =  Generate_element_Sij(obs(m),gA(iens))
    enddo
   enddo

! Construct ensemble of measurements D=d+E
   do j=1,nrens
   do m=1,nrobs
      D(m,j)=obs(m)%d
   enddo
   enddo

! Compute innovation D_prime=D-HA
   D=D-S

! Compute mean(HA) 
   meanS=0.0
   do j=1,nrens
   do m=1,nrobs
      meanS(m)=meanS(m)+S(m,j)
   enddo
   enddo
   meanS=(1.0/float(nrens))*meanS

! Compute HA_prime=HA-mean(HA)  Heb
   do j=1,nrens
      S(:,j)=S(:,j)-meanS(:)
   enddo

end subroutine prep_4_EnOI
end module m_prep_4_EnOI

