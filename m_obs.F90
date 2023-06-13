module m_obs
   use ice_kinds_mod
   use mod_measurement

   implicit none
   integer, save :: ininrobs,lskip
   type(inimeasurement4),allocatable,save :: iniobs(:),sit_var(:),sit_var_cs2(:),sit_var_smos(:),sit_var_landy(:)
   type(measurement),allocatable,save :: obs(:), tmp_obs(:)
   real*8,allocatable,save :: gE(:)
   real*8,allocatable,save :: gathgE(:)
   real,allocatable,save :: E(:,:)
   real,allocatable,save :: tmp_DD(:,:), tmp_S(:,:)
end module m_obs
