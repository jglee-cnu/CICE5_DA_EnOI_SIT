## CICE5_DA_EnOI_SIT
This repository contains the source codes for assimilating the satellite-based sea ice thickness using CICE5 (Community Ice CodE) model 

The CICE5 sea ice model is available via (https://github.com/CICE-Consortium)
The source code in this repo is partly based on the code created by Geir Evensen (https://enkf.nersc.no/Code/)

## Process of sea ice data assimilation in CICE5
This code is compiled together with source code for CICE5 to configure, and data assimilation can be performed within the CICE5 driver source code by calling the subroutine in m_main_enoi.F90.

```

  do iblk = 1, nblocks
  if (mod(real(istep),secday/dt) .eq. 0) then
  call main_enoi(iblk) 
  endif;enddo

```

The command lines above should be inserted after dynamic/thermodynamic processes in the driver code of CICE5 (subroutine ice_step in CICE_RunMod.F90).

## A brief description of each code

### m_main_enoi.F90
: This code is used for coupling the data assimilation system with CICE5.
The background error perturbations, background fields, and analysis fields are dealt with in this code.

### EnOI.F90
: This code processes the observations and model information needed for data assimilation and includes analysis step.

### m_active_obs.F90
: This module code activates the information of observation to be assimilated.

### m_prep_4_EnOI.F90
: This module code deals with observation increments and model error information in the observation space.

### m_local_analysis.F90
: This module code contains the analysis step based on the local analysis method.

### m_Generate_element_Sij.F90  &  m_modstate_point
: These module codes are used to get the model information in the observation space.

### m_spherdist.F90  &  spherical_dists.F90 
: They are used for computing the distance between geographical positions.

### m_read_obs.F90
: This module code is used for reading the information of satellite observation data.

### m_obs.F90  &  mod_measurement.F90
: They are used to save information related to the observations.

### mod_dimensions.F90  &  mod_states.F90
: They are used to save model information.

### common_mpi.F90
: This code is used for conducting the parallelization based on MPI (Message Passing Interface).

## References
### Lee J-G and Ham Y-G (2022) Satellite-Based Data Assimilation System for the Initialization of Arctic Sea Ice Concentration and Thickness Using CICE5. Front. Clim. 4:797733. doi: 10.3389/fclim.2022.797733
: This paper explains about the data assimilation method, generation of background error covariance, observation information. The quality of the sea ice reanalysis from data assimilation is also investigated.

## Requirement 
* intel compiler
* LAPACK — Linear Algebra PACKage (version 3.7.0 is used in this work)
* MPI parallelization program
