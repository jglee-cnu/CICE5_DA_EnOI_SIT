# CICE5_DA_EnOI_SIT
This repository contains the source codes for sea ice data assimilation using CICE5 model 

They are source codes for assimilating satellite-derived sea ice thickness data for the CICE5 (Community Ice CodE) model.

A brief description of each code is provided below.

The source code is partly based on the code created by Geir Evensen (https://enkf.nersc.no/Code/)

* m_main_enoi.F90
: This code is used for coupling the data assimilation system with CICE5.
The background error perturbations, background fields, and analysis fields are dealt with in this code.

* EnOI.F90
: This code processes the observations and model information needed for data assimilation and includes analysis step.

* m_active_obs.F90
: This module code activates the information of observation to be assimilated.

* m_prep_4_EnOI.F90
: This module code deals with observation increments and model error information in the observation space.

* m_local_analysis.F90
: This module code contains the analysis step based on the local analysis method.

* m_Generate_element_Sij.F90  &  m_modstate_point
: These module codes are used to get the model information in the observation space.

* m_spherdist.F90  &  spherical_dists.F90 
: They are used for computing the distance between geographical positions.

* m_read_obs.F90
: This module code is used for reading the information of satellite observation data.

* m_obs.F90  &  mod_measurement.F90
: They are used to save information related to the observations.

* mod_dimensions.F90  &  mod_states.F90
: They are used to save model information.

* common_mpi.F90
: This code is used for conducting the parallelization based on MPI (Message Passing Interface).
