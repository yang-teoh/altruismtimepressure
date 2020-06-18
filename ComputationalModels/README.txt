The Gaze-Informed Attentional Drift Diffusion Model

This folder contains files for the running the computational model outlined in "Attentional priorities drive effects of time pressure on altruistic choice" by Y.Y. Teoh, Z. Yao, W. Cunningham & C. Hutcherson. 

:Version: 1.0
:Authors: Yi Yang Teoh 
:Copyright (c) 2020 Yi Yang Teoh
:License: GNU-GPLv3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License 
along with this program (license.txt).  If not, see <https://www.gnu.org/licenses/>

Python Files: 

giADDM.py
Custom python functions implementing & modifying the ADDM (Krabjich et al,2010) to use real-time gaze data to amplify attribute weights in drift.

simul_giADDM.py
Custom python script containing code to simulate fixations & choice data for one subject.

evolveMCMC.py
Custom python functions adapting the deMCMC fitting procedures described in Holmes & Trueblood (2018) with modifications to include migrating chains described in (Hu & Tsui, 2010; Turner et al., 2013).
See: https://doi.org/10.1016/j.cogpsych.2015.11.002 & https://osf.io/dkmbk/ for Holmes & Trueblood (2018) 

Fit_giADDM.py
Custom python script containing fitting procedures for the giADDM one subject. [Fitting time for Simulated Subject 1 on full data: ~ 1600 core-hours, ~ 400hours x 4 cores (standard dual-core (4 logical cores) computer)]


inspect_giADDM.py
Custom python script to analyze fitted MCMC parameter estimates: retrieves MCMC parameters, transformed giADDM parameters, convergence statistics, model predictions/simulations for each trial


Fit_multiDDM.py
Custom python script containing fitting procedures for the MultiDDM one subject. [Fitting time for Simulated Subject 1 on full data: ~ 1600 core-hours, ~ 400hours x 4 cores (standard dual-core (4 logical cores) computer)]


inspect_multiDDM.py
Custom python script to analyze fitted MCMC parameter estimates: retrieves MCMC parameters, transformed giADDM parameters, convergence statistics, model predictions/simulations for each trial

Directories:
Analysis	
## Models: biasedDDM; giADDM (full); giADDM_ngb (fixed genbias = 0)
Output from Fit_giADDM.py for 50 subjects
Output from inspect_giADDM.py [summary]

	
