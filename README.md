# networks_eco_evo

This repository contains simulated data and associated code for the manuscript 'Evolution reverses the effect of network structure on metapopulation persistence.'

**Authors:**   
Lisa C. McManus<sup>1,2</sup>, Edward W. Tekwa<sup>1</sup>, Daniel E. Schindler<sup>3</sup>, Timothy E. Walsworth<sup>4</sup>, Madhavi A. Colton<sup>5</sup>, Michael M. Webster<sup>6</sup>, Timothy E. Essington<sup>3</sup>, Daniel L. Forrest<sup>1</sup>, Stephen R. Palumbi<sup>7</sup>, Peter J. Mumby<sup>8</sup> and Malin L. Pinsky<sup>1</sup>

**Affiliations:**  
<sup>1</sup>Department of Ecology, Evolution, and Natural Resources, Rutgers University, New Brunswick, NJ, USA.  
<sup>2</sup>Hawaiʻi Institute of Marine Biology, University of Hawaiʻi at Manoa, Kaneʻohe, HI 96744, USA.  
<sup>3</sup>School of Aquatic and Fishery Sciences, university of Washington, Seattle, WA, USA.  
<sup>4</sup>Department of Watershed Sciences, Utah State University, Logan, UT, USA.  
<sup>5</sup>Coral Reef Alliance, Oakland, CA, USA.  
<sup>6</sup>Department of Environmental Studies, New York University, 285 Mercer St., New York, NY 10003, USA.  
<sup>7</sup>Department of Biology, Hopkins Marine Station, Stanford University, Pacific Grove, CA, USA.  
<sup>8</sup>Marine Spatial Ecology Laboratory, School of Biological Sciences, The University of Queensland, St Lucia, Queensland, Australia.  

**Notes:**
To regenerate figures in the main text, use the Jupyter notebook 'Figures for McManus et al. (2021) Ecology.ipnyb.' The data to generate these figures should be placed in the 'Data' directory. Some files are already there, but most of the data need to be downloaded from here: https://doi.org/10.6084/m9.figshare.14245733.v1.  
We have also included the scripts used for calculating numerical solutions on a high-performance computing system. Note that you may need to modify the syntax and install packages depending on the particular HPC that you are using. These data were generated on the Rutgers Amarel cluster.

A description for each of the HPC scripts is as follows:

**Random network simulations** (in HPC_scripts/random):  
1. functions.py: contains all functions called in the numerical solver
2. generate_random_matrix.py: additional functions that allow generation of random network configurations  
3. random_batch_param.py: contains the parameter values and number of iterations  
4. routine_random_network.py: the actual routine to run simulations under constant and increasing temperatures (produces data)
5. submit_random_batch.sh: file to submit 'routine_random_network.py' to the HPC
6. random_summary_files.py: calculate mean final values generated from the routine_random_network.py runs (produces summary data used in figure generation)
7. submit_random_summary_files.sh: file to submit 'random_summary_files.py' to the HPC
8. random_summary_files_TS.py: calculate mean time series generated from the routine_random_network.py runs (produces summary data used in figure generation)
9. submit_random_summary_files.sh: file to submit 'random_summary_files_TS.py' to the HPC


**Regular network simulations** (in HPC_scripts/random):  
1. functions.py: contains all functions called in the numerical solver
2. D_reg.npy: connectivity/dispersal matrix for regular network simulations
3. regular_batch_param.py: contains the parameter values and number of iterations  
4. routine_regular_network.py: the actual routine to run simulations under constant and increasing temperatures (produces data)
5. submit_regular_batch.sh: file to submit 'routine_regular_network.py' to the HPC
6. regular_summary_files.py: calculate mean final values generated from the routine_regular_network.py runs (produces summary data used in figure generation)
7. submit_regular_summary_files.sh: file to submit 'regular_summary_files.py' to the HPC
8. regular_summary_files_TS.py: calculate mean time series generated from the routine_regular_network.py runs (produces summary data used in figure generation)
9. submit_regular_summary_files.sh: file to submit 'regular_summary_files_TS.py' to the HPC

