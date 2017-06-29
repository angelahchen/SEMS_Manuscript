# SEMS_Manuscript

Sixteen simulation settings were analyzed using Stepwise Epistatic Model Selection (SEMS), Joint Linkage Analysis (JL), and Generalized Linear Model (GLM).

All the raw data and results are in an external cloud repository. You can access them through this link: https://uofi.box.com/s/krdkzkui1tzoh3jvov5a2cvsh7rrwl3z 

There are five folders in this repository: 
Simulation: contains phenotypic and genotypic data organized by setting 
Analysis: contains analysis results organized by method 
Summary_Statistics: contains summary statistics tables organized by setting 
Code: contains the seven R scripts on Github 
Read_Me: contains this read.me file

“Paths” in the following description refer to where specific files were stored in this repository.  

First, genotypic values were randomly selected from the NAM dataset (path to file: ~/Simulation/Input/) to simulate phenotypic values for each setting (path to simulations for each setting: ~/Simulation/Output/). Simulate_QTNs.R is applied. 

Then, all sixteen simulation settings were analyzed using SEMS (path to raw results: ~/Analysis/SEMS), JL (path to raw results: ~/Analysis/JL), and GLM (path to results after correcting for multiple testing: ~/Analysis/GLM_after_FDR_at_0.5). 

Note: Create_FDR_corrected_results_GLM.R corrects for multiple testing from raw GLM results. 

Finally, summary statistics tables for each analysis were created by the following R scripts:
* For SEMS: 
Parameters_for_summary_statistics_SEMS.R 
Generate_summary_statistics_SEMS.R 
* For JL:
Parameters_for_summary_statistics_JL.R
Generate_summary_statistics_JL.R
* For GLM: 
Parameters_for_summary_statistics_GLM.R
Generate_summary_statistics_GLM.R 
