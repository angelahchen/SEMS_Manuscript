
############# Parameters for the setting ###################
  simulated.QTL.info.dir = "/Simulation/Output/Setting_1/"; # where the simulated QTNs are 
  summary.dir = "/Summary_Statistics/Setting_1/GLM/"; # where the summary statistics tables go
  simulation.results.dir = "/Analysis/GLM_after_FDR_at_0.5/"; # where the analysis results (after applying FDR) file is 
  file.name.prefix= "Setting_1_GLM_FDR_adjusted" # file name of the analysis results (after applying FDR)
  number.of.simulated.additive.QTN = 1; # when there's no additive QTN, enter 1
  number.of.simulated.epistatic.QTN.pairs = 4;
  setting.number = 1;
  include.additive.QTL = FALSE;
  include.epistatic.QTL = TRUE;

############# Global parameters ###################
  setwd()
  home.dir = getwd() 
  
  #Set up some global variables
  flanking.bp.dist.from.QTN = 2000000 # Adjust this to create different window sizes
  number.of.traits.per.setting = 100
  chm.to.analyze = c(1,2,3,4,5,6,7,8,9,10)
  number.of.simulations = 100
  size.of.maize.chromosomes = c(301331005, 237007151, 232096209, 241402927, 217748528,
                                169132432, 176664055, 175775017, 156591811, 150168908) #This will make the Manhattan plots look good
  
################Create statistics and plots #####################
  setwd("~/Code")
  source("Generate_summary_statistics_GLM.R")
  assess.the.results()
  