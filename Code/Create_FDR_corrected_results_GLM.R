source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
library(multtest)

############# Parameters for the setting ###################
GLM_results = "Setting_1_GLM"; # file name of the analysis results (before applying FDR)


### Post-process the GLM results 
#####################################################
setwd("~/Analysis/GLM/") # where the analysis results file is (original file before applying FDR)
PWI<- read.table(paste(GLM_results,".txt",sep=""), header = TRUE) 

#Parse results by traits (100 traits totally; 100*1106=110600 phenotype values)
PWI.list <- lapply(1:100, function(i){
  return(PWI[((i-1)*1106+1):((i-1)*1106+1106),])
})

###Create a modified version of the BH_FDR function 
#####################################################
"GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure"<-
  function(PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH"){
    #Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
    #Output: PWIP, number.of.significant.SNPs
    #Authors: Alex Lipka and Zhiwu Zhang 
    # Last update: May 5, 2011 
    ##############################################################################################
    #Make sure that your compouter has the latest version of Bioconductor (the "Biobase" package) and multtest
    
    # if(is.null(PWI))
    # {
    #   PWIP=NULL
    #   number.of.significant.SNPs = 0
    # }
    
    # if(!is.null(PWI))
    # {  
    
    #library(multtest)
    
    #  if(dim(PWI)[1] == 1){
    #    PWIP <- cbind(PWI, PWI[4])
    #    colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
    #  }
    
    #if(dim(PWI)[1] > 1){  # number of row 
    #mt.rawp2adjp Performs the Simes procedure.  The output should be two columns, Left column: originial p-value
    #Right column: Simes corrected p-value
    
    res <- mt.rawp2adjp(PWI[,6], proc=FDR.Procedure, alpha=FDR.Rate, na.rm = FALSE)
    
    #This command should order the p-values in the order of the SNPs in the data set
    adjp <- res$adjp[order(res$index), ]
    
    #round(adjp[1:7,],4)
    #Logical statment: 0, if Ho is not rejected; 1, if  Ho is rejected, by the Simes corrected p-value
    #  temp <- mt.reject(adjp[,2], FDR.Rate)
    
    #Lists all number of SNPs that were rejected by the BY procedure
    #temp$r
    
    #Attach the FDR adjusted p-values to AS_Results
    
    PWIP <- cbind(PWI, adjp[,2])
    
    #Sort these data by lowest to highest FDR adjusted p-value
    PWIP <- PWIP[order(PWIP[,14]),]
    
    colnames(PWIP)[14] <- "FDR_Adjusted_P-values"
    #  number.of.significant.SNPs = temp$r
    #}
    #print("GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure accomplished successfully!")
    # }  
    #return(list(PWIP=PWIP, number.of.significant.SNPs = number.of.significant.SNPs))
    return(PWIP)
  }#GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure ends here

###Execute the BH-FDR function for GLM 
#####################################################

##Markers are adjusted within each trait (100 traits totally)
##PWIP[[1]] stores FDR.adjusted p values that are smaller than 0.05 from the first trait
PWIP<- lapply(1:100, function(i){
  tmp_PWIP<-GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.list[[i]], FDR.Rate = 0.05, FDR.Procedure = "BH")
  return(tmp_PWIP[which(tmp_PWIP[,14]<0.05),])
})

# PWIP<- lapply(1:100, function(i){
#   tmp_PWIP<-GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.list[[i]], FDR.Rate = 0.05, FDR.Procedure = "BH")
#   return(tmp_PWIP)
# })
###Post-processing GLM results 
#####################################################

## This concatenate each data frame. Using Reduce() to "rbind" is faster than using a for loop!
GLM <- Reduce(rbind, PWIP) 

## To check whether all data frames are merged together correctly
#nrow=c()
#for (i in 1:100){
#  nrow[i] <- nrow(PWIP[[i]])
#}
#sum(nrow)

################Create FDR corrected results ##################### 
## Write a table for the GLM_results  
#write.table(GLM, file=paste("FDR_adjusted_", GLM_results, ".txt", sep = ""), sep="", eol="\n", na="NA", row.names = TRUE, col.names = TRUE) 
setwd("~/Analysis/GLM_after_FDR_at_0.5")
write.csv(GLM, file=paste(GLM_results, "_FDR_adjusted",".csv", sep = ""),
          eol = "\n", na = "NA", row.names = TRUE)

## Note: the generated csv file contains 1 extra column [,1]. You need to delete this column before processing "Parameters_for_summary_statistics_GLM.R"

# #######