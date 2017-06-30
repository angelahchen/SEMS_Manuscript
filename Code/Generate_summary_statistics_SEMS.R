# #####################################################

######################################
################Begin the function
assess.the.results <- function(){
  #Read in the information on thesimulated QTN
  #setwd(paste(home.dir, simulated.QTL.info.dir, sep = ""))
  the.simulated.additive.QTN <- read.table(paste("Genotypic.information.for.",number.of.simulated.additive.QTN,".Additive.QTN.txt", sep = ""), head = TRUE)
  the.simulated.epistatic.QTN <- read.table(paste("Genotypic.information.for.",number.of.simulated.epistatic.QTN.pairs,".Epistatic.QTN.txt", sep = ""), head = TRUE)
  the.results <- read.csv(paste(file.name.prefix, ".csv",sep = ""), head = TRUE)
  #the.results <- the.results[-1,]
  #For each heritability
  #setwd(paste(home.dir, simulation.results.dir, sep = ""))
  count <- 0
  
  
  
  #Remove the "mean", "family", and "error" terms from the output
  the.results <- the.results[-which((as.character(the.results[,2])=="mean")|
                                      (as.character(the.results[,2])=="family")|
                                      (as.character(the.results[,2])=="Error")|
                                      (as.character(the.results[,2])=="NaN")|
                                      is.na(the.results[,2]) ), ]
  if(nrow(the.results) != 0){
    #Make graph 3 #NOTE: Turn this into a function
    create.graph.3(the.results = the.results, the.simulated.additive.QTN  =  the.simulated.additive.QTN,
                   the.simulated.epistatic.QTN = the.simulated.epistatic.QTN, setting.number = setting.number,
                   include.additive.QTL = include.additive.QTL, include.epistatic.QTL = include.epistatic.QTL,
                   summary.dir = summary.dir)  
    
    the.summary.statistics <- create.summary.statistics(the.results = the.results, the.simulated.additive.QTN  =  the.simulated.additive.QTN,
                                                        the.simulated.epistatic.QTN = the.simulated.epistatic.QTN, setting.number = setting.number,
                                                        include.additive.QTL = include.additive.QTL, include.epistatic.QTL = include.epistatic.QTL,
                                                        flanking.bp.dist.from.QTN = flanking.bp.dist.from.QTN, number.of.traits.per.setting =number.of.traits.per.setting) 
    
    print(the.summary.statistics)
    
    write.table(the.summary.statistics$additive.QTN.signals.detected, paste("Setting.",setting.number,".Flank.Dist.",flanking.bp.dist.from.QTN,
                                                                            ".Add.Signals.txt", sep = " "), quote = FALSE, sep = "\t", 
                row.names = FALSE,col.names = TRUE)
    write.table(the.summary.statistics$epistatic.QTN.signals.detected, paste("Setting.",setting.number,".Flank.Dist.",flanking.bp.dist.from.QTN,
                                                                             ".Epi.Signals.txt", sep = " "), quote = FALSE, sep = "\t", 
                row.names = FALSE,col.names = TRUE)
    write.table(the.summary.statistics$proportion.of.remaining.additive.signals, paste("Setting.",setting.number,".Flank.Dist.",flanking.bp.dist.from.QTN,
                                                                                       ".Add.False.Signals.txt", sep = " "), quote = FALSE, sep = "\t", 
                row.names = FALSE,col.names = TRUE)
    write.table(the.summary.statistics$proportion.of.remaining.epistatic.signals, paste("Setting.",setting.number,".Flank.Dist.",flanking.bp.dist.from.QTN,
                                                                                        ".Epi.False.Signals.txt", sep = " "), quote = FALSE, sep = "\t", 
                row.names = FALSE,col.names = TRUE)
    
  }
  
}#end assess.the.results()


######################################
create.graph.3 <- function(the.results = the.results, the.simulated.additive.QTN  =  the.simulated.additive.QTN,
                           the.simulated.epistatic.QTN = the.simulated.epistatic.QTN, setting.number = setting.number,
                           include.additive.QTL = NA, include.epistatic.QTL = NA, summary.dir = NULL){
  
  #Use either table() to create a frequency table for each of the markers that appear
  the.proportions.of.markers <- as.data.frame(table(as.character(the.results[,2])))
  the.proportions.of.markers[,1] <- as.character(the.proportions.of.markers[,1])
  
  
  #Create a matrix consisting of: Peak Marker Name, Chromosome, Position
  #For loop through each marker (i.e., rows of "the.proportions.of.markers")
  count <- 1
  the.chromosomes <- NULL
  the.positions <- NULL
  the.line.color <- NULL
  the.lty <- NULL
  the.pch <- NULL
  the.proportion <- NULL
  for(i in 1:nrow(the.proportions.of.markers)){
    #If it is just one marker:
    if(grepl("\\+", the.proportions.of.markers[i,1])){
      
      #Parse out the chromosomes, positions, line color = purple, lty = count, pch = count
      this.two.chromosomes.and.positions <- unlist(strsplit(as.character(the.proportions.of.markers[i,1]),"\\+"))
      this.chromosome <- as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[1]),"_"))[1],4,5))
      this.chromosome.2 <- as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[2]),"_"))[1],4,5))
      
      this.position <- as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[1]),"_"))[2],1,25))
      this.position.2 <- as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[2]),"_"))[2],1,25))  
      
      these.line.colors <- c("green", "green")
      these.lty <- c(count,count)
      these.pch <- c(15, 15)
      
      the.chromosomes <- c(the.chromosomes, this.chromosome, this.chromosome.2)
      the.positions <- c(the.positions, this.position, this.position.2)
      the.line.color <- c(the.line.color, these.line.colors)
      the.lty <- c(the.lty, these.lty)
      the.pch <- c(the.pch, these.pch)
      
      #Make sure you add the proportion of times it is added to the model
      the.proportion <- c(the.proportion, the.proportions.of.markers[i,2], the.proportions.of.markers[i,2])
      count <- count + 1
      
      #Remove all of the objects created here
      rm(this.two.chromosomes.and.positions) 
      rm(this.chromosome) 
      rm(this.chromosome.2) 
      
      rm(this.position) 
      rm(this.position.2)  
      
      rm(these.line.colors)
      rm(these.lty)
      rm(these.pch) 
      
    }else{#If it is an pair of markers:
      #Parse out the chromosome, position, line color = black, lty = 1, pch = 1
      this.chromosome <- as.numeric(substr(unlist(strsplit(as.character(the.proportions.of.markers[i,1]),"_"))[1],4,5))
      this.position <- as.numeric(substr(unlist(strsplit(as.character(the.proportions.of.markers[i,1]),"_"))[2],1,25))
      this.line.color <- "orange"
      this.lty <- 1
      this.pch <- 17
      
      the.chromosomes <- c(the.chromosomes, this.chromosome)
      the.positions <- c(the.positions, this.position)
      the.line.color <- c(the.line.color, this.line.color)
      the.lty <- c(the.lty, this.lty)
      the.pch <- c(the.pch, this.pch)
      
      #Make sure you add the proportion of times it is added to the model
      the.proportion <- c(the.proportion, the.proportions.of.markers[i,2])
      
      #Remove all of the objects created here
      rm(this.chromosome) 
      
      rm(this.position) 
      
      rm(this.line.color)
      rm(this.lty)
      rm(this.pch)
      
    }
    
  }#End for(i in 1:nrow(the.proportions.of.markers))
  
  #Create the input file for the Manhattan plot
  the.input.file.for.Manhattan.plot <- cbind(the.proportion, the.chromosomes, the.positions)
  
  
  #Create the additive QTN
  additive.QTN.chr <- NULL
  additive.QTN.pos <- NULL
  for(i in 1:nrow(the.simulated.additive.QTN)){
    this.chromosome <- as.numeric(substr(unlist(strsplit(as.character(the.simulated.additive.QTN[i,1]),"_"))[1],4,5))
    this.position <- as.numeric(substr(unlist(strsplit(as.character(the.simulated.additive.QTN[i,1]),"_"))[2],1,25))
    additive.QTN.chr <- c(additive.QTN.chr, this.chromosome)
    additive.QTN.pos <- c(additive.QTN.pos, this.position)
    rm(this.chromosome)
    rm(this.position)
  }
  
  the.additive.QTN <- cbind(additive.QTN.chr, additive.QTN.pos)
  
  color.vector.QTN <- NULL
  lty.vector.QTN <- NULL
  for(i in 1:nrow(the.additive.QTN)){
    QTN.row <- c(NA, the.additive.QTN[i,1], the.additive.QTN[i,2])
    if(i == 1){
      QTN.matrix <- QTN.row
    }else{
      QTN.matrix <- rbind(QTN.matrix, QTN.row)
    }
    if(include.additive.QTL){
      color.vector.QTN <- c(color.vector.QTN, "black")
      lty.vector.QTN <- c(lty.vector.QTN, 1)
    }
  }#end for(i in 1:number.of.simluated.QTL)
  if(!is.null(dim(QTN.matrix))){
    colnames(QTN.matrix) <- colnames(the.input.file.for.Manhattan.plot)
  }else{
    names(QTN.matrix) <- colnames(the.input.file.for.Manhattan.plot)
  }
  if(include.additive.QTL){
    data <- rbind(the.input.file.for.Manhattan.plot, QTN.matrix)
  }else{
    data <- the.input.file.for.Manhattan.plot
  }
  
  
  
  #Create the epistatic QTN
  epistatic.QTN.chr <- NULL
  epistatic.QTN.pos <- NULL
  for(i in 1:nrow(the.simulated.epistatic.QTN)){
    #Old code
    this.chromosome <- as.numeric(substr(unlist(strsplit(as.character(the.simulated.epistatic.QTN[i,1]),"_"))[1],4,5))
    this.position <- as.numeric(substr(unlist(strsplit(as.character(the.simulated.epistatic.QTN[i,1]),"_"))[2],1,25))
    epistatic.QTN.chr <- c(epistatic.QTN.chr, this.chromosome)
    epistatic.QTN.pos <- c(epistatic.QTN.pos, this.position)
    rm(this.chromosome)
    rm(this.position)
  }
  
  the.epistatic.QTN <- cbind(epistatic.QTN.chr, epistatic.QTN.pos)
  
  
  count.mod <- 1
  for(i in 1:nrow(the.epistatic.QTN)){
    QTN.row <- c(NA, the.epistatic.QTN[i,1], the.epistatic.QTN[i,2])
    if(i == 1){
      QTN.matrix.epi <- QTN.row
    }else{
      QTN.matrix.epi <- rbind(QTN.matrix.epi, QTN.row)
    }
    if(include.epistatic.QTL){ 
      color.vector.QTN <- c(color.vector.QTN, "purple")
      lty.vector.QTN <- c(lty.vector.QTN, count.mod)
    }
    if((i %% 2) == 0) count.mod <- count.mod + 1
  }#end for(i in 1:number.of.simluated.QTL)
  colnames(QTN.matrix.epi) <- colnames(the.input.file.for.Manhattan.plot)
  if(include.epistatic.QTL){
    data <- rbind(data, QTN.matrix.epi) # I changed the name "the.input.file.for.Manhattan.plot" to "data" so that I will not have to change too much code
  }
  
  
  
  #Add the approximate size of each chromosome to the data input file
  the.chromosomes.end <- cbind(rep(NA, length(chm.to.analyze)), seq(1:length(chm.to.analyze)), size.of.maize.chromosomes)
  colnames(the.chromosomes.end) <- colnames(data)
  data <- rbind(data, the.chromosomes.end)
  
  #Add some additional rows ofa bp of "1" for each chromosome. This should help with getting the chromosome number to be at the center between the chromsomes on
  # the X-axis of the Manhattan plot
  the.chromosomes.start <- cbind(rep(NA, length(chm.to.analyze)), seq(1:length(chm.to.analyze)), rep(1, length(chm.to.analyze)))
  colnames(the.chromosomes.start) <- colnames(data)
  data <- rbind(data, the.chromosomes.start)
  
  
  #Use the NAM GWAS results plot code to make an informative Manhattan Plot
  lastbase <- 0
  ticks <- NULL
  
  #change base position to accumulatives
  for (j in chm.to.analyze)
  {
    index=(data[,2]==j)
    ticks <- c(ticks, lastbase+mean(data[index,3]))
    data[index,3]=data[index,3]+lastbase
    lastbase=max(data[index,3])
  }
  
  
  ticks.end <- NULL
  for(j in chm.to.analyze[-length(chm.to.analyze)]){
    data.subset <- data[which(data[,2]==j),]
    ticks.end <- c(ticks.end, max(data.subset[,3]))
  }
  
  #Set the limits of the y (RMIP) 
  
  
  #y.lim <- ceiling(max(data[which(!is.na(data[,1])),1]))/number.of.simulations
  y.lim <- 1
  x <- data[-which(is.na(data[,1])) , 3]
  y <- data[-which(is.na(data[,1])),1]/number.of.simulations
  the.QTN <- data[which(rownames(data)!=""), 3]
  
  
  
  #setwd(paste(home.dir, summary.dir, sep = ""))  
  pdf(paste("Manhattan_Plot_",file.name.prefix,i,"_", number.of.simulated.additive.QTN, "_Add_QTN_",number.of.simulated.epistatic.QTN.pairs,
            "_Epi_QTN_Pairs.pdf", sep = ""), width = 10)
  par(mar = c(5,5,5,5))
  
  plot(x,y, pch = the.pch, col = the.line.color, ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
  axis(1, at=ticks,cex.axis=1.5,labels=chm.to.analyze,tick=F)
  mtext("Chromosome", side=1,line=3, cex = 1.5)
  axis(2, at=c(0,0.25, 0.5, 0.75, 1),cex.axis=1.5,labels=c("0.0","0.25","0.5","0.75","1.0"),tick=F)
  mtext("Proportion Selected", side=2,line=3, cex = 1.5)
  #par(new = TRUE)
  #plot(x,z, type="h",col="darkgreen",xaxt="n",yaxt="n",xlab="",ylab="", lwd = 2, ylim = c(0,z.lim))
  #axis(4, cex.axis = 1.5)
  #mtext("-log base 10 P-values", side=4,line=3, cex = 1.5)
  abline(v = ticks.end, col = "lightgrey", lwd = 1.5, cex = 1.0)
  abline(v = the.QTN, col = color.vector.QTN, lwd = 2, cex = 1.0, lty = lty.vector.QTN)
  points(x,y, pch = the.pch, col = the.line.color)
  #legend("topleft", "QTN", col = "black", lty = 4,, lwd = 2, cex=1.0, bty = "n")
  box()
  
  dev.off()
  
  
  #setwd(paste(home.dir, simulation.results.dir, sep = ""))
}#End create.graph.3

######Begin function
create.positive.rate.graph <- function(results.set.1 = NA, results.set.2 = NA, label.1 = NA, label.2 = NA){
  setwd(paste(home.dir, summary.dir, sep = ""))  
  pdf(paste("True_Positive_Rate_", number.of.simulated.QTN,"_QTN.pdf", sep = ""), width = 14)
  count <- 1
  for(i in heritability.vector){
    if(i == 0) next;
    par(mar = c(5,5,4,3), lab = c(8,5,7))
    x <- seq(1:number.of.simulated.QTN)
    y <- results.set.1$the.true.positive.rate.matrix[count,]
    plot(y~x,type="l", cex = 1.5, lty = 1, lwd = 1.5, ylim=c(0,1), 
         xlim = c(1,(number.of.simulated.QTN)), col = "red",
         axes = FALSE, xlab = "", ylab = "")
    axis(1, at=seq(1,number.of.simulated.QTN),cex.axis=1.5,labels=seq(1,number.of.simulated.QTN),tick=F)
    mtext("QTL.Number", side=1,line=3, lwd = 1.5, cex = 1.5)
    axis(2, at=c(0,0.25, 0.5, 0.75, 1),cex.axis=1.5,labels=c("0.0","0.25","0.5","0.75","1.0"),tick=F)
    mtext("Proportion Selected", side=2,line=3, lwd = 1.5, cex = 1.5)
    y <- results.set.2$the.true.positive.rate.matrix[count,]
    lines(y~x, type="l", cex = 1.5, col = "blue", lty = 2, lwd = 1.5)
    title(paste("Heritability = ", as.character(i), sep = ""))
    legend("topleft", c(paste(label.1, sep = ""), paste(label.2, sep = "")), col = c("red", "blue"), 
           lty = c(1,2), lwd = c(1.5,1.5), cex=1.0, bty = "n")
    box()
    count <- count + 1
  }#end for(i in heritability.vector)
  
  dev.off()
}#end  create.positive.rate.graph()

#######End function

create.false.positive.rate.graph <- function(results.set.1 = NA, results.set.2 = NA, label.1 = NA, label.2 = NA){
  setwd(paste(home.dir, summary.dir, sep = ""))  
  pdf(paste("False_Positive_Rate_", number.of.simulated.QTN,"_QTN.pdf", sep = ""), width = 14)
  par(mar = c(5,5,4,3), lab = c(8,5,7))
  x <- heritability.vector
  y <- as.vector(results.set.1$the.false.positive.rate)
  plot(y~x,type="l", cex = 1.5, lty = 1, lwd = 1.5, ylim=c(0,1), 
       xlim = c(heritability.vector[1],heritability.vector[length(heritability.vector)]), col = "red",
       axes = FALSE, xlab = "", ylab = "")
  axis(1, at=heritability.vector,cex.axis=1.5,labels=heritability.vector,tick=F)
  mtext("Heritability", side=1,line=3, lwd = 1.5, cex = 1.5)
  axis(2, at=c(0,0.25, 0.5, 0.75, 1),cex.axis=1.5,labels=c("0.0","0.25","0.5","0.75","1.0"),tick=F)
  mtext("Proportion Selected", side=2,line=3, lwd = 1.5, cex = 1.5)
  y <- as.vector(results.set.2$the.false.positive.rate)
  lines(y~x, type="l", cex = 1.5, col = "blue", lty = 2, lwd = 1.5)
  abline(h = 0.05, col = "grey", lwd = 1.0, lty = 3)
  title("False positive rate")
  legend("topleft", c(paste(label.1, sep = ""), paste(label.2, sep = "")), col = c("red", "blue"), 
         lty = c(1,2), lwd = c(1.5,1.5), cex=1.0, bty = "n")
  box()
  dev.off()
  
}#end create.false.positive.rate.graph()



##################################################################
##################################################################
############## Begin the function
create.summary.statistics <- function(the.results = the.results, the.simulated.additive.QTN  =  the.simulated.additive.QTN,
                                      the.simulated.epistatic.QTN = the.simulated.epistatic.QTN, setting.number = setting.number,
                                      include.additive.QTL = include.additive.QTL, include.epistatic.QTL = include.epistatic.QTL,
                                      flanking.bp.dist.from.QTN = flanking.bp.dist.from.QTN, number.of.traits.per.setting =number.of.traits.per.setting ){
  #Initialize some of the output files
  additive.QTN.signals.detected <- NULL
  epistatic.QTN.signals.detected <- NULL
  
  
  #Subdivide the results into i.) Additive and ii.) Epistatic signals
  the.additive.results <- the.results[-which(grepl("\\+", the.results[,2])),]
  #the.additive.results <- the.results
  the.add.chrom <- NULL
  the.add.pos <- NULL
  for(i in 1:nrow(the.additive.results)){
    the.add.chrom <- c(the.add.chrom, as.numeric(substr(unlist(strsplit(as.character(the.additive.results[i,2]),"_"))[1],4,5)) )
    the.add.pos <- c(the.add.pos, as.numeric(substr(unlist(strsplit(as.character(the.additive.results[i,2]),"_"))[2],1,25))) 
  } 
  
  
  the.epistatic.results <- the.results[which(grepl("\\+", the.results[,2])),]
  the.epi.chrom.1 <- NULL
  the.epi.chrom.2 <- NULL
  the.epi.pos.1 <- NULL
  the.epi.pos.2 <- NULL
  for(i in 1:nrow(the.epistatic.results)){
    this.two.chromosomes.and.positions <- unlist(strsplit(as.character(the.epistatic.results[i,2]),"\\+"))
    the.epi.chrom.1 <- c(the.epi.chrom.1, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[1]),"_"))[1],4,5)) )
    the.epi.chrom.2 <- c(the.epi.chrom.2, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[2]),"_"))[1],4,5)) )
    
    the.epi.pos.1 <- c(the.epi.pos.1, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[1]),"_"))[2],1,25)) )
    the.epi.pos.2 <- c(the.epi.pos.2, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[2]),"_"))[2],1,25)) )  
    
  } 
  
  
  if(include.additive.QTL){ 
    proportion.of.true.add.positives <- rep(0, nrow(the.simulated.additive.QTN))
    proportion.of.misspec.epi.signals <- rep(0, nrow(the.simulated.additive.QTN))
    #For loop through each additive QTN
    for(i in 1:nrow(the.simulated.additive.QTN)){
      
      #Parse out the additive markers that are within +/- 2 Mb of that QTN
      
      logical.statement.chr <- the.add.chrom == the.simulated.additive.QTN[i,3]
      logical.statement.leq.pos.plus.distance <- the.add.pos < (the.simulated.additive.QTN[i,4] + flanking.bp.dist.from.QTN)
      logical.statement.geq.pos.minus.distance <- the.add.pos > (the.simulated.additive.QTN[i,4] - flanking.bp.dist.from.QTN)
      
      the.add.results.for.this.add.QTN <-  the.additive.results[which(logical.statement.chr & logical.statement.leq.pos.plus.distance &  logical.statement.geq.pos.minus.distance)  , ]
      
      
      if(nrow(the.add.results.for.this.add.QTN) != 0){ 
        
        #Remove all additive additive signals corresponding to this QTN
        logical.statement.true.positive <- paste(the.additive.results[,1],the.additive.results[,2],sep = "") %in% paste(the.add.results.for.this.add.QTN[,1], the.add.results.for.this.add.QTN[,2], sep = "")
        the.additive.results <- the.additive.results[which(!logical.statement.true.positive),]
        
        the.add.chrom <- the.add.chrom[which(!logical.statement.true.positive)]
        the.add.pos <- the.add.pos[which(!logical.statement.true.positive)]
        
        # If collinear markers are present in one model, count it only once
        if(nrow(the.add.results.for.this.add.QTN) > 1){
          vector.of.collinear.markers <- rep(TRUE, nrow(the.add.results.for.this.add.QTN))
          for(j in 2:nrow(the.add.results.for.this.add.QTN)){
            if(the.add.results.for.this.add.QTN[j,1] == the.add.results.for.this.add.QTN[(j-1),1]) vector.of.collinear.markers[j] = FALSE
          }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
          the.add.results.for.this.add.QTN <- the.add.results.for.this.add.QTN[which(vector.of.collinear.markers),]
          rm(vector.of.collinear.markers)
        }#end if(nrow(the.add.results.for.this.add.QTN) > 1)
        #Record the proportion of such signals that are detected; this will
        # be the "true positive rate"
        proportion.of.true.add.positives[i] <- nrow(the.add.results.for.this.add.QTN)/number.of.traits.per.setting
        
      }#end if(nrow(the.add.results.for.this.add.QTN) != 0)
      
      
      #Copy the number of epistatic signals that have one at least one marker within the 
      # +/- 1 Mb region surrounding the QTN
      logical.statement.chr.1 <- (the.epi.chrom.1  == the.simulated.additive.QTN[i,3]) 
      logical.statement.chr.2 <- (the.epi.chrom.2  == the.simulated.additive.QTN[i,3]) 
      
      
      logical.statement.leq.pos.plus.distance.1 <- (the.epi.pos.1 < (the.simulated.additive.QTN[i,4] + flanking.bp.dist.from.QTN) ) 
      logical.statement.leq.pos.plus.distance.2 <- (the.epi.pos.2 < (the.simulated.additive.QTN[i,4] + flanking.bp.dist.from.QTN) )
      
      logical.statement.geq.pos.minus.distance.1 <- (the.epi.pos.1  > (the.simulated.additive.QTN[i,4] - flanking.bp.dist.from.QTN) ) 
      logical.statement.geq.pos.minus.distance.2 <- (the.epi.pos.2  > (the.simulated.additive.QTN[i,4] - flanking.bp.dist.from.QTN) )
      
      #Record the proportion of such signals; this will be the proportion of times
      # in which the given additive QTN is misspecified as an epistatic signal
      the.epi.results.for.this.add.QTN <-  the.epistatic.results[which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1)  |
                                                                         (logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  ), ]

      if(length(which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1)  |
                      (logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  )) >0){
          the.epistatic.results <-  the.epistatic.results[-which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1)  |
                                                                             (logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  ), ] 
       

      
      
          the.epi.chrom.1 <-  the.epi.chrom.1[-which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1)  |
                                                                   (logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  ) ]  
          the.epi.chrom.2 <-  the.epi.chrom.2[-which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1)  |
                                                       (logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  ) ] 
          
          the.epi.pos.1 <-  the.epi.pos.1[-which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1)  |
                                                       (logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  ) ]    
          the.epi.pos.2 <-  the.epi.pos.2[-which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1)  |
                                                   (logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  ) ]    
      
      }#end if(length(which((logical.statement.chr.1 & logical.statement.leq.pos.plus.distance.1 &  logical.statement.geq.pos.minus.distance.1) |(logical.statement.chr.2 & logical.statement.leq.pos.plus.distance.2 &  logical.statement.geq.pos.minus.distance.2)  )) >0)
      
      if(nrow(the.epi.results.for.this.add.QTN) != 0){ 
        
        # If collinear markers are present in one model, count it only once
        if(nrow(the.epi.results.for.this.add.QTN) > 1){
          vector.of.collinear.markers <- rep(TRUE, nrow(the.epi.results.for.this.add.QTN))
          for(j in 2:nrow(the.epi.results.for.this.add.QTN)){
            if(the.epi.results.for.this.add.QTN[j,1] == the.epi.results.for.this.add.QTN[(j-1),1]) vector.of.collinear.markers[j] = FALSE
          }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
          the.epi.results.for.this.add.QTN<- the.epi.results.for.this.add.QTN[which(vector.of.collinear.markers),]
          rm(vector.of.collinear.markers)
        }#end if(nrow(the.epi.results.for.this.add.QTN) > 1)
        #Record the proportion of such signals that are detected; this will
        # be the "true positive rate"
        proportion.of.misspec.epi.signals[i] <- nrow(the.epi.results.for.this.add.QTN)/number.of.traits.per.setting
        
        
      }#end if(nrow(the.add.results.for.this.add.QTN) != 0)
      
      
      
    }#End for(i in 1:nrow(the.simulated.additive.QTN))
    
    #Calculate the proportion of remaining additive markers; these are the false positive additive markers
    the.remaining.additive.results <- the.additive.results

    #Calculate the proportion of remaining epistatic markers; these are the epistatic markers that are not in the vicinity of any additive QTN
    the.epistatic.signals.not.vicinity.of.additive.QTN <- the.epistatic.results  
    
    
     if(nrow(the.additive.results) != 0){
      if(nrow(the.additive.results) > 1){
        vector.of.collinear.markers <- rep(TRUE, nrow(the.additive.results))
        for(j in 2:nrow(the.additive.results)){
          if(the.additive.results[j,1] == the.additive.results[(j-1),1]) vector.of.collinear.markers[j] = FALSE
        }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
        the.additive.results <- the.additive.results[which(vector.of.collinear.markers),]
        rm(vector.of.collinear.markers)
      }#end if(nrow(the.additive.results) > 1)
      proportion.of.remaining.additive.signals <- nrow(the.additive.results)/number.of.traits.per.setting
    }else{
      proportion.of.remaining.additive.signals <- NULL
    }
 
    
    
    #Regenerate the list of additive signals
    the.additive.results <- the.results[-which(grepl("\\+", the.results[,2])),]
    #the.additive.results <- the.results
    the.add.chrom <- NULL
    the.add.pos <- NULL
    for(i in 1:nrow(the.additive.results)){
      the.add.chrom <- c(the.add.chrom, as.numeric(substr(unlist(strsplit(as.character(the.additive.results[i,2]),"_"))[1],4,5)) )
      the.add.pos <- c(the.add.pos, as.numeric(substr(unlist(strsplit(as.character(the.additive.results[i,2]),"_"))[2],1,25))) 
    } 
    
    
    #Regenerate the epistatic signals
    the.epistatic.results <- the.results[which(grepl("\\+", the.results[,2])),]
    the.epi.chrom.1 <- NULL
    the.epi.chrom.2 <- NULL
    the.epi.pos.1 <- NULL
    the.epi.pos.2 <- NULL
    for(i in 1:nrow(the.epistatic.results)){
      this.two.chromosomes.and.positions <- unlist(strsplit(as.character(the.epistatic.results[i,2]),"\\+"))
      the.epi.chrom.1 <- c(the.epi.chrom.1, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[1]),"_"))[1],4,5)) )
      the.epi.chrom.2 <- c(the.epi.chrom.2, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[2]),"_"))[1],4,5)) )
      
      the.epi.pos.1 <- c(the.epi.pos.1, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[1]),"_"))[2],1,25)) )
      the.epi.pos.2 <- c(the.epi.pos.2, as.numeric(substr(unlist(strsplit(as.character(this.two.chromosomes.and.positions[2]),"_"))[2],1,25)) )  
      
    } 
    
    
    additive.QTN.signals.detected <- rbind(c("True_Positives", proportion.of.true.add.positives), c("Misspec_as_Epistatic", proportion.of.misspec.epi.signals))
    colnames(additive.QTN.signals.detected) <- c("Type", as.character(the.simulated.additive.QTN[,1]))
    
  }else{ 
    proportion.of.remaining.additive.signals <- NULL
    the.remaining.additive.results <- the.additive.results
    the.epistatic.signals.not.vicinity.of.additive.QTN <- the.epistatic.results
  }#end if(include.additive.QTL)
  
  
  
  if(include.epistatic.QTL){
    proportion.of.true.epi.positives.both.QTN <- NULL
    proportion.of.true.epi.positives.one.QTN <- NULL
    proportion.of.misspec.add.signals <- NULL
    
    
    #For loop through each epistatic QTN
    for(i in 1: nrow(the.simulated.epistatic.QTN)){ 
      #Skip over the even number QTNs; as each pairs of rows in the simulated.epistaic.QTN consitute one simulated
      # epistatic effect
      if(i%%2 == 0) next;
      
      #Parse out the marker pairs that are respectively within +/- 2 Mb of both QTN pairs
      # If collinear markers are present in the model, count it only once
      
      logical.statement.chr.Marker1.QTN1 <- the.epi.chrom.1 == the.simulated.epistatic.QTN[i,3]
      logical.statement.chr.Marker1.QTN2 <- the.epi.chrom.1 == the.simulated.epistatic.QTN[(i+1),3]
      
      logical.statement.leq.pos.plus.distance.Marker1.QTN1 <-  the.epi.pos.1 < (the.simulated.epistatic.QTN[i,4] + flanking.bp.dist.from.QTN)
      logical.statement.leq.pos.plus.distance.Marker1.QTN2 <-  the.epi.pos.1 < (the.simulated.epistatic.QTN[(i+1),4] + flanking.bp.dist.from.QTN)    
      
      logical.statement.geq.pos.minus.distance.Marker1.QTN1 <- the.epi.pos.1 > (the.simulated.epistatic.QTN[i,4] - flanking.bp.dist.from.QTN)
      logical.statement.geq.pos.minus.distance.Marker1.QTN2 <- the.epi.pos.1 > (the.simulated.epistatic.QTN[(i+1),4] - flanking.bp.dist.from.QTN)
      
      marker1.is.near.QTN1 <- (logical.statement.chr.Marker1.QTN1) & (logical.statement.leq.pos.plus.distance.Marker1.QTN1) & 
        (logical.statement.geq.pos.minus.distance.Marker1.QTN1)
      
      marker1.is.near.QTN2 <- (logical.statement.chr.Marker1.QTN2) & (logical.statement.leq.pos.plus.distance.Marker1.QTN2) & 
        (logical.statement.geq.pos.minus.distance.Marker1.QTN2)
      
      logical.statement.chr.Marker2.QTN1 <- the.epi.chrom.2 == the.simulated.epistatic.QTN[i,3]
      logical.statement.chr.Marker2.QTN2 <- the.epi.chrom.2 == the.simulated.epistatic.QTN[(i+1),3]
      
      logical.statement.leq.pos.plus.distance.Marker2.QTN1 <-  the.epi.pos.2 < (the.simulated.epistatic.QTN[i,4] + flanking.bp.dist.from.QTN)
      logical.statement.leq.pos.plus.distance.Marker2.QTN2 <-  the.epi.pos.2 < (the.simulated.epistatic.QTN[(i+1),4] + flanking.bp.dist.from.QTN)    
      
      logical.statement.geq.pos.minus.distance.Marker2.QTN1 <- the.epi.pos.2 > (the.simulated.epistatic.QTN[i,4] - flanking.bp.dist.from.QTN)
      logical.statement.geq.pos.minus.distance.Marker2.QTN2 <- the.epi.pos.2 > (the.simulated.epistatic.QTN[(i+1),4] - flanking.bp.dist.from.QTN)  
      
      
      marker2.is.near.QTN1 <- (logical.statement.chr.Marker2.QTN1) & (logical.statement.leq.pos.plus.distance.Marker2.QTN1) & 
        (logical.statement.geq.pos.minus.distance.Marker2.QTN1)
      
      marker2.is.near.QTN2 <- (logical.statement.chr.Marker2.QTN2) & (logical.statement.leq.pos.plus.distance.Marker2.QTN2) & 
        (logical.statement.geq.pos.minus.distance.Marker2.QTN2)
      
      #Record the proportion of such signals that are detected; this will be the "true positive rate"
      the.epi.results.for.both.of.these.epi.QTN <- the.epistatic.results[which(( (marker1.is.near.QTN1) & (marker2.is.near.QTN2))  | ((marker1.is.near.QTN2) & (marker2.is.near.QTN1))), ] 
      
      if(nrow(the.epi.results.for.both.of.these.epi.QTN) != 0){ 
        
        #Remove all epistatic signals corresponding to identifying both QTN
        logical.statement.true.positive.both.QTN <- paste(the.epistatic.results[,1],the.epistatic.results[,2],sep = "") %in% 
          paste(the.epi.results.for.both.of.these.epi.QTN[,1], the.epi.results.for.both.of.these.epi.QTN[,2], sep = "")
        
        
        the.epistatic.results <- the.epistatic.results[which(!logical.statement.true.positive.both.QTN),]
        
        the.epi.chrom.1 <- the.epi.chrom.1[which(!logical.statement.true.positive.both.QTN)]
        the.epi.chrom.2 <- the.epi.chrom.2[which(!logical.statement.true.positive.both.QTN)]
        the.epi.pos.1 <- the.epi.pos.1[which(!logical.statement.true.positive.both.QTN)]
        the.epi.pos.2 <- the.epi.pos.2[which(!logical.statement.true.positive.both.QTN)]
        
        
        # If collinear markers are present in one model, count it only once
        if(nrow(the.epi.results.for.both.of.these.epi.QTN) > 1){
          vector.of.collinear.markers <- rep(TRUE, nrow(the.epi.results.for.both.of.these.epi.QTN))
          for(j in 2:nrow(the.epi.results.for.both.of.these.epi.QTN)){
            if(the.epi.results.for.both.of.these.epi.QTN[j,1] == the.epi.results.for.both.of.these.epi.QTN[(j-1),1]) vector.of.collinear.markers[j] = FALSE
          }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
          the.epi.results.for.both.of.these.epi.QTN <- the.epi.results.for.both.of.these.epi.QTN[which(vector.of.collinear.markers),]
          rm(vector.of.collinear.markers)
        }#end if(nrow(the.epi.results.for.both.of.these.epi.QTN) > 1)
        
        #Record the proportion of such signals that are detected; this will
        # be the "true positive rate"
        proportion.of.true.epi.positives.both.QTN <- c(proportion.of.true.epi.positives.both.QTN,
                                                       (nrow(the.epi.results.for.both.of.these.epi.QTN)/number.of.traits.per.setting) )
        
      }else{
        proportion.of.true.epi.positives.both.QTN <- c(proportion.of.true.epi.positives.both.QTN, 0)      
      }#end if(nrow(the.epi.results.for.both.of.these.epi.QTN) != 0)
      
      #Reintroduce the logical statments; these will need to get updated because we deleted some columns
      logical.statement.chr.Marker1.QTN1 <- the.epi.chrom.1 == the.simulated.epistatic.QTN[i,3]
      logical.statement.chr.Marker1.QTN2 <- the.epi.chrom.1 == the.simulated.epistatic.QTN[(i+1),3]
      
      logical.statement.leq.pos.plus.distance.Marker1.QTN1 <-  the.epi.pos.1 < (the.simulated.epistatic.QTN[i,4] + flanking.bp.dist.from.QTN)
      logical.statement.leq.pos.plus.distance.Marker1.QTN2 <-  the.epi.pos.1 < (the.simulated.epistatic.QTN[(i+1),4] + flanking.bp.dist.from.QTN)    
      
      logical.statement.geq.pos.minus.distance.Marker1.QTN1 <- the.epi.pos.1 > (the.simulated.epistatic.QTN[i,4] - flanking.bp.dist.from.QTN)
      logical.statement.geq.pos.minus.distance.Marker1.QTN2 <- the.epi.pos.1 > (the.simulated.epistatic.QTN[(i+1),4] - flanking.bp.dist.from.QTN)
      
      marker1.is.near.QTN1 <- (logical.statement.chr.Marker1.QTN1) & (logical.statement.leq.pos.plus.distance.Marker1.QTN1) & 
        (logical.statement.geq.pos.minus.distance.Marker1.QTN1)
      
      marker1.is.near.QTN2 <- (logical.statement.chr.Marker1.QTN2) & (logical.statement.leq.pos.plus.distance.Marker1.QTN2) & 
        (logical.statement.geq.pos.minus.distance.Marker1.QTN2)
      
      logical.statement.chr.Marker2.QTN1 <- the.epi.chrom.2 == the.simulated.epistatic.QTN[i,3]
      logical.statement.chr.Marker2.QTN2 <- the.epi.chrom.2 == the.simulated.epistatic.QTN[(i+1),3]
      
      logical.statement.leq.pos.plus.distance.Marker2.QTN1 <-  the.epi.pos.2 < (the.simulated.epistatic.QTN[i,4] + flanking.bp.dist.from.QTN)
      logical.statement.leq.pos.plus.distance.Marker2.QTN2 <-  the.epi.pos.2 < (the.simulated.epistatic.QTN[(i+1),4] + flanking.bp.dist.from.QTN)    
      
      logical.statement.geq.pos.minus.distance.Marker2.QTN1 <- the.epi.pos.2 > (the.simulated.epistatic.QTN[i,4] - flanking.bp.dist.from.QTN)
      logical.statement.geq.pos.minus.distance.Marker2.QTN2 <- the.epi.pos.2 > (the.simulated.epistatic.QTN[(i+1),4] - flanking.bp.dist.from.QTN)  
      
      
      marker2.is.near.QTN1 <- (logical.statement.chr.Marker2.QTN1) & (logical.statement.leq.pos.plus.distance.Marker2.QTN1) & 
        (logical.statement.geq.pos.minus.distance.Marker2.QTN1)
      
      marker2.is.near.QTN2 <- (logical.statement.chr.Marker2.QTN2) & (logical.statement.leq.pos.plus.distance.Marker2.QTN2) & 
        (logical.statement.geq.pos.minus.distance.Marker2.QTN2)
      
      
      
      
      #Parse out the marker pairs that are respectively within +/- 2 Mb of only one of the two QTN pairs
      the.epi.results.for.one.of.these.epi.QTN <- the.epistatic.results[which(((marker1.is.near.QTN1) & (!marker2.is.near.QTN2)) |                                                                           
                                                                                ((marker1.is.near.QTN2) & (!marker2.is.near.QTN1)) |                                                                              
                                                                                ((marker2.is.near.QTN1) & (!marker1.is.near.QTN2)) |                                                                              
                                                                                ((marker2.is.near.QTN2) & (!marker1.is.near.QTN1))), ] 
      
      #This following loop will look to see if there are any marker pairs in "the.epi.results.for.one.of.these.epi.QTN" are actually within +/- 2 Mb
      # of other pairs of epistatic QTN. This can happen when, for example, two QTNs have overlapping +/- 2 Mb windows.
      
      
      
      if(nrow(the.epi.results.for.one.of.these.epi.QTN) != 0){ 
        logical.statement.true.positive.one.epi.QTN <- paste(the.epistatic.results[,1],the.epistatic.results[,2],sep = "") %in% 
          paste(the.epi.results.for.one.of.these.epi.QTN[,1], the.epi.results.for.one.of.these.epi.QTN[,2], sep = "")
        
        the.epistatic.results <- the.epistatic.results[which(!logical.statement.true.positive.one.epi.QTN),]
        
        the.epi.chrom.1 <- the.epi.chrom.1[which(!logical.statement.true.positive.one.epi.QTN)]
        the.epi.chrom.2 <- the.epi.chrom.2[which(!logical.statement.true.positive.one.epi.QTN)]
        the.epi.pos.1 <- the.epi.pos.1[which(!logical.statement.true.positive.one.epi.QTN)]
        the.epi.pos.2 <- the.epi.pos.2[which(!logical.statement.true.positive.one.epi.QTN)]
        
        
        # If collinear markers are present in the model, count it only once
        vector.of.collinear.markers <- rep(TRUE, nrow(the.epi.results.for.one.of.these.epi.QTN))
        if(nrow(the.epi.results.for.one.of.these.epi.QTN) > 1){
          for(j in 2:nrow(the.epi.results.for.one.of.these.epi.QTN)){
            if(the.epi.results.for.one.of.these.epi.QTN[j,1] == the.epi.results.for.one.of.these.epi.QTN[(j-1),1]) vector.of.collinear.markers[j] = FALSE
          }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
          the.epi.results.for.one.of.these.epi.QTN <- the.epi.results.for.one.of.these.epi.QTN[which(vector.of.collinear.markers),]
          rm(vector.of.collinear.markers)
        }#end if(nrow(the.epi.results.for.one.of.these.epi.QTN)
        proportion.of.true.epi.positives.one.QTN <- c(proportion.of.true.epi.positives.one.QTN,
                                                      (nrow(the.epi.results.for.one.of.these.epi.QTN)/number.of.traits.per.setting) )
      }else{
        proportion.of.true.epi.positives.one.QTN <- c(proportion.of.true.epi.positives.one.QTN, 0)
      }#end if(nrow(the.epi.results.for.one.of.these.epi.QTN) != 0)
      
      
      #Copy the number of additive signals that are within +/- 2 Mb of the QTN pairs
      
      logical.statement.epi.chr.1 <- (the.add.chrom  == the.simulated.epistatic.QTN[i,3]) 
      logical.statement.epi.chr.2 <- (the.add.chrom  == the.simulated.epistatic.QTN[(i+1),3]) 
      
      logical.statement.epi.leq.pos.plus.distance.1 <- (the.add.pos < (the.simulated.epistatic.QTN[i,4] + flanking.bp.dist.from.QTN) ) 
      logical.statement.epi.leq.pos.plus.distance.2 <- (the.add.pos < (the.simulated.epistatic.QTN[(i+1),4] + flanking.bp.dist.from.QTN) )
      
      logical.statement.epi.geq.pos.minus.distance.1 <- (the.add.pos  > (the.simulated.epistatic.QTN[i,4] - flanking.bp.dist.from.QTN) ) 
      logical.statement.epi.geq.pos.minus.distance.2 <- (the.add.pos  > (the.simulated.epistatic.QTN[(i+1),4] - flanking.bp.dist.from.QTN) )
      
      
      the.add.results.for.these.epi.QTN <-  the.additive.results[which((logical.statement.epi.chr.1 & logical.statement.epi.leq.pos.plus.distance.1 &  logical.statement.epi.geq.pos.minus.distance.1)  |
                                                                         (logical.statement.epi.chr.2 & logical.statement.epi.leq.pos.plus.distance.2 &  logical.statement.epi.geq.pos.minus.distance.2)  ), ]
      
      if(nrow(the.add.results.for.these.epi.QTN) != 0){ 
        the.additive.results <- the.additive.results[-which((logical.statement.epi.chr.1 & logical.statement.epi.leq.pos.plus.distance.1 &  logical.statement.epi.geq.pos.minus.distance.1)  |
                                                             (logical.statement.epi.chr.2 & logical.statement.epi.leq.pos.plus.distance.2 &  logical.statement.epi.geq.pos.minus.distance.2)  ), ]
        the.add.chrom <- the.add.chrom[-which((logical.statement.epi.chr.1 & logical.statement.epi.leq.pos.plus.distance.1 &  logical.statement.epi.geq.pos.minus.distance.1)  |
                                                (logical.statement.epi.chr.2 & logical.statement.epi.leq.pos.plus.distance.2 &  logical.statement.epi.geq.pos.minus.distance.2)  )]
        the.add.pos <- the.add.pos[-which((logical.statement.epi.chr.1 & logical.statement.epi.leq.pos.plus.distance.1 &  logical.statement.epi.geq.pos.minus.distance.1)  |
                                            (logical.statement.epi.chr.2 & logical.statement.epi.leq.pos.plus.distance.2 &  logical.statement.epi.geq.pos.minus.distance.2)  )]
        
        
        
        
        # If collinear markers are present in one model, count it only once
        if(nrow(the.add.results.for.these.epi.QTN) > 1){
          vector.of.collinear.markers <- rep(TRUE, nrow(the.add.results.for.these.epi.QTN))
          for(j in 2:nrow(the.add.results.for.these.epi.QTN)){
            if(the.add.results.for.these.epi.QTN[j,1] == the.add.results.for.these.epi.QTN[(j-1),1]) vector.of.collinear.markers[j] = FALSE
          }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
          the.add.results.for.these.epi.QTN <- the.add.results.for.these.epi.QTN[which(vector.of.collinear.markers),]
          rm(vector.of.collinear.markers)
        }#end if(nrow(the.add.results.for.these.epi.QTN) > 1)
        #Record the proportion of such signals that are detected; this will
        # be the "true positive rate"
        proportion.of.misspec.add.signals <-c(proportion.of.misspec.add.signals, (nrow(the.add.results.for.these.epi.QTN)/number.of.traits.per.setting))
        
        
      }else{
        proportion.of.misspec.add.signals <-c(proportion.of.misspec.add.signals, 0)
      }#end if(nrow(the.add.results.for.these.epi.QTN) != 0)
      
      
    }#End for(i in 1: (nrow(the.simulated.epistatic.QTN)/2)) 
    
    remaining.additive.signals.not.misspec.epi <- the.additive.results
    the.remaining.epistatic.results <- the.epistatic.results
    
    epistatic.QTN.signals.detected <- rbind( c("True_Positives_Both_QTN", proportion.of.true.epi.positives.both.QTN),
                                             c("True_Positives_One_QTN", proportion.of.true.epi.positives.one.QTN),
                                             c("Misspec_as_Additive", proportion.of.misspec.add.signals))
    
    the.colnames.this.object <- "Type"
    for(j in 1:nrow(the.simulated.epistatic.QTN)){
      if(j%%2 == 0) next;
      the.colnames.this.object <- c(the.colnames.this.object, paste(as.character(the.simulated.epistatic.QTN[j,1]), "_and_", 
                                                                    as.character(the.simulated.epistatic.QTN[(j+1),1]),sep = ""))
    }#end for(j in 1:nrow(the.simulated.epistatic.QTN))
    
    colnames(epistatic.QTN.signals.detected) <- the.colnames.this.object
    
  }else{
    remaining.additive.signals.not.misspec.epi <- the.additive.results
    the.remaining.epistatic.results <- the.epistatic.results
  } #end if(include.epistatic.QTL)
  
  #Obtain all markers in common between "the.remaining.additive.results" and "remaining.additive.signals.not.misspec.epi"
  if(include.additive.QTL | include.epistatic.QTL){
    the.remaining.additive.results <- data.frame(as.character(paste(the.remaining.additive.results[,1], the.remaining.additive.results[,2], sep = "")), the.remaining.additive.results)
    remaining.additive.signals.not.misspec.epi <- data.frame(as.character(paste(remaining.additive.signals.not.misspec.epi[,1], remaining.additive.signals.not.misspec.epi[,2], sep = "")), 
                                                             remaining.additive.signals.not.misspec.epi) 
    
    remaining.additive.signals.outside.of.any.QTN <- remaining.additive.signals.not.misspec.epi[remaining.additive.signals.not.misspec.epi[,1] %in% intersect(remaining.additive.signals.not.misspec.epi[,1], the.remaining.additive.results[,1]),]
    remaining.additive.signals.outside.of.any.QTN <- remaining.additive.signals.outside.of.any.QTN[,-1]
  

    if(nrow(remaining.additive.signals.outside.of.any.QTN) > 1){
        vector.of.collinear.markers <- rep(TRUE, nrow(remaining.additive.signals.outside.of.any.QTN))
        for(j in 2:nrow(remaining.additive.signals.outside.of.any.QTN)){
          if(remaining.additive.signals.outside.of.any.QTN[j,1] == remaining.additive.signals.outside.of.any.QTN[(j-1),1]) vector.of.collinear.markers[j] = FALSE
        }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
        remaining.additive.signals.outside.of.any.QTN <- remaining.additive.signals.outside.of.any.QTN[which(vector.of.collinear.markers),]
        rm(vector.of.collinear.markers)
      }#end if(nrow(the.additive.results) > 1)
      proportion.of.remaining.additive.signals <- nrow(remaining.additive.signals.outside.of.any.QTN)/number.of.traits.per.setting
      
      the.epistatic.signals.not.vicinity.of.additive.QTN <- data.frame(as.character(paste(the.epistatic.signals.not.vicinity.of.additive.QTN[,1], the.epistatic.signals.not.vicinity.of.additive.QTN[,2], sep = "")), 
                                                                       the.epistatic.signals.not.vicinity.of.additive.QTN)
      the.remaining.epistatic.results <- data.frame(as.character(paste(the.remaining.epistatic.results[,1], the.remaining.epistatic.results[,2], sep = "")), 
                                                    the.remaining.epistatic.results) 
      
      remaining.epistatic.signals.outside.of.any.QTN <- the.epistatic.signals.not.vicinity.of.additive.QTN[the.epistatic.signals.not.vicinity.of.additive.QTN[,1] %in% intersect(the.epistatic.signals.not.vicinity.of.additive.QTN[,1], the.remaining.epistatic.results[,1]),]
 
      
      remaining.epistatic.signals.outside.of.any.QTN <- remaining.epistatic.signals.outside.of.any.QTN[,-1]
      
    if(nrow(remaining.epistatic.signals.outside.of.any.QTN) > 1){
        vector.of.collinear.markers <- rep(TRUE, nrow(remaining.epistatic.signals.outside.of.any.QTN))
        for(j in 2:nrow(remaining.epistatic.signals.outside.of.any.QTN)){
          if(remaining.epistatic.signals.outside.of.any.QTN[j,1] == remaining.epistatic.signals.outside.of.any.QTN[(j-1),1]) vector.of.collinear.markers[j] = FALSE
        }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
        remaining.epistatic.signals.outside.of.any.QTN <- remaining.epistatic.signals.outside.of.any.QTN[which(vector.of.collinear.markers),]
        rm(vector.of.collinear.markers)
        proportion.of.remaining.epistatic.signals <- nrow(remaining.epistatic.signals.outside.of.any.QTN )/number.of.traits.per.setting
      }else{
      proportion.of.remaining.epistatic.signals <- nrow(remaining.epistatic.signals.outside.of.any.QTN)/number.of.traits.per.setting
     }#end if(nrow(remaining.additive.signals.outside.of.any.QTN) != 0)

    }else{
    #Calculate the proportion of remaining epistatic markers; these are the false positive epistatic markers
    if(nrow(the.epistatic.results) > 1){
      vector.of.collinear.markers <- rep(TRUE, nrow(the.epistatic.results))
      for(j in 2:nrow(the.epistatic.results)){
        if(the.epistatic.results[j,1] == the.epistatic.results[(j-1),1]) vector.of.collinear.markers[j] = FALSE
      }# end for(j in 1:nrow(the.add.results.for.this.add.QTN))){
      the.epistatic.results <- the.epistatic.results[which(vector.of.collinear.markers),]
      rm(vector.of.collinear.markers)
  }#end if(nrow(the.epistatic.results) > 1)

  
  proportion.of.remaining.epistatic.signals <- nrow(the.epistatic.results)/number.of.traits.per.setting
  }#end if(include.additive.QTL | include.epistatic.QTL)
  return(list(additive.QTN.signals.detected=additive.QTN.signals.detected, epistatic.QTN.signals.detected=epistatic.QTN.signals.detected,
              proportion.of.remaining.additive.signals=proportion.of.remaining.additive.signals, 
              proportion.of.remaining.epistatic.signals=proportion.of.remaining.epistatic.signals))
  
  
}#end create.summary.statistics


