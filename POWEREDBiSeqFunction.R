### Functions needed for TOOLNAME

## Dorothea Seiler Vellame
## ds420@exeter.ac.uk


## RRBSMatrixMaker ####################################
## Make a matrix of methylation and coverage values from RRBS data

##  INPUT: filePath - file path to cov files 
##         coordinatesOfUniqueCharactersInSampleName - the coordinates of the unique section of sample name
##                                                     so that column names can be shorter

## OUTPUT: RRBS matrix containing chromosome_startPosition, methylation 
##         value per sample and coverage value per sample. 

RRBSMatrixMaker<- function(filePath, coordinatesOfUniqueCharactersInSampleName = "all"){
  setwd(filePath)
  temp=list.files(pattern="*.cov")
  
  ## read data into a list
  if (coordinatesOfUniqueCharactersInSampleName[1] != "all"){
    sampleList<-list()
    for (i in temp){
      sampleList[[i]]<-assign(substr(i,coordinatesOfUniqueCharactersInSampleName[1],coordinatesOfUniqueCharactersInSampleName[2]),
                              as.data.frame(read.table(i,header=FALSE,sep="\t",
                                                       stringsAsFactors = FALSE,quote=""))) 
    }
    
    names<-substr(temp,coordinatesOfUniqueCharactersInSampleName[1],coordinatesOfUniqueCharactersInSampleName[2])
  }else{
    sampleList<-list()
    for (i in temp){
      sampleList[[i]]<-assign(i,
                              as.data.frame(read.table(i,header=FALSE,sep="\t",
                                                       stringsAsFactors = FALSE,quote=""))) 
    }
    
    names<-temp
  }
  numberOfSamples<-length(names)
  
  ## For each sample, format so that you have 3 cols: chr_startPosition, %m, cov
  # To use, loop through samples and names together so that they match
  RRBSColSelect<- function(RRBS, name){
    RRBStemp<-as.data.frame(matrix(ncol=3,nrow=nrow(RRBS)))
    RRBStemp[,1]<-paste(RRBS[,1],RRBS[,2],sep="_")
    RRBStemp[,2]<-RRBS[,4]
    RRBStemp[,3]<-RRBS[,5]+RRBS[,6]
    colnames(RRBStemp)<-c(paste("chr_start",name,sep="_"),paste("%m",name,sep="_"),
                          paste("cov",name,sep="_"))
    return(RRBStemp)
  }
  
  # loop through the sampleList applying ColSelect
  for(i in 1:numberOfSamples){
    sampleList[[i]]<-RRBSColSelect(sampleList[[i]],names[i])
  }
  
  ## create unique list of all sites
  # combine all site names into one column
  allSites<-c()
  for(i in 1:numberOfSamples){
    allSites<-c(allSites,sampleList[[i]][,1])
  }
  
  # keep only unique probes
  allSites<-unique(allSites)
  
  # create RRBSmatrix
  RRBSmatrix<-as.data.frame(matrix(ncol=(1+2*numberOfSamples),
                                   nrow=length(allSites)))
  RRBSmatrix[,1]<-allSites
  colnames(RRBSmatrix)[1]<-"chr_start"
  colnames(RRBSmatrix)[2:(numberOfSamples+1)]<-paste(names,"m",sep="_")
  colnames(RRBSmatrix)[(numberOfSamples+2):(2*numberOfSamples+1)]<-paste(names,"cov",sep="_")
  
  # match samples to probes
  for(i in 1:numberOfSamples){
    k<-match(sampleList[[i]][,1],RRBSmatrix[,1])
    RRBSmatrix[k,i+1]<-sampleList[[i]][,2]
    RRBSmatrix[k,i+1+numberOfSamples]<-sampleList[[i]][,3]
  }
  return(RRBSmatrix)
}



## parameterChecksFunction ############################
## Used within POWEREDBiSeq
## calculate negative binomial parameters of RRBS for data simulation

## INPUT: RRBS data coverage columns
##        assumptionCheck=TRUE means you check if var>mean holds (nbinomial assumption)

## OUTPUT: parameter values, p and r (two sets)

parameterChecksFunction = function(RRBSCov, assumptionCheck=FALSE){
  
  # convert 0s to NAs 
  convertZeroToNA = function(RRBSCol){
    RRBSCol[which(RRBSCol == 0)]<-NA
    return(RRBSCol)
  }
  
  RRBScov = apply(RRBSCov,2,convertZeroToNA)
  
  # calculate the mean and variance of the data and store in var
  var<-matrix(ncol=3,nrow=nrow(RRBSCov))
  colnames(var)<-c("sd","var","mean")
  # use the apply function
  var[,1]<-apply(RRBSCov,1,sd,na.rm=TRUE)
  var[,2]<-(var[,1])^2
  var[,3]<-apply(RRBSCov,1,mean,na.rm=TRUE)
  
  # create a, b, c for the quadratic formula
  medianVar<-median(var[,2],na.rm=TRUE)      # take the median values of var and mean
  medianMean<-median(var[,3],na.rm=TRUE)
  
  # hardcode quadratic solver 
  # a=var b=mu-2var c=var-mu (to calculate p)
  a=medianVar
  b=medianMean-2*medianVar
  c=medianVar-medianMean
  
  p1=(-b+sqrt(b^2-4*a*c))/(2*a)
  p2=(-b-sqrt(b^2-4*a*c))/(2*a)
  r1= medianMean*(1-p1)/p1
  r2= medianMean*(1-p2)/p2
  
  # do assumption check if assumptionCheck==TRUE
  if (assumptionCheck==TRUE){
    var<-cbind(var,matrix(ncol=1,nrow=nrow(RRBSCov)))
    colnames(var)[4]<-"var - mean"
    var[,4]<-var[,2]-var[,3]
    print(paste("Number of sites where var<mean = ",sum(var[,4]<0,na.rm=TRUE)))
  }
  return(c(p1,r1,p2,r2))
}



## POWEREDBiSeq #######################################
## A script to get the power at a given rd threshold

##  INPUT: rd - the read depth threshold to be used
##         RRBSMatrix made using the RRBSMatrix script
##         meanDiff - the difference in DNA methylation that you want to be able to detect between the two groups
##                    Between 0 - 1
##         pathToFunctions - the path to the functions folder 
##         pheno - factor vector of the two groups. If not input it will be assumed that the first half will be compared
##                 to the second (given that data is missing at random this should be representative)
##         nSampleNeeded - Number of samples that will be in each group as a minimum. Must be >2 for t-test to work

## OUTPUT: optimumRD and matrix showing:
##         rd - read depth at simulation 
##         power - power in simulation
##         bonferroniP - suggested bonferroni correction to use
##         nCpGTestedProp - the proportion of the data that will be kept given the read depth threshold
##         meanDiffPowerCheck - matrix showing the power you have detecting smaller meanDiffs given rdOptimum

POWEREDBiSeq = function(rd,
                        RRBSMatrix, 
                        meanDiff, 
                        pheno = FALSE, 
                        nSampleNeeded = 2){
  
  ## Check that scale of RRBS data and mean difference are proportions not percentages
  meanDiff = meanDiff*100
  nSamples = (ncol(RRBSMatrix)-1)/2 
  covCols = (nSamples+2):ncol(RRBSMatrix)
  
  ## set nPerm parameters for each subsetting (cannot be larger than number of rows in data)
  ifelse(nrow(RRBSMatrix) <40000,  
         optimalSearchNPerm = nrow(RRBSMatrix),  
         optimalSearchNPerm = 40000)
  ifelse(nrow(RRBSMatrix) <100000, 
         nPermPrior = nrow(RRBSMatrix),
         nPermPrior = 100000) ## must be bigger or equal to nPermR
  ifelse(nrow(RRBSMatrix) <60000, {
    nPermCpGRD = nrow(RRBSMatrix)
    nPermR = nrow(RRBSMatrix)
    nPermrdTrue = nrow(RRBSMatrix)
  },{
    nPermCpGRD = 60000
    nPermR = 60000
    nPermrdTrue = 60000
  })

  
  if(nSampleNeeded < 2){
    return("nSampleNeeded must be greater than 2")          
  }
  
  if(pheno[1] != FALSE){ 
    if (class(pheno) != "factor"){
      pheno = as.factor(pheno)
    }
    if (length(levels(pheno)) > 2){
      return("Pheno array given has too many levels - must have 2")
    }
  }
  
  # if n samples is odd nSamplesPerGroup1 and 2 will be different 
  if(pheno[1] != FALSE){ 
    nSamplesPerGroup1 = table(pheno)[1]
    nSamplesPerGroup2 = table(pheno)[2]
  }else{
    nSamplesPerGroup1 = floor(nSamples/2)
    nSamplesPerGroup2 = ceiling(nSamples/2)
    pheno = as.factor(c(rep("case", each = nSamplesPerGroup1),rep("control", each = nSamplesPerGroup2)))
  }
  
  if(nSampleNeeded > min(nSamplesPerGroup1, nSamplesPerGroup2)){
    return("nSamplesNeeded must be smaller than the sample size of each group")
  }
  
  
  ## functions ######
  # Each group must contain at least nSampleNeeded samples 
  filteredProbesCalc = function(covMatrixRow, rd){
    CpGRD1 = sum(covMatrixRow[pheno == levels(pheno)[1]] > rd, na.rm = TRUE)
    CpGRD2 = sum(covMatrixRow[pheno == levels(pheno)[2]] > rd, na.rm = TRUE)
    CpGRD = ifelse(CpGRD1 > nSampleNeeded & CpGRD2 > nSampleNeeded, T, F)
    return(CpGRD)
  }
  
  ## calculate proportions of DNAm in range
  priorCalc = function(RRBSm){
    
    if(is.null(ncol(RRBSm))){
      RRBSm = cbind(RRBSm, RRBSm)
    }
    
    sumCalc = function(RRBSmRow){
      g1 = sum(RRBSmRow < 5, na.rm = T)
      g3 = sum(RRBSmRow > 95, na.rm = T)
      tot = sum(!is.na(RRBSmRow))
      return(c(g1, g3, tot))
    }
    
    propDat = colSums(t(apply(RRBSm, 1, sumCalc)))
    
    priors = c(propDat[1]/propDat[3],
               (propDat[3]-(propDat[1]+propDat[2]))/propDat[3],
               propDat[2]/propDat[3])
    
    return(priors)
  }
  
  ## A function that returns a DNAm matrix filtered for rd
  rdFiltMatrix = function(nSampleInt, RRBSFull, nSamples, rd, cols = "meth"){
    RRBSCol = RRBSFull[,c(nSampleInt+1, nSampleInt+1+nSamples)]
    
    if(cols == "meth"){
      RRBSCol[which(RRBSCol[,2] <= rd), 1] = NA
      return(RRBSCol[,1])
    }
    RRBSCol[which(RRBSCol[,2] <= rd), 2] = NA
    return(RRBSCol[,2])
    
  }
  
  muDistrib = list(function(){runif(1, min = 0, max = 5)},
                   function(){runif(1, min = 5, max = 95)},
                   function(){runif(1, min = 95, max = 100)})
  
  getR = function(RRBSCov, nPermR){
    ## remove empty rows
    RRBSCov = RRBSCov[apply(RRBSCov,1,function(x){sum(is.na(x))})<125,]
    
    # get r parameter for coverage simulation (subset if remaining bigger than subset)
    ifelse(nrow(RRBSCov) > nPermR,
           parameterChecksFunction(RRBSCov[sample(1:nrow(RRBSCov),nPermR),])[4],
           parameterChecksFunction(RRBSCov)[4])
  }
  
  ## function to simulate the pvalues comparing groups
  simFunc = function(meanDiff, nSamplesPerGroup1, nSamplesPerGroup2, rd, r, nSampleNeeded, prior){
    mu1 = muDistrib[[sample(3, size = 1, prob = prior)]]()
    mu2 = mu1 - meanDiff
    mu2[mu1 - meanDiff < 0] = mu1[mu1 - meanDiff < 0] + meanDiff
    
    ## sample read depth from nbinomial distribution with given mean
    if(rd < rdTrue){ 
      rdSim = rdTrue
    }else{
      rdSim = rd
    }
    
    rd1 = rnbinom(nSampleNeeded, r, mu = rdSim)
    rd2 = rnbinom(nSampleNeeded, r, mu = rdSim)
    
    while(any(rd1<rd)){
      rd1[rd1<rd] = rnbinom(sum(rd1<rd), r, mu = rdSim) 
    }
    
    while(any(rd2<rd)){
      rd2[rd2<rd] = rnbinom(sum(rd2<rd), r, mu = rdSim) 
    }
    
    if(nSampleNeeded<nSamplesPerGroup1){
      rd1new = rnbinom(nSamplesPerGroup1 - nSampleNeeded, r, mu = rdSim)
      rd1new = rd1new[rd1new>rd]
      rd1 = c(rd1, rd1new)
    }
    
    if(nSampleNeeded<nSamplesPerGroup1){
      rd2new = rnbinom(nSamplesPerGroup1 - nSampleNeeded, r, mu = rdSim)
      rd2new = rd2new[rd2new>rd]   
      rd2 = c(rd2, rd2new)
    }
    
    ## use binomial distibution to sample from read depth number of methylated reads and calc DNAm value
    obs1<-rep(NA, length(rd1))
    for(i in 1:length(rd1)){
      obs1[i]<-rbinom(1, rd1[i], mu1/100)/rd1[i]
    }
    obs2<-rep(NA, length(rd2))
    for(i in 1:length(rd2)){
      obs2[i]<-rbinom(1, rd2[i], mu2/100)/rd2[i]
    }
    
    # t-test wont work if data are the same, assign p value larger than cut off rather than actually testing
    sim1T<-try(t.test(obs1, obs2)$p.value,silent=TRUE)
    if (class(sim1T)=="try-error"){
      sim1T<-1   # large so as not to be taken as significant
    }
    return(sim1T)
  }
  
  rForRD = getR(RRBSMatrix[sample(1:nrow(RRBSMatrix),nPermPrior), covCols], nPermR)
  
  rdTrue = mean(as.matrix(RRBSMatrix[sample(1:nrow(RRBSMatrix),nPermrdTrue),covCols]), na.rm = T)
  
  ## rdPowerFunction to get power results
  rdPower = function(rdIN,
                     meanDiffIN = meanDiff,
                     nSamplesPerGroup1IN = nSamplesPerGroup1,
                     nsamplesPerGroup2IN = nsamplesPerGroup2,
                     phenoIN = pheno,
                     nSamplesIN = nSamples,
                     RRBSMatrixIN = RRBSMatrix,
                     covColsIN = covCols,
                     nPermPriorIN = nPermPrior,
                     nPermCpGRDIN = nPermCpGRD,
                     nPermRIN = nPermR,
                     nPermrdTrueIN = nPermrdTrue,
                     optimalSearchNPermIN = optimalSearchNPerm,
                     nSampleNeededIN = nSampleNeeded,
                     rForRDIN = rForRD,
                     rdTrueIN = rdTrue){
    
    optimisationResults = matrix(nrow = 1, ncol = 5)
    colnames(optimisationResults) = c("rd", "power", "bonferroniP", "nCpGTestedProp", "nCpGTested")
    optimisationResults[1, 1] = rdIN
    
    ## calculate the proportion of CpGs that will be tested for the Bonferroni correction
    CpGRD = apply(RRBSMatrixIN[sample(1:nrow(RRBSMatrixIN),nPermCpGRDIN),covColsIN], 1, filteredProbesCalc, rdIN)
    optimisationResults[1, 4] = sum(CpGRD > 0)/nPermCpGRDIN
    optimisationResults[1, 5] = floor(optimisationResults[1, 4]*nrow(RRBSMatrixIN))
    optimisationResults[1, 3] = 0.05/optimisationResults[1, 5]
    
    ## calculate methylation priors
    priorForFilt = priorCalc(sapply(1:nSamplesIN, rdFiltMatrix, RRBSMatrixIN[sample(1:nrow(RRBSMatrixIN),nPermPriorIN),], nSamplesIN, rd = rdIN)) 
    
    
    ## simulate p values when comparing between groups
    sim1T = replicate(optimalSearchNPermIN, simFunc(meanDiffIN, nSamplesPerGroup1IN, nSamplesPerGroup2IN, rdTrue, 
                                                    rForRD, nSampleNeededIN, 
                                                    prior = priorForFilt))
    
    ## bootstrap sim1T so that length(sim1T) = optimisationResults[1, 5]
    if(length(sim1T) > optimisationResults[1, 5]){
      sim = sim1T[sample(1:optimisationResults[1, 5])]
    }else{
      sim = sim1T[sample(1:length(sim1T), optimisationResults[1, 5], replace = T)]
    }
    
    
    optimisationResults[1, 2] = sum(sim < optimisationResults[1, 3])/optimisationResults[1, 5]*100
    
    return(as.data.frame(optimisationResults))
  }
  
  powerOptim = data.frame(t(sapply(rep(rd,5), rdPower)))
  
  
  powerOut = data.frame(rd = rd, 
                        power = mean(unlist(powerOptim$power)), 
                        minPower = min(unlist(powerOptim$power)),
                        maxPower = max(unlist(powerOptim$power)),
                        bonferroniP = mean(unlist(powerOptim$bonferroniP)),
                        nCpGTestedProp = mean(unlist(powerOptim$nCpGTestedProp)),
                        nCpGTested = floor(mean(unlist(powerOptim$nCpGTested))))
  
  return(powerOut)
}
