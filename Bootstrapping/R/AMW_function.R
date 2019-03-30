#NOTES

#The script contains 2 functions

#The first function, "AMWBootstrap" is to carry bootstrapping using the AMW-BF method.
#The function requires 5 inputs (although 2 of those inputs are solely for testing purposes and will probably not be required in the future)

#The 5 inputs required for ODPBootstrap are "data","prior",dataDeter1","varying","deterministic"
#INPUTS:
#"data". This requires an input of triangle of claims
#"prior". Thi is the prior estimate for the BF method
#"dataDeter1". This requires a triangle of resampled Pearson residual. We use this solely for testing purposes to ensure the R model matches the R script
#"varying". This requires a 'TRUE' or 'FALSE' input. If set at 'TRUE', the model will perform the calculations for the varying ODP. Otherwise, it will perform the calculations for the constant ODP.
#"determinsitic". This requires a 'TRUE' or 'FALSE' input. When set at 'TRUE', this tells the model to use "dataDeter1" input. Again, this is solely for testing purposes.If set at FALSE, the proper resampling is done.   


#The second function is "AMW_AIB" is to carry the Actuary in the Box
#The function requires 7 inputs (although 3 of those inputs are solely for testing purposes)
#INPUTS:
#"data1". This requires an input of triangle of claims
#"data". This requires a fully developed bootstrapped cumulative claims. Essentially, the output of "ODPBootstrap" will the input for "ODP_AIB"
#"prior". Prior estimate for the BF method
#"dataDeter2".This requires a projected triangle of resampled Pearson Residuals. This is solely for testing purposes
#"dataDeter3". This is a 10 random variables between 0 and 1 used to fix the variability in the Pseudo Prior Estimate#
#"varying". This requires a 'TRUE' or 'FALSE' input. If set at 'TRUE', the model will perform the calculations for the varying ODP. Otherwise, it will perform the calculations for the constant ODP.
#"determinsitic". This requires a 'TRUE' or 'FALSE' input. When set at 'TRUE', this tells the model to use "dataDeter1" input. Again, this is solely for testing purposes.If set at FALSE, the proper resampling is done.   

#FOr the program to work, we need to download/load the Library ChainLadder

#The script is strutured as follows:

#Part 1 - Code for ODP Bootsrap
#Part 2 - Code for Actuary in the Box
#Part 3 - Code to run 1 simulation. The idea is to match the results in the excel
#Part 4 - Running the code for several simulations

#######################################################################################################




AMWBootstrap <- function(data=claims, prior=dataPrior, dataDeter1 = dataResample1, varying=FALSE, deterministic=TRUE){
  
  
  n <- nrow(data)
  
  #Start by cumulating data
  cum <- incr2cum(data)
  
  #Develop the triangle to Ultimate using the Mack Chain Ladder function and obtain the relevant development factors
  mack <- MackChainLadder(cum, est.sigma="Mack")
  f<- mack$f
  
  
  #Next, obain the cumulative development factors
  cumF <- matrix(data=NA, nrow=1, ncol=n)
  cumF[n] <- f[n]
  for (i in (n-1):1){
    cumF[i] <- cumF[i+1]*f[i]
  }
  
  
  #We can then calculate the reserve from the prior estimate
  reverse <- rev(cumF)
  Reserve <- matrix(data=NA, nrow=n, ncol=1)
  Reserve[1] <- 0
  for (k in 2:n){
    Reserve[k]<-prior[k,1]*(1-(1/reverse[k]))
  }
  
  #We can also calculate the Ultimate from the prior
  Ultimate <- matrix(data=NA, nrow=n, ncol=1)
  for(k in 1:n){
    Ultimate[k] <- cum[k,n+1-k]+Reserve[k,1]
  }
  
  #Next, we calculate the fitted cumulative triangle
  fcc <- matrix(data=NA, nrow=n, ncol=n)
  for(j in n:1){
    for (i in n:1){
      if(i+j>n+1){
        fcc[i,j]=NA
      }
      else if(i+j==n+1){
        fcc[i,j]=cum[i,j]
      }
      else{
        fcc[i,j]=fcc[i,j+1]/f[j]
      }
    }
  }
  
  #We then calculate the fitted incremental claims
  fic <- cum2incr(fcc)
  
  #Next, we calculate the pearson residuals
  pr <- (data-fic)/sqrt(fic)
  
  #number of elements
  noe <- length(data)^2-sum(is.na(data))
  
  #number of origin periods
  noo <- length(data)
  
  #number of parameters
  nop <- noo - 1
  
  #degrees of freedom
  dof <- noe - noo - nop
  
  #pearson scale paramter
  psp <- sum(pr^2,na.rm=TRUE)/dof
  
  #Calculate pspAdj. This is required for the varying version
  sumProduct <- matrix(data=NA, nrow=1,ncol=n)
  pspAdj <- matrix(data=NA, nrow=1,ncol=n)
  for(i in 1:ncol(cum)){
    sumProduct[i] <- sum(pr[,i]^2,na.rm=TRUE)
    pspAdj[i] <- (sumProduct[i]/(nrow(cum)+1-i))*noe/dof
  }
  pspAdj[ncol(cum)] <- min(pspAdj[ncol(cum)-1],pspAdj[ncol(cum)-2])
  pspAdj
  
  
  #We then obtain the standardised pearson residuals for constant and varying which we arrange as a list...
  if(varying==TRUE){
    spr <- matrix(data=NA, nrow=nrow(cum),ncol=ncol(cum))
    for (j in 1:ncol(cum)){
      for (i in 1:nrow(cum)){
        spr[i,j] <- sqrt(noe/dof)*pr[i,j]/sqrt(pspAdj[j])  
      }
    }  
    spr[1,n]=NA
    spr[n,1]=NA
    spr <- as.vector(as.matrix(spr))
    spr <- spr[is.finite(spr)]
  }
  else{
    spr <- sqrt(noe/dof)*pr/sqrt(psp)
    spr[1,n]=NA
    spr[nrow(cum),1]=NA
    spr <- as.vector(as.matrix(spr))
    spr <- spr[is.finite(spr)]  
  }
  
  #Calculate the resampled pearson residuals (we provide 2 versions)
  
  #The 1st version version applies if we are running the deterministic version
  #In this case, the program will run resample residuals as provided in the excel "deterministicResamplePearsonResidual_AMW_constant"
  #In the 2nd version, actual resampling is done and the model will be stochastic
  #Note that the 1st version is only for testing purposes
  
  k <- length(spr)
  
  rpr <- matrix(data=NA, nrow=n,ncol=n) 
  
  #1st version
  if(deterministic==TRUE){
    rpr = dataDeter1  
  }
  
  #2nd version
  else{
    rpr <- matrix(data=NA, nrow=n,ncol=n)
    for(j in 1:n){
      for(i in 1:row(cum)){
        rpr[i,j]<-is.na(cum[i,j])#changes cum matrix to true and false
        if(rpr[i,j]==TRUE){
          rpr[i,j]=0
        }
        else{
          rpr[i,j]=spr[min(round(runif(1)*k)+1,k)]
        }
      }
    }
   
  }
  
  #We calculate the Psuedo Incremental data for constant and varying....
  if(varying==TRUE){
    pic <- matrix(data=NA, ncol=ncol(cum),nrow = nrow(cum))
    for (j in 1:ncol(cum)){
      for (i in 1:nrow(cum)){
        pic[i,j] <- rpr[i,j]*sqrt(pspAdj[j]*fic[i,j])+fic[i,j]  
      }
    } 
  }
  else{
  pic <- rpr*sqrt(psp*fic)+fic
  }
  #and the Pseudo Cumulative Claims data
  pcc <- incr2cum(pic)
  
  pcc_CL <- MackChainLadder(pcc)
  full_pcc <- pcc_CL$FullTriangle
  
  return(full_pcc)
  
}

#################################################################################################################

AMW_AIB <- function(data1=claims, data=BootstrappedData, prior=dataPrior, dataDeter2=dataResample2, dataDeter3 = rand, varying=FALSE,deterministic=TRUE){


  n <- nrow(data1)
  
  #Start by cumulating data
  cum <- incr2cum(data1)
  
  #Develop the triangle to Ultimate using the Mack Chain Ladder function and obtain the relevant development factors
  mack <- MackChainLadder(cum, est.sigma="Mack")
  f<- mack$f
  
  
  #Next, obain the cumulative development factors
  cumF <- matrix(data=NA, nrow=1, ncol=n)
  cumF[n] <- f[n]
  for (i in (n-1):1){
    cumF[i] <- cumF[i+1]*f[i]
  }
  
  
  #We can then calculate the reserve from the prior estimate
  reverse <- rev(cumF)
  Reserve <- matrix(data=NA, nrow=n, ncol=1)
  Reserve[1] <- 0
  for (k in 2:n){
    Reserve[k]<-prior[k,1]*(1-(1/reverse[k]))
  }
  
  #We can also calculate the Ultimate from the prior
  Ultimate <- matrix(data=NA, nrow=n, ncol=1)
  for(k in 1:n){
    Ultimate[k] <- cum[k,n+1-k]+Reserve[k,1]
  }
  
  #Next, we calculate the fitted cumulative triangle
  fcc <- matrix(data=NA, nrow=n, ncol=n)
  for(j in n:1){
    for (i in n:1){
      if(i+j>n+1){
        fcc[i,j]=NA
      }
      else if(i+j==n+1){
        fcc[i,j]=cum[i,j]
      }
      else{
        fcc[i,j]=fcc[i,j+1]/f[j]
      }
    }
  }
  
  #We then calculate the fitted incremental claims
  fic <- cum2incr(fcc)
  
  #Next, we calculate the pearson residuals
  pr <- (data1-fic)/sqrt(fic)
  
  #number of elements
  noe <- length(data1)^2-sum(is.na(data1))
  
  #number of origin periods
  noo <- length(data1)
  
  #number of parameters
  nop <- noo - 1
  
  #degrees of freedom
  dof <- noe - noo - nop
  
  #pearson scale paramter
  psp <- sum(pr^2,na.rm=TRUE)/dof
  
  
  #Calculate pspAdj. This is required for the varying version
  sumProduct <- matrix(data=NA, nrow=1,ncol=n)
  pspAdj <- matrix(data=NA, nrow=1,ncol=n)
  for(i in 1:ncol(cum)){
    sumProduct[i] <- sum(pr[,i]^2,na.rm=TRUE)
    pspAdj[i] <- (sumProduct[i]/(nrow(cum)+1-i))*noe/dof
  }
  pspAdj[ncol(cum)] <- min(pspAdj[ncol(cum)-1],pspAdj[ncol(cum)-2])
  pspAdj
  
  
  #We then obtain the standardised pearson residuals for constant and varying which we arrange as a list...
  if(varying==TRUE){
    spr <- matrix(data=NA, nrow=nrow(cum),ncol=ncol(cum))
    for (j in 1:ncol(cum)){
      for (i in 1:nrow(cum)){
        spr[i,j] <- sqrt(noe/dof)*pr[i,j]/sqrt(pspAdj[j])  
      }
    }  
    spr[1,n]=NA
    spr[n,1]=NA
    spr <- as.vector(as.matrix(spr))
    spr <- spr[is.finite(spr)]
  }
  else{
    spr <- sqrt(noe/dof)*pr/sqrt(psp)
    spr[1,n]=NA
    spr[nrow(cum),1]=NA
    spr <- as.vector(as.matrix(spr))
    spr <- spr[is.finite(spr)]  
  }
  
  
  #Calculate the development factors of the Bootstrapped data....
  newf <- sapply(1:(n-1),function(i){
    sum(data[c(1:(n-i+1)),i+1])/sum(data[c(1:(n-i+1)),i])
  }
  )
 
  
  #...and calculate the cumulative development factors
  newcumF <- matrix(data=NA, nrow=1, ncol=n)
  newcumF[n] <- 1
  for (i in (n-1):1){
    newcumF[i] <- newcumF[i+1]*newf[i]
  }
  newcumF
  
  newquo <- matrix(data=NA, nrow=1, ncol=n)
  newquo[1] <- 1/newcumF[1]
  for (i in 2:n){
    newquo[i] = (1/newcumF[i]) - (1/newcumF[i-1])
  }
  
  
  
  
  rrp <- matrix(data=NA, nrow=n,ncol=n)
  if(deterministic==TRUE){
    rrp <- dataDeter2
  }
  
  else{
  #resampled for projection
  rrp[,1]=0
  for(j in 2:ncol(cum)){
    for(i in 1:nrow(cum)){
      rrp[i,j]<-is.na(cum[i,j])
      if(rrp[i,j]==FALSE){
        rrp[i,j]=0
      }
      else{
        rrp[i,j]=spr[min(round(runif(1)*k)+1,k)]
      }
    }
  }
  }
  
  
  #pseudo prior estimate
  ppe <- matrix(data=NA, nrow=n, ncol <- 1)
  if(deterministic==TRUE){
    for(i in 1:nrow(cum)){
      ppe[i] <-  qnorm(dataDeter3[i,1],prior[i,1],prior[i,1]*0.05) 
    } 
  }
    else{
      for(i in 1:nrow(cum)){
        ppe[i] <-  qnorm(runif(1),prior[i,1],prior[i,1]*0.05) 
      }  
    }
  
  
  #We then calculate the projected incremental claim
  pjic <- matrix (data=NA, nrow=n, ncol=n)
  if(varying==TRUE){
    pjic[,1]=0
    for(j in 2:ncol(cum)){
      for(i in 1:nrow(cum)){
        pjic[i,j]<-is.na(cum[i,j])
        if(pjic[i,j]==FALSE){
          pjic[i,j]=0
        }
        else{
          pjic[i,j]=ppe[i]*newquo[j]+rrp[i,j]*sqrt(abs(pspAdj[j]*newquo[j]*ppe[i]))
        }
      }
    }
  }
  else{
    pjic[,1]=0
    for(j in 2:n){
      for(i in 1:n){
        pjic[i,j]<-is.na(cum[i,j])
        if(pjic[i,j]==FALSE){
          pjic[i,j]=0
        }
        else{
          pjic[i,j]=ppe[i]*newquo[j]+rrp[i,j]*sqrt(abs(psp*newquo[j]*ppe[i]))
        }
      }
    }    
  }
  
  #We then calculate reserve calculations from the projected incremental claims
  reserve <- apply(pjic,1,sum)
  
  
  #Next, we do the one year projections
  oneYear <- matrix(data=NA, nrow=nrow(cum), ncol=ncol(cum))
  for(j in 1:ncol(cum)){
    for(i in 1:nrow(cum)){
      if(pjic[i,j]==0){
        oneYear[i,j] <- cum[i,j] 
      }
      else{
        oneYear[i,j] <- oneYear[i,j-1]+pjic[i,j]
      }
      
    }
  }
  
  for(j in 1:ncol(cum)){
    for(i in 1:nrow(cum)){
      if(i+j<=ncol(cum)+2){
        oneYear[i,j]=oneYear[i,j]
      }
      else{
        oneYear[i,j]=NA
      }
    } 
  }
  
  #and obtain the one year development factors
  oneYearf <- sapply(1:(n-1),function(i){
    sum(oneYear[c(1:(n-i+1)),i+1])/sum(oneYear[c(1:(n-i+1)),i])
  }
  )
  
  #We also obtain the cumulative one year development factors
  oneYearcumF <- matrix(data=NA, nrow=1, ncol=length(oneYearf))
  oneYearcumF[length(oneYearf)] <- oneYearf[length(oneYearf)]
  for (i in (length(oneYearf)-1):1){
    oneYearcumF[i] <- oneYearcumF[i+1]*oneYearf[i]
  }
  
  #We obtain the one-year reserves...
  reverse <- rev(oneYearcumF)
  oneYearReserve <- matrix(data=NA, nrow=nrow(cum), ncol=1)
  oneYearReserve[1] <- 0
  for (k in 2:n){
    oneYearReserve[k]<-ppe[k,1]*(1-(1/reverse[k-1]))
  }
  oneYearReserve
  
  #...and the full one year triangle
  fullOneYear <- matrix(data=NA, nrow=n, ncol=1)
  fullOneYear[1] <- oneYear[1,n]
  for(k in 2:n){
    fullOneYear[k] <- oneYear[k,n+2-k]+oneYearReserve[k]
  }
  fullOneYear
  
  #Finally, we calculate the CDR
  CDR <- Ultimate - fullOneYear
  
  return(oneYearReserve)
    
    }  

###############################################################################################################
#Part 3 - Code to run 1 simulation. The idea is to match the results in the excel

#These are where my files are located on my PC. These will need to be changed to the new location on user's computer for program to run. All files are provided in the folder attached

claims <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\DataBF.csv", header=FALSE)
dataPrior <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\dataPrior_BF.csv",header=FALSE)
dataResample1 = read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\DataRPR_BF.csv",header=FALSE)
dataResample2<- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\DataRRP_BF.csv",header=FALSE)
prior <- dataPrior
rand <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\randBF.csv",header=FALSE)

Boot <- AMWBootstrap(data=claims, prior=dataPrior, dataDeter1 = dataResample1, varying=FALSE, deterministic=TRUE)
AMW_AIB(data1=claims, data=Boot, prior=dataPrior, dataDeter2=dataResample2, dataDeter3 = rand, varying=FALSE,deterministic=TRUE)

Boot <- AMWBootstrap(data=claims, prior=dataPrior, dataDeter1 = dataResample1, varying=TRUE, deterministic=TRUE)
AMW_AIB(data1=claims, data=Boot, prior=dataPrior, dataDeter2=dataResample2, dataDeter3 = rand, varying=TRUE,deterministic=TRUE)

################################################################################################################
#Part 4 - Stochastic Simulations
sim<-10
var<-FALSE
deter<- FALSE

#These are where my files are located on my PC. These will need to be changed to the new location on user's computer for program to run. All files are provided in the folder attached
claims <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\DataBF.csv", header=FALSE)
dataPrior <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\dataPrior_BF.csv",header=FALSE)
dataResample1 = read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\DataRPR_BF.csv",header=FALSE)
dataResample2<- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\DataRRP_BF.csv",header=FALSE)
prior <- dataPrior
rand <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\randBF.csv",header=FALSE)

totalReserve <- matrix(data=NA, nrow=1,ncol=sim)
for(i in 1:sim){
  Boot <- ODPBootstrap(data=claims, dataDeter1 = dataResample1, varying=var, deterministic=deter)
  reserve <- AMW_AIB(data1=claims, data=Boot, prior=dataPrior, dataDeter2=dataResample2, dataDeter3 = rand, varying=var,deterministic=deter)

  totalReserve[1,i] <- sum(reserve)
}
totalReserve



