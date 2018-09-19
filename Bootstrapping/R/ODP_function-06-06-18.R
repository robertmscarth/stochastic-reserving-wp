#NOTES

#The script contains 2 functions

#The first function, "ODPBootstrap" is to carry bootstrapping using the ODP method.
#The function requires 5 inputs (although 2 of those inputs are solely for testing purposes and will probably not be required in the future)

#The 4 inputs required for ODPBootstrap are "data","dataDeter1","varying","deterministic"
#INPUTS:
#"data". This requires an input of triangle of claims
#"dataDeter1". This requires a triangle of resampled Pearson residual. We use this solely for testing purposes to ensure the R model matches the R script
#"varying". This requires a 'TRUE' or 'FALSE' input. If set at 'TRUE', the model will perform the calculations for the varying ODP. Otherwise, it will perform the calculations for the constant ODP.
#"determinsitic". This requires a 'TRUE' or 'FALSE' input. When set at 'TRUE', this tells the model to use "dataDeter1" input. Again, this is solely for testing purposes.If set at FALSE, the proper resampling is done.   


#The second function is "ODP_AIB" is to carry the Actuary in the Box
#The function requires 5 inputs (although 2 of those inputs are solely for testing purposes)
#INPUTS:
#"data1". This requires an input of triangle of claims
#"data". This requires a fully developed bootstrapped cumulative claims. Essentially, the output of "ODPBootstrap" will the input for "ODP_AIB"
#"dataDeter2".This requires a projected triangle of resampled Pearson Residuals. This is solely for testing purposes
#"varying". This requires a 'TRUE' or 'FALSE' input. If set at 'TRUE', the model will perform the calculations for the varying ODP. Otherwise, it will perform the calculations for the constant ODP.
#"determinsitic". This requires a 'TRUE' or 'FALSE' input. When set at 'TRUE', this tells the model to use "dataDeter1" input. Again, this is solely for testing purposes.If set at FALSE, the proper resampling is done.   

#FOr the program to work, we need to download/load the Library ChainLadder

#The script is strutured as follows:

#Part 1 - Code for ODP Bootsrap
#Part 2 - Code for Actuary in the Box
#Part 3 - Code to run 1 simulation. The idea is to match the results in the excel
#Part 4 - Running the code for several simulations

######################################################################################################

#Part 1 - Code for ODP Bootsrap

ODPBootstrap <- function(data=claims, dataDeter1 = Mackdata1, varying=FALSE, deterministic=TRUE){
  
  
  n <- nrow(data)
  
  #Start by cumulating data
  cum <- incr2cum(data)
  
  #Develop the triangle to ultimate using the Basic Chain Ladder Method and calculates the development factors 
  mack <- MackChainLadder(cum, est.sigma="Mack")
  f<- mack$f
  
  #Calculate the fitted cumulative claims based on development factors f
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
  
  #Decumulate the fitted cumulative claims
  fic <- cum2incr(fcc)
  
  #Calculate the Pearson Residuals
  pr <- (data-fic)/sqrt(fic)
  
  #number of elements
  noe <- length(data)^2-sum(is.na(data))
  
  
  #number of origin periods & parameters
  if (varying==TRUE){
    noo <- length(data)
    nop <- noo -1
  }
  
  else{
    noo <- length(data)-1
    nop <- noo
  }
  
  
  #degrees of freedom
  dof <- noe - noo - nop
  
  #pearson scale paramter
  psp <- sum(pr^2,na.rm=TRUE)/dof
  
  #We also calculate the PSP for the varying case. This will be used only if varying is set to TRUE.
  sumProduct <- matrix(data=NA, nrow=1,ncol=n)
  pspAdj <- matrix(data=NA, nrow=1,ncol=n)
  
  for(i in 1:n){
    sumProduct[i] <- sum(pr[,i]^2,na.rm=TRUE)
    pspAdj[i] <- (sumProduct[i]/(n+1-i))*noe/dof
  }
  pspAdj[ncol(cum)] <- min(pspAdj[ncol(cum)-1],pspAdj[ncol(cum)-2])
  
  #Calculate standardised pearson residuals for the constant and varying case
  if (varying==TRUE){
    spr <- matrix(data=NA, nrow=nrow(cum),ncol=ncol(cum))
    for (j in 1:ncol(cum)){
      for (i in 1:nrow(cum)){
        spr[i,j] <- sqrt(noe/dof)*pr[i,j]/sqrt(pspAdj[j])  
      }
    }
  }
  
  else{
    spr <- sqrt(nop/dof)*pr/sqrt(psp)
    spr[1,n]=NA
    spr[n,1]=NA
    spr <- as.vector(as.matrix(spr))
    spr <- spr[is.finite(spr)]  
  }
  
  
  
  #Calculate the resampled pearson residuals (we provide 2 versions)
  
  #The 1st version version applies if we are running the deterministic version
  #In this case, the program will run resample residuals as provided in the excel "deterministicResamplePearsonResidual_ODP_constant"
  #In the 2nd version, actual resampling is done and the model will be stochastic
  #Note that the 1st version is only for testing purposes
   
  k <- length(spr)
  
  rpr <- matrix(data=NA, nrow=nrow(cum),ncol=ncol(cum))
  
  #1st version  
  if(deterministic==TRUE){
  rpr <- dataDeter1 #Mackdata1 is a fixed resampled pearson residual
  }
  
  #2nd version
  else{
  
  for(j in 1:n){
    for(i in 1:n){
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
  
  
  #Calculate the Psuedo Incremental Claims for the Constant and Varying Case
  if (varying==TRUE){
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
  
  
  #Calculate the Pseudo Cumulative Claims
  pcc <- incr2cum(pic)
  
  #Develops the Full PCC to ultimate for a full triangle
  pcc_CL <- MackChainLadder(pcc)
  full_pcc <- pcc_CL$FullTriangle
  
  return(full_pcc)
  }
  

######################################################################################################

#Part 2 - Code for Actuary in the Box

ODP_AIB <- function(data1=claims, data=BootstrappedData,dataDeter2=dataResample2,varying=FALSE,deterministic=TRUE){

  n <- nrow(data1)

  #Start by cumulating data
  cum <- incr2cum(data1)
  
  #Develop the triangle to ultimate using the Basic Chain Ladder Method and calculates the development factors 
  mack <- MackChainLadder(cum, est.sigma="Mack")
  f<- mack$f
  
  #Calculate the fitted cumulative claims based on development factors f
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
  
  #Decumulate the fitted cumulative claims
  fic <- cum2incr(fcc)
  
  #Calculate the Pearson Residuals
  pr <- (data1-fic)/sqrt(fic)  

  #number of elements
  noe <- length(data1)^2-sum(is.na(data1))
  
  
  #number of origin periods & parameters
  if (varying==TRUE){
    noo <- length(data1)
    nop <- noo -1
  }
  
  else{
    noo <- length(data1)-1
    nop <- noo
  }
  
  
  #degrees of freedom
  dof <- noe - noo - nop
  
  #pearson scale paramter
  psp <- sum(pr^2,na.rm=TRUE)/dof
  
  
  #We also calculate the PSP for the varying case. This will be used only if varying is set to TRUE.
  sumProduct <- matrix(data=NA, nrow=1,ncol=n)
  pspAdj <- matrix(data=NA, nrow=1,ncol=n)
  
  for(i in 1:n){
    sumProduct[i] <- sum(pr[,i]^2,na.rm=TRUE)
    pspAdj[i] <- (sumProduct[i]/(n+1-i))*noe/dof
  }
  pspAdj[n] <- min(pspAdj[n-1],pspAdj[n-2])
  
  #Calculate standardised pearson residuals for the constant and varying case
  if (varying==TRUE){
    spr <- matrix(data=NA, nrow=n,ncol=n)
    for (j in 1:n){
      for (i in 1:n){
        spr[i,j] <- sqrt(noe/dof)*pr[i,j]/sqrt(pspAdj[j])  
      }
    }
  }
  
  else{
    spr <- sqrt(nop/dof)*pr/sqrt(psp)
    spr[1,n]=NA
    spr[n,1]=NA
    spr <- as.vector(as.matrix(spr))
    spr <- spr[is.finite(spr)]  
  }
  
  #data <- Boot
  #data
  
  
  
  #Calculate the development factors of the Bootstrapped data
  newf <- sapply(1:(n-1),function(i){
    sum(data[c(1:(n-i+1)),i+1])/sum(data[c(1:(n-i+1)),i])
  }
  )
  
  
  
  #Calculate the cumulative development factors for the Bootstrapped data
  newcumF <- matrix(data=NA, nrow=1, ncol=n)
  newcumF[n] <- 1
  for (i in (n-1):1){
    newcumF[i] <- newcumF[i+1]*newf[i]
  }
  
 
  #Calculate 1/newCumF
  newquo <- matrix(data=NA, nrow=1, ncol=n)
  newquo[1] <- 1/newcumF[1]
  for (i in 2:n){
    newquo[i] = (1/newcumF[i]) - (1/newcumF[i-1])
  }
  
  #Again, 2 versions are given here. 1st version is provided for testing purpose
  #The 1st version applies if we are running the deterministic version
  #In this case, the program will run resample residuals as provided in the excel "deterministicResidualForProjection_ODP_constant"
  #In the 2nd version, actual resampling is done and the model will be stochastic
  #Note that the 1st version is only for testing purposes
  
  rrp <- matrix(data=NA, nrow=n,ncol=n)
  
  
  if(deterministic==TRUE){
    rrp <- dataDeter2
  }
  
  else{
    
    rrp[,1]=0
    for(j in 2:n){
      for(i in 1:n){
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

  #Calculate projected incremental claim for constant and varying
  if(varying==TRUE){
  pjic <- matrix (data=NA, nrow=n, ncol=n)
  pjic[,1]=0
  for(j in 2:n){
    for(i in 1:n){
      pjic[i,j]<-is.na(cum[i,j])
      if(pjic[i,j]==FALSE){
        pjic[i,j]=0
      }
      else{
        
          pjic[i,j]=data[i,n]*newquo[j]+rrp[i,j]*sqrt(abs(pspAdj[j]*newquo[j]*data[i,n]))
      }
    }
  }
  }
  
  else{
    pjic <- matrix (data=NA, nrow=n, ncol=n)
    pjic[,1]=0
    for(j in 2:n){
      for(i in 1:n){
        pjic[i,j]<-is.na(cum[i,j])
        if(pjic[i,j]==FALSE){
          pjic[i,j]=0
        }
        else{
          pjic[i,j]=data[i,n]*newquo[j]+rrp[i,j]*sqrt(abs(psp*newquo[j]*data[i,n]))
        }
      }
    } 
  }
  
  
  
  #reserve calculations
  reserve <- apply(pjic,1,sum)
  
  
  #Calculate one-year projections
  oneYear <- matrix(data=NA, nrow=n, ncol=n)
  for(j in 1:n){
    for(i in 1:n){
      if(pjic[i,j]==0){
        oneYear[i,j] <- cum[i,j] 
      }
      else{
        oneYear[i,j] <- oneYear[i,j-1]+pjic[i,j]
      }
      
    }
  }
  
  for(j in 1:n){
    for(i in 1:n){
      if(i+j<=n+2){
        oneYear[i,j]=oneYear[i,j]
      }
      else{
        oneYear[i,j]=NA
      }
    } 
  }
  
  
  #Calculate One-year development factors
  oneYearf <- sapply(1:(n-1),function(i){
    sum(oneYear[c(1:(n-i+1)),i+1])/sum(oneYear[c(1:(n-i+1)),i])
  }
  )
  
  
  #Calculate One year fully developed triangle
  fullOneYear <- matrix(data=NA, nrow=n, ncol=n)
  fullOneYear <- oneYear
  for(k in 3:n){
    fullOneYear[(n-k+3):n, k] <- fullOneYear[(n-k+3):n,k-1]*oneYearf[k-1]
  }
  
  #Calculate the one Year Ultimate
  oneYearUltimate <- fullOneYear[,n]
  
  
  #Calculate the oneYearReserve
  oneYearReserve <- matrix(data=NA, nrow=n, ncol=1)
  oneYearReserve[1] <- 0
  oneYearReserve[2] <- 0
  for (k in 3:n){
    oneYearReserve[k]<-oneYearUltimate[k]-fullOneYear[k,n-k+2]
  }
  
  
  mackUltimate <- mack$FullTriangle[,n]
  
  #Calculate the CDR
  CDR <- mackUltimate - oneYearUltimate
  
  return(oneYearReserve)
  
  }

#########################################################################################################

#Part 3 - Code to run 1 simulation. The idea is to match the results in the excel

#These are where my files are located on my PC. These will need to be changed to the new location on user's computer for program to run. All files are provided in the folder attached
dataCL <- read.csv("dataCL.csv", header=FALSE)
dataResample1 <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\deterministicResamplePearsonResidual_ODP_constant.csv",header=FALSE)
dataResample2 <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\deterministicResidualForProjection_ODP_constant.csv",header=FALSE)


deter <- TRUE # We want to match to the results from the excel, so deter is set to TRUE

#We first output the constant ODP, so varying has been set to FALSE
Boot <- ODPBootstrap(data=dataCL, dataDeter1 = dataResample1, varying=FALSE, deterministic=deter)
Boot
reser <- ODP_AIB(data1=dataCL, data=Boot, dataDeter2 = dataResample2, varying=FALSE, deterministic=deter)
reser

#Next, we output the varying ODP, so varying is now set to TRUE
Boot <- ODPBootstrap(data=dataCL, dataDeter1 = dataResample1, varying=TRUE, deterministic=deter)
Boot
reser <- ODP_AIB(data1=dataCL,data=Boot, dataDeter2 = dataResample2, varying=TRUE, deterministic=deter)
reser

################################################################################################################
#Part 4 - Stochastic Simulations
sim<-10
var<-FALSE
deter<- FALSE

#These are where my files are located on my PC. These will need to be changed to the new location on user's computer for program to run. All files are provided in the folder attached
dataCL <- read.csv("dataCL.csv", header=FALSE)
dataResample1 <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\deterministicResamplePearsonResidual_ODP_constant.csv",header=FALSE)
dataResample2 <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\deterministicResidualForProjection_ODP_constant.csv",header=FALSE)


totalReserve <- matrix(data=NA, nrow=1,ncol=sim)
for(i in 1:sim){
Boot <- ODPBootstrap(data=dataCL, dataDeter1 = dataResample1, varying=var, deterministic=deter)
reserve <- ODP_AIB(data1=dataCL,data=Boot, dataDeter2 = dataResample2, varying=var, deterministic=deter)
totalReserve[1,i] <- sum(reserve)
}
totalReserve
