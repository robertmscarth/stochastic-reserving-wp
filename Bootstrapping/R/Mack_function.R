#NOTES

#The script contains 2 functions

#The first function, "MackBootstrap" is to carry bootstrapping using the Mack method.
#The function requires 1 input which is the input of triangle of claims


#The second function is "Mack_AIB" is to carry the Actuary in the Box
#The function requires 2 inputs
#INPUTS:
#"data1". This requires an input of triangle of claims
#"data". This requires a fully developed bootstrapped cumulative claims. Essentially, the output of "ODPBootstrap" will the input for "ODP_AIB"

#FOr the program to work, we need to download/load the Library ChainLadder

#The script is strutured as follows:

#Part 1 - Code for ODP Bootsrap
#Part 2 - Code for Actuary in the Box
#Part 3 - Code to run simulations


######################################################################################################

#Part 1 - Code for Mack Bootsrap

MackBootstrap <- function(data=claims){
  
  
  n <- nrow(data)
  
  #Start by cumulating data
  cum <- incr2cum(data)
  
  
  mack <- MackChainLadder(cum, est.sigma="Mack")
  mack$FullTriangle
  
  df <- mack$f
  
  
  #Observed Development Factors
  odf <- matrix(data=NA, nrow=(n-1), ncol=(n-1))
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      odf[i,j] = cum[i,j+1]/cum[i,j]    
    }
  }
  
  
  #Pearson Residuals
  pr = matrix(data=NA, nrow=(n-1), ncol=(n-1))
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      pr[i,j] = sqrt(cum[i,j])*(odf[i,j]-df[j])
    }
  }
  
  
  #Square Sum
  prr = matrix(data=NA, nrow=(n-1), ncol=(n-1))
  ss = matrix(data=NA, nrow=1, ncol=(n-1))
  for(j in 1:(n-1)){
    for(i in 1:(n-1)){
      prr[i,j] = pr[i,j]^2
    }
    ss[j] = sum(prr[,j],na.rm=TRUE)
  }
  
  
  #Bias Correct
  bs <- 1
  
  #sum square divided by count
  #Bias Correct
  ssc=matrix(data=NA, nrow=1, ncol=(n-1))
  bsc=matrix(data=NA, nrow=1, ncol=(n-1))
  for (j in 1:(n-1)){
    ssc[j] = ss[j]/(n-j)
    if(j<(n-1)){
      bsc[j] = ssc[j]*(n-j)/(n-(j+1))  
    }
    else{
      bsc[j] = min(bsc[j-1],bsc[j-2])
    }
  }
  
  
  #standardised pearson residuals
  spr = matrix(data=NA, nrow=(n-1), ncol=(n-2))
  for(i in 1:(n-1)){
    for(j in 1:(n-2)){
      spr[i,j]=(sqrt((n-j)/(n-(j+1)))*pr[i,j])/sqrt(bsc[j])
    }
  }
  spr
  
  
  spr <- as.vector(as.matrix(spr))
  spr <- spr[is.finite(spr)]
  spr
  
  
  k <- length(spr)
  
  
  
  #resampled pearson residuals
  rpr <- matrix(data=NA, nrow=(n-1),ncol=(n-1))
  for(j in 1:(n-1)){
    for(i in 1:(n-1)){
      rpr[i,j]<-is.na(cum[i,j+1])#changes cum matrix to true and false
      if(rpr[i,j]==TRUE){
        rpr[i,j]=0
      }
      else{
        rpr[i,j]=spr[min(round(runif(1)*k)+1,k)]
      }
    }
  }
  
  
  
  #Pseudo Development Factors
  pdf <- matrix(data=NA, nrow=(n-1),ncol=(n-1))
  for(j in 1:(n-1)){
    for(i in 1:(n-1)){
      pdf[i,j]=mack$f[j]+rpr[i,j]*sqrt((bsc[j]/cum[i,j]))
      if(i+j>n){
        pdf[i,j]=NA
      }
    }
  }
  
  
  #To calculate the new cumulative triangle
  newCum <- matrix(data=NA, nrow=(n-1),ncol=(n-1))
  for(j in 1:(n-1)){
    for(i in 1:(n-1)){
      newCum[i,j]<-pdf[i,j]*cum[i,j]
    }
  }
  newCum=cbind(c(rep(NA,times=nrow(newCum))),newCum)
  newCum=rbind(newCum,c(rep(NA,times=ncol(newCum))))
  newCum
  
  return(newCum)
  
    }
  

######################################################################################################

#Part 2 - Code for Actuary in the Box

Mack_AIB <- function(data1=claims, data=BootstrappedData){

  n <- nrow(data1)

  #Cumulate the data
  cum <- incr2cum(data1)
  
  
  #Observed Development Factors
  odf <- matrix(data=NA, nrow=(n-1), ncol=(n-1))
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      odf[i,j] = cum[i,j+1]/cum[i,j]    
    }
  }
  
  
  #Pearson Residuals
  pr = matrix(data=NA, nrow=(n-1), ncol=(n-1))
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      pr[i,j] = sqrt(cum[i,j])*(odf[i,j]-df[j])
    }
  }
  
  
  #Square Sum
  prr = matrix(data=NA, nrow=(n-1), ncol=(n-1))
  ss = matrix(data=NA, nrow=1, ncol=(n-1))
  for(j in 1:(n-1)){
    for(i in 1:(n-1)){
      prr[i,j] = pr[i,j]^2
    }
    ss[j] = sum(prr[,j],na.rm=TRUE)
  }
  
  
  #Bias Correct
  bs <- 1
  
  
  #sum square divided by count
  #Bias Correct
  ssc=matrix(data=NA, nrow=1, ncol=(n-1))
  bsc=matrix(data=NA, nrow=1, ncol=(n-1))
  for (j in 1:(n-1)){
    ssc[j] = ss[j]/(n-j)
    if(j<(n-1)){
      bsc[j] = ssc[j]*(n-j)/(n-(j+1))  
    }
    else{
      bsc[j] = min(bsc[j-1],bsc[j-2])
    }
  }
  
  
  #standardised pearson residuals
  spr = matrix(data=NA, nrow=(n-1), ncol=(n-2))
  for(i in 1:(n-1)){
    for(j in 1:(n-2)){
      spr[i,j]=(sqrt((n-j)/(n-(j+1)))*pr[i,j])/sqrt(bsc[j])
    }
  }
  spr <- as.vector(as.matrix(spr))
  spr <- spr[is.finite(spr)]
  
  
  k <- length(spr)
  
  #To calculate the new development triangles
  newf <- sapply(1:(n-1),
                 function(i){
                   sum(data[c(1:(n-i)),i+1])/sum(cum[c(1:(n-i)),i])
                 }
  )
  
  
  #Resampled residuals for projection
  rrp <- matrix(data=NA, nrow=n,ncol=n)
  rrp[,1]=0
  for(j in 2:n){
    for(i in 1:n){
      rrp[i,j]<-is.na(data[i,j])
      if(rrp[i,j]==FALSE){
        rrp[i,j]=0
      }
      else{
        rrp[i,j]=spr[min(round(runif(1)*k)+1,k)]
      }
    }
  }
  
  #pseudo cumulative data 
  pcc <- matrix(data=NA, nrow=n,ncol=n)
  for(j in 1:n){
    for(i in 1:n){
      if(rrp[i,j]==0){
        pcc[i,j] <- cum[i,j] 
      }
      else{
        pcc[i,j] <- pcc[i,j-1]*newf[j-1]+sqrt(bsc[j-1]*pcc[i,j-1])*rrp[i,j]
      }
    }
  }
  for(j in 1:n){
    for(i in 1:n){
      if(i+j<=n){
        pcc[i,j]=NA
      }
      else{
        pcc[i,j]=pcc[i,j]
      }
    } 
  }
  
  
  
  #Calculate the ultimo reserves
  reserve <- matrix(data=NA, nrow=n , ncol=1)
  reserve[1]=0
  for(i in 2:n){
    reserve[i]=pcc[i,n]-pcc[i,n-i+1]
  }
  
  
  #Calculate the one year reserves
  oneYear <- matrix(data=NA, nrow=n, ncol=n)
  for(j in 1:n){
    for(i in 1:n){
      if(rrp[i,j]==0){
        oneYear[i,j] <- cum[i,j] 
      }
      else{
        oneYear[i,j] <- oneYear[i,j-1]*newf[j-1]+sqrt(bsc[j-1]*oneYear[i,j-1])*rrp[i,j]
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
  
  
  #Calculate the one developemnt factors
  oneYearf <- sapply(1:(n-1),function(i){
    sum(oneYear[c(1:(n-i+1)),i+1])/sum(oneYear[c(1:(n-i+1)),i])
  }
  )
  
  #Obtain the full one year triangle
  fullOneYear <- matrix(data=NA, nrow=n, ncol=n)
  fullOneYear <- oneYear
  for(k in 3:n){
    fullOneYear[(n-k+3):n, k] <- fullOneYear[(n-k+3):n,k-1]*oneYearf[k-1]
  }
  
  
  oneYearUltimate <- fullOneYear[,n]
  
  #Calculate the one year CDR
  oneYearReserve <- matrix(data=NA, nrow=n, ncol=1)
  oneYearReserve[1] <- 0
  oneYearReserve[2] <- 0
  for (k in 3:n){
    oneYearReserve[k]<-oneYearUltimate[k]-fullOneYear[k,n-k+2]
  }
  
  #Calculate the Mack Ultimate
  mackUltimate <- mack$FullTriangle[,n]
  
  #Calculate the CDR
  CDR <- mackUltimate - oneYearUltimate
  
  return(list(reserve=reserve,oneYearUltimate=oneYearUltimate,oneYearReserve=oneYearReserve,CDR=CDR))
  
  }

################################################################################################################
#Part 3 - Stochastic simulations using my Bootstrapping function
sim<-1000

#These are where my files are located on my PC. These will need to be changed to the new location on user's computer for program to run. All files are provided in the folder attached
claims <- read.csv("C:\\Users\\Aniketh\\Dropbox\\Stochastic Reserving\\Mack Method\\dataCL.csv", header=FALSE)



totalReserve <- matrix(data=NA, nrow=1,ncol=sim)
set.seed(12345)
for(i in 1:sim){
Boot <- MackBootstrap(data=dataCL)
AIB <- Mack_AIB(data1=dataCL,data=Boot)
totalReserve[i] <- sum(AIB$oneYearReserve)
}
totalReserve


#####################################################################################################################
#Statistics

avg <- mean(totalReserve)
std <- sd(totalReserve)
cv  <- std/avg
