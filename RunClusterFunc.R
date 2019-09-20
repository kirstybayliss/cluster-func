################################################################################
#' This code is used for the analysis of spatio-temporal earthquake data as described in "Probabilistic identification of earthquake clusters using rescaled nearest #' neighbour distance networks", Kirsty Bayliss, Mark Naylor and Ian G Main, 2019. 
#' Geophysical Journal International, Volume 217, Issue 1, April 2019, Pages 487–503, https://doi.org/10.1093/gji/ggz034
#'
#' This includes: * code to identify space-time-magnitude nearest neighbours in the manner of Zaliapin and Ben-Zion, 2013 (https://doi.org/10.1002/jgrb.50179)
#'                * code to fit a Weibull mixture model to the nearest neighbour distances using MCMC
#'                * code to construct probabilistic cluster networks of individual earthquake clusters
#'
#' The individual functions are clearly documented at the beginning, things get messier near the end (sorry! I will clean this up someday!).
#'
#' Please get in touch with any problems/questions.
#' Kirsty.Bayliss@ed.ac.uk :)
##################################################################################


### Need gtools for MCMC/ClusterConstructor functions, the rest are for plotting
library(gtools)
library(igraph)
library(intergraph)
library(GGally)
library(ggnetwork)
library(ggplot2)

############################################# Functions for constructing nearest neighbour networks
#' Spherical distance between events
#'
#' Calculates spherical distance between two events using the cosine rule, earth radius 6371km.
#' Necessary for EarthquakeNearestNeighbour calculations.
#'
#' @param lat1 latitude of 1st event
#' @param lon1 longitude of 1st event
#' @param lat2 latitude of 2nd event
#' @param lon2 longitude of 2nd event
#'
#' @return Distance Spherical distance in km
#' @export


DistanceCalc <- function(lat1, lon1, lat2, lon2){
  la1 <- (lat1*pi)/180
  la2 <- (lat2*pi)/180
  lo1 <- (lon1*pi)/180
  lo2 <- (lon2*pi)/180
  lodiff <- abs(lo1-lo2)
  ladiff <- abs(la1-la2)
  if (ladiff==0.0 & lodiff == 0.0){
    Distance=0.0001
  }
  else{
    Ang <- acos((sin(la1)*sin(la2))+(cos(la1)*cos(la2)*cos(lodiff)))
    Distance <- Ang*6371
  }
  return(Distance)

}




#' Find earthquake nearest neighbours
#'
#' Identify a potential parent event for each event in the catalogue which is the event which is closest in a space-time magnitude sense.
#' Based on the nearest neighbour distance metric of Zaliapin and Ben-Zion (2013).
#' Takes in a dataframe which will need specifically named columns outlined below.
#'
#' @param EQcat  a dataframe containing columns "long", "lat" for event locations, "mag" for magnitude and a "datetime" column in POSIX format.
#' @return EQCat with added columns for event ID, nearest neighbour (NN), nearest neighbour distance (NND), rescaled (IED) and regular (Dist) distance between parent and child, rescaled (IET) and regular time (Time) in years and the magnitude of the identified parent event (parMag)
#' @export
#'


EarthquakeNearestNeighbour=function(EQcat){
  ##### function to calculate nearest neighbours in an earthquake catalogue. must include columns called "long", "lat" for location
  ### "mag" for magnitude and "datetime" which should be POSIXct.
  ### returns catalogue with columns for event ID,  nearest neighbour (NN), nearest neighbour distance (NND), rescaled distance(IED)
  ### and rescaled time (IET). UPDATED to also output non-recaled time and distance and parent event magnitude (Time, Dist and parMag respectively)
  EQcat$ID <- vector(mode="numeric", length=length(EQcat$lat))
  EQcat$NN <- vector(mode="numeric", length=length(EQcat$lat))
  EQcat$NND <- rep(1, length(EQcat$lat))
  EQcat$IED <- vector(mode="numeric", length=length(EQcat$lat))
  EQcat$Dist <- vector(mode="numeric", length=length(EQcat$lat))
  EQcat$IET <- vector(mode="numeric", length=length(EQcat$lat))
  EQcat$Time <- vector(mode="numeric", length=length(EQcat$lat))
  EQcat$parMag <- vector(mode="numeric", length=length(EQcat$lat))
  #EQcat$NND2 <- rep(1, length(EQcat$lat))

  EQcat$ID[1] <- 1
  for ( i in 2: length(EQcat$lat)){
    EQcat$ID[i] <- i
    prev = i-1
    #print(i)
    for (j in 1:prev){
      M <- EQcat$mag[j]
      TDiff<- difftime(as.POSIXct(EQcat$datetime[i], origin="1970-01-01"), as.POSIXct(EQcat$datetime[j], origin="1970-01-01"), units="days")
      DT <- as.numeric(TDiff)/365.25
      if (is.na(DT) == TRUE){
        print ("no time difference case1")
        print(i)
        print(j)
        DT=0.0000001
      }
      if (is.null(DT) == TRUE){
        print ("no time difference case2")
        print(i)
        print(j)
        DT=0.0000001
      }
      if (DT ==0) {
        print ("no time difference case3")
        print(i)
        print(j)
        DT=0.0000001
      }
      Dist <- DistanceCalc(EQcat$lat[i], EQcat$long[i], EQcat$lat[j], EQcat$long[j])
      if (EQcat$lat[i] == EQcat$lat[j] & EQcat$long[i]==EQcat$long[j]) {
        Dist = 0.0001
        print ("same location")
        print(i)
        print(j)

      }
      DR <- Dist**(1.6)
      ND <- DT*DR*(10**(-1*M))
      #print($NND[i])
      if (ND < EQcat$NND[i]){
        EQcat$NND[i] <- ND
        EQcat$NN[i] <- EQcat$ID[j]
        EQcat$IET[i] <- DT*(10**(-0.5*1*M))
        EQcat$Time[i] <- DT
        EQcat$Dist[i] <- Dist
        EQcat$IED[i] <- (Dist**1.6)*(10**(-0.5*1*M))
        EQcat$parMag[i] <- M
        #SCalMc$NND2[i] <- SCalMc$IET[i]*SCalMc$IED[i]

      }

    }
  }
  EQcat$NND[1] <- max(EQcat$NND)
  return(EQcat)
}

#' Weibull density function
#'
#' Alternative parameterisation of the Weibull distribution used in Weibmcmc() function using shape parameter a and scale parameter theta.
#'
#' R defines the Weibull distribution with shape parmeter a and scale parameter b as f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a)
#' However for ease of modelling we define the Weibull distribution  with shape parameter a and scale parameter theta as f(x) = (a*theta)*(x^(a-1))*exp(-theta*(x^a))
#' where theta = b^(-a) or b = theta^(-1/a)
#'
#' @param a shape parameter
#' @param theta scale parameter
#' @keywords Weibull
#' @export
#' @examples
#' dreweib(x, a, theta)

dreweib= function(x, a, theta){
  (a*theta)*(x**(a-1))*exp(-(x**a)*theta)
}



#' Fit Weibull mixture model with MCMC
#'
#' This function uses Gibbs sampling and a Metropolis step to fit a mixture of two Weibull distributions to the input data.
#' Uses a gamma conjugate prior for the scale parameter theta and a gamma prior for shape parameter a updated with Metropolis step of size ceta.
#' Will require list of component values for each component (ie 2 values for each hyperparameter)
#' Default parameter values have been found to work for several different catalogues but may need changing if NND is scaled differently.
#' For full details see Bayliss et al, 2019, "Probabilistic identification of earthquake clusters using rescaled nearest neighbour distance networks", GJI.
#' Modified from Gamma mixture model from Mohammadi et al, 2013. "Using mixture of gamma distributions for Bayesian analysis in an M/G/1 queue with optional second  #` service", Computational statistics.
#'
#' @import gtools
#' @param x data to fit mixture to
#' @param iter required number of MCMC iterations, adjust until confident of parameter convergence - default 100 000
#' @param burn number of iterations to remove from start of chain, used to remove events at beginning of chain which will have signifcant variation - default 50 000
#' @param alpha hyperparameter for gamma prior for shape parameter a - default c(5,5)
#' @param beta hyperparameter for gamma prior for shape parameter a - default c(5,5)
#' @param mu hyperparameter for gamma prior for scale parameter theta - default c(5,5)
#' @param nu hyperparameter for gamma prior for scale parameter theta - default c(0.01, 0.01)
#' @param thin set thinning of samples to reduce correlations in chain - default 2
#' @param ceta size of step for Metropolis jump - default c(0.0004, 0.0004)
#' @param a starting values for a paramter - default c(0.8, 0.4)
#' @param theta starting values for theta parameter - default c(100, 900)
#' @param pa prior for mixture component weightings - default c(0.5, 0.5)
#' @return MCMC chains in dataframe pa, theta and a
#' @keywords mixture
#' @export
#' @examples
#' output=weibmcmc(x1,iter, burn, alpha, beta, mu, nu, thin, ceta, a, theta, pa)
#'

weibmcmc=function(x,iter=100000,burn=50000, alpha=c(5,5), beta=c(5,5), mu=c(5,5), nu=c(0.01, 0.01), thin=2, ceta=c(0.0004, 0.0004), a=c(0.8, 0.4), theta=c(100, 900), pa=c(0.5, 0.5)){

  k=2
  aaa=0
  pasample=asample=thetasample=list()
  acceptrate=c(0,0)
  total=c(0,0)


  for(r in 1:iter){

    ## run mcmc algorithm for z, nn, pa, beta, alpha respectively.
    z=matrix(NA,k,length(x))
    for (i in 1:length(x)){

      gam= pa*dreweib(x[i], a, theta)
      check=rmultinom(1,size=1,prob=(gam)/sum(gam))
      if (dim(check)[1] != k) {print("Error in check. Iteration:"); print(i); print(dim(check))}
      z[,i]=check
    }
    nn=vector()
    anew=vector()
    #ll1=vector()
    for (i in 1:k){
      for (g in 1:thin){
        nn[i]=sum(z[i,])
        #Multinomial*x vals => x vals zero where x is not in distro i
        pr=(x*z[i,])
        # only keep elements of pr where x is non-zero
        xpr=pr[pr>0]


        # gamma conj prior for scale param theta
        theta[i]=rgamma(1,mu[i] + nn[i],nu[i]+ sum(xpr**a[i]))


        #Metropolis step, keeping anew positive valued
        step= rnorm(1, 0, ceta[i])
        anew[i]=a[i] + step

        if (anew[i]<0){
          repeat
          {
            anew[i]= a[i] + rnorm(1, 0, ceta[i])
            if (anew[i] >0)
              break
          }
        }


        # Sum of all x in distro i, for calculating acceptance ratio. Broken up in case of weird errors during run.
        sumthing=(sum(log(xpr)))
        prodbit=(beta[i] - sumthing)*(a[i] - anew[i])
        sumbit= theta[i]*(sum(xpr**a[i])- sum(xpr**anew[i]))

        logadiff = (log(anew[i]/a[i]))
        accept = exp(((nn[i]+alpha[i]-1)*logadiff)+( prodbit + sumbit))
        accept=accept+ (log(pnorm((anew[i]/ceta[i]), 0, ceta[i])) - log(pnorm((a[i]/ceta[i]), 0, ceta[i])))
        total[i]=total[i]+1
        if (is.na(accept)) {
          print("accept is NA:"); print(accept)
          a[i]=a[i]

        }

        else if(is.finite(accept) == TRUE){
          if (accept >= (runif(1))) {
            a[i] = anew[i]
            #print("Accepted new a")
            #print(a)
            #print(theta)
            if (r>burn) {acceptrate[i]=acceptrate[i]+1; total[i]=total[i]+1}

          }
          else {a[i]=a[i]}
        }
        else {a[i]= a[i]}
      }
    }

    if (r%%100 == 0){
      print(r)
    }

    pa=rdirichlet(1,1+nn)
    #print(pa)
    if (r>burn){
      aaa = aaa+1
      pasample[[(aaa)]]=pa
      thetasample[[aaa]]=theta
      asample[[aaa]]=a


    }
  }
  return(list(aaa=aaa,pa=pasample,theta=thetasample,a=asample, acceptrate=acceptrate, total=total))
}

#' Weibull cumulative distribution at given x given Weibull parameters
#'
#' Uses R stats pweibull function to calculate distribution function of x for given parameters, scaled by weighting in mixture. Used by MCMCFitsPlot.
#'
#' @param x point(s) of interest
#' @param a Weibull shape parameter
#' @param lambda Weibull scale parameter as in pweibull function. Use rscale to convert from preferred parameterisation
#' @param pa weighting for mixture model
#'
#' @return distribution function of x
#' @export

pweibmix <- function(x, a, lambda, pa){r <- pa*pweibull(x, a, lambda)}


#' Cumulative distribution for total mixture
#'
#' Combines cumulative distributions for both components to get total cumulative distribution for the mixture model
#'
#' @param x data
#' @param shape1 Weibull shape parameter of first component
#' @param  scale1 Weibull scale parameter (pweibull parameterisation - use rscale from preferred parameterisation) of first component
#' @param shape2 Weibull shape parameter of second mixture component
#' @param scale2 Weibull scale parameter of second mixture component
#' @param prob weighting for first component of mixture - second component calculated from 1-prob
#'
#' @return r total cumulative distribution function for the mixture at given x
#'
#' @export


cweibullmix <- function(x,shape1,scale1,shape2,scale2,prob) {
  r <- prob*pweibull(x,shape1,scale1)+(1-prob)*pweibull(x,shape2,scale2)
}

#' Convert between Weibull parameterisations
#'
#' Convert from standard r Weibull parameteristaion for scale to our preferred parametrisation.
#'
#' R defines the Weibull distribution with shape parmeter a and scale parameter b as f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a).
#' However for ease of modelling we define the Weibull distribution with shape parameter a and scale parameter theta as f(x) = (a*theta)*(x^(a-1))*exp(-theta*(x^a))
#' where theta = b^(-a) or b = theta^(-1/a)
#' This function converts from theta paremeterisation to b so that r pweibull functions can be used.
#'
#' See also dreweib function
#'
#' @param a shape parameter a of Weibull distribution
#' @param theta scale parameter theta
#'
#' @return b scale parameter b for use with R Weibull parameterisations
#'
#' @export
rscale= function (a, theta) {theta**(-1/a)}


#' Plots resulting mixture model from MCMC fits
#'
#' Plots histogram and then a sample of fits from the resulting MCMC mixture. Fits are sampled randomly from the chains.
#'
#'
#' @param x1 data to plot
#' @param a1s Weibull shape parameter outputs for background component
#' @param a2s Weibull shape parameter outputs for clustered component
#' @param theta1s Weibull scale parameter outputs for background component
#' @param theta2s Weibull scale parameter for clustered component
#' @param pa1s weighting output for background component of Weibull mixture
#' @param pa2s weighting output for clustered component
#' @param N number of random samples to plot
#' @param breaks number of bin breaks in histogram, default 100
#' @param joint Fit of complete mixture, default=FALSE
#'
#' @return plots graph
#'
#' @export


MCMCFitsPlot <- function(x1, a1s, a2s, theta1s, theta2s, pa1s, pa2s, N, breaks=100, xlim=c(-15, 0), ylim=c(0, 0.4), joint=FALSE){

  tmp=hist(log10(x1),freq=FALSE,breaks=breaks, ylim=ylim, main="", xlab="log10(NND)", xlim=xlim)
  binWidth = tmp$breaks[2]-tmp$breaks[1]
  par(new=TRUE)

  for (i in (1:N)){
    j <- runif(1, 1, length(a1s))
    #print(j)
    a1i=a1s[j]
    a2i=a2s[j]
    pa1i=pa1s[j]
    pa2i=pa2s[j]
    theta1i=theta1s[j]
    theta2i= theta2s[j]
    lambda1i=rscale(a1i, theta1i)
    lambda2i=rscale(a2i, theta2i)
    true <- cweibullmix(10^(tmp$breaks),a1i,lambda1i,a2i,lambda2i, pa1i)
    n=length(true)
    trueDiff <- true[2:n]-true[1:n-1]
    Normalisationi = sum(trueDiff*10^(tmp$breaks[2:n]-tmp$breaks[1:n-1]))

    jset1i <- pweibmix(10^(tmp$breaks), a1i, lambda1i, pa1i)
    n1i = length(jset1i)
    jsetdiff1i <- jset1i[2:n1i]-jset1i[1:n1i-1]
    points(tmp$mids, jsetdiff1i*10^(tmp$breaks[2:n1i]-tmp$breaks[1:n1i-1])/Normalisationi/binWidth, type="l", lwd=1, col=rgb(.6, .6, .6, .1), ylim=ylim, xlim=xlim)
    par(new=TRUE)

    jset2i <- pweibmix(10^(tmp$breaks), a2i, lambda2i, pa2i)
    n2i = length(jset2i)
    jsetdiff2i <- jset2i[2:n2i]-jset2i[1:n2i-1]
    points(tmp$mids, jsetdiff2i*10^(tmp$breaks[2:n2i]-tmp$breaks[1:n2i-1])/Normalisationi/binWidth, type="l",lwd=1, col=rgb(.6, .6, .6, .1), ylim=ylim, xlim=xlim)
    par(new=TRUE)

    if(joint == TRUE){
      JointSet <- jsetdiff1i+jsetdiff2i
      points(tmp$mids, JointSet*10^(tmp$breaks[2:n2i]-tmp$breaks[1:n2i-1])/Normalisationi/binWidth, type="l",lwd=1, col="black", ylim=ylim, xlim=xlim)
      par(new=TRUE)
    }
  }
}

#' Identify clusters from nearest neighbours
#'
#' Takes in MCMC mixture parameter chains and an earthquake catalogue (preferably output of EarthquakeNearestNeighbour) and builds clusters based on nearest neighbour events.
#' Calculates for each event a probability that it belongs in the clustered distribution using Weibull mixture fit and compares this value with a random variable to decide which
#' links to retain and which to cut. It then links sequentially backwards through the catalogue to create earthquake cluster trees.
#' Repeats this for N realisations to quantify link uncertainty.
#' Will save N dataframes of cluster realisations to specified location.
#' See paper for full details (Bayliss et al, 2019)
#'
#' @param a1s Weibull shape parameter outputs for background component
#' @param a2s Weibull shape parameter outputs for clustered component
#' @param theta1s Weibull scale parameter outputs for background component
#' @param theta2s Weibull scale parameter for clustered component
#' @param pa1s weighting output for background component of Weibull mixture
#' @param pa2s weighting output for clustered component
#' @param EQcat An earthquake catalogue, preferably output of EarthquakeNearestNeighbour() as requires most of these columns ("ID", "NN", "NND", "lat", "long", "mag", "datetime").
#' @param N number of realisations. Should be less than chain length to avoid disaster!
#' @param SaveLoc specify location to save output
#' @param CatName name for chosen catalogue
#'
#' @return CatNameClusters csv file containing information for each identified cluster, including cluster size, head node, average leaf depth and information on the mainshock.
#' @return CatName csv file containing earthquake catalogue for specified realisation. Will have added root (headnode) column and probability of clustering.
#'
#' @export


ClusterConstructor=function( a1, a2, theta1, theta2, pa1, pa2, EQcat, N, SaveLoc, CatName){

  for (k in 1: N)
  {
    print(k)
    x1 <- EQcat$NND
    print(length(x1))

    r = round(runif(1, 1, length(a1)))
    print(r)
    a1i=a1[r]
    a2i=a2[r]
    pa1i=pa1[r]
    pa2i=pa2[r]
    theta1i=theta1[r]
    theta2i= theta2[r]

   ## Some set-up
   #EQcat$PC <- NA
   ClusterSet <- as.data.frame(matrix(0, ncol = 10, nrow = 1))
   #colnames(ClusterSet) <- c("ldcount", "NLeaves", "HN", "ClustSize", "AvgLD", "HMag", "HLon", "HLat", "HTime", "MMag")
   len  <- length(EQcat$NND)

   PC <- vector(mode="numeric", length = length(x1))

   ## Calculate probability event is clustered given distribution parameters
   for (n in 1:len){
     PC[n] = (pa1i*dreweib(x1[n], a1i, theta1i))/((pa1i*(dreweib(x1[n], a1i, theta1i)))+((1-pa1i)*dreweib(x1[n], a2i, theta2i)))
   }
   #pr <- as.numeric(unlist(pr))


   EQcat$PC <- 1-PC
   print ("starting realisation")

   ## Identify if event is linked to nearest neighbour by comparing cluster probability to stochastic variable
   for (i in 1:length(EQcat$ID)){

     if (EQcat$PC[i] > runif(1,0,1)){

       EQcat$parent[i] = EQcat$NN[i]
     }
     else EQcat$parent[i]= NA
   }

   print("Parents calculated")
   EQcat$ld <- vector(mode="numeric", length=length(EQcat$ID))
   EQcat$root <- NA


   ## For each event, link parents back to root event
   for (i in 1:length(EQcat$parent)){
     p=EQcat$parent[i]
     if (is.na(p) == TRUE){
       EQcat$ld[i] == 0
     }
     else {
       while (is.na(p) == FALSE){
         if (p == 0){
           break
         }
         else {
           EQcat$root[i] <- p
           ploc=which(EQcat$ID == p)
           EQcat$ld[i]=EQcat$ld[i]+1
           p=EQcat$parent[ploc]
         }
       }
     }
   }

   EQcat$parent[is.na(EQcat$parent)] <- 0
   EQcat$NChildren <- vector(mode="numeric", length=length(EQcat$ID))

   EQcat$NChildren = 0

   for (i in 1:length(EQcat$ID)){
     EQcat$NChildren[i] = sum(EQcat$parent == EQcat$ID[i])
   }

   print(EQcat$ld)
   print(length(which(is.na(EQcat$parent))))
   EQcat$parent[EQcat$parent == 0] <- NA
   print("Calculated leaf depths")
   Root_occur <- data.frame(table(EQcat$root))
   plot(Root_occur$Freq, pch=16)
   EQcat$avgld <- vector(mode="numeric", length=length(EQcat$root))

   HeadNodes <- unique(EQcat$root)
   HeadNodes <- HeadNodes[-1]
   Cluster <- as.data.frame(matrix(0, ncol = 21, nrow = length(HeadNodes)))
   colnames(Cluster) <- c("ldcount", "NLeaves", "HN", "ClustSize", "AvgLD", "HMag", "HLon", "HLat", "HTime", "MMag", "longestBranch", "Mag2", "MMagTime", "M2Time", "Num24", "NBranches", "MSid", "NFS", "MSOffspring", "MeanNND", "MedianNND")

   print(length(HeadNodes))
   for (i in 1:length(HeadNodes)){
     HN <- HeadNodes[i]
     Cluster$HN[i] <- HN
     Cluster$longestBranch[i] <- 0
     clust <- subset(EQcat, EQcat$root == HN | EQcat$ID == HN)

     lC <- length(clust$mag)

  
     MMagidx <- which.max(clust$mag)
     Cluster$MSid[i] <- clust$ID[MMagidx]
     Cluster$MSOffspring[i] <- clust$NChildren[MMagidx]
     Cluster$MMagTime[i] <- clust$datetime[MMagidx]
     Cluster$MeanNND[i] <- mean(clust$NND)
     Cluster$MedianNND[i] <- median(clust$NND)
    #print(Cluster$MMagTime[i])
     PreMSClust <- subset(clust, clust$datetime < clust$datetime[MMagidx])
     Cluster$NFS[i] <- length(PreMSClust$datetime)
    #print(Cluster$NFS[i])
     ASClust <- subset(clust, clust$datetime > clust$datetime[MMagidx])
    #print(length(ASClust$datetime))
     if (length(ASClust$datetime) > 1){
       Cluster$Mag2[i] <- max(ASClust$mag)
      #print(Cluster$Mag2[i])
      ### Getting errors here - fix!
       M2idx <- which.max(ASClust$mag)
      #print(M2idx)
      #print(ASClust[M2idx])
       Cluster$M2Time[i] <- ASClust$datetime[M2idx]
     }
     else {
       #print("MS at end of chain")
       Cluster$Mag2[i] <- NA
       Cluster$M2Time[i] <- NA
     }
    #Cluster$Mag2[i] <- sort(clust$mag, partial=lC-1)[lC-1]


     Cluster$NBranches[i] <- sum(clust$NChildren >= 2)
     len <- length(clust$long)
     if ( len != 0){
       ## If there's more than one event
       for (j in 1:len){
         ## For each event in cluster
         Node=clust$ID[j]
         Check <- Node %in% EQcat$parent
         ## Check if event is a parent
         ## Want to identify leaves, which are not parents
         if (Check == FALSE){
           ## If event is a leaf, add it's leafdepth count to cluster total
           Cluster$ldcount[i]=Cluster$ldcount[i] + clust$ld[j]
           if (clust$ld[j] > Cluster$longestBranch[i]){
             Cluster$longestBranch[i] <- clust$ld[j]
           }
           if (Cluster$ldcount[i] != 0){
             ## If total ldcount is not zero, add one to leaf count. Not sure why I put this in as ldcount should never be 0...
             Cluster$NLeaves[i] = Cluster$NLeaves[i] +1
           }
         }
       }
     }
     else{
       Cluster$ldcount[i] = 0
       Cluster$NLeaves[i] = 0
     }


     EQcat$avgld[i] <- (Cluster$ldcount[i])/(Cluster$NLeaves[i])
     Cluster$AvgLD[i] <- (Cluster$ldcount[i])/(Cluster$NLeaves[i])
   }

   for (i in 1:length(Cluster$HN)){
     T <- Cluster$HN[i]
     clust <- subset(EQcat, EQcat$root == Cluster$HN[i] | EQcat$ID == Cluster$HN[i])
     Cluster$ClustSize[i] <- length(clust$long)
     Cluster$HMag[i] <- EQcat$mag[T]
     Cluster$HLon[i] <- EQcat$long[T]
     Cluster$HLat[i] <- EQcat$lat[T]
     Cluster$HTime[i] <- EQcat$datetime[T]
     if (length(clust$lat) != 0){
       Cluster$MMag[i] <- max(clust$mag)}
     else { Cluster$MMag[i] <- Cluster$HMag[i]}
   }

   print("clusters identified")
   myfile <- file.path(SaveLoc, paste0(CatName,"_", k, ".csv"))
   print(myfile)
   write.table(EQcat, file = myfile, sep = ",", row.names = FALSE, col.names = FALSE,
               quote = FALSE, append = FALSE)

   myfile <- file.path(SaveLoc, paste0(CatName, "Clusters_", k, ".csv"))
   write.table(Cluster, file = myfile, sep = ",", row.names = FALSE, col.names = FALSE,
               quote = FALSE, append = FALSE)

   k=k+1
  }
}
################################################################################
# Assumes InputCat is your data in R dataframe, with columns long, lat, mag and datetime (as POSIXct)
InputCat <- read_delim("InputCat.csv", delim=",")
colnames(InputCat) <- c("long", "lat", "mag", "datetime")
Cat <- EarthquakeNearestNeighbour(InputCat)

# Change this if you want to save all of the (many) output dataframes in some other folder!
wd <- getwd()

## Weibull mixture priors - I have found that these work for all datasets I have used so far, BUT if there is different/no rescaling and NNDs are not on range 10^-15 - 10^1 (or there abouts), you might need to change these so that you are working with a realistic prior for the given data (this took me ages to properly understand - essentially you need to be sure that your prior is broad enough that reasonable values of the parameter posteriors are within the prior range)
a=c(0.8, 0.4)
theta=c(100, 900)
## parameters for a 
alpha=c(5,5)
beta=c(5, 5)
## parameters for theta 
mu=c(5, 5)
nu=c(0.01, 0.01)

## Run MCMC fit
k=2; pa=c(0.5,0.5)
# When choosing ceta, remember this is for sampling a values
ceta=c(0.0004, 0.0004)
thin=2
## How long this will take depends on the sixe of x and the number of iterations. If x < 10000, 100000 iterations should take a few hours but your milage may vary
## For big datasets, this can be a lot longer.
iter=100000
burn=50000

output=weibmcmc(Cat$NND, iter, burn, alpha, beta, mu, nu, thin, ceta, a, theta, pa)

k1=output$k; allk1=output$allk; pa1=output$pa; a1=output$a; theta1=output$theta; pa=output$pa

aout = as.numeric(unlist(output$a))
thetaout = as.numeric(unlist(output$theta))
pa=as.numeric(unlist(output$pa))

## Seperate MCMC chains
a1s <- aout[seq(1, length(aout), 2)]
a2s <- aout[seq(2, length(aout), 2)]
theta1s <- thetaout[seq(1, length(thetaout), 2)]
theta2s <- thetaout[seq(2, length(thetaout), 2)]
pa1s <- pa[seq(1, length(aout), 2)]
pa2s <- pa[seq(2, length(aout), 2)]

# Use these to see how fits look visually
# MCMCFitsPlot(Cat$NND, a1s, a2s, theta1s, theta2s, pa1s, pa2s, 100, breaks=100)

# I recommend changing the CatName to something actually useful
ClusterConstructor(a1s, a2s, theta1s, theta2s, pa1s, pa2s, Cat, N=100, SaveLoc= wd, CatName = "Cat")


########## Things get super messy from here on...

#### Read back in all N realisations of cluster and catalogue data - I cheat and seperate outputs into two folders (Cluster for the cluster outputs and Out for the catalogue outputs). I should really add this into writing step earlier, but the user would have to set up the files manually anyway.
## Generally this works for my usual set-up of N=100 for catalogues of around 100 000 events. It will work with larger catalogues (up to 500 000) BUT will potentially crash R at some point, normally when trying to pick out an individual cluster. I'm working on a less messy way to do this!

setwd(paste0(wd,"/Cluster"))
temp <- list.files(pattern="*.csv")

for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i], header=FALSE, col.names= c("ldcount", "NLeaves", "HN", "ClustSize", "AvgLD", "HMag", "HLon", "HLat", "HTime", "MMag", "longestBranch", "Mag2", "MMagTime", "M2Time", "Num24", "NBranches", "MSid", "NFS","MSOffspring", "MeanNND", "MedianNND")))

setwd(paste0(wd,"/Out"))
temp2 <- list.files(pattern="*.csv")
for (i in 1:length(temp2)) assign(temp2[i], read.csv(temp2[i], header=FALSE , col.names=c("long", "lat", "Mag", "datetime", "ID", "NN", "NND", "IED", "Dist", "IET", "Time", "parMag", "PC", "parent", "ld", "root", "NChildren", "avgld")))

### There is almost definitely a neater/better way to do this, but I have not implemented it yet.
AvgLD <- list()
MMag <- list()
ClustSize <- list()
NLeaves <- list()
LDCount <- list()
HNLats <- list()
HNLons <- list()
HNTime <- list()
Mag2 <- list()
LB <- list()
HN <- list()
MMagT <- list()
M2T <- list()
Num24 <- list()
NBranches <- list()
MSid <- list()
NFS <- list()
NOffspring <- list()
MeanNND <- list()
MedianNND <- list()

for ( i in 1:length(temp2)){
  DF <- get(paste0(CatName, "Clusters_", i, ".csv"))
  AvgLD <- append(AvgLD, DF$AvgLD)
  MMag <- append(MMag, DF$MMag)
  ClustSize <- append(ClustSize, DF$ClustSize)
  NLeaves <- append(NLeaves, DF$NLeaves)
  LDCount <- append(LDCount, DF$ldcount)
  HNLats <- append(HNLats, DF$HLat)
  HNLons <- append(HNLons, DF$HLon)
  HNTime <- append(HNTime, DF$HTime)
  Mag2 <- append(Mag2, DF$Mag2)
  LB <- append(LB, DF$longestBranch)
  HN <- append(HN, DF$HN)
  MMagT <- append(MMagT, DF$MMagTime)
  M2T <- append(M2T, DF$M2Time)
  Num24 <- append(Num24, DF$Num24)
  NBranches <- append(NBranches, DF$NBranches)
  MSid <- append(MSid, DF$MSid)
  NFS <- append(NFS, DF$NFS)
  NOffspring <- append(NOffspring, DF$MSOffspring)
  MeanNND <- append(MeanNND, DF$MeanNND)
  MedianNND <- append(MedianNND, DF$MedianNND)
  
}

LDCount <- as.numeric(LDCount)
MMag <- as.numeric(unlist(MMag))
AvgLD <- as.numeric(unlist(AvgLD))
ClustSize <- as.numeric(unlist(ClustSize))
NLeaves <- as.numeric(unlist(NLeaves))
HNLons <- as.numeric(unlist(HNLons))
HNLats <- as.numeric(unlist(HNLats))
HNTime <- as.numeric(unlist(HNTime))
Mag2 <- as.numeric(unlist(Mag2))
LB <- as.numeric(unlist(LB))
HN <- as.numeric(unlist(HN))
MMagT <- as.numeric(MMagT)
M2T <- as.numeric(M2T)
Num24 <- as.numeric(Num24)
NBranches <- as.numeric(NBranches)
MSid <- as.numeric(MSid)
NFS <- as.numeric(unlist(NFS))
NOffspring <- as.numeric(unlist(NOffspring))
MedianNND <- as.numeric(unlist(MedianNND))
MeanNND <- as.numeric(unlist(MeanNND))

# Make dataframe LS with information on all clusters over 100 realisations - obviously this can be massive so beware!
LS <- data.frame(MMag, AvgLD, LDCount, ClustSize, NLeaves, HNLons, HNLats, HNTime, Mag2, LB, HN, MMagT, M2T, Num24, NBranches, MSid, NFS, NOffspring, MeanNND, MedianNND)

LDPal <- colorRampPalette(c("orange", "green", "black"))
LDCol <- LDPal(10)[as.numeric(cut(LS$AvgLD, breaks=10))]

### Basic plots to look at whole distribution
plot(LS$MMag , LS$AvgLD, pch=16, col=adjustcolor(LDCol, alpha=0.2), xlab="Mainshock magnitude", ylab = "Average leaf depth")
plot(LS$MMag , log10(LS$MedianNND), pch=16, col=adjustcolor(LDCol, alpha=0.2), xlab="Mainshock magnitude", ylab = "log10(Median NND)")
plot(LS$MMag , LS$NFS, pch=16, col=adjustcolor(LDCol, alpha=0.2), xlab="Mainshock magnitude", ylab = "number of foreshocks")

###########################################
### Code to look at individual clusters

### My favourite method for this is identifying clusters by their mainshock ID, which is the ID of the largest event in the cluster. 
## If you want a cluster containing an individual event that is NOT the mainshock this will involve a bit of detective work, and you made find it easier to check the
## root of the event in question over a few realisations to get the headnode (confusing terminology, but basically just the root event), and identify the largest 
## event in the cluster from there.
## Using mainshock is potentially problematic for more swarm-type clustering, but it was the best method I could come up with at the time.
MSID <- 8740
savename="MS8740"

Cluster1 <- subset(LS, LS$MSid == MSID)
hist(Cluster1$AvgLD, breaks=10, main=savename, xlab="Average Leaf Depths")

### See all possible headnodes (roots) of cluster
HN <- table(Cluster1$HN) 
HNfreq <- as.numeric(HN)
HNs <- as.numeric(names(HN))
hist(Cluster1$HN, main=savename, xlab="Headnode")

### I don't remember why I called this MaxLD, but basically it finds all events that should be in your cluster given the possible roots.
# This is important because the headnode/root is more likely to change over many realisations than the mainshock, for both mainshock-aftershock and swarms.

MaxLD <- list() 
for (i in HNs){
  rh <- i
  print(rh)
  for (i in 1:length(temp2)){
    DF <- get(paste0(CatName, "_", i, ".csv"))
    LDSet <- subset(DF, DF$root == rh)
    MaxLD <- rbind(MaxLD, LDSet)
  } 
  #DF <- get(paste0(CatName, "_1.csv"))
  #HNInf <- subset(DF,DF$ID == rh )
  #MaxLD <- rbind(MaxLD, HNInf)
}

### link certainty - how many times is a link kept?
## You'll need all those other libraries for this I'm afraid.

df <- data.frame(A=MaxLD$ID, B=MaxLD$NN)
df.g <- graph.data.frame(d=df, directed=TRUE)
E(df.g)$weight <- 0.01
gSimp <- simplify(df.g, edge.attr.comb = list(weight="sum"))

cert <- E(gSimp)$weight
G50 <- sum(cert > 50)
fractG50 <- G50/length(cert)
h1 <-hist(cert, breaks=100, xlab="certainty in link", main="")

#### Find details of all events in cluster
## Use MaxLD IDs of events in cluster to get event details
## Remove any old IDs2 kicking around.
rm(IDs2)

#IDs2 <- unique(MaxLD$ID)
IDs2 <- append(unique(MaxLD$ID), HNs[1])
for (i in HNs){
  j = which(IDs2 == i)
  print(j)
}
Len <- length(IDs2)

Mag <- vector(mode="numeric", length=Len)
Time <- vector(mode="numeric", length=Len)
Lat <- vector(mode="numeric", length=Len)
Lon <- vector(mode="numeric", length=Len)
LD <- vector(mode="numeric", length=Len)
# ClC can be used if there is some other/existing cluster allocation info
#ClC <- vector(mode="numeric", length=Len)
#Depth <- vector(mode="numeric", length=Len)
Weights <- vector(mode="numeric", length=Len)

MLDDat <- as.data.frame(table(MaxLD$ID))

for (i in 1:Len){
  k <- which(MaxLD$ID == IDs2[i])
  k1 <- k[1]
  print(k1)
  #Time[i] <- as.numeric(as.POSIXct(MaxLD$Time[k1]))
  Time[i] <- MaxLD$datetime[k1]
  Lat[i] <- MaxLD$lat[k1]
  Lon[i] <- MaxLD$lon[k1]
  Mag[i] <- MaxLD$Mag[k1]
  LD[i] <- MaxLD$ld[k1]
  WID <- which(MLDDat$Var1 == IDs2[i])
  Weights[i] <- MLDDat$Freq[WID]
  #ClC[i] <- MaxLD$ClAlloc[k1]
  #Depth[i] <- MaxLD$depth[k1]
}

CatSet <- get(paste0(CatName,"_1.csv"))
Lat[Len] <- CatSet$lat[which(CatSet$ID==HNs[1])]
Lon[Len] <- CatSet$lon[which(CatSet$ID==HNs[1])]
Time[Len] <- CatSet$datetime[which(CatSet$ID==HNs[1])]
Mag[Len] <- CatSet$Mag[which(CatSet$ID==HNs[1])] 
Weights[Len] <- E(gSimp)$weight[1]

LD2 <- LD[-Len] 
MagR <- Mag-1.5

##### Cumulative time plots
T2 <- sort(Time, decreasing=F)

MS <- which.max(Mag)
MST <- as.numeric(Time[MS])
Tdiff <- T2 - MST
TDiff <- Tdiff/(24*60*60)
T24 <- MST + 24*60*60
Set1Day <- subset(T2, T2 < T24)

NPal<- colorRampPalette(c('pink','red', 'gray', 'black'))

## Universal Pallete
NCol <- NPal(10)[as.numeric(cut(TDiff, breaks=c(-5000, -1,  0,  1, 2, 3, 10, 20, 50, 365.25, 10000 )) )]
events <- rep(1, length(TDiff))
par(xpd=FALSE)
plot(TDiff, cumsum(events), xlab="Days since Mainshock", ylab="Cumulative number of earthquakes", pch=16, col=NCol,  xlim=c(-5, 10))
lines(TDiff, cumsum(events))

### Plot black vertical line for events above mag 5
for (i in 1:length(IDs2)){
  #print(MaxLD$Mag[i])
  if (Mag[i] >= 6){
    print(Mag[i])
    print(i)
    TLoc <- Time[i]
    L1 <- which(T2 == TLoc)
    abline(v=TDiff[L1])
  }
}

### Plot red line at most productive event
t1 <- table(MaxLD$parent) 
mostChildren <- which.max(t1/100)
ID <- as.numeric(names(mostChildren))
TimeMC <- Time[which(IDs2==ID)[1]]
abline(v=TDiff[which(T2 == TimeMC)], col="red")

###### Network graph set-up
netg <- asNetwork(gSimp)
w1 <- as.numeric(E(gSimp)$weight)
lPal <- colorRampPalette(c('yellow', 'green', 'blue', 'black'))
ECol <- lPal(4)[as.numeric(cut(w1, breaks=c(-0.01, 0.25, 0.50, 0.75, 1.01)))]

MagClass <- vector(mode="numeric", length=length(Mag))
MagAlpha <- vector(mode="numeric", length=length(Mag))
#Weights <- c(E(gSimp)$weight[1], MLDDat$Freq)
WeightCol <- lPal(4)[as.numeric(cut(Weights, breaks=c(0, 0.25, 0.50, 0.75, 1)))]

for (i in 1:length(Mag)){
  #print(Mag[i])
  if (Mag[i] >= 6){
    MagClass[i] <- 8
    MagAlpha[i] <- 1
  }
  else if (Mag[i] >= 5){
    MagClass[i] <- 7
    MagAlpha[i] <- 0.7
  }
  else if (Mag[i] >= 4){
    MagClass[i] <- 17
    MagAlpha[i] <- 0.4
  }
  else if (Mag[i] >= 3){
    MagClass[i] <- 16
    MagAlpha[i] <- 0.1
  }
  else {MagClass[i] <- 15}
}

TDiff2 <- (Time-MST)/(24*60*60)
max(TDiff2)
which.max(TDiff2)
NPal<- colorRampPalette(c('pink','red', 'gray', 'black'))
NCol <- NPal(10)[as.numeric(cut(TDiff2, breaks=c(-13000, -1,  0,  1, 2, 3, 10, 20, 50, 365.25, 20000 )) )]

### Plot cluster network
ggnet2(netg,node.size=Mag-1, node.color=NCol, shape=MagClass, edge.label=LD2, edge.color=ECol)+guides(color = FALSE, size = FALSE)

######################## Cumulative events with symbols #############
FixedCumEvents <- as.data.frame(cbind(TDiff2, MagClass, MagR))
SortCumEvents <- FixedCumEvents[order(FixedCumEvents$TDiff2),]
events2 <- rep(1, length(TDiff2))
NPal<- colorRampPalette(c('pink','red', 'gray', 'black'))
NCol2 <- NPal(10)[as.numeric(cut(SortCumEvents$TDiff2, breaks=c(-5000, -1,  0,  1, 2, 3, 10, 20, 50, 100, 10000 )) )]
Size <- vector(mode = "numeric", length=length(SortCumEvents$MagR))
MC <- vector(mode = "numeric", length=length(SortCumEvents$MagR))
for(i in 1:length(SortCumEvents$MagR)){
  if (SortCumEvents$MagR[i] >= 4.5){
    MC[i] <- 8
    Size[i] <- 4
  }
  else if (SortCumEvents$MagR[i] >= 3.5){
    MC[i] <- 7
    Size[i] <- 3
  }
  else if (SortCumEvents$MagR[i] >= 2.5){
    MC[i] <- 17 
    Size[i]  <- 2
  }
  else 
  {MC[i] <- 16
  Size[i] <- 1}
}

plot(SortCumEvents$TDiff2, cumsum(events2), xlab="Days since Mainshock", ylab="Cumulative number of earthquakes", pch=MC, cex= Size, col=NCol2, xlim=c(-5, 10))
lines(SortCumEvents$TDiff2, cumsum(events2))
legend("topleft",legend=c( "M6+", "M5+", "M4+", "M3+"), pch=c(8, 7, 17, 16))

## Plot red vertical line for M7+ events
for (i in 1:length(SortCumEvents$TDiff2)){
  #print(MaxLD$Mag[i])
  if (SortCumEvents$MagR[i] >= 5.5){
    #print(Mag[i])
    #print(i)
    TLoc <- SortCumEvents$TDiff2[i]
    abline(v=TLoc, col="red")
  }
}

