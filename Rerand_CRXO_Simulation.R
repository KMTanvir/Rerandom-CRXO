######################################################################################################
# Re-randomized Cluster Cross Over Design
######################################################################################################

# Required library
library("lme4")

# Function to generate a re-randomized cluster crossover design matrix
RCRXO_Desmat <- function(Ts, clusters) {
  # Ts: Number of time periods
  # clusters: number of clusters
  
  # Initialize matrix
  CRXO <- matrix(data=0, ncol=Ts, nrow=clusters)  
  
  # Rerandomize treatment allocation
  for (i in 1:clusters) {
    CRXO[i, ] <- sample(rep(c(0,1), length.out=Ts))  
  }
  
  return(CRXO)
}


# Power Simulation for a Re-randomized Cluster Crossover (ReRand-CRXO) Trial
ReRandCRXO <- function(nrep, Ts, clusters, m, TimeEffsInd, Treat_effect, ICC, CAC){
  # nrep: number of replications
  # m: number of observations per cluster per period
  # TimeEffsInd: 0 = shared time effects, 1 = separate time effects
  # Treat_effect: treatment effect
  # ICC: intra-cluster correlation
  # CAC: cluster autocorrelation
  
  # Generate the re-randomized cluster crossover (CRXO) design matrix
  DesMatrix <- RCRXO_Desmat(Ts, clusters)
  
  
  # Time effects matrix creation
  if (TimeEffsInd == 0) {
    # Shared time effects across clusters
    TimeEffects <- 0.1 * matrix(data = seq(1:ncol(DesMatrix)), 
                                nrow = nrow(DesMatrix), ncol = ncol(DesMatrix), byrow = TRUE)
  } 
  else if (TimeEffsInd == 1) {
    # Unique time effects for each cluster
    TimeEffects <- matrix(data = rnorm(n = ncol(DesMatrix) * clusters, mean = 0, sd = 0.1), 
                          nrow = clusters, ncol = ncol(DesMatrix))
  }
  
  # Simulating model results
  output <- replicate(nrep, rerand_result_cont(DesMatrix, clusters, m, TimeEffects, Treat_effect, ICC, CAC))
  
  rerand_results <- NULL
  
  #Bias
  rerand_results[1] <- mean(output[1,]-Treat_effect)
  
  #Empirical SE
  rerand_results[2] <- sqrt(sum((output[1,]-mean(output[1,]))^2)/(nrep-1))
  
  #MSE
  rerand_results[3] <- (sum((output[1,]-Treat_effect)^2)/(nrep))
  
  #Average model SE
  rerand_results[4] <- sqrt(sum((output[2,])^2)/(nrep))
  
  #Rejection percentage
  rerand_results[5] <-  sum(abs(output[1,])/output[2,]>1.96)/nrep    # Hypothesis test
  
  return(rerand_results)
}

rerand_result_cont <- function(DesMatrix, clusters, m, Teffs, Treat_effect, ICC, CAC){
  # Teffs: Time effects matrix
  
  #Simulate the data:
  simulated_dataset <- Rerand_CRXO_dataset(DesMatrix, clusters, m, Teffs, Treat_effect, ICC, CAC)
  
  #Fit a model including shared time effects only
  fit_T <- lmer(Y ~ treat_vec + time_vec + (1|cluster_vec), simulated_dataset)
  
  return(c(fixef(fit_T)[2], sqrt(vcov(fit_T)[2,2])))
}

# A function to generate one dataset from a longitudinal CRT with design matrix given by DesMatrix
Rerand_CRXO_dataset <- function(DesMatrix, clusters, m, Teffs, Treat_effect, ICC, CAC){
  
  # Assume total variance of 1 (for simplified interpretation)
  # sigma_eps2: error variance
  # sigmaA2: variance of cluster random effect
  
  sigmaA2 = ICC  # All of ICC is assigned to between-cluster variance
  sigma_eps2 = 1 - sigmaA2  # Remaining variance goes to individual error
  
  # n_period = total number of periods
  n_period = ncol(Teffs)  # number of period
  
  # Generate a vector for the design and a vector for the time effects
  
  fulldesmat <- DesMatrix[rep(1:nrow(DesMatrix), 1), ] # no sorting here
  treat_vec <- rep(as.vector(t(fulldesmat)), each = m)
  
  fulltimemat <- Teffs[rep(1:nrow(Teffs), 1), ]
  Time_eff <- rep(as.vector(t(fulltimemat)), each = m)
  
  #Error terms:
  epsi = rnorm(clusters * n_period * m, mean=0, sd=sqrt(sigma_eps2)) # epsilon (error term)
  
  #Cluster random effects:
  clus_rand_eff <- rnorm(clusters,mean=0, sd=sqrt(sigmaA2))
  clus_rand_effect <- rep(clus_rand_eff, each=n_period*m) #one for each participant in each cluster
  
  clusterVi <- rep(seq(1 : clusters), each = n_period * m)
  cluster_vec = factor(clusterVi)
  
  timeVi <- rep(seq(1:n_period), each=m)
  timeVi <- rep(timeVi, times=clusters)
  time_vec = factor(timeVi)
  
  #Put everything together to get the outcomes:
  Y = Treat_effect + clus_rand_effect + Time_eff + epsi
  full_data = data.frame(Y, cluster_vec, time_vec, treat_vec)
  
  return(full_data)  
}

# Example of running the function
# ReRandCRXO(100, 6, 5, 10, 1, 0.15, 0.1, 0.95)

