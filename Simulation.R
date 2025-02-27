##############################################################
# Simulation at different settings (RCRXO)
##############################################################

# Define parameter values
ts <- c(4,5,6)  # Number of time periods
clusters <- c(10)  # Number of clusters
m <- c(15)  # Participants per cluster per period
icc <- c(0.01, 0.05, 0.1)  # Intra-cluster correlation
cacl <- c(1, 0.95, 0.75)  # Cluster autocorrelation
TimeEffsInd <- c(0)  # Time effects indicator
Treat_effect <- c(0, 0.2)  # Treatment effect

nrep <- 1000

# Generate all possible parameter combinations
simulation_params <- expand.grid(
  Time = ts, 
  Clusters = clusters, 
  Participants = m, 
  ICC = icc, 
  CAC = cacl, 
  TimeEffect = TimeEffsInd, 
  TreatmentEffect = Treat_effect
)

num_output_metrics <- 5  # Adjust this based on the actual number of outputs from ReRandCRXO
simulation_results <- matrix(NA, nrow = nrow(simulation_params), ncol = num_output_metrics) 

# Loop through each parameter combination and run the simulation
set.seed(3728)
for (i in 1:nrow(simulation_params)) {
  
  # Extract parameters from the current row of simulation_params
  params <- simulation_params[i, ]
  
  # Run the simulation function with extracted parameters
  simulation_results[i, ] <- ReRandCRXO(nrep, 
                                        params$Time, params$Clusters, params$Participants, 
                                        params$TimeEffect, params$TreatmentEffect,
                                        params$ICC, params$CAC)
}

# Convert results to a data frame and merge with parameter values
simulation_results_RCRXO <- cbind(simulation_params, as.data.frame(simulation_results))

# Assign column names for results
colnames(simulation_results_RCRXO)[(ncol(simulation_params) + 1):ncol(simulation_results_RCRXO)] <- 
  c("Bias", "ESE", "ASE", "MSE", "Alpha")

simulation_results_RCRXO

write.csv(simulation_results_RCRXO, file = "simulation_results_RCRXO.csv", row.names = FALSE)


##############################################################
# Simulation at different settings (CRXO)
##############################################################

# Define parameter values
ts <- c(4)  # Number of time periods
clusters <- c(10)  # Number of clusters
m <- c(15)  # Participants per cluster per period
icc <- c(0.01, 0.05, 0.1)  # Intra-cluster correlation
cacl <- c(1, 0.95, 0.75)  # Cluster autocorrelation
TimeEffsInd <- c(0)  # Time effects indicator
Treat_effect <- c(0, 0.2)  # Treatment effect

nrep <- 1000

# Generate all possible parameter combinations
simulation_params <- expand.grid(
  Time = ts, 
  Clusters = clusters, 
  Participants = m, 
  ICC = icc, 
  CAC = cacl, 
  TimeEffect = TimeEffsInd, 
  TreatmentEffect = Treat_effect
)

num_output_metrics <- 5  # Adjust this based on the actual number of outputs
simulation_results <- matrix(NA, nrow = nrow(simulation_params), ncol = num_output_metrics) 

# Loop through each parameter combination and run the simulation
set.seed(3728)
for (i in 1:nrow(simulation_params)) {
  
  # Extract parameters from the current row of simulation_params
  params <- simulation_params[i, ]
  
  # Run the simulation function with extracted parameters
  simulation_results[i, ] <- RandCRXO(nrep, 
                                        params$Time, params$Clusters, params$Participants, 
                                        params$TimeEffect, params$TreatmentEffect,
                                        params$ICC, params$CAC)
}

# Convert results to a data frame and merge with parameter values
simulation_results_CRXO <- cbind(simulation_params, as.data.frame(simulation_results))

# Assign column names for results
colnames(simulation_results_CRXO)[(ncol(simulation_params) + 1):ncol(simulation_results_RCRXO)] <- 
  c("Bias", "ESE", "ASE", "MSE", "Alpha")

simulation_results_CRXO

write.csv(simulation_results_CRXO, file = "simulation_results_CRXO.csv", row.names = FALSE)


##############################################################
# Simulation at different settings (MCRXO)
##############################################################

# Define parameter values
ts <- c(4, 6, 8)  # Number of time periods
clusters <- c(10)  # Number of clusters
m <- c(15)  # Participants per cluster per period
icc <- c(0.01, 0.05, 0.1)  # Intra-cluster correlation
cacl <- c(1, 0.95, 0.75)  # Cluster autocorrelation
TimeEffsInd <- c(0)  # Time effects indicator
Treat_effect <- c(0, 0.2)  # Treatment effect

nrep <- 1000

# Generate all possible parameter combinations
simulation_params <- expand.grid(
  Time = ts, 
  Clusters = clusters, 
  Participants = m, 
  ICC = icc, 
  CAC = cacl, 
  TimeEffect = TimeEffsInd, 
  TreatmentEffect = Treat_effect
)

num_output_metrics <- 5  # Adjust this based on the actual number of outputs
simulation_results <- matrix(NA, nrow = nrow(simulation_params), ncol = num_output_metrics) 

# Loop through each parameter combination and run the simulation
set.seed(3728)
for (i in 1:nrow(simulation_params)) {
  
  # Extract parameters from the current row of simulation_params
  params <- simulation_params[i, ]
  
  # Run the simulation function with extracted parameters
  simulation_results[i, ] <- MCRXO(nrep, 
                                      params$Time, params$Clusters, params$Participants, 
                                      params$TimeEffect, params$TreatmentEffect,
                                      params$ICC, params$CAC)
}

# Convert results to a data frame and merge with parameter values
simulation_results_MCRXO <- cbind(simulation_params, as.data.frame(simulation_results))

# Assign column names for results
colnames(simulation_results_MCRXO)[(ncol(simulation_params) + 1) : ncol(simulation_results_MCRXO)] <- c("Bias", "ESE", "ASE", "MSE", "Alpha")

simulation_results_MCRXO

write.csv(simulation_results_MCRXO, file = "simulation_results_MCRXO.csv", row.names = FALSE)


##############################################################
# Simulating average number of crossover in RCRXO
#############################################################

Ts <- 6  # for 4,5,6
clusters <- 10
n_sim <- 1000
total_shifts <- NULL

set.seed(123)
for (i in 1:n_sim) {
  design_matrix <- RCRXO_Desmat(Ts, clusters)
  # For each cluster (each row), count the shifts
  shifts_per_cluster <- apply(design_matrix, 1, count_shifts)
  # Store the total shift count for this simulation
  total_shifts[i] <- sum(count_shifts(design_matrix))
}

mean(total_shifts)/clusters

