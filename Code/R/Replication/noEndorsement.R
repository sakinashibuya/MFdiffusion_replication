rm(list = ls())

# Load packages

packages <- c(
  "R.matlab",
  "tidyverse",
  "igraph",
  "pracma",
  "ggplot2",
  "plotly",
  "dqrng",
  "Matrix"
)

not_installed <- !packages %in% installed.packages()
if (any(not_installed)) install.packages(packages[not_installed])
lapply(packages,require,character.only = TRUE)

# Global setting
user <- "Sakina"
if (user == "Hiro"){
  project_path <- "/Users/mizuhirosuzuki/Dropbox/MFdiffusion_replication/"
}
if (user == "Sakina"){
  project_path <- "/home/sakina/Github/MFdiffusion_replication"
}

# random seed
dqRNGkind("Xoroshiro128+")
dqset.seed(453779)
# Which moment conditions to use
case <- 1
# Whether first step or second step (0 = first step, 1 = secodn step)
twoStepOptimal <- 0
# Number of simulations
S <- 75
# Time span for one period
timeVector <- 'trimesters'
# Model type (1 -> qN = qP, 3 -> qN \ne qP)
modelType <- 1
# Whether bootstrap step or not (0 = No, 1 = Yes)
bootstrap <- 0

# --------------------------

# Village indices used for analyses ----------------
vills <- c(1:4, 6, 9, 12, 15, 19:21, 23:25, 29, 31:33, 36, 
           39, 42, 43, 45:48, 50:52, 55, 57, 59:60, 62, 
           64:65, 67:68, 70:73, 75)
num_vills <- length(vills)

# Select moments

if (case == 1){
  m <- 5
} else if (case == 2){
  m <- 3
} else if (case == 3){
  m <- 3
} else if (case == 4){
  m <- 3
}

# Select time vector and number of repetitions per trial --------------
# Months
TMonths <- c(31, 35, 15, 35, 13, 2, 32, 5, 
             31, 35, 31, 29, 19, 22, 25, 25, 
             23, 23, 24, 25, 26, 24, 17, 16, 
             17, 13, 19, 20, 20, 19, 22, 14, 
             12, 15, 10, 19, 18, 18, 19, 19, 19, 17, 17)

if (timeVector == 'months'){
    t_period <- TMonths + 1
} else if (timeVector == 'quarters'){
    t_period <- ceiling(TMonths / 3) + 1 # Quarters have 3 months in them
} else if (timeVector == 'trimesters'){
    t_period <- ceiling(TMonths / 4) + 1 # Trimesters have 4 months in them
}

# Load data ------------------------
# Adjacency matrix
X <- readMat(file.path(project_path, "datav4.0/Matlab Replication/India Networks/adjacencymatrix.mat"))
X <- X$X

# Load other matrices
leaders <- vector("list", length = num_vills)
TakeUp <- vector("list", length = num_vills)
EmpRate <- vector("list", length = num_vills)
inGiant <- vector("list", length = num_vills)
hermits <- vector("list", length = num_vills)
Z <- vector("list", length = num_vills)
TakingLeaders <- vector("list", length = num_vills)
ZLeaders <- vector("list", length = num_vills)
Outcome <- c()
Covars <- tibble()
Sec <- vector("list", length = num_vills)

for (i in seq_along(vills)){
  # Giant-component vectors (ie. vectors of indicators of whether belonging to giant components in each village)
  # (contained in a list)
  tempinGiant <- read_tsv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/inGiant", as.character(vills[i]), ".csv")), col_names = FALSE)
  tempinGiant <- as.logical(pull(tempinGiant))
  inGiant[[i]] <- tempinGiant
  
  # Leaders vectors (contained in a list)
  templeaders <- read_tsv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/HHhasALeader", as.character(vills[i]), ".csv")), col_names = FALSE)
  templeaders_all <- as.logical(pull(templeaders[,2]))
  templeaders <- as.logical(pull(templeaders[tempinGiant,2]))
  leaders[[i]] <- templeaders
  
  # Take-up vectors (contained in a list)
  tempTakeUp <- read_tsv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/MF", as.character(vills[i]), ".csv")), col_names = FALSE)
  EmpRate[[i]] <- mean(pull(tempTakeUp[!templeaders_all,]))
  tempTakeUp <- as.logical(pull(tempTakeUp[tempinGiant,]))
  TakeUp[[i]] <- tempTakeUp

  # Hermits (isolated HHs)
  d <- rowSums(X[[i]][[1]]) # number of neighbors
  hermits[[i]] <- (d == 0)
  
  # Covariates (only used ones, since W is used as a weight matrix later)
  tempZ <- read_tsv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/hhcovariates", as.character(vills[i]), ".csv")), col_names = FALSE)
  tempZ <- tempZ[tempinGiant,1:6]
  Z[[i]] <- tempZ
  
  # Leader statistics
  TakingLeaders[[i]] <- tempTakeUp[templeaders]
  ZLeaders[[i]] <- tempZ[templeaders,]
  Outcome <- c(Outcome, tempTakeUp[templeaders])
  Covars <- bind_rows(Covars, tempZ[templeaders,])
  
  # Second neighbors
  tempSec <- (X[[i]][[1]] %*% X[[i]][[1]] > 0)
  diag(tempSec) <- 0
  Sec[[i]] <- (tempSec - X[[i]][[1]] > 0)
  
}

# Select parameter grid --------------------

if (modelType == 1){
  qN <- c(seq(0, 0.01, 0.001), seq(0.05, 1, 0.05))
} else if (modelType == 3){
  qN <- c(seq(0, 0.01, 0.001), seq(0.05, 1, 0.05))
  qP <- c(seq(0, 0.1, 0.005), seq(0.15, 1, 0.05))
}

# Logistic fit to get coefficients for covariates and the constant ----------------------
leader_df <- as_tibble(Outcome) %>%
  bind_cols(Covars)
glm_res <- glm(value ~ ., data = leader_df, family = "binomial") 
Betas <- glm_res$coefficients

# Calculate network statistics (netstats) for each village ------------
netstats <- list()

for (i in seq_along(vills)){
  
  X_graph <- graph_from_adjacency_matrix(X[[i]][[1]], mode = "undirected")
  D <- distances(X_graph)
  
  minDistFromLeaders <- apply(D[, which(leaders[[i]])], 1, min)
  avgDistFromLeaders <- apply(D[, which(leaders[[i]])], 1, mean)
  
  infected <- TakeUp[[i]]
  
  if (dot(as.numeric(infected), as.numeric(leaders[[i]])) > 0){
    minDistInfectedLeaders <- apply(as.matrix(D[, which(leaders[[i]] & infected)]), 1, min)
  } else {
    minDistInfectedLeaders <- rep(0, nrow(leaders[[i]]))
  }
  
  if (dot(1 - as.numeric(infected), as.numeric(leaders[[i]])) > 0){
    minDistNonInfectedLeaders <- apply(D[, which(leaders[[i]] & !infected)], 1, min)
  } else {
    minDistNonInfectedLeaders <- rep(0, nrow(leaders[[i]]))
  }
  
  # variable of an indicator if neighboring leader is infected:
  # minimum distance to infected leaders is 1 & minimum distance to non-infected leaders is 0 or > 1
  ### If a HH has both a neighboring infected leader and a neighboring non-infected leader, then
  ### "neighborOfInfected" will be 0, but should it be 1?
  neighborOfInfected <- (((minDistInfectedLeaders == 1) - (minDistNonInfectedLeaders == 1)) > 0) 
  
  # variable of an indicator if neighboring leader is not infected:
  # minimum distance to infected leaders is 0 or >1 & minimum distance to non-infected leaders is 1 
  neighborOfNonInfected <- (((minDistInfectedLeaders == 1) - (minDistNonInfectedLeaders == 1)) < 0) 
  
  network_degree <- rowSums(X[[i]][[1]])
  
  netstats_j <- list(
    minDistFromLeaders = minDistFromLeaders,
    avgDistFromLeaders = avgDistFromLeaders,
    minDistInfectedLeaders = minDistInfectedLeaders,
    minDistNonInfectedLeaders = minDistNonInfectedLeaders,
    neighborOfInfected = neighborOfInfected,
    neighborOfNonInfected = neighborOfNonInfected,
    network_degree = network_degree
  )

  netstats[[i]] <- netstats_j
  
}

# Define a function to calculate moments (moments) ------------
moments <- function(X, leaders, netstats, infected, Sec, j, case){
  N <- nrow(X)
  
  network_degree <- netstats['network_degree'][[1]]
  neighborOfInfected <- netstats['neighborOfInfected'][[1]]
  neighborOfNonInfected <- netstats['neighborOfNonInfected'][[1]]


  if (case == 1){
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums(outer(rep(1, N), infected) * X) # Number of infected neighbors
    if (sum(infectedNeighbors == 0 & network_degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & infected == 1 & network_degree > 0)) / sum(infectedNeighbors == 0 & network_degree > 0)
    } else if (sum(infectedNeighbors == 0 & network_degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Fraction of individuals that are infected in the neighborhood of infected leaders stats_1 == 0
    if (sum(neighborOfInfected) > 0){
      stats_2 <- sum(infected * neighborOfInfected) / sum(neighborOfInfected)
    } else {
      stats_2 <- 0
    }

    # 3. Fraction of individuals that are infected in the neighborhood of non-infected leaders
    if (sum(neighborOfInfected) > 0){
      stats_3 <- sum(infected * neighborOfNonInfected) / sum(neighborOfNonInfected)
    } else {
      stats_3 <- 0
    }

    # 4. Covariance of individuals taking with share of neighbors taking
    NonHermits = (network_degree > 0)
    ShareofTakingNeighbors = infectedNeighbors[NonHermits] / network_degree[NonHermits]
    NonHermitTakers = infected[NonHermits]
    stats_4 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermits)

    # 5. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermits] / network_degree[NonHermits]
    stats_5 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermits)

    return(c(stats_1, stats_2, stats_3, stats_4, stats_5))

  } else if (case == 2){
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums(outer(rep(1, N), infected) * X) # Number of infected neighbors

    if (sum(infectedNeighbors == 0 & network_degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & infected == 1 & network_degree > 0)) / sum(infectedNeighbors == 0 & network_degree > 0)
    } else if (sum(infectedNeighbors == 0 & network_degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Covariance of individuals taking with share of neighbors taking
    NonHermits = (network_degree > 0)
    ShareofTakingNeighbors = infectedNeighbors[NonHermits] / network_degree[NonHermits]
    NonHermitTakers = infected[NonHermits]
    stats_2 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermits)

    # 3. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermits] / network_degree[NonHermits]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermits)

    return(c(stats_1, stats_2, stats_3))

  } else if (case == 3){
    # same as case 2, but purged of leader injection points.
    leaderTrue = (leaders > 0) # a variable that denotes whether a node is either a leader

    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums((outer(rep(1, N), infected)) * X) # Number of infected neighbors

    if (sum(infectedNeighbors == 0 & network_degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & (leaderTrue == 0) & infected == 1 & network_degree > 0)) / sum(infectedNeighbors == 0 & network_degree > 0)
    } else if (sum(infectedNeighbors == 0 & (leaderTrue == 0) & network_degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Covariance of individuals taking with share of neighbors taking
    NonHermits = (network_degree > 0)
    NonHermitsNonLeaders = (NonHermits & (1 - leaderTrue)) # not isolates, not leaders

    ShareofTakingNeighbors = infectedNeighbors[NonHermitsNonLeaders] / network_degree[NonHermitsNonLeaders]
    NonHermitTakers = infected[NonHermitsNonLeaders]
    stats_2 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermitsNonLeaders)

    # 3. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermitsNonLeaders] / network_degree[NonHermitsNonLeaders]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermitsNonLeaders)

    return(c(stats_1, stats_2, stats_3))

  } else if (case == 4){
    # same as case 3, but purged of ALL leader nodes.
    leaderTrue = (leaders > 0) # a variable that denotes whether a node is either a leader
    
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums(outer(rep(1, N), infected) * X %*% (1 - leaderTrue)) # Number of infected neighbors

    if (sum(infectedNeighbors == 0 & network_degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & (leaderTrue == 0) & infected == 1 & network_degree > 0)) / sum(infectedNeighbors == 0 & network_degree > 0)
    } else if (sum(infectedNeighbors == 0 & (leaderTrue == 0) & network_degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Covariance of individuals taking with share of neighbors taking
    NonHermits = (network_degree > 0)
    NonHermitsNonLeaders = (NonHermits & (1 - leaderTrue)) # not isolates, not leaders

    ShareofTakingNeighbors = infectedNeighbors[NonHermitsNonLeaders] / network_degree[NonHermitsNonLeaders]
    NonHermitTakers = infected[NonHermitsNonLeaders]
    stats_2 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermitsNonLeaders)

    # 3. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)) %*% (1 - leaderTrue))
    ShareofSecond = infectedSecond[NonHermitsNonLeaders] / network_degree[NonHermitsNonLeaders]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermitsNonLeaders)

    return(c(stats_1, stats_2, stats_3))

  }

}

# test
moments(X[[1]][[1]], leaders[[1]], netstats[[1]], TakeUp[[1]], Sec[[1]], 1, 1)

# Define a function to simulate diffusion of MF (diffusion_model) ---------------------

diffusion_model <- function(parms, Z, Betas, X, leaders, j, t_period, EmpRate){
  
  qN <- parms[1] # Probability non-taker transmits information
  qP <- parms[2] # Probsbility that a just-informed-taker transmits information
  N <- nrow(X) # Number of households
  
  infected <- rep(FALSE, N) # Nobody has been infected yet.
  infectedbefore <- rep(FALSE, N) # Nobody has been infected yet.
  contagiousbefore <- rep(FALSE, N) # People who were contagious before
  contagious <- leaders # Newly informed/contagious.
  dynamicInfection <- rep(0, t_period) # Will be a vector that tracks the infection rate for the number of periods it takes place
  
  x <- matrix(dqrunif(N * t_period), N, t_period)
  t <- 1
  for (t in seq(t_period)){
    qNt <- qN
    qPt <- qP
    
    # Step 1: Take-up decision based on newly informed
    LOGITprob <- 1 / (1 + exp(- cbind(rep(1, N), as.matrix(Z)) %*% Betas))
    infected <- ((!contagiousbefore & contagious & as.vector(x[,t] < LOGITprob)) | infected)
    s1 <- sum(infected)
    s2 <- sum(infectedbefore)
    infectedbefore <- (infectedbefore | infected)
    contagiousbefore <- (contagious | contagiousbefore)
    C <- sum(contagious)
    
    # Step 2: Information flows
    transmitPROB <- (contagious & infected) * qPt + (contagious & !infected) * qNt
    contagionlikelihood <- X[contagious,] * outer(transmitPROB[contagious], rep(1, N))
    
    # Step 3
    # start_time <- Sys.time()
    contagious <- ((colSums(contagionlikelihood > matrix(dqrunif(C * N), C, N)) > 0) | contagiousbefore)
    # contagious <- ((t(contagionlikelihood > matrix(dqrunif(C * N), C, N)) %*% rep(1, C) > 0) | contagiousbefore)
    dynamicInfection[t] <- sum(infectedbefore) / N
    # end_time <- Sys.time()
    # print(end_time - start_time)

  }
  
  return(list(infectedbefore, dynamicInfection, contagious))
  
}

# test
system.time(diffusion_model(c(0.3, 0.1), Z[[1]], Betas, X[[1]][[1]], leaders[[1]], 1, t_period[1], EmpRate[[1]]))

# Define a function to compute the deviation of the empirical moments from the simulated ones (divergence_model) ---------------------------
# Model 1 = q
# Model 3 = qN, qP
# This function relies on diffusion_model as the transmission process and moments

divergence_model <- function(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, t_period, EmpRate, case){

  # Parameters
  G <- length(X)

  # Computation of the vector of divergences across all the moments
  EmpiricalMoments <- matrix(0, G, m)
  MeanSimulatedMoments <- matrix(0, G, m)
  D <- matrix(0, G, m)
  TimeSim <- matrix(0, G, S)

  for (g in seq(G)){
    # Compute moments - G x m object
    start_time <- Sys.time()
    EmpiricalMoments[g,] <- moments(X[[g]][[1]], leaders[[g]], netstats[[g]], TakeUp[[g]], Sec[[g]], g, case)
    end_time <- Sys.time()
    print(end_time - start_time)
    
    # Compute simulated moments
    SimulatedMoments <- matrix(0, S, m)
    for (s in seq(S)){
      infectedSIM <- diffusion_model(theta, Z[[g]], Betas, X[[g]][[1]], leaders[[g]], g, t_period[g], EmpRate[[g]])
      SimulatedMoments[s,] <- moments(X[[g]][[1]], leaders[[g]], netstats[[g]], as.vector(infectedSIM[[1]]), Sec[[g]], g, case)
    }
    
    # Compute the mean simulated moment - a G x m object
    MeanSimulatedMoments[g,] <- colMeans(SimulatedMoments)
    D[g,] <- MeanSimulatedMoments[g,] - EmpiricalMoments[g,]
  }
  #return(D)
  div_list <- list(D, EmpiricalMoments, MeanSimulatedMoments)
  return(div_list)
}

# test
#system.time(divergence_model(X, Z, Betas, leaders, TakeUp, Sec, c(0.3, 0.1), 5, 75, t_period, EmpRate, 1))

# Running the model

if (modelType == 1){ # Case where qN = qP
  
  #D <- array(rep(0, num_vills * m * length(qN)), dim = c(num_vills, m, length(qN)))
  div_output <- list()
  
  for (i in seq(length(qN))){
    #print(i)
    theta <- c(qN[i], qN[i])
    #D[,,i] <- divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, t_period, EmpRate, case)
    div_output[[i]] <- divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, t_period, EmpRate, case)
  }
} else if (modelType == 3){ # Case where qN \ne qP
  D <- array(rep(0, num_vills * m * length(qN) * length(qP)), dim = c(num_vills, m, length(qN), length(qP)))
  for (i in seq(length(qN))){
    for (j in seq(length(qP))){
      # print(i)
      theta <- c(qN[i], qP[j])
      D[,,i,j] <- divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, t_period, EmpRate, case)
    }
  }
}

# Save the output ------------------
# file_name <- paste0('data_model_', as.character(modelType), '_mom_', as.character(case), '_', timeVector, '.RData')
# save(D, file = file.path(project_path, 'Rdata', file_name))

# Save the div_output
file_name <- paste0('data_model_', as.character(modelType), '_div_output_', as.character(case), '_', timeVector, '.RData')
save(div_output, file = file.path(project_path, 'Rdata', file_name))

# Running the aggregator ------------------

file_name <- paste0('data_model_', as.character(modelType), '_mom_', as.character(case), '_', timeVector, '.RData')
load(file = file.path('Rdata', file_name))

if (bootstrap == 0){
  B <- 1
} else if (bootstrap == 1){
  B <- 1000
}
Q <- zeros(B, 1)

# Two step optimal weights
if (twoStepOptimal == 1){
  if (modelType == 1){
    qN_info <- 0.1
    qP_info <- 0.1
    theta <- c(qN_info, qP_info)
    D <- divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, t_period, EmpRate, case);
    A <- (t(D) %*% D) / num_vills
    W <- inv(A)
  } else if (modelType == 3){
    # Load the first-step estimates
    file_name <- paste0('param_est_', as.character(modelType), '_mom_', as.character(case), '_', timeVector, '_', as.character(0), '.RData')
    load(param_est)
    
    qN_info <- param_est[1,1]
    qP_info <- param_est[1,2]
    theta <- c(qN_info, qP_info)
    D <- divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, t_period, EmpRate, case);
    A <- (t(D) %*% D) / num_vills
    W <- inv(A)
  }
} else if (twoStepOptimal == 0){
  W <- eye(m)
}

# Aggregate -------------

# weights for bootstrap
wt <- zeros(B, num_vills)
for (b in seq(B)){

  # Generate weights b
  if (bootstrap == 1){
    wt[b,] <- rexp(num_vills)
    wt[b,] <- wt[b,] / mean(wt[b,])
  } else if (bootstrap == 0){
    wt[b,] <- rep(1 / num_vills, num_vills)
  }

  # For each model, generate the criterion function value for this
  # bootstrap run

  # Info model
  if (modelType == 1){
    momFunc <- array(rep(0, length(qN) * B * m), dim = c(length(qN), B, m))
    Qa <- array(rep(0, length(qN) * B), dim = c(length(qN), B))
    param_est <- zeros(B, 1)
    for (i in seq(length(qN))){
      # Compute the moment function
      momFunc[i,b,] <- wt[b,] %*% D[,,i] / num_vills
      # Criterion function
      Qa[i,b] <- t(momFunc[i,b,]) %*% W %*% momFunc[i,b,]
    }
    param_est[b] <- qN[which.min(Qa[,b])]
  } else if (modelType == 3){
    momFunc <- array(rep(0, length(qN) * length(qP) * B * m), dim = c(length(qN), length(qP), B, m))
    Qa <- array(rep(0, length(qN) * length(qP) * B), dim = c(length(qN), length(qP), B))
    param_est <- zeros(B, 2)
    for (i in seq(length(qN))){
      for (j in seq(length(qP))){
        # Compute the moment function
        momFunc[i,j,b,] <- wt[b,] %*% D[,,i,j] / num_vills
        # Criterion function
        Qa[i,j,b] <- t(momFunc[i,j,b,]) %*% W %*% momFunc[i,j,b,]
      }
    }
    min_ind <- which(Qa[,,b] == min(Qa[,,b]), arr.ind = TRUE)
    param_est[b,1] <- qN[min_ind[1]]
    param_est[b,2] <- qP[min_ind[2]]
  }
}

file_name <- paste0('param_est_', as.character(modelType), '_mom_', as.character(case), '_', timeVector, '_', as.character(twoStepOptimal), '.RData')
save(param_est, file = file.path('Rdata', file_name))

# plot_ly(x = qP, y = qN, z = Qa[,,1], type = "heatmap")
# plot_ly(x = qP, y = qN, z = Qa[,,1], type = "surface")
# plot_ly(x = qP, y = qN, z = log(Qa[,,1]), type = "contour")

if (modelType == 3 & bootstrap == 0){
  plot_df <- expand_grid(qP, qN) %>%
    mutate(Qa = c(Qa[,,1]))
  
  filename <- paste0("Figures/", 'data_model_', as.character(modelType), '_mom_', as.character(case), '_', timeVector, '_', as.character(twoStepOptimal), '.pdf')
  pdf(file = filename)
  ggplot(plot_df, aes(x = qP, y = qN, z = log(Qa))) + 
    geom_contour_filled() +
    geom_point(aes(x = param_est[1,2], y = param_est[1,1]), color = 'red', shape = 3) +
    scale_fill_viridis_d(name = 'log(criterion function)') 
  dev.off()
}


## breadthdistRAL -------------------------
#
#CIJ <- X[[1]][[1]]
#dummies <- leaders[[1]]
#
#breadthdistRAL <- function(CIJ, dummies){
#  N <- nrow(CIJ)
#
#  D <- matrix(rep(0, N * N), N, N)
#  for (i in seq_along(dummies)){
#    if (dummies[i] != 0){
#      D[, i] <- breadth(CIJ, i)
#    }
#  }
#
#  # replace zeros with 'Inf's
#  D[is.infinite(D)] <- 999999
#
#  # construct R
#  R <- as.double(D[D != 999999]) # ???
#
#  return(list(R, D))
#}
#
## breadth -------------------------
#
#source <- 2
#
#breadth <- function(CIJ, source){
#  N <- nrow(CIJ)
#
#  # colors: white, gray, black
#  white <- 0
#  gray <- 1
#  black <- 2
#
#  # initialize colors
#  color <- t(rep(0, N))
#  # initialize distances
#  distance <- Inf * t(rep(1, N))
#  # initialize branches
#  branch <- t(rep(0, N))
#
#  # start on vertex 'source'
#  color[source] <- gray
#  distance[source] <- 0
#  branch[source] <- -1
#  Q <- source
#
#  # keep going until the entire graph is explored
#  while (length(Q) > 0){
#    u <- Q[1]
#    ns <- (CIJ[u,] > 0)
#    for (v in seq(ns)){
#      if (distance[v] == 0){
#        distance[v] <- distance[u] + 1
#      }
#      if (color[v] == white){
#        color[v] <- gray
#        distance[v] <- distance[u] + 1
#        branch[v] <- u
#        Q <- c(Q, v)
#      }
#    }
#    Q <- Q[2:length(Q)] # ??
#    color[u] <- black
#  }
#  
#  return(list(distance, branch))
#
#}



