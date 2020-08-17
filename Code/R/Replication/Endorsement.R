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
  "Matrix",
  "microbenchmark",
  "Rfast",
  "RcppArmadillo"
)

pacman::p_load(packages, character.only = TRUE)

# Global setting
user <- "Hiro"
if (user == "Hiro"){
  project_path <- "/Users/mizuhirosuzuki/Dropbox/MFdiffusion_replication/"
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
# Model type (2 -> qN = qP, 4 -> qN \ne qP)
modelType <- 2
# Whether bootstrap step or not (0 = No, 1 = Yes)
bootstrap <- 0

relative = 0

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
OmegaE <- vector("list", length = num_vills)
OmegaD <- vector("list", length = num_vills)
OmegaN <- vector("list", length = num_vills)
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
  
  # Omega matrices (contained in a list)
  if (relative == 0){
    tempOmegaE <- read_csv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/Omega_abs", as.character(vills[i]), ".csv")), col_names = FALSE)
    tempOmegaD <- read_csv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/DOmega_abs", as.character(vills[i]), ".csv")), col_names = FALSE)
    tempOmegaN <- read_csv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/NOmega_abs", as.character(vills[i]), ".csv")), col_names = FALSE)
  } else if (relative == 1){
    tempOmegaE <- read_csv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/Omega_rel", as.character(vills[i]), ".csv")), col_names = FALSE)
    tempOmegaD <- read_csv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/DOmega_rel", as.character(vills[i]), ".csv")), col_names = FALSE)
    tempOmegaN <- read_csv(file.path(project_path, paste0("datav4.0/Matlab Replication/India Networks/NOmega_rel", as.character(vills[i]), ".csv")), col_names = FALSE)
  }
  tempOmegaE <- tempOmegaE[tempinGiant, tempinGiant]
  tempOmegaD <- tempOmegaD[tempinGiant, tempinGiant]
  tempOmegaN <- tempOmegaN[tempinGiant, tempinGiant]
  OmegaE[[i]] <- as.matrix(tempOmegaE)
  OmegaD[[i]] <- as.matrix(tempOmegaD)
  OmegaN[[i]] <- as.matrix(tempOmegaN)
  
  # Hermits (isolated HHs)
  d <- rowsums(X[[i]][[1]]) # number of neighbors
  hermits[[i]] <- (d == 0)
  
  # Covariates (only variables that will be used in the analysis below are loaded, since W is used as a weight matrix later)
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
lambda <- c(seq(-1, -0.3, 0.1), seq(-0.25, 0.3, 0.05), seq(0.4, 1, 0.1))

if (modelType == 2){
  qN <- c(seq(0, 0.5, 0.05), seq(0.6, 1, 0.1))
} else if (modelType == 4){
  qN <- c(seq(0, 0.5, 0.05), seq(0.6, 1, 0.1))
  qP <- c(seq(0, 0.5, 0.05), seq(0.6, 1, 0.1))
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
  
  network_degree <- rowsums(X[[i]][[1]])
  
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
    infectedNeighbors <- rowsums(outer(rep(1, N), infected) * X) # Number of infected neighbors
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
    infectedSecond = rowsums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermits] / network_degree[NonHermits]
    stats_5 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermits)

    return(c(stats_1, stats_2, stats_3, stats_4, stats_5))

  } else if (case == 2){
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowsums(outer(rep(1, N), infected) * X) # Number of infected neighbors

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
    infectedSecond = rowsums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermits] / network_degree[NonHermits]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermits)

    return(c(stats_1, stats_2, stats_3))

  } else if (case == 3){
    # same as case 2, but purged of leader injection points.
    leaderTrue = (leaders > 0) # a variable that denotes whether a node is either a leader

    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowsums((outer(rep(1, N), infected)) * X) # Number of infected neighbors

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
    infectedSecond = rowsums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermitsNonLeaders] / network_degree[NonHermitsNonLeaders]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermitsNonLeaders)

    return(c(stats_1, stats_2, stats_3))

  } else if (case == 4){
    # same as case 3, but purged of ALL leader nodes.
    leaderTrue = (leaders > 0) # a variable that denotes whether a node is either a leader
    
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowsums(outer(rep(1, N), infected) * X %*% (1 - leaderTrue)) # Number of infected neighbors

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
    infectedSecond = rowsums(Sec * outer(infected, rep(1, N)) %*% (1 - leaderTrue))
    ShareofSecond = infectedSecond[NonHermitsNonLeaders] / network_degree[NonHermitsNonLeaders]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermitsNonLeaders)

    return(c(stats_1, stats_2, stats_3))

  }

}

# test
moments(X[[1]][[1]], leaders[[1]], netstats[[1]], TakeUp[[1]], Sec[[1]], 1, 1)

# Define a function to simulate diffusion of MF (endorsement_model) ---------------------

endorsement_model <- function(parms, Z, Betas, X, leaders, OmegaE, OmegaD, OmegaN, j, t_period, EmpRate){
  
  qN <- parms[1] # Probability non-taker transmits information
  qP <- parms[2] # Probability that a just-informed-taker transmits information
  lambda <- parms[3] # Coefficient in logit governing the probability of being a taker
  N <- nrow(X) # Number of households
  
  infectedbefore_list <- list()
  dynamicInfection_list <- list()
  contagiousbefore_list <- list()
  x <- matrix(dqrunif(N * t_period), N, t_period)
  Omega_list <- list(OmegaE, OmegaD, OmegaN)
  
  for (k in seq(3)){
    
    infected <- rep(FALSE, N) # Nobody has been infected yet.
    infectedbefore <- rep(FALSE, N) # Nobody has been infected yet.
    contagiousbefore <- rep(FALSE, N) # People who were contagious before
    contagious <- leaders # Newly informed/contagious.
    transmissionHist <- matrix(rep(FALSE, N * N), nrow = N, ncol = N) 
    dynamicInfection <- rep(0, t_period) # Will be a vector that tracks the infection rate 
                                         # for the number of periods it takes place
    
    Omega <- Omega_list[[k]]
    
    for (t in seq(t_period)){
      
      # Step 1: Take-up decision based on newly informed
      z <- outer(infectedbefore, rep(1, N))
      # start_time <- Sys.time()
      # regressor <- diag(Omega %*% (transmissionHist * z)) / diag(Omega %*% transmissionHist)
      # regressor <- colSums(t(Omega) * (transmissionHist * z)) / colSums(t(Omega) * transmissionHist)
      regressor <- colsums(t(Omega) * (transmissionHist * z)) / colsums(t(Omega) * transmissionHist)
      # end_time <- Sys.time()
      # print(end_time - start_time)
      regressor[is.na(regressor) | is.infinite(regressor)] <- 0
      LOGITprob <- 1 / (1 + exp(- cbind(rep(1, N), as.matrix(Z)) %*% Betas) - regressor * lambda)
      
      infected <- ((!contagiousbefore & contagious & as.vector(x[,t] < LOGITprob)) | infected)
      s1 <- sum(infected)
      s2 <- sum(infectedbefore)
      infectedbefore <- (infectedbefore | infected)
      contagiousbefore <- (contagious | contagiousbefore)
      C <- sum(contagious)
      
      # Step 2: Information flows
      transmitPROB <- (contagious & infected) * qP + (contagious & !infected) * qN
      contagionlikelihood <- X * outer(transmitPROB, rep(1, N))
      
      # Step 3
      t0 <- matrix(dqrunif(N * N), nrow = N, ncol = N)
      t1 <- (contagionlikelihood > t0) # a full transmission matrix
      
      # zero out stuff because one can only transmit if contagious, and one
      # cannot be transmitted to unless they were not contagious before
      t1[!contagious,] <- FALSE
      transmissionHist <- (transmissionHist | t1)
      
      t2 <- t1[contagious, !contagiousbefore] # which contagious households transmit to previously non-contagious
      contagious[!contagiousbefore] <- (sum(t2) > 0)
      dynamicInfection[t] <- sum(infectedbefore) / N
  
    } 
    
    infectedbefore_list[[k]] <- infectedbefore
    dynamicInfection_list[[k]] <- dynamicInfection
    contagiousbefore_list[[k]] <- contagiousbefore
    
  }
   
  return(list(infectedbefore_list, dynamicInfection_list, contagiousbefore_list))
  
}

# test
system.time(endorsement_model(c(0.3, 0.1, 0.2), Z[[2]], Betas, X[[2]][[1]], leaders[[2]], OmegaE[[2]], OmegaD[[2]], OmegaN[[2]], 2, t_period[2], EmpRate[[2]]))

# Define a function to compute the deviation of the empirical moments from the simulated ones (divergence_endorsement_model) ---------------------------
# Model 2 = q
# Model 4 = qN, qP
# This function relies on endorsement_model as the transmission process and moments

divergence_endorsement_model <- function(X, Z, Betas, leaders, OmegaE, OmegaD, OmegaN, TakeUp, Sec, theta, m, S, t_period, EmpRate, case){

  # Parameters
  G <- length(X)

  # Computation of the vector of divergences across all the moments
  EmpiricalMoments <- matrix(0, G, m)
  MeanSimulatedMomentsE <- matrix(0, G, m)
  MeanSimulatedMomentsD <- matrix(0, G, m)
  MeanSimulatedMomentsN <- matrix(0, G, m)
  
  DE <- matrix(0, G, m)
  DD <- matrix(0, G, m)
  DN <- matrix(0, G, m)

  for (g in seq(G)){
    
    print(g)
    # Compute moments - G x m object
    EmpiricalMoments[g,] <- moments(X[[g]][[1]], leaders[[g]], netstats[[g]], TakeUp[[g]], Sec[[g]], g, case)
    
    # Compute simulated moments
    SimulatedMomentsE <- matrix(0, S, m)
    SimulatedMomentsD <- matrix(0, S, m)
    SimulatedMomentsN <- matrix(0, S, m)
    for (s in seq(S)){
      endorsement_model_output <- endorsement_model(theta, Z[[g]], Betas, X[[g]][[1]], leaders[[g]], OmegaE[[g]], OmegaD[[g]], OmegaN[[g]], g, t_period[g], EmpRate[[g]])
      infectedbefore_output <- endorsement_model_output[[1]]
      SimulatedMomentsE[s,] <- moments(X[[g]][[1]], leaders[[g]], netstats[[g]], as.vector(infectedbefore_output[[1]]), Sec[[g]], g, case)
      SimulatedMomentsD[s,] <- moments(X[[g]][[1]], leaders[[g]], netstats[[g]], as.vector(infectedbefore_output[[2]]), Sec[[g]], g, case)
      SimulatedMomentsN[s,] <- moments(X[[g]][[1]], leaders[[g]], netstats[[g]], as.vector(infectedbefore_output[[3]]), Sec[[g]], g, case)
    }
    
    # Compute the mean simulated moment - a G x m object
    MeanSimulatedMomentsE[g,] <- colMeans(SimulatedMomentsE)
    MeanSimulatedMomentsD[g,] <- colMeans(SimulatedMomentsD)
    MeanSimulatedMomentsN[g,] <- colMeans(SimulatedMomentsN)
    
    DE[g,] <- MeanSimulatedMomentsE[g,] - EmpiricalMoments[g,]
    DD[g,] <- MeanSimulatedMomentsD[g,] - EmpiricalMoments[g,]
    DN[g,] <- MeanSimulatedMomentsN[g,] - EmpiricalMoments[g,]
  }
  
  return(list(DE, DD, DN))

}

# test
system.time(divergence_endorsement_model(X, Z, Betas, leaders, OmegaE, OmegaD, OmegaN, TakeUp, Sec, c(0.3, 0.1, 0.1), 5, 75, t_period, EmpRate, 1))

# Running the model

if (modelType == 2){ # Case where qN = qP
  DE <- array(rep(0, num_vills * m * length(qN) * length(lambda)), dim = c(num_vills, m, length(qN), length(lambda)))
  DD <- array(rep(0, num_vills * m * length(qN) * length(lambda)), dim = c(num_vills, m, length(qN), length(lambda)))
  DN <- array(rep(0, num_vills * m * length(qN) * length(lambda)), dim = c(num_vills, m, length(qN), length(lambda)))
  
  for (i in seq(length(qN))){
    for (j in seq(length(lambda))){
      print(i)
      theta <- c(qN[i], qN[i], lambda[j])
      model_output <- divergence_endorsement_model(X, Z, Betas, leaders, OmegaE, OmegaD, OmegaN, TakeUp, Sec, theta, m, S, t_period, EmpRate, case)
      DE[,,i,j] <- model_output[[1]]
      DD[,,i,j] <- model_output[[2]]
      DN[,,i,j] <- model_output[[3]]
    }
  }
} else if (modelType == 4){ # Case where qN \ne qP
  DE <- array(rep(0, num_vills * m * length(qN) * length(qP) * length(lambda)), dim = c(num_vills, m, length(qN), length(qP), length(lambda)))
  DD <- array(rep(0, num_vills * m * length(qN) * length(qP) * length(lambda)), dim = c(num_vills, m, length(qN), length(qP), length(lambda)))
  DN <- array(rep(0, num_vills * m * length(qN) * length(qP) * length(lambda)), dim = c(num_vills, m, length(qN), length(qP), length(lambda)))
  for (i in seq(length(qN))){
    for (j in seq(length(qP))){
      for (k in seq(length(lambda))){
        # print(i)
        theta <- c(qN[i], qP[j], lambda[k])
        model_output <- divergence_endorsement_model(X, Z, Betas, leaders, OmegaE, OmegaD, OmegaN, TakeUp, Sec, theta, m, S, t_period, EmpRate, case)
        DE[,,i,j,k] <- model_output[[1]]
        DD[,,i,j,k] <- model_output[[2]]
        DN[,,i,j,k] <- model_output[[3]]
      }
    }
  }
}

D_list <- list(DE, DD, DN)

# Save the output ------------------
file_name <- paste0('data_model_', as.character(modelType), '_mom_', as.character(case), '_', timeVector, '.RData')
save(D_list, file = file.path('Rdata', file_name))

# Running the aggregator ------------------

file_name <- paste0('data_model_', as.character(modelType), '_mom_', as.character(case), '_', timeVector, '.RData')
load(file = file.path('Rdata', file_name))
DE <- D_list[[1]]
DD <- D_list[[2]]
DN <- D_list[[3]]

if (bootstrap == 0){
  B <- 1
} else if (bootstrap == 1){
  B <- 1000
}


if (modelType == 2){
  DDTotal <- array(rep(0, length(qN) * length(lambda) * G * m), dim = c(length(qN), length(lambda), G, m))
  DETotal <- array(rep(0, length(qN) * length(lambda) * G * m), dim = c(length(qN), length(lambda), G, m))
  DNTotal <- array(rep(0, length(qN) * length(lambda) * G * m), dim = c(length(qN), length(lambda), G, m))
  
  for (i in seq(length(qN))){
    for (k in seq(length(lambda))){
      DDTotal[i,k,,] <- DD[i,k]
      DETotal[i,k,,] <- DE[i,k]
      DNTotal[i,k,,] <- DN[i,k]
    }
  }
} else if (modelType == 4){
  DDTotal <- array(rep(0, length(qN) * length(qP) * length(lambda) * G * m), dim = c(length(qN), length(qP), length(lambda), G, m))
  DETotal <- array(rep(0, length(qN) * length(qP) * length(lambda) * G * m), dim = c(length(qN), length(qP), length(lambda), G, m))
  DNTotal <- array(rep(0, length(qN) * length(qP) * length(lambda) * G * m), dim = c(length(qN), length(qP), length(lambda), G, m))
  
  for (i in seq(length(qN))){
    for (k in seq(length(lambda))){
      DDTotal[i,j,k,,] <- DD[i,j,k]
      DETotal[i,j,k,,] <- DE[i,j,k]
      DNTotal[i,j,k,,] <- DN[i,j,k]
    }
  }
}

# Two step optimal weights
if (twoStepOptimal == 1){
  # For E matrix
  #theta <- c(qN_E, qP_E, lambda_E)
  theta <- c(0.1, 0.1, 0.3)
  E_out <- divergence_endorsement_model(
    X, Z, Betas, leaders, OmegaE, OmegaD, OmegaN, TakeUp, 
    Sec, theta, m, S, t_period, EmpRate, case
    )
  DE <- E_out[[1]]
  AE <- (t(DE) %*% DE) / num_vills
  WE <- inv(AE)
  
  # For D matrix
  #theta <- c(qN_D, qP_D, lambda_D)
  theta <- c(0.1, 0.1, 0.3)
  D_out <- divergence_endorsement_model(
    X, Z, Betas, leaders, OmegaE, OmegaD, OmegaN, TakeUp, 
    Sec, theta, m, S, t_period, EmpRate, case
    )
  DD <- D_out[[1]]
  AD <- (t(DD) %*% DD) / num_vills
  WD <- inv(AD)
  
  # For N matrix
  #theta <- c(qN_N, qP_N, lambda_N)
  theta <- c(0.1, 0.1, 0.3)
  N_out <- divergence_endorsement_model(
    X, Z, Betas, leaders, OmegaE, OmegaD, OmegaN, TakeUp, 
    Sec, theta, m, S, t_period, EmpRate, case
    )
  DN <- N_out[[1]]
  AN <- (t(DN) %*% DN) / num_vills
  WN <- inv(AN)
} else if (twoStepOptimal == 0){
  WE <- eye(m)
  WD <- eye(m)
  WN <- eye(m)
}

# Pre-allocation ---------
QEndE <- rep(0, B)
QEndD <- rep(0, B)
QEndN <- rep(0, B)
TestEndEEndD <- rep(0, B)
TestEndDEndN <- rep(0, B)
TestEndEEndN <- rep(0, B)
importantparmsE <- c()
importantparmsD <- c()
importantparmsN <- c()
valE <- c()
valD <- c()
valN <- c()

goopy

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

  # Endorsement model
  if (modelType == 2){
    momFuncE <- array(rep(0, length(qN) * length(lambda) * B * m), dim = c(length(qN), length(lambda), B, m))
    momFuncD <- array(rep(0, length(qN) * length(lambda) * B * m), dim = c(length(qN), length(lambda), B, m))
    momFuncN <- array(rep(0, length(qN) * length(lambda) * B * m), dim = c(length(qN), length(lambda), B, m))
    
    for (i in seq(length(qN))){
      for (k in seq(length(lambda))){
        # Compute the moment function
        momFuncE[i,k,b,] <- wt[b,] %*% t(DETotal[i,k,,]) / G
        momFuncD[i,k,b,] <- wt[b,] %*% t(DDTotal[i,k,,]) / G
        momFuncN[i,k,b,] <- wt[b,] %*% t(DNTotal[i,k,,]) / G
        
        # Criterion function
        QEndorseE[i,k,b] <- t(momFuncE[i,k,b,]) %*% WE %*% momFuncE[i,k,b,]
        QEndorseD[i,k,b] <- t(momFuncD[i,k,b,]) %*% WD %*% momFuncD[i,k,b,]
        QEndorseN[i,k,b] <- t(momFuncN[i,k,b,]) %*% WN %*% momFuncN[i,k,b,]
      }
    }
    
    
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
    
    
    
  } else if (modelType == 4){
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



