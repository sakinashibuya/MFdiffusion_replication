rm(list = ls())

packages <- c(
  "R.matlab",
  "tidyverse"
)

not_installed <- !packages %in% installed.packages()
if (any(not_installed)) install.packages(packages[not_installed])
lapply(packages,require,character.only=TRUE)

# Global setting
user <- "Hiro"
if (user == "Hiro"){
  project_path <- "/Users/mizuhirosuzuki/Dropbox/MFdiffusion_replication/"
}

# --------------------------

# Village indices used for analyses
vills <- c(1:4, 6, 9, 12, 15, 19:21, 23:25, 29, 31:33, 36, 
           39, 42, 43, 45:48, 50:52, 55, 57, 59:60, 62, 
           64:65, 67:68, 70:73, 75)
num_vills <- length(vills)

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

# Logistic fit to get coefficients for covariates and the constant
leader_df <- as_tibble(Outcome) %>%
  bind_cols(Covars)
glm_res <- glm(value ~ ., data = leader_df, family = "binomial") 
Betas <- glm_res$coefficients

# Divergence_Model ---------------------------
# This computes the deviation of the empirical moments from the simulated ones for
# Model 1 = q
# Model 3 = qN, qP
# This function relies on diffusion_model as the transmission process and moments

divergence_model <- function(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, T, EmpRate, version){

  # Parameters
  G <- length(X)

  # Computation of the vector of divergences across all the moments
  EmpiricalMoments <- matrix(0, G, m)
  MeanSimulatedMoments <- matrix(0, G, m)
  D <- matrix(0, G, m)
  TimeSim <- matrix(0, G, S)

  for (g in seq(G)){
    # Compute moments - G x m object
    EmpiricalMoments[g,] <- moments(X[[g]][[1]], leaders[[g]], TakeUp[[g]], Sec[[g]], g, version)
    
    # Compute simulated moments
    SimulatedMoments <- matrix(0, S, m)
    for (s in seq(S)){
      infectedSIM <- diffusion_model(theta, Z[[g]], Betas, X[[g]][[1]], leaders[[g]], g, T[g], EmpRate[[g]])
      SimulatedMoments[s,] <- moments(X[[g]][[1]], leaders[[g]], infectedSIM, Sec[[g]], g, version)
    }
    
    # Compute the mean simulated moment - a G x m object
    MeanSimulatedMoments[g,] <- colMeans(SimulatedMoments)
    D[g,] <- MeanSimulatefMoments[g,] - EmpiricalMoments[g,]
  }
  
  return(D)

}

# diffusion_model ---------------------

diffusion_model <- function(parms, Z, Betas, X, leaders, j, T, EmpRate){
  
  qN <- parms[1] # Probability non-taker transmits information
  qP <- parms[2] # Probsbility that a just-informed-taker transmits information
  N <- nrow(X) # Number of households
  
  infected <- rep(FALSE, N) # Nobody has been infected yet.
  infectedbefore <- rep(FALSE, N) # Nobody has been infected yet.
  contagiousbefore <- rep(FALSE, N) # People who were contagious before
  contagious <- leaders # Newly informed/contagious.
  dynamicInfection <- rep(0, T) # Will be a vector that tracks the infection rate for the number of periods it takes place
  
  x <- matrix(runif(N * T), N, T)
  t <- 1
  for (t in seq(T)){
    print(t)
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
    contagious <- ((t(contagionlikelihood > matrix(runif(C * N), N, C)) %*% rep(1, C) > 0) | contagiousbefore)
    dynamicInfection[t] <- sum(infectedbefore) / N
    
  }
  
  return(list(infectedbefore, dynamicInfection, contagious))
  
}

# moments ------------
moments <- function(X, leaders, infected, Sec, j, version){
  N <- nrow(X)

  #####

  #####

  if (case == 1){
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums(matrix(1,N,1) %*% t(infected) * X) # Number of infected neighbors

    if (sum(infectedNeighbors == 0 & netstats[j].degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & infected == 1 & netstats[j].degree > 0)) / sum(infectedNeighbors == 0 & netstats[j].degree > 0)
    } else if (sum(infectedNeighbors == 0 & netstats[j].degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Fraction of individuals that are infected in the neighborhood of infected leaders stats_1 == 0
    if (sum(netstas[j].neighborOfInfected) > 0){
      stats_2 <- sum(infected * netstats[j].neighborOfInfected) / sum(netstats[j].neighborOfInfected)
    } else {
      stats_2 <- 0
    }

    # 3. Fraction of individuals that are infected in the neighborhood of non-infected leaders
    if (sum(netstats[j].neighborOfInfected) > 0){
      stats_3 <- sum(infected * netstats[j].neighborOfNonInfected) / sum(netstats[j].neighborOfNonInfected)
    } else {
      stats_3 <- 0
    }

    # 4. Covariance of individuals taking with share of neighbors taking
    NonHermits = (netstats[j].degree > 0)
    ShareofTakingNeighbors = infectedNeighbors[NonHermits] / netstats[j].degree[NonHermits]
    NonHermitTakers = infected[NonHermits]
    stats_4 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermits)

    # 5. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermits] / netstats[j].degree[NonHermits]
    stats_5 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermits)

    return(c(stats_1, stats_2, stats_3, stats_4, stats_5))

  } else if (case == 2){
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums((ones(N, 1) * t(infected)) * X) # Number of infected neighbors

    if (sum(infectedNeighbors == 0 & netstats[j].degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & infected == 1 & netstats[j].degree > 0)) / sum(infectedNeighbors == 0 & netstats[j].degree > 0)
    } else if (sum(infectedNeighbors == 0 & netstats[j].degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Covariance of individuals taking with share of neighbors taking
    NonHermits = (netstats[j].degree > 0)
    ShareofTakingNeighbors = infectedNeighbors[NonHermits] / netstats[j].degree[NonHermits]
    NonHermitTakers = infected[NonHermits]
    stats_2 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermits)

    # 3. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermits] / netstats[j].degree[NonHermits]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermits)

    return(c(stats_1, stats_2, stats_3))

  } else if (case == 3){
    # same as case 2, but purged ofleader injection points.
    leaderTrue = (leaders > 0) # a variable that denotes whether a node is either a leader

    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums((ones(N, 1) * t(infected)) * X) # Number of infected neighbors

    if (sum(infectedNeighbors == 0 & netstats[j].degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & (leaderTrue == 0) & infected == 1 & netstats[j].degree > 0)) / sum(infectedNeighbors == 0 & netstats[j].degree > 0)
    } else if (sum(infectedNeighbors == 0 & (leaderTrue == 0) & netstats[j].degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Covariance of individuals taking with share of neighbors taking
    NonHermits = (netstats[j].degree > 0)
    NonHermitsNonLeaders = (NonHermits & (1 - leaderTrue)) # not isolates, not leaders

    ShareofTakingNeighbors = infectedNeighbors[NonHermitsNonLeaders] / netstats[j].degree[NonHermitsNonLeaders]
    NonHermitTakers = infected[NonHermitsNonLeaders]
    stats_2 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermitsNonLeaders)

    # 3. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)))
    ShareofSecond = infectedSecond[NonHermitsNonLeaders] / netstats[j].degree[NonHermitsNonLeaders]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermitsNonLeaders)

    return(c(stats_1, stats_2, stats_3))

  } else if (case == 4){
    # same as case 3, but purged of ALL leader nodes.
    leaderTrue = (leaders > 0) # a variable that denotes whether a node is either a leader
        
    # 1. Fraction of nodes that have no taking neighbors but are takers themselves
    infectedNeighbors <- rowSums((ones(N, 1) * t(infected)) * X %*% (1 - leaderTrue)) # Number of infected neighbors

    if (sum(infectedNeighbors == 0 & netstats[j].degree > 0) > 0){
      stats_1 <- sum((infectedNeighbors == 0 & (leaderTrue == 0) & infected == 1 & netstats[j].degree > 0)) / sum(infectedNeighbors == 0 & netstats[j].degree > 0)
    } else if (sum(infectedNeighbors == 0 & (leaderTrue == 0) & netstats[j].degree > 0) == 0){
      stats_1 <- 0
    }

    # 2. Covariance of individuals taking with share of neighbors taking
    NonHermits = (netstats[j].degree > 0)
    NonHermitsNonLeaders = (NonHermits & (1 - leaderTrue)) # not isolates, not leaders

    ShareofTakingNeighbors = infectedNeighbors[NonHermitsNonLeaders] / netstats[j].degree[NonHermitsNonLeaders]
    NonHermitTakers = infected[NonHermitsNonLeaders]
    stats_2 <- sum(NonHermitTakers * ShareofTakingNeighbors) / sum(NonHermitsNonLeaders)

    # 3. Covariance of individuals taking with share of second neighbors taking
    infectedSecond = rowSums(Sec * outer(infected, rep(1, N)) %*% (1 - leaderTrue))
    ShareofSecond = infectedSecond[NonHermitsNonLeaders] / netstats[j].degree[NonHermitsNonLeaders]
    stats_3 <- sum(NonHermitTakers * ShareofSecond) / sum(NonHermitsNonLeaders)

    return(c(stats_1, stats_2, stats_3))

  }

}

# breadthdistRAL -------------------------

CIJ <- X[[1]][[1]]
dummies <- leaders[[1]]

breadthdistRAL <- function(CIJ, dummies){
  N <- nrow(CIJ)

  D <- matrix(rep(0, N * N), N, N)
  for (i in seq_along(dummies)){
    if (dummies[i] != 0){
      D[, i] <- breadth(CIJ, i)
    }
  }

  # replace zeros with 'Inf's
  D[is.infinite(D)] <- 999999

  # construct R
  R <- as.double(D[D != 999999]) # ???

  return(list(R, D))
}

# breadth -------------------------

source <- 2

breadth <- function(CIJ, source){
  N <- nrow(CIJ)

  # colors: white, gray, black
  white <- 0
  gray <- 1
  black <- 2

  # initialize colors
  color <- t(rep(0, N))
  # initialize distances
  distance <- Inf * t(rep(1, N))
  # initialize branches
  branch <- t(rep(0, N))

  # start on vertex 'source'
  color[source] <- gray
  distance[source] <- 0
  branch[source] <- -1
  Q <- source

  # keep going until the entire graph is explored
  while (length(Q) > 0){
    u <- Q[1]
    ns <- (CIJ[u,] > 0)
    for (v in seq(ns)){
      if (distance[v] == 0){
        distance[v] <- distance[u] + 1
      }
      if (color[v] == white){
        color[v] <- gray
        distance[v] <- distance[u] + 1
        branch[v] <- u
        Q <- c(Q, v)
      }
    }
    Q <- Q[2:length(Q)] # ??
    color[u] <- black
  }
  
  return(list(distance, branch))

}



