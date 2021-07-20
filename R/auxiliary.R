#library(Rcpp)
#library(RcppArmadillo)
#library(doParallel)
#registerDoParallel()

#Rcpp::sourceCpp("./src/fast_dist_calc.cpp")

#' A Hamming Distance Vector Calculator
#'
#' Computes all pairwise Hamming distances for a binary matrix X.
#'
#' Dependencies: hamming_bitwise from fast_dist_calc.cpp
#' @param X The N x P binary matrix
#' @return A length \eqn{{N \choose 2}} vector of pairwise Hamming distances
#' @examples
#' X <- matrix(nrow = 5, ncol = 10, rbinom(50, 1, 0.5))
#' getHammingDistance(X)
#'
getHammingDistance <- function(X) {
  # Check that each entry of X is 0 or 1
  if (!all(X == 1 | X == 0)) {
    stop("X is non-binary. Check that all entries of X are 0 or 1.")
  }

  # The number of columns of X is 1 (edge case)
  if (is.null(dim(X))) {
    return(as.vector(dist(X, "manhattan")))
  }
  # The number of columns of X is >1 and X is high-dim
  if (dim(X)[1] > 100 & dim(X)[2] > 64) {
    return(hamming_bitwise(t(X)))
  }
  # X is not too large
  return(as.vector(dist(X, "manhattan")))
}

#' A \eqn{l_p^p} Distance Vector Calculator
#'
#' Computes all pairwise \eqn{l_p^p} distances for a real matrix X,
#' for a specified choice of Minkowski norm exponent p.
#'
#' Dependencies: lp_distance from fast_dist_calc.cpp
#' @param X The N x P real matrix
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return A length \eqn{{N \choose 2}} vector of pairwise \eqn{l_p^p} distances
#' @examples
#' X <- matrix(nrow = 5, ncol = 10, rnorm(50))
#' getLpDistance(X, p = 2)
#'
getLpDistance <- function(X, p) {
  # The number of columns of X is 1 (edge case)
  if (is.null(dim(X))) {
    return(as.vector(dist(X, p = p)^p))
  }

  # The number of columns of X is >1 and X is high-dim
  if (dim(X)[1] > 100 & dim(X)[2] > 64) {
    return(lp_distance(t(X), p))
  }

  # X is not too large
  return(as.vector(dist(X, p = p)^p))
  #Rfast::Dist(X, method = "euclidean", vector = TRUE)
}

#' V Statistic for Binary Matrices
#'
#' Computes V statistic for a binary matrix X, as defined in
#' Aw, Spence and Song (2021+).
#'
#' Dependencies: getHammingDistance
#' @param X The N x P binary matrix
#' @return V(X), the variance of the pairwise Hamming distance between samples
#' @examples
#' X <- matrix(nrow = 5, ncol = 10, rbinom(50, 1, 0.5))
#' getBinVStat(X)
#'
getBinVStat <- function(X) {
  # The number of columns of X is 1 (edge case)
  if (is.null(dim(X))) {
    return(var(getHammingDistance(X)))
  }
  # V is constant when X has only two rows (edge case)
  if (dim(X)[1] == 2) {
    return(0)
  }
  # All other cases
  return(var(getHammingDistance(X)) / dim(X)[2])
}

#' V Statistic for Real Matrices
#'
#' Computes V statistic for a real matrix X,
#' where V(X) = scaled variance of \eqn{l_p^p} distances between the
#' row samples of X.
#'
#' Dependencies: getLpDistance
#' @param X The N x P real matrix
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}s
#' @return V(X), the variance of the pairwise \eqn{l_p^p} distance between samples
#' @examples
#' X <- matrix(nrow = 5, ncol = 10, rnorm(50))
#' getRealVStat(X, p = 2)
#'
getRealVStat <- function(X, p) {
  # The number of columns of X is 1 (edge case)
  if (is.null(dim(X))) {
    return(var(getLpDistance(X, p)))
  }
  # V is constant when X has only two rows (edge case)
  if (dim(X)[1] == 2) {
    return(0)
  }
  # All other cases
  return(var(getLpDistance(X, p)) / dim(X)[2])
}

#' Resampling V Statistic (Version 1)
#'
#' Generates a new array X' under the permutation null and then
#' returns the V statistic computed for X'.
#'
#' This is Version 1, which takes in the block labels. It is suitable in
#' the most general setting, where the features are grouped by labels.
#' Given original X and a list denoting labels of each feature,
#' independently permutes the rows within each block of X and returns resulting V.
#' If block labels are not specified, then features are assumed independent, which
#' is to say that block_labels is set to 1:ncol(X).
#'
#' Dependencies: getBinVStat, getRealVStat
#' @param X The N x P binary or real matrix
#' @param block_labels A vector of length P, whose pth component indicates the block membership of feature p
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return V(X'), where X' is a resampled by permutation of entries blockwise
#' @examples
#' X <- matrix(nrow = 5, ncol = 10, rnorm(50)) # real matrix example
#' naiveBlockPermute1(X, block_labels = c(1,1,2,2,3,3,4,4,5,5), p = 2) # use Euclidean distance
#'
#' X <- matrix(nrow = 5, ncol = 10, rbinom(50, 1, 0.5)) # binary matrix example
#' naiveBlockPermute1(X, block_labels = c(1,1,2,2,3,3,4,4,5,5))
#'
naiveBlockPermute1 <- function(X,
                               block_labels,
                               p = 2) {
  # Check block labels (done in the downstream function)
  # Count number of independent blocks
  num_blocks <- max(block_labels)

  # Generate permutations for each block
  perms <- do.call(cbind, lapply(1:num_blocks, function(x) sample(1:dim(X)[1], dim(X)[1])))

  # For each feature apply the permutation from that feature's block
  X_permute <- do.call(cbind, lapply(1:dim(X)[2], function(i) X[perms[, block_labels[i]], i]))

  # Compute the test statistic V for permuted data
  if (all(X == 1 | X == 0)) {
    to_return <- getBinVStat(X_permute)
  } else {
    to_return <- getRealVStat(X_permute, p = p)
  }

  # Return
  return(to_return)
}

#' Resampling V Statistic (Version 2)
#'
#' Generates a new array X' under the permutation null and then
#' returns the V statistic computed for X'.
#'
#' This is Version 2, which takes in the block boundaries. It is suitable
#' for use when the features are already arranged such that the block
#' memberships are determined by index delimiters. Given original X and
#' a list denoting labels of each feature, independently permutes the rows
#' within each block of X and returns resulting V. If block labels are not specified,
#' then features are assumed independent, which is to say that block_labels is set to 1:ncol(X).
#'
#' Dependencies: getBinVStat, getRealVStat
#' @param X The N x P binary or real matrix
#' @param block_boundaries A vector of length at most P, whose entries indicate positions at which to demarcate blocks
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return V(X'), where X' is a resampled by permutation of entries blockwise
#' @examples
#' X <- matrix(nrow = 5, ncol = 10, rnorm(50)) # real matrix example
#' naiveBlockPermute2(X, block_boundaries = c(4,7,9), p = 2) # use Euclidean distance
#'
#' X <- matrix(nrow = 5, ncol = 10, rbinom(50, 1, 0.5)) # binary matrix example
#' naiveBlockPermute2(X, block_boundaries = c(4,7,9))
#'
naiveBlockPermute2 <- function(X,
                               block_boundaries,
                               p = 2) {
  # Save the number of blocks
  num_blocks <- length(block_boundaries)
  # Pemuted version of X to be returned
  X_permute <- matrix(0,nrow = dim(X)[1], ncol = dim(X)[2])
  counter <- 1
  # Permute within each block
  while (counter < num_blocks) {
    X_permute[,block_boundaries[counter]:(block_boundaries[counter + 1] - 1)] <-
      as.matrix(X[sample(nrow(X)),block_boundaries[counter]:(block_boundaries[counter + 1] - 1)])
    counter <- counter + 1
  }
  X_permute[,block_boundaries[num_blocks]:dim(X)[2]] <-
    as.matrix(X[sample(nrow(X)),block_boundaries[num_blocks]:dim(X)[2]])

  # Compute and return the test statistic V for permuted data
  if (all(X == 1 | X == 0)) {
    to_return <- getBinVStat(X_permute)
  } else {
    to_return <- getRealVStat(X_permute, p = p)
  }

  # Return
  return(to_return)
}

#' Map from Indices to Label Pairs
#'
#' Builds a map from indexes to pairs of labels. This is
#' for caching distances, to avoid recomputing Hamming distances
#' especially when dealing with high-dimensional (large P) matrices.
#'
#' Dependencies: None
#' @param N Sample size, i.e., nrow(X)
#' @return N x N matrix whose entries record the index
#' corresponding  to the pair of labels (indexed by the matrix dims)
#' @noRd
#'
buildForward <- function(N) {
  forward <- matrix(0, ncol = N, nrow = N)
  idx <- 1
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      forward[i,j] <- idx
      forward[j,i] <- idx
      idx <- idx + 1
    }
  }
  return(forward)
}

#' Map from Label Pairs to Indices
#'
#' Builds a map from pairs of labels to indexes. This is
#' for caching distances, to avoid recomputing Hamming distances
#' especially when dealing with high-dimensional (large P) matrices.
#'
#' Dependencies: None
#' @param N Sample size, i.e., nrow(X)
#' @return N x N matrix whose entries record the index
#' corresponding  to the pair of labels (indexed by the matrix dims)
#' @noRd
#'
buildReverse <- function(N) {
  reverse <- matrix(0,ncol = 2, nrow = choose(N,2))
  idx <- 1
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      reverse[idx,1] <- i
      reverse[idx,2] <- j
      idx <- idx + 1
    }
  }
  return(reverse)
}

#' Permutation by Caching Distances
#'
#' What do you do when you have to compute pairwise distances many times, and those
#' damn distances take a long time to compute? Answer: You cache the distances and
#' permute the underlying sample labels!
#'
#' This function permutes the distances (Hamming, \eqn{l_p^p}, etc.) within blocks.
#' Permutations respect the fact that we are actually permuting the
#' underlying labels. Arguments forward and reverse should be
#' precomputed using buildForward and buildReverse.
#'
#' Dependencies: buildForward, buildReverse
#' @param dists \eqn{{N \choose 2}} by B matrix, with each column
#' containing the distances (ex: Hamming, l_p) for the block
#' @param forward N x N matrix mapping the pairs of sample labels
#' to index of the \eqn{{N \choose 2}}-length vector
#' @param reverse \eqn{{N \choose 2}} x 2 matrix mapping the index to
#' pairs of sample labels
#' @return A matrix with same dimensions as dists containing
#' the block-permuted pairwise distances
#'
cachePermute <- function(dists, forward, reverse) {
  # Permuted version of dists to be returned
  new_dists <- matrix(0, nrow = dim(dists)[1], ncol = dim(dists)[2])
  # Iterate over blocks
  for (i in 1:dim(dists)[2]) {
    # Generate a permutation of the N labels
    perm <- sample(dim(forward)[1])
    for (j in 1:dim(dists)[1]) {
      # Get labels for this index
      n1 <- reverse[j,1]; n2 <- reverse[j,2]
      # Assign the permuted labels to them
      n1 <- perm[n1]; n2 <- perm[n2]
      # Find index for permuted pair of labels
      new_j <- forward[n1,n2]
      # Assign the value to output
      new_dists[j,i] <- dists[new_j,i]
    }
  }
  return(new_dists)
}

#' Resampling Many V Statistics (Version 1)
#'
#' Generates a block permutation distribution of V statistic.
#' Precomputes distances and some indexing arrays to quickly
#' generate samples from the block permutation distribution of the V
#' statistic of X.
#'
#' This version is with block labels specified.
#'
#' Dependencies: buildForward, buildReverse, cachePermute, getHammingDistance, getLpDistance
#' @param X The binary or real matrix on which to perform
#' permutation resampling
#' @param block_labels Length P vector recording the block label of each feature
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @param nruns The resampling number (use at least 1000)
#' @return A vector of resampled values of the V statistic
#'
cacheBlockPermute1 <- function(X,
                               block_labels,
                               nruns,
                               p = 2) {
  # The number of blocks is B
  B <- max(block_labels)

  # Create the {N choose 2} x B matrix of Hamming distances
  dists <- matrix(0, nrow = choose(dim(X)[1],2), ncol = B)
  if (all(X == 1 | X == 0)) {
    # compute Hamming distance because X is binary
    for (b in 1:B) {
      dists[,b] <- getHammingDistance(X[,which(block_labels == b)])
    }
  } else {
    # else compute l_p^p distance because X is non-binary
    for (b in 1:B) {
      dists[,b] <- getLpDistance(X[,which(block_labels == b)], p)
    }
  }

  # Cache indexing arrays
  forward <- buildForward(dim(X)[1])
  reverse <- buildReverse(dim(X)[1])

  # Write a local function for computing V_vec from permuting
  newVLocal <- function(dists,forward,reverse, p) {
    partials <- cachePermute(dists, forward, reverse)
    return(var(rowSums(partials)) / dim(X)[2])
  }

  # Compute vector of V statistics
  to_return <- foreach(i=1:nruns, .combine = c) %dopar% newVLocal(dists, forward, reverse, p = 2)
  # to_return <- c()
  # for (r in 1:nruns) {
  #     partials <- cachePermuteHamming(ham_dists, forward, reverse)
  #     to_return <- c(to_return,
  #         var(rowSums(partials)) / dim(X)[2])
  # }
  return(to_return)
}

#' Resampling Many V Statistics (Version 2)
#'
#' Generates a block permutation distribution of V statistic.
#' Precomputes distances and some indexing arrays to quickly
#' generate samples from the block permutation distribution of the V
#' statistic of X.
#'
#' This version is with block boundaries specified.
#'
#' Dependencies: buildForward, buildReverse, cachePermute, getHammingDistance, getLpDistance
#' @param X The binary or real matrix on which to perform
#' permutation resampling
#' @param block_boundaries Vector denoting the positions where a new
#' block of non-independent features starts
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @param nruns The resampling number (use at least 1000)
#' @return A vector of resampled values of the V statistic
#'
cacheBlockPermute2 <- function(X,
                               block_boundaries,
                               nruns,
                               p = 2) {
  # The number of blocks is B
  B <- length(block_boundaries)

  # Create the {N choose 2} x B matrix of Hamming distances
  dists <- matrix(0, nrow = choose(dim(X)[1],2), ncol = B)
  if (all(X == 1 | X == 0)) {
    # compute Hamming distance because X is binary
    for (b in 1:(B-1)) {
      dists[,b] <- getHammingDistance(X[,block_boundaries[b]:(block_boundaries[b+1]-1)])
    }
    dists[,B] <- getHammingDistance(X[,block_boundaries[B]:dim(X)[2]])

  } else {
    # else compute l_p^p distance because X is non-binary
    for (b in 1:(B-1)) {
      dists[,b] <- getLpDistance(X[,block_boundaries[b]:(block_boundaries[b+1]-1)], p)
    }
    dists[,B] <- getLpDistance(X[,block_boundaries[B]:dim(X)[2]], p)
  }

  # Cache indexing arrays
  forward <- buildForward(dim(X)[1])
  reverse <- buildReverse(dim(X)[1])

  # Write a local function for computing V_vec from permuting
  newVLocal <- function(dists, forward, reverse) {
    partials <- cachePermute(dists, forward, reverse)
    return(var(rowSums(partials)) / dim(X)[2])
  }

  # Compute vector of V statistics
  to_return <- foreach(i=1:nruns, .combine = c) %dopar% newVLocal(dists, forward, reverse)
  # to_return <- c()
  # for (r in 1:nruns) {
  #     partials <- cachePermuteHamming(ham_dists, forward, reverse)
  #     to_return <- c(to_return,
  #         var(rowSums(partials)) / dim(X)[2])
  # }
  return(to_return)
}

#' p-value Computation for Test of Exchangeability Using Distance Data
#' 
#' Generates a block permutation p-value.
#' 
#' Generates a block permutation distribution of V statistic by storing
#' the provided list of distance data as an \eqn{{N\choose2} \times B} array,
#' and then permuting the underlying indices of each individual to generate 
#' resampled \eqn{{N\choose2} \times B} arrays. The observed V statistic is 
#' also computed from the distance data.   
#'
#' Each element of dist_list should be a \eqn{N\times N} distance matrix.
#' 
#' Dependencies: buildForward, buildReverse, cachePermute, foreach/%dopar% (from doParallel)
#' @param dist_list The list (length \eqn{B}) of pairwise distance data. 
#' Each element in list should be either a distance matrix or a table recording
#' pairwise distances. 
#' @param nruns The resampling number (use at least 1000)
#' @return The p-value obtained from comparing the empirical tail cdf of the observed 
#' V statistic computed from distance data. 
#' 
distDataPermute <- function(dist_list,
                            nruns) {
  # Get number of samples
  N <- unique(sapply(dist_list, function(x) dim(x)[1]))
  
  # Get length of distance list
  B <- length(dist_list)
  
  # Construct a (N choose 2) x B array of independent distances
  all_dist_matrix <- matrix(0, nrow = choose(N,2), ncol = B)
  # Obtain all pairwise distances for each feature
  for (b in 1:B) { # fill entries of all_dist_matrix
    all_dist_matrix[,b] <- dist_list[[b]][t(combn(colnames(dist_list[[b]]), 2))] # combing trick
  }
  
  # Compute observed V statistic from distance data
  V_obs <- var(rowSums(all_dist_matrix)) / B
  
  ##if (B < 50 & largeP) {
  #  message(date(), ": No. independent features whose distances are provided is too few for asymptotic test to be reliable.")
  #} 
  
  ## In all other cases, largeP is false, so permutation test is performed.
  #message(date(), ": Performing approximately exact resampling test with ", nruns, "permutations.")
  
  # Cache indexing arrays
  forward <- buildForward(N)
  reverse <- buildReverse(N)
  
  # Write a local function for computing V_vec from permuting
  newVLocal <- function(dists, forward, reverse) {
    partials <- cachePermute(dists, forward, reverse)
    return(var(rowSums(partials)) / B)
  }
  
  # Compute vector of V statistics
  V_vec <- foreach(i=1:nruns, .combine = c) %dopar% newVLocal(all_dist_matrix, forward, reverse)
  
  # Return
  return(mean(V_vec > V_obs)) # strictly greater than for conservativeness
}

#' p-value Computation for Test of Exchangeability with Block Dependencies
#'
#' Generates a block permutation p-value. Uses a heuristic to
#' decide whether to use distance caching or simple block permutations.
#'
#' Dependencies: buildForward, buildReverse, cachePermute,
#' cacheBlockPermute1, cacheBlockPermute2, getHammingDistance,
#' getLpDistance, naiveBlockPermute1, naiveBlockPermute2
#' @param X The binary or real matrix on which to perform
#' permutation resampling
#' @param block_boundaries Vector denoting the positions where a new
#' block of non-independent features starts. Default is NULL.
#' @param block_labels Length P vector recording the block label of each feature.
#' Default is NULL.
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @param nruns The resampling number (use at least 1000)
#' @return The block permutation p-value
#'
blockPermute <- function(X,
                         block_boundaries = NULL,
                         block_labels = NULL,
                         nruns,
                         p = 2) {
  # Check that exactly one of block_labels or block_boundaries is NULL
  if (!is.null(block_labels) & !is.null(block_boundaries)) {
    stop("Both block labels and block boundaries were specified. Exactly one of these should be specified.")
  }
  if (is.null(block_labels) & is.null(block_boundaries)) {
    stop("Neither block labels nor block boundaries was specified. Exactly one of these should be specified.")
  }

  # Observed value of V statistic from data array
  if (all(X == 1 | X == 0)) {
    token <- "binary"
    V_obs <- getBinVStat(X)
  } else {
    V_obs <- getRealVStat(X, p)
  }

  # Estimate speedup from caching (= B/P ratio)
  if (!is.null(block_labels)) {
    speedup_ratio <- max(block_labels) / dim(X)[2]
  } else {
    speedup_ratio <- length(block_boundaries) / dim(X)[2]
  }

  # If speedup is minimal, perform simple block permutations
  if (speedup_ratio > 0.01) {
    #V_vec <- sapply(1:nruns, function(x) {naiveBlockPermute(X, block_boundaries)})
    if (!is.null(block_labels)) {
      V_vec <- foreach(i=1:nruns, .combine = c) %dopar% naiveBlockPermute1(X, block_labels, p)
    } else {
      V_vec <- foreach(i=1:nruns, .combine = c) %dopar% naiveBlockPermute2(X, block_boundaries, p)
    }
  } else {
    # Else speedup is worth the overhead in caching
    if (!is.null(block_labels)) {
      V_vec <- cacheBlockPermute1(X, block_labels, nruns, p)
    } else {
      V_vec <- cacheBlockPermute2(X, block_boundaries, nruns, p)
    }
  }
  return(mean(V_vec > V_obs)) # strictly greater than for conservativeness
}

#' Tail Probability for Chi Square Convolution Random Variable
#'
#' Computes \eqn{P(X > val)} where \eqn{X = w_1 Y + w_2 Z}, where
#' \eqn{Y} is chi square distributed with \eqn{d_1} degrees of freedom,
#' \eqn{Z} is chi square distributed with \eqn{d_2} degrees of freedom,
#' and \eqn{w_1} and \eqn{w_2} are weights with \eqn{w_2} assumed positive.
#' The probability is computed using numerical integration of the
#' densities of the two chi square distributions. (Method: trapezoidal rule)
#'
#' This is used in the large P asymptotics of the permutation test.
#'
#' Dependencies: None
#' @param val observed statistic
#' @param w1 weight of first chi square rv
#' @param w2 weight of second chi square rv, assumed positive
#' @param d1 degrees of freedom of first chi square rv
#' @param d2 degrees of freedom of second chi square rv
#' @return 1 - CDF = P(X > val)
#'
weightedChi2P <- function(val, w1, w2, d1, d2){
  # Setup grid with uniform mesh size 1/10000
  grid <- seq(0, val, length.out=10000)
  probs <- pchisq((val - grid)/w1, df=d1, lower.tail=FALSE) * dchisq(grid/w2, df = d2) / w2
  # (f(x_k) + f(x_{k-1}))/2
  probs <- (probs[-1] + probs[-length(probs)]) * 0.5
  # diff(grid) = Delta(x_k) = mesh size
  return(sum(diff(grid)*probs) + pchisq(val/w2, df=d2, lower.tail=FALSE))
}

#' Covariance Computations Between Pairs of Distances (Independent Case)
#'
#' Computes covariance matrix entries and associated alpha, beta
#' and gamma quantities defined in Aw, Spence and Song (2021),
#' assuming the P features of the dataset X are independent.
#'
#' This is used in the large P asymptotics of the permutation test.
#'
#' Dependencies: buildReverse, getLpDistance
#' @param X The binary or real matrix
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return The three distinct entries of cov matrix, (alpha, beta, gamma)
#'
getCov <- function(X,
                   p = 2) {
  # Count number of features and samples
  N <- dim(X)[1]
  P <- dim(X)[2]

  # Compute for BINARY MATRIX
  if (all(X == 0 | X == 1)) {
    # Count column sums
    col_sums <- colSums(X)
    # Compute alpha, beta and gamma analytically
    alpha <- sum(sapply(col_sums, function(x) { x * (N-x) / choose(N,2) * (1 - x * (N-x) / choose(N,2))})) / P
    beta <- sum(sapply(col_sums, function(x) { x * (N-x) / (N*(N-1)) - (x * (N-x) / choose(N,2))^2 })) / P
    gamma <- sum(sapply(col_sums,
                        function(x) { 4 * x * (N-x) * (x-1) * (N-x-1) / (N*(N-1)*(N-2)*(N-3)) - (x * (N-x) / choose(N,2))^2 })) / P
  } else {
    # Compute for REAL MATRIX
    # Create the {N choose 2} x P matrix of distances
    dists <- matrix(0, nrow = choose(N,2), ncol = P)
    for (b in 1:P) {
      dists[,b] <- getLpDistance(X[,b], p)
    }
    # Compute average distance per block
    dist_ave <- colSums(dists) / choose(N,2)

    # Compute variance of d(X_1,X_2)
    alpha_vec <- colSums(dists^2) / choose(N,2) - dist_ave^2 # get alpha for each block
    alpha <- mean(alpha_vec) # average across all B blocks

    # Create an array for easy computation of other covariances
    dist_mat <- array(0, dim = c(N,N,P))
    reverse <- buildReverse(N) # build reverse map
    for (b in 1:P) {
      for (i in 1:choose(N,2)) {
        idx_pair <- reverse[i,]
        dist_mat[idx_pair[1],idx_pair[2],b] <- dists[i,b]
        dist_mat[idx_pair[2],idx_pair[1],b] <- dists[i,b]
      }
    }

    # Compute E[d(X_1,X_2) * d(X_1,X_3)]
    beta_vec <- c()
    for (b in 1:P) {
      beta_vec <- c(beta_vec, sum((rowSums(dist_mat[,,b]) - dist_mat[,,b]) * dist_mat[,,b]) / N / (N-1) / (N-2))
    }
    # Compute covariance of d(X_1,X_2) and d(X_1,X_3)
    beta_vec <- beta_vec - dist_ave^2
    beta <- mean(beta_vec)

    # Compute E[d(X_1,X_2) * d(X_3,X_4)]
    gamma_vec <- c()
    for (b in 1:P) {
      # Some arithmetic reasoning needed here: note R uses column-major filling in of entries
      nu <- matrix(dist_ave[b] * N * (N-1), nrow = N, ncol = N) # compute sum_i sum_j d(X_i,X_j)
      nu_row <- matrix(2 * rowSums(dist_mat[,,b]), nrow = N, ncol = N)
      nu_col <- t(matrix(2 * colSums(dist_mat[,,b]), nrow = N, ncol = N))
      gamma <- sum((nu - nu_row - nu_col + 2 * dist_mat[,,b]) * dist_mat[,,b]) / (N * (N-1) * (N-2) * (N-3))
      gamma_vec <- c(gamma_vec, gamma)
    }
    # Compute covariance of d(X_1,X_2) and d(X_3,X_4)
    gamma_vec <- gamma_vec -  dist_ave^2
    gamma <- mean(gamma_vec)
  }

  # Create and annotate the output vector to be returned
  to_return <- c(alpha, beta, gamma)
  names(to_return) <- c("alpha", "beta", "gamma")
  return(to_return)
}

#' Get Chi Square Weights
#'
#' Computes weights for the asymptotic random variable
#' from the alpha, beta and gamma computed of data array X.
#'
#' This is used in the large P asymptotics of the permutation test.
#'
#' Dependencies: None
#' @param alpha covariance matrix entry computed from getCov
#' @param beta covariance matrix entry computed from getCov
#' @param gamma covariance matrix entry computed from getCov
#' @param N The sample size, i.e., nrow(X) where X is the original dataset
#' @return The weights (w1, w2)
#'
getChi2Weights <- function(alpha, beta, gamma, N) {
  # Compute the weights
  w1 <- (alpha + (N-4) * beta - (N-3)*gamma) / choose(N,2)
  w2 <- (alpha - 2*beta + gamma) / choose(N,2)
  # Create and annotate the output vector to be returned
  to_return <- c(w1,w2)
  names(to_return) <- c("w1","w2")
  return(to_return)
}

#' Approximate p-value for Test of Exchangeability (Assuming Large P)
#'
#' Computes the large P asymptotic p-value for dataset X,
#' assuming its P features are independent.
#'
#' This is the large P asymptotics of the permutation test.
#'
#' Dependencies: getBinVStat, getRealVStat, getChi2Weights, weightedChi2P, getCov
#' @param X The binary or real matrix on which to perform test of exchangeability
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return The asymptotic p-value
#'
indLargeP <- function(X,
                      p = 2) {
  N <- dim(X)[1]
  greeks <- getCov(X, p)
  weights <- getChi2Weights(greeks[1],greeks[2],greeks[3],N)
  d1 <- N - 1
  d2 <- choose(N-1,2) - 1
  if (all(X == 0 | X == 1)) {
    V_obs <- getBinVStat(X)
  } else {
    V_obs <- getRealVStat(X, p)
  }
  to_return <- as.numeric(weightedChi2P(V_obs,weights[1],weights[2],d1,d2))
  return(to_return)
}

#' Approximate p-value for Test of Exchangeability (Assuming Large N and P)
#'
#' Computes the large (N,P) asymptotic p-value for dataset X,
#' assuming its P features are independent
#'
#' This is the large N and large P asymptotics of the permutation test.
#'
#' Dependencies: getBinVStat, getRealVStat, getCov, getChi2Weights
#' @param X The binary or real matrix on which to perform test of exchangeability
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return The asymptotic p-value
#'
indGaussian <- function(X,
                        p = 2) {
  N <- dim(X)[1]
  greeks <- getCov(X, p)
  weights <- getChi2Weights(greeks[1],greeks[2],greeks[3],N)
  d1 <- N - 1
  d2 <- choose(N-1,2) - 1
  if (all(X == 0 | X == 1)) {
    V_obs <- getBinVStat(X)
  } else {
    V_obs <- getRealVStat(X, p)
  }
  to_return <- 1 - pnorm(V_obs,
                         mean = weights[1]*d1 + weights[2]*d2,
                         sd = sqrt(weights[1]^2 * 2 * d1 + weights[2]^2 * 2 * d2))
  return(to_return)
}

#' Covariance Computations Between Pairs of Distances (Block Dependencies Case)
#'
#' Computes covariance matrix entries and associated alpha, beta
#' and gamma quantities defined in Aw, Spence and Song (2021),
#' for partitionable features that are grouped into blocks. Uses
#' precomputation to compute the unique entries of the asymptotic
#' covariance matrix of the pairwise Hamming distances in O(N^2) time.
#'
#' This is used in the large P asymptotics of the permutation test.
#'
#' Dependencies: buildReverse, getHammingDistance, getLpDistance
#' @param X The binary or real matrix
#' @param block_boundaries Vector denoting the positions where a new
#' block of non-independent features starts.
#' @param block_labels Length P vector recording the block label of each feature.
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return The three distinct entries of cov matrix, (alpha, beta, gamma)
#'
getBlockCov <- function(X,
                        block_boundaries,
                        block_labels,
                        p = 2) {
  # Dimensions of the problem
  N <- dim(X)[1]
  P <- dim(X)[2]
  if (!is.null(block_boundaries)) {
    B <- length(block_boundaries)
  } else {
    B <- max(block_labels)
  }

  # Create the {N choose 2} x B matrix of distances
  dists <- matrix(0, nrow = choose(N,2), ncol = B)

  # If BINARY MATRIX, fill with Hamming distances
  if (all(X == 0 | X == 1)) {
    if (!is.null(block_boundaries)) {
      for (b in 1:(B-1)) {
        dists[,b] <- getHammingDistance(X[,block_boundaries[b]:(block_boundaries[b+1]-1)])
      }
      dists[,B] <- getHammingDistance(X[,block_boundaries[B]:dim(X)[2]])
    } else {
      for (b in 1:B) {
        dists[,b] <- getHammingDistance(X[,which(block_labels == b)])
      }
    }
  } else {
    # Or else, fill with Lp distances
    if (!is.null(block_boundaries)) {
      for (b in 1:(B-1)) {
        dists[,b] <- getLpDistance(X[,block_boundaries[b]:(block_boundaries[b+1]-1)], p)
      }
      dists[,B] <- getLpDistance(X[,block_boundaries[B]:dim(X)[2]], p)
    } else {
      for (b in 1:B) {
        dists[,b] <- getLpDistance(X[,which(block_labels == b)], p)
      }
    }
  }

  # Compute average distance per block
  dist_ave <- colSums(dists) / choose(N,2)

  # Compute variance of d(X_1,X_2)
  alpha_vec <- colSums(dists^2) / choose(N,2) - dist_ave^2 # get alpha for each block
  alpha <- mean(alpha_vec) # average across all B blocks

  # Create an array for easy computation of other covariances
  dist_mat <- array(0, dim = c(N,N,B))
  reverse <- buildReverse(N) # build reverse map
  for (b in 1:B) {
    for (i in 1:choose(N,2)) {
      idx_pair <- reverse[i,]
      dist_mat[idx_pair[1],idx_pair[2],b] <- dists[i,b]
      dist_mat[idx_pair[2],idx_pair[1],b] <- dists[i,b]
    }
  }

  # Compute E[d(X_1,X_2) * d(X_1,X_3)]
  beta_vec <- c()
  for (b in 1:B) {
    beta_vec <- c(beta_vec, sum((rowSums(dist_mat[,,b]) - dist_mat[,,b]) * dist_mat[,,b]) / N / (N-1) / (N-2))
  }
  # Compute covariance of d(X_1,X_2) and d(X_1,X_3)
  beta_vec <- beta_vec - dist_ave^2
  beta <- mean(beta_vec)

  # Compute E[d(X_1,X_2) * d(X_3,X_4)]
  gamma_vec <- c()
  for (b in 1:B) {
    # Some arithmetic reasoning needed here: note R uses column-major filling in of entries
    nu <- matrix(dist_ave[b] * N * (N-1), nrow = N, ncol = N) # compute sum_i sum_j d(X_i,X_j)
    nu_row <- matrix(2 * rowSums(dist_mat[,,b]), nrow = N, ncol = N)
    nu_col <- t(matrix(2 * colSums(dist_mat[,,b]), nrow = N, ncol = N))
    gamma <- sum((nu - nu_row - nu_col + 2 * dist_mat[,,b]) * dist_mat[,,b]) / (N * (N-1) * (N-2) * (N-3))
    gamma_vec <- c(gamma_vec, gamma)
  }
  # Compute covariance of d(X_1,X_2) and d(X_3,X_4)
  gamma_vec <- gamma_vec -  dist_ave^2
  gamma <- mean(gamma_vec)
  # Create and annotate the output vector to be returned
  to_return <- c(alpha,beta,gamma) * B / P
  names(to_return) <- c("alpha","beta","gamma")
  return(to_return)
}

#' Approximate p-value for Test of Exchangeability (Assuming Large P with Block Dependencies)
#'
#' Computes the large P asymptotic p-value for dataset X,
#' assuming its P features are independent within specified blocks.
#'
#' This is the large P asymptotics of the permutation test.
#'
#' Dependencies: getBinVStat, getRealVStat, getChi2Weights, weightedChi2P, getBlockCov
#' @param X The binary or real matrix on which to perform test of exchangeability
#' @param block_boundaries Vector denoting the positions where a new
#' block of non-independent features starts.
#' @param block_labels Length P vector recording the block label of each feature.
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return The asymptotic p-value
#'
blockLargeP <- function(X,
                        block_boundaries,
                        block_labels,
                        p = 2) {
  N <- dim(X)[1]
  greeks <- getBlockCov(X,
                        block_boundaries,
                        block_labels,
                        p)
  weights <- getChi2Weights(greeks[1],greeks[2],greeks[3],N)
  d1 <- N - 1
  d2 <- choose(N-1,2) - 1
  if (all(X == 0 | X == 1)) {
    V_obs <- getBinVStat(X)
  } else {
    V_obs <- getRealVStat(X, p)
  }
  to_return <- as.numeric(weightedChi2P(V_obs,weights[1],weights[2],d1,d2))
  return(to_return)
}

#' Approximate p-value for Test of Exchangeability (Assuming Large N and P with Block Dependencies)
#'
#' Computes the large (N,P) asymptotic p-value for dataset X,
#' assuming its P features are independent within specified blocks.
#'
#' This is the large N and large P asymptotics of the permutation test.
#'
#' Dependencies: getBinVStat, getRealVStat, getBlockCov, getChi2Weights
#' @param X The binary or real matrix on which to perform test of exchangeability
#' @param block_boundaries Vector denoting the positions where a new
#' block of non-independent features starts.
#' @param block_labels Length P vector recording the block label of each feature.
#' @param p The power p of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}
#' @return The asymptotic p-value
#'
blockGaussian <- function(X,
                          block_boundaries,
                          block_labels,
                          p) {
  N <- dim(X)[1]
  greeks <- getBlockCov(X, block_boundaries, block_labels, p)
  weights <- getChi2Weights(greeks[1],greeks[2],greeks[3],N)
  d1 <- N - 1
  d2 <- choose(N-1,2) - 1
  if (all(X == 0 | X == 1)) {
    V_obs <- getBinVStat(X)
  } else {
    V_obs <- getRealVStat(X, p)
  }
  to_return <- 1 - pnorm(V_obs,
                         mean = weights[1]*d1 + weights[2]*d2,
                         sd = sqrt(weights[1]^2 * 2 * d1 + weights[2]^2 * 2 * d2))
  return(to_return)
}

#' Asymptotic p-value of Exchangeability Using Distance Data
#' 
#' Generates an asymptotic p-value.
#' 
#' Generates a weighted convolution of chi-squares distribution of V statistic 
#' by storing the provided list of distance data as an \eqn{{N\choose2} \times B} array,
#' and then using large-P theory to generate the asymptotic null distribution 
#' against which the p-value of observed V statistic is computed. 
#' 
#' Each element of dist_list should be a \eqn{N\times N} distance matrix.
#' 
#' Dependencies: buildReverse, getChi2Weights, weightedChi2P
#' @param dist_list The list (length \eqn{B}) of pairwise distance data. 
#' Each element in list should be either a distance matrix or a table recording
#' pairwise distances. 
#' @return The asymptotic p-value obtained from the weighted convolution of chi-squares
#' distribution.
#' 
distDataLargeP <- function(dist_list) {
  # Check that all matrices have the same dimension 
  assertthat::assert_that(var(sapply(dist_list, function(x) dim(x)[1])) == 0,
                          msg = "Not all matrices have the same dimension. Check distance matrices provided.")
  
  # Get number of samples
  N <- unique(sapply(dist_list, function(x) dim(x)[1]))
  
  # Get length of distance list
  B <- length(dist_list)
  
  # Construct a (N choose 2) x B array of independent distances
  all_dist_matrix <- matrix(0, nrow = choose(N,2), ncol = B)
  
  # Obtain all pairwise distances for each feature
  for (b in 1:B) { # fill entries of all_dist_matrix
    all_dist_matrix[,b] <- dist_list[[b]][t(combn(colnames(dist_list[[b]]), 2))] # combing trick
  }
  
  # Compute observed V statistic from distance data
  V_obs <- var(rowSums(all_dist_matrix)) / B
  
  # Compute average distance per block
  dist_ave <- colSums(all_dist_matrix) / choose(N,2)
  
  # Compute variance of d(X_1,X_2)
  alpha_vec <- colSums(all_dist_matrix^2) / choose(N,2) - dist_ave^2 # get alpha for each block
  alpha <- mean(alpha_vec) # average across all B blocks
  
  # Create a tensor for easy computation of other covariances
  dist_tensor <- array(0, dim = c(N,N,B))
  reverse <- buildReverse(N) # build reverse map
  for (b in 1:B) {
    for (i in 1:choose(N,2)) {
      idx_pair <- reverse[i,]
      dist_tensor[idx_pair[1], idx_pair[2],b] <- all_dist_matrix[i,b]
      dist_tensor[idx_pair[2], idx_pair[1],b] <- all_dist_matrix[i,b]
    }
  }
  
  # Compute E[d(X_1,X_2) * d(X_1,X_3)]
  beta_vec <- c()
  for (b in 1:B) {
    beta_vec <- c(beta_vec, sum((rowSums(dist_tensor[,,b]) - dist_tensor[,,b]) * dist_tensor[,,b]) / N / (N-1) / (N-2))
  }
  
  # Compute covariance of d(X_1,X_2) and d(X_1,X_3)
  beta_vec <- beta_vec - dist_ave^2
  beta <- mean(beta_vec)
  
  # Compute E[d(X_1,X_2) * d(X_3,X_4)]
  gamma_vec <- c()
  for (b in 1:B) {
    # Some arithmetic reasoning needed here: note R uses column-major filling in of entries
    nu <- matrix(dist_ave[b] * N * (N-1), nrow = N, ncol = N) # compute sum_i sum_j d(X_i,X_j)
    nu_row <- matrix(2 * rowSums(dist_tensor[,,b]), nrow = N, ncol = N)
    nu_col <- t(matrix(2 * colSums(dist_tensor[,,b]), nrow = N, ncol = N))
    gamma <- sum((nu - nu_row - nu_col + 2 * dist_tensor[,,b]) * dist_tensor[,,b]) / (N * (N-1) * (N-2) * (N-3))
    gamma_vec <- c(gamma_vec, gamma)
  }
  
  # Compute covariance of d(X_1,X_2) and d(X_3,X_4)
  gamma_vec <- gamma_vec -  dist_ave^2
  gamma <- mean(gamma_vec)
  
  # Get weights for chi-square distribution
  weights <- getChi2Weights(alpha, beta, gamma, N)
  
  # Compute degrees of freedom
  d1 <- N - 1
  d2 <- choose(N-1,2) - 1
  
  to_return <- as.numeric(weightedChi2P(V_obs, weights[1], weights[2], d1, d2))
  
  # Return
  return(to_return)
}
