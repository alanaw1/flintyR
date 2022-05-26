#' A Non-parametric Test for Exchangeability and Homogeneity
#'
#' Computes the p-value of a multivariate dataset \eqn{\mathbf{X}}, which
#' informs the user if the sample is exchangeable at a given
#' significance level, while simultaneously accounting for feature
#' dependencies. See Aw, Spence and Song (2021) for details.
#'
#' Automatically detects if dataset is binary, and runs the Hamming
#' distance version of test if so. Otherwise, computes the squared
#' Euclidean distance between samples and evaluates whether the
#' variance of Euclidean distances, \eqn{V}, is atypically large under the
#' null hypothesis of exchangeability. Note the user may tweak the
#' choice of power \eqn{p} if they prefer an \eqn{l_p^p} distance other than Euclidean.
#'
#' Under the hood, the variance statistic, \eqn{V}, is computed efficiently.
#' Moreover, the user can specify their choice of block permutations,
#' large \eqn{P} asymptotics, or large \eqn{P} and large \eqn{N} asymptotics. The latter two
#' return reasonably accurate p-values for moderately large dimensionalities.
#'
#' User recommendations: When the number of independent blocks \eqn{B} or number of
#' independent features \eqn{P} is at least 50, it is safe to use large \eqn{P} asymptotics.
#' If \eqn{P} or \eqn{B} is small, however, stick with permutations.
#'
#' Dependencies: All functions in auxiliary.R
#' @param X The binary or real matrix on which to perform test of exchangeability.
#' @param block_boundaries Vector denoting the positions where a new
#' block of non-independent features starts. Default is NULL.
#' @param block_labels Length \eqn{P} vector recording the block label of each feature.
#' Default is NULL.
#' @param largeP Boolean indicating whether to use large \eqn{P} asymptotics. Default is FALSE.
#' @param largeN Boolean indicating whether to use large \eqn{N} asymptotics. Default is FALSE.
#' @param nruns Resampling number for exact test. Default is 5000.
#' @param type Either an unbiased estimate of (`'unbiased'`, default), or valid, but biased estimate of, (`'valid'`) p-value (see Hemerik and Goeman, 2018), or both (`'both'`). Default is `'unbiased'`.
#' @param p The power \eqn{p} of \eqn{l_p^p}, i.e., \eqn{||x||_p^p = (x_1^p+...x_n^p)}. Default is 2.
#' @return The p-value to be used to test the null hypothesis of exchangeability.
#' @export
#' @importFrom doParallel %dopar%
#' @examples
#' # Example 1 (get p-value of small matrix with independent features using exact test)
#' suppressWarnings(require(doParallel))
#' # registerDoParallel(cores = 2)
#'
#' X1 <- matrix(nrow = 5, ncol = 10, rbinom(50, 1, 0.5)) # binary matrix, small
#' getPValue(X1) # perform exact test with 5000 permutations
#'
#' # should be larger than 0.05
#'
#' # Example 2 (get p-value of high-dim matrix with independent features using asymptotic test)
#' X2 <- matrix(nrow = 10, ncol = 1000, rnorm(1e4)) # real matrix, large enough
#' getPValue(X2, p = 2, largeP = TRUE) # very fast
#'
#' # should be larger than 0.05
#' # getPValue(X2, p = 2) # slower, do not run (Output: 0.5764)
#'
#' # Example 3 (get p-value of high-dim matrix with partitionable features using exact test)
#'
#' X3 <- matrix(nrow = 10, ncol = 1000, rbinom(1e4, 1, 0.5))
#' getPValue(X3, block_labels = rep(c(1,2,3,4,5), 200))
#'
#' # Warning message: # there are features that have zero variation (i.e., all 0s or 1s)
#' # In getPValue(X3, block_labels = rep(c(1, 2, 3, 4, 5), 200)) :
#' # There exist columns with all ones or all zeros for binary X.
#'
#' # Example 4 (get p-value of high-dim matrix with partitionable features using asymptotic test)
#'
#' ## This elaborate example generates binarized versions of time series data.
#'
#'# Helper function to binarize a marker
#'# by converting z-scores to {0,1} based on
#'# standard normal quantiles
#'binarizeMarker <- function(x, freq, ploidy) {
#'  if (ploidy == 1) {
#'    return((x > qnorm(1-freq)) + 0)
#'  } else if (ploidy == 2) {
#'    if (x <= qnorm((1-freq)^2)) {
#'      return(0)
#'    } else if (x <= qnorm(1-freq^2)) {
#'      return(1)
#'    } else return(2)
#'  } else {
#'    cat("Specify valid ploidy number, 1 or 2")
#'  }
#'}
#'
#' getAutoRegArray <- function(B, N, maf_l = 0.38, maf_u = 0.5, rho = 0.5, ploid = 1) {
#' # get minor allele frequencies by sampling from uniform
#' mafs <- runif(B, min = maf_l, max = maf_u)
#' # get AR array
#' ar_array <- t(replicate(N, arima.sim(n = B, list(ar=rho))))
#' # theoretical column variance
#' column_var <- 1/(1-rho^2)
#' # rescale so that variance per marker is 1
#' ar_array <- ar_array / sqrt(column_var)
#' # rescale each column of AR array
#' for (b in 1:B) {
#'   ar_array[,b] <- sapply(ar_array[,b],
#'                          binarizeMarker,
#'                          freq = mafs[b],
#'                          ploidy = ploid)
#' }
#' return(ar_array)
#' }
#'
#' ## Function to generate the data array with desired number of samples
#' getExHaplotypes <- function(N) {
#'   array <- do.call("cbind",
#'                    lapply(1:50, function(x) {getAutoRegArray(N, B = 20)}))
#'   return(array)
#' }
#'
#' ##  Generate data and run test
#' X4 <- getExHaplotypes(10)
#' getPValue(X4, block_boundaries = seq(from = 1, to = 1000, by = 25), largeP = TRUE)
#'
#' # stopImplicitCluster()
#'
getPValue <- function(X,
                      block_boundaries = NULL,
                      block_labels = NULL,
                      largeP = FALSE,
                      largeN = FALSE,
                      nruns = 5000,
                      type = 'unbiased',
                      p = 2) {
  # [!] Check data matrix X is a matrix
  assertthat::assert_that(is.matrix(X),
                          msg = "X must be a 2D matrix, see ?getPValue for examples. Try as.matrix(X).")

  # [!] Warning if there exist non-varying features for binary matrix
  if (all(X == 0 | X == 1)) {
    col_sums <- colSums(X)
    if (any(col_sums == 0) | any(col_sums == dim(X)[1])) {
      warning("There exist columns with all ones or all zeros for binary X.")
    }
  }

  # [!] Check that type is well-defined
  assertthat::assert_that((type == "unbiased" | type == "valid" | type == "both"),
                          msg = "Please specify a valid type (`unbiased`, `valid`, or `both`)")

  # [!] Check classes of largeP, largeN and nruns
  assertthat::assert_that(is.logical(largeP),
                          msg = "largeP must be either TRUE or FALSE.")
  assertthat::assert_that(is.logical(largeN),
                          msg = "largeN must be either TRUE or FALSE.")
  assertthat::assert_that(isTRUE(all.equal(nruns, as.integer(nruns))) & is.numeric(nruns) & nruns >= 1,
                          msg = "nruns must be a positive integer.")

  # [!] Check that large N only asymptotics are not input by user
  assertthat::assert_that(largeP || !largeN,
                          msg = "Large P and large N asymptotics not yet implemented. Pick another version.")

  # [!] Check that not both block_boundaries and block_labels are specified
  assertthat::assert_that(is.null(block_boundaries) || is.null(block_labels),
                          msg = "Block boundaries and block labels are specified simultaneously. Specify at most one.")

  # CASE 1: Assume features are independent if neither block labels nor block boundaries is specified
  if (is.null(block_boundaries) & is.null(block_labels)) {
    block_boundaries = 1:dim(X)[2]
  }

  # CASE 2: Block boundaries specified -- perform block permutations with block boundaries
  if (!is.null(block_boundaries) & is.null(block_labels)) {
    # [!] Check block boundary sequence is valid
    assertthat::assert_that(max(block_boundaries) <= dim(X)[2],
                            msg = "Block boundaries exceed number of columns. Check block boundaries.")

    # [!] Check block boundary sequence is monotone increasing
    assertthat::assert_that(all(block_boundaries == cummax(block_boundaries)),
                            msg = "Block boundaries not monotone increasing. Check block boundaries.")

    # [!] Check block boundary sequence has no repeats
    assertthat::assert_that(!any(diff(block_boundaries) == 0),
                            msg = "Block boundaries have repeats. Check block boundaries.")
    # Add 1 to the vector of positions for downstream permutation
    # functionality if it's not there already
    if (block_boundaries[1] != 1) {
      block_boundaries <- c(1, block_boundaries)
    }
  }

  # CASE 3: Block labels specified -- perform block permutations with block labels
  if (!is.null(block_labels)) {
    # [!] Check block label sequence is valid
    assertthat::assert_that(max(block_labels) <= dim(X)[2],
                            msg = "Number of blocks exceeds number of columns. Check block labels.")

    # [!] Check for non-conventional labeling
    assertthat::assert_that(max(block_labels) == length(unique(block_labels)),
                            msg = "Block labels are not from 1 to B. Please relabel blocks using this convention.")

    # EDGE CASE: a single block is specified -- sample should automatically be exchangeable
    if (max(block_labels) == 1) {
      cat("All blocks are labeled 1, i.e., no independent sets of features detected, so samples are assumed exchangeable.\n")
      return(1)
    }
  }

  # Perform the test corresponding to largeP and largeN switches
  if (largeP & largeN & (identical(block_boundaries, 1:dim(X)[2]) | identical(block_labels, 1:dim(X)[2]))) {
    return(indGaussian(X, p)) # large P large N ind features
  }
  if (largeP & largeN & !(identical(block_boundaries, 1:dim(X)[2]) | identical(block_labels, 1:dim(X)[2]))) {
    return(blockGaussian(X, block_boundaries, block_labels, p)) # large P large N block
  }
  if (largeP & !largeN & (identical(block_boundaries, 1:dim(X)[2]) | identical(block_labels, 1:dim(X)[2]))) {
    return(indLargeP(X, p)) # large P ind features
  }
  if (largeP & !largeN & !(identical(block_boundaries, 1:dim(X)[2]) | identical(block_labels, 1:dim(X)[2]))) {
    return(blockLargeP(X, block_boundaries, block_labels, p)) # large P block
  }
  if (!largeP & !largeN) {
    return(blockPermute(X, block_boundaries, block_labels, nruns, type, p))
  }
}


#' A Non-parametric Test for Exchangeability and Homogeneity (Distance List Version)
#'
#' Computes the p-value of a multivariate dataset, which
#' informs the user if the sample is exchangeable at a given
#' significance level, while simultaneously accounting for feature
#' dependencies. See Aw, Spence and Song (2021) for details.
#'
#' This version takes in a list of distance matrices recording
#' pairwise distances between individuals across \eqn{B} independent features.
#'
#' Dependencies: distDataLargeP and distDataPermute from auxiliary.R
#' @param dist_list The list of distances.
#' @param largeP Boolean indicating whether to use large \eqn{P} asymptotics. Default is FALSE.
#' @param nruns Resampling number for exact test. Default is 1000.
#' @param type Either an unbiased estimate of (`'unbiased'`, default), or valid, but biased estimate of, (`'valid'`) p-value (see Hemerik and Goeman, 2018), or both (`'both'`). Default is `'unbiased'`.
#' @return The p-value to be used to test the null hypothesis of exchangeability.
#' @export
#' @importFrom doParallel %dopar%
#'
distDataPValue <- function(dist_list,
                           largeP = FALSE,
                           nruns = 1000,
                           type = 'unbiased') {

  # [!] Check that dist_list is a list
  assertthat::assert_that(is.list(dist_list),
                          msg = "Input is not a valid list of distance matrices.")

  # [!] Check that each element of dist_list is a distance matrix
  assertthat::assert_that(all(sapply(dist_list, function(x) is.matrix(x))),
                          msg = "Not all elements of dist_list is a distance matrix. Check list provided.")

  # [!] Check all distance matrices have the same dimension
  assertthat::assert_that(var(sapply(dist_list, function(x) dim(x)[1])) == 0,
                          msg = "Not all matrices have the same dimension. Check distance matrices provided.")

  # [!] Check all distance matrices are square
  assertthat::assert_that(all(sapply(dist_list, function(x) {dim(x)[1] == dim(x)[2]})),
                          msg = "Not all matrices are square. Check distance matrices provided.")

  # [!] Check that largeP and no. independent features are consistent
  if (length(dist_list) < 50 & largeP) {
    stop("Too few independent distance matrices for chi-square approximation to be valid, please set largeP = FALSE.")
  }

  # Get p-value
  if (largeP) {
    return(distDataLargeP(dist_list))
  } else {
    return(distDataPermute(dist_list, nruns, type))
  }
}
