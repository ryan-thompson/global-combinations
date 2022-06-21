# Description: Contains functions that implement forecast combination methods from the paper
# "Global combinations of expert forecasts"
# Author: Ryan Thompson

#--------------------------------------------------------------------------------------------------#
# cv.combine function for cross-validated forecast combination
#--------------------------------------------------------------------------------------------------#

# Performs leave-one-out cross-validation for the globalisation and shrinkage parameters of the 
# combine function

#--- Input ---#
# errors: a list of forecast error matrices where each column corresponds to a forecaster; each 
# matrix should have the same number of columns

#--- Output ---#
# A matrix of cross-validated weights where the kth column is a vector of forecast combination 
# weights for the kth task

cv.combine <- \(errors, ...) {
  
  # Save data dimensions
  n <- min(sapply(errors, nrow))
  m <- length(errors)
  
  # Fit weights to full data
  fit <- combine(errors, ...)
  ncv <- length(fit)
  
  # Determine best gamma
  if (ncv != 1) {
    
    # Cross-validate gamma
    cv.folds <- lapply(1:n, \(i) {
      
      # Set up training and validation sets
      train <- lapply(errors, \(x) x[- i, , drop = F])
      valid <- lapply(errors, \(x) x[i, , drop = F])
      
      # Remove any forecasters unavailable in validation set
      keep <- !apply(do.call(rbind, valid), 2, anyNA)
      train <- lapply(train, \(x) x[, keep, drop = F])
      valid <- lapply(valid, \(x) x[keep, drop = F])
      
      # Fit on training set and evaluate on validation set
      fit.i <- combine(train, ...)
      sapply(1:ncv, \(j) sapply(1:m, \(k) (valid[[k]] %*% fit.i[[j]][, k]) ^ 2))
      
    })
    
    cv <- Reduce('+', cv.folds)
    sapply(1:m, \(k) fit[[which.min(cv[k, ])]][, k])

  } else {
    
    fit[[1]]
    
  }
  
}

#--------------------------------------------------------------------------------------------------#
# combine function for forecast combination
#--------------------------------------------------------------------------------------------------#

# Fits the forecast combination weights for a grid of globalisation and shrinkage parameters. The
# globalisation parameter is gamma and the shrinkage parameter is lambda

#--- Input ---#
# errors: a list of forecast error matrices where each column corresponds to a forecaster; each 
# matrix should have the same number of columns
# scheme: the type of forecast combination scheme to apply; 'equal' for equal weights, 'optimal' for
# optimal weights, and so on
# group: a vector of integers indicating the grouping of the tasks; by default all tasks are in the 
# same group
# sigma: an optional list of covariance matrices; if supplied errors will be ignored otherwise the 
# covariance matrices will be computed from errors
# scale: a logical indicating whether to scale the tasks to put them on the same footing
# ngamma: the number of globalisation parameters to evaluate
# gamma.min: the minimum value for gamma
# gamma.max: the maximum value for gamma
# gamma: an optional sequence of globalisation parameters
# nlambda: the number of shrinkage parameters to evaluate
# lambda.min: the minimum value for lambda
# lambda.max: the maximum value for lambda
# lambda: an optional sequence of shrinkage parameters
# correct: a logical indicating whether to perform a positive definite correction on any nonpositive 
# definite covariance matrices

#--- Output ---#
# A list of matrices of weights; the kth column of each matrix is a vector of forecast combination 
# weights for the kth task

combine <- \(errors, scheme = c('equal', 'optimal', 'optimal convex', 'optimal equal'),
             group = rep(1, length(sigma)), sigma = NULL, scale = T,
             ngamma = 10, gamma.min = 1e-3, gamma.max = 1e3, gamma = NULL,
             nlambda = 10, lambda.min = 1e-3, lambda.max = 1e3, lambda = NULL,
             correct = F, ...) {
  
  # Check arguments
  scheme <- match.arg(scheme)
  
  # Estimate covariance matrices
  if (is.null(sigma)) sigma <- lapply(1:length(errors), \(k) covariance(errors[[k]], correct))
  
  # Save dimensions
  p <- ncol(sigma[[1]])
  m <- length(sigma)
  
  # Construct gamma sequence
  if (is.null(gamma)) {
    if (scheme == 'equal') {
      gamma <- 0
    } else {
      gamma <- exp(seq(log(gamma.max), log(gamma.min), length.out = ngamma))
    }
  } else {
    ngamma <- length(gamma)
  }
  
  # Construct lambda sequence
  if (is.null(lambda)) {
    if (scheme == 'equal') {
      lambda <- 0
    } else {
      lambda <- exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda))
    }
  } else {
    nlambda <- length(lambda)
  }
  
  # Estimate scaling factor
  scale.factor <- matrix(1, m, nlambda)
  if (scale & scheme != 'equal') {
    fit <- optimise(sigma, scheme, group, 0, lambda, scale.factor, ...)
    scale.factor <- sapply(1:nlambda, \(i) sapply(1:m, \(k) fit[[i]][, k] %*%
                                                   (sigma[[k]] + diag(lambda[i], p)) %*%
                                                   fit[[i]][, k]))
  }
  
  # Perform fit
  optimise(sigma, scheme, group, gamma, lambda, scale.factor, ...)
  
}

optimise <- \(sigma, scheme, group, gamma, lambda, scale.factor,
              params = list(OutputFlag = 0, NonConvex = 2, Threads = 1, TimeLimit = 60)) {

  # Save dimensions
  p <- ncol(sigma[[1]])
  m <- length(sigma)
  g <- length(unique(group))
  
  # Set up space for results
  weights <- list()
  pm <- p * m
  pg <- p * g
  nvar <- pm + pm + m + pg
  
  for (i1 in 1:length(gamma)) {
    for (i2 in 1:length(lambda)) {
      
      # Scale sigma
      sigma.scaled <- lapply(1:m, \(k) (sigma[[k]] + diag(lambda[i2], p)) / scale.factor[k, i2])
      
      model <- list()
      
      # Set up variable schemes
      model$vtype <- c(
        rep('C', pm), # Weight vectors (continuous)
        rep('B', pm), # Selection vectors (binary)
        rep('I', m),  # Numbers of selected forecasts (integer)
        rep('C', pg)  # Auxiliary weight vector (continuous)
      )
      model$lb <- c(
        rep(ifelse(scheme == 'equal', 1 / p, ifelse(scheme == 'optimal', - Inf, 0)), pm),
        rep(0, pm), 
        rep(1, m), 
        rep(- Inf, pg)
      )
      model$ub <- c(
        rep(ifelse(scheme == 'equal', 1 / p, ifelse(scheme == 'optimal', Inf, 1)), pm),
        rep(1, pm), 
        rep(p, m), 
        rep(Inf, pg)
      )
      
      # Set up objective matrix
      model$Q <- Matrix::bdiag(
        Matrix::bdiag(sigma.scaled) + Matrix::Diagonal(pm, gamma[i1]),
        Matrix::Matrix(0, pm, pm),
        Matrix::Matrix(0, m, m),
        Matrix::Diagonal(pg, rep(table(group) * gamma[i1], each = p))
      )
      for (k in 1:m) model$Q[cbind(1:p + p * (k - 1), 1:p + pm + pm + m + p * (group[k] - 1))] <- 
        - 2 * gamma[i1]
      
      # Set up linear constraints
      model$A <- Matrix::sparseMatrix(rep(1:m, each = p), 1:pm, x = 1, dims = c(m, nvar))
      model$sense <- rep('=', each = m)
      model$rhs <- rep(1, each = m)
      
      # Set up bilinear constraints
      if (scheme == 'optimal equal') {
        model$quadcon <- list()
        for (j in 1:p) {
          for (k in 1:m) {
            model$quadcon[[j + p * (k - 1)]] <- list()
            model$quadcon[[j + p * (k - 1)]]$Qc <- 
              Matrix::sparseMatrix(2 * pm + k, (k - 1) * p + j, x = 1, dims = c(nvar, nvar))
            model$quadcon[[j + p * (k - 1)]]$q <- 
              Matrix::sparseVector(- 1,  pm + (k - 1) * p + j, nvar)
            model$quadcon[[j + p * (k - 1)]]$sense <- '='
            model$quadcon[[j + p * (k - 1)]]$rhs <- 0
          }
        }
      }
      
      # Use warm start if available
      if (length(weights) >= 1) model$start <- solved$x
      
      # Save weight vectors in p × m matrix
      solved <- gurobi::gurobi(model, params)
      weights[[length(weights) + 1]] <- matrix(solved$x[1:pm], p, m)
      
    }
    
  }
  
  weights

}

#--------------------------------------------------------------------------------------------------#
# covariance function for covariance estimation
#--------------------------------------------------------------------------------------------------#

# x: a matrix of mean zero variables
# correct: a logical indicating whether to perform a positive definite correction

covariance <- \(x, correct = F) {
  p <- ncol(x)
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:i) {
      xij <- na.omit(cbind(x[, i], x[, j]))
      if (length(xij) == 0) stop('Too many missing values')
      sigma[i, j] <- t(xij[, 1]) %*% xij[, 2] / nrow(xij)
    }
  }
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
  if (correct) as.matrix(Matrix::nearPD(sigma)$mat) else sigma
}
