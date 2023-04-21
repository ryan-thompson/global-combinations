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
    if (m != 1) {
      w <- sapply(1:m, \(k) fit[[which.min(cv[k, ])]][, k])
      attributes(w)$gamma <- sapply(1:m, \(k) attributes(fit[[which.min(cv[k, ])]])$gamma)
      attributes(w)$lambda <- sapply(1:m, \(k) attributes(fit[[which.min(cv[k, ])]])$lambda)
      w
    } else {
      fit[[which.min(cv)]]
    }

  } else {

    fit[[1]]

  }

}

#--------------------------------------------------------------------------------------------------#
# combine function for forecast combination
#-------------------------------------------------------------------------------------------------#

# Fits the forecast combination weights for a grid of globalisation and shrinkage parameters. The
# globalisation parameter is gamma and the shrinkage parameter is lambda

#--- Input ---#
# errors: a list of forecast error matrices where each column corresponds to a forecaster; each
# matrix should have the same number of columns
# scheme: the type of forecast combination scheme to apply; 'scheme = equal' for equal weights,
# 'scheme = optimal' for optimal weights, and so on
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
# q: type of globalisation penalty; `q = 1` for l1 penalty or `q =2 ` for l2^2 penalty
# correct: a logical indicating whether to perform a positive definite correction on any nonpositive
# definite covariance matrices

#--- Output ---#
# A list of matrices of weights; the kth column of each matrix is a vector of forecast combination
# weights for the kth task

combine <- \(errors, scheme = c('equal', 'optimal', 'optimal convex', 'optimal equal'),
             group = rep(1, length(sigma)), sigma = NULL, scale = T,
             ngamma = 10, gamma.min = 1e-3, gamma.max = 1e3, gamma = NULL,
             nlambda = 10, lambda.min = 1e-3, lambda.max = 1e3, lambda = NULL,
             q = 2, correct = F, ...) {

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
  scale.factor <- replicate(nlambda, rep(1, m), simplify = FALSE)
  if (scale & scheme != 'equal') {
    fit <- optimise(sigma, scheme, group, 0, lambda, scale.factor, q, ...)
    scale.factor <- lapply(1:nlambda, \(i) sapply(1:m, \(k) fit[[i]][, k] %*% 
                            sigma[[k]]%*% fit[[i]][, k] + lambda[i] * sum(abs(fit[[i]][, k]) ^ q)))
  }

  # Perform fit
  optimise(sigma, scheme, group, gamma, lambda, scale.factor, q, ...)

}

optimise <- \(sigma, scheme, group, gamma, lambda, scale.factor, q,
              params = list(OutputFlag = 1, NonConvex = 2, Threads = 1, TimeLimit = 60)) {

  # Save dimensions
  p <- ncol(sigma[[1]])
  m <- length(sigma)
  g <- length(unique(group))

  # Set up space for results
  weights <- list()
  pm <- p * m
  pg <- p * g
  pmg <- p * m * g
  nvar <- pm + pm + m + pg + pm + pmg + pmg

  for (i1 in 1:length(gamma)) {
    for (i2 in 1:length(lambda)) {

      model <- list()

      # Set up variable schemes
      model$vtype <- c(
        rep('C', pm),  # Weight vectors (continuous)
        rep('B', pm),  # Selection vectors (binary)
        rep('I', m),   # Numbers of selected forecasts (integer)
        rep('C', pg),  # Auxiliary weight vectors (continuous)
        rep('C', pm),  # Absolute value of weight vectors (continuous)
        rep('C', pmg), # Difference of weight vectors (continuous)
        rep('C', pmg)  # Absolute value of difference vectors (continuous)
      )
      model$lb <- c(
        rep(ifelse(scheme == 'equal', 1 / p, ifelse(scheme == 'optimal', - Inf, 0)), pm),
        rep(0, pm),
        rep(1, m),
        rep(- Inf, pg),
        rep(0, pm),
        rep(- Inf, pmg),
        rep(0, pmg)
      )
      model$ub <- c(
        rep(ifelse(scheme == 'equal', 1 / p, ifelse(scheme == 'optimal', Inf, 1)), pm),
        rep(1, pm),
        rep(p, m),
        rep(Inf, pg),
        rep(Inf, pm),
        rep(Inf, pmg),
        rep(Inf, pmg)
      )

      # Set up objective matrix
      sigma.scaled <- lapply(1:m, \(k) (sigma[[k]] + diag(ifelse(q == 2, lambda[i2], 0), p)) /
                               scale.factor[[i2]][k])
      global.penalty.mat <- Matrix::Matrix(0, pmg, pmg)
      if (q == 2) for (l in 1:g)
        global.penalty.mat[replicate(2, ((l - 1) * pm + 1:pm)[rep(group == l, each = p)])] <-
        gamma[i1]
      model$Q <- Matrix::bdiag(
        Matrix::bdiag(sigma.scaled),
        Matrix::Matrix(0, pm, pm),
        Matrix::Matrix(0, m, m),
        Matrix::Matrix(0, pg, pg),
        Matrix::Matrix(0, pm, pm),
        global.penalty.mat,
        Matrix::Matrix(0, pmg, pmg)
      )

      # Set up objective vector
      global.penalty.vec <- rep(0, pmg)
      if (q == 1) for (l in 1:g)
        global.penalty.vec[((l - 1) * pm + 1:pm)[rep(group == l, each = p)]] <- gamma[i1]
      shrink.penalty.vec <- rep(0, pm)
      if (q == 1) for (k in 1:m)
        shrink.penalty.vec[p * (k - 1) + 1:p] <- lambda[i2] * scale.factor[[i2]][k]
      model$obj <- c(rep(0, pm + pm + m + pg), shrink.penalty.vec, rep(0, pmg), global.penalty.vec)

      # Set up linear constraints
      # Sum to one constraint
      model$A <- Matrix::sparseMatrix(rep(1:m, each = p), 1:pm, x = 1, dims = c(m, nvar))
      # Difference of weight vectors constraint
      for (l in 1:g)
        model$A <- rbind(
          model$A,
          cbind(
            Matrix::Diagonal(pm),
            Matrix::Matrix(0, pm, pm),
            Matrix::Matrix(0, pm, m),
            Matrix::sparseMatrix(1:pm, (l - 1) * p + rep(1:p, m), x = - 1, dims = c(pm, pg)),
            Matrix::Matrix(0, pm, pm),
            Matrix::sparseMatrix(1:pm, (l - 1) * pm + 1:pm, x = - 1, dims = c(pm, pmg)),
            Matrix::Matrix(0, pm, pmg)
          )
        )
      model$sense <- c(rep('=', m), rep('=', pmg))
      model$rhs <- c(rep(1, m), rep(0, pmg))

      # Set up absolute value constraints
      model$genconabs <- list()
      for (j in 1:pm) {
        # Absolute values for weight vectors
        model$genconabs[[j]] <- list()
        model$genconabs[[j]]$resvar <- pm + pm + m + pg + j
        model$genconabs[[j]]$argvar <- j
        # Absolute values for difference of weight vectors
        model$genconabs[[pm + j]] <- list()
        model$genconabs[[pm + j]]$resvar <- pm + pm + m + pg + pm + pmg + j
        model$genconabs[[pm + j]]$argvar <- pm + pm + m + pg + pm + j
      }

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
      weights.matrix <- matrix(solved$x[1:pm], p, m)
      attributes(weights.matrix)$gamma <- gamma[i1]
      attributes(weights.matrix)$lambda <- lambda[i2]
      weights[[length(weights) + 1]] <- weights.matrix

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