

# This program is based largely on the methodology in:
#   Khitatrakun, Surachai, Gordon B T Mermin, and Norton Francis. “Incorporating State Analysis into the 
#   Tax Policy Center’s Microsimulation Model: Documentation and Methodology.” Working Paper, March 2016.
#   https://www.taxpolicycenter.org/sites/default/files/alfresco/publication-pdfs/2000697-Incorporating-State-Analysis-into-the-TPCs-Microsimulation-Model.pdf.


calc_delta <- function(hh_weights, beta, xmatrix){
  beta_x <- exp(beta %*% t(xmatrix))
  log(hh_weights / colSums(beta_x))
}


calc_weights <- function(beta, delta, xmatrix){
  # calculate all weights
  beta_x <- beta %*% t(xmatrix)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1, function(mat) mat + delta) 
  exp(beta_xd)
}


poisson_weights2 <- function(hh_weights, targets, xmatrix, x_scale=NULL, step_scale=NULL, maxiter=200){
  # function to solve for state weights for each household, so that:
  #   (a) the weights for each state follow a poisson distribution
  #   (b) the sum of weighted values for each state, for each household characteristic, equals or comes close to
  #       the targets, and
  #   (c) the sum of state weights for each household equals the total national weight for that household

  #. definitions ----
  # s prefix means a variable is the scaled version of an input variable (multiplied by a scale factor)
  # i prefix means a variable is an initial value set before iteration
  # e prefix means a variable is an estimate developed during iteration
  # f prefix means a variable is the final value, after iteration

  # hh_weights:
  #   vector of national weights (total weight for each household) -- one per household

  # targets:
  #   matrix of targets
  #     1 row per state
  #     1 column per characteristic

  # xmatrix:
  #   matrix of household characteristics (e.g., wages, interest, dividends)
  #     1 row per household
  #     1 column per characteristic


  # scale parameters - to scale, we MULTIPLY by these values:
  #   x_scale:      for scaling xmatrix and targets
  #   step_scale:   for scaling the computed step size


  #. initialization ----
  t1 <- proc.time()

  #.  . define indexes ----
  n_states <- nrow(targets) # i.e., s
  n_hh <- nrow(xmatrix) # number of households
  n_characteristics <- ncol(xmatrix)

  #.  . compute scaling factors ----
  # to scale, we MULTIPLY by these values
  if(is.null(x_scale)) x_scale <- 1000 / sum(xmatrix) # seems to be a good default value -- scaled x will sum to 1000
  if(is.null(step_scale)) step_scale <- n_hh

  #.  . scale the problem ----
  s_xmatrix <- xmatrix * x_scale
  s_targets <- targets * x_scale

  #.  . compute inverse of xprime-x matrix ----
  #   When we calculate search direction in the loop, we will need the inverse of a derivatives matrix that
  #   involves the xprime-x matrix. Because the inverse of matrix that has been multiplied by a non-zero scalar equals
  #   the inverse of the scalar multiplied by inverse of the matrix, we can calculate the inverse of xprime-x once
  #   before entering the loop. Within the loop we will multiply by the inverse of the relevant scalar (computed in the loop).
  xpx <- t(s_xmatrix) %*% s_xmatrix
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible

  #.  . define initial values ----
  i_beta <- matrix(0, nrow=n_states, ncol=n_characteristics) # tpc uses 0 as beta starting point
  i_delta <- calc_delta(hh_weights, i_beta, s_xmatrix) # tpc uses initial delta based on initial beta

  i_hh_state_weights <- calc_weights(i_beta, i_delta, s_xmatrix)

  i_dist <- s_targets - t(i_hh_state_weights) %*% s_xmatrix # distance from targets using initial values

  i_sse <- sum(i_dist^2) # initial sum of squared distances

  #.  . set estimated values used in loop to intial values, before entering loop for first time ----
  e_beta <- i_beta
  e_delta <- i_delta

  # begin iterations ----
  t2 <- proc.time()
  for(i in 1:maxiter){
    # iterate, searching for best beta and delta coefficients for the poisson model

    e_hh_state_weights <- calc_weights(e_beta, e_delta, s_xmatrix) # estimated weights for each household for each state
    e_hh_weights <- rowSums(e_hh_state_weights) # estimated total weights for each households
    e_state_weights <- colSums(e_hh_state_weights) # estimated total weights for each state

    e_targets <- t(e_hh_state_weights) %*% s_xmatrix
    dist <- s_targets - e_targets # distance from scaled targets
    dist_u <- dist / x_scale
    sse <- sum(dist^2) # sum of squared errors
    # print(sse)
    # sse <- 1; sseu <- 1
    sseu <- sse
    sseu <- sse * (x_scale^-2) # unscaled sse
    maxdist <- max(abs(rowSums(dist_u)))

    # if(sseu < 1e-6 | (sse < 1e-10 & sseu < 1e-2)) {
    #   print(sprintf("DONE at iteration %i: scaled sse: %.5e, unscaled sse: %.5e", i, sse, sseu))
    #   break
    # }

    if(maxdist < 500) {
      print(sprintf("DONE at iteration %i: max error: %.5e", i, maxdist))
      break
    }

    if(i <= 10 | (i %% 10)==0) {
      print(sprintf("iteration %i: scaled sse: %.5e, unscaled sse: %.5e max error: %.5e", i, sse, sseu, maxdist))
    }

    step <- (1 / e_state_weights) * dist %*% invxpx * step_scale # get the step (search direction)

    e_beta <- e_beta + step
    e_delta <- calc_delta(e_hh_weights, e_beta, s_xmatrix)
  }
  t3 <- proc.time()
  # after we exit the loop, e_beta and e_delta have the desired coefficients

  # post-loop: calculate weights and return results ----
  f_beta <- e_beta
  f_delta <- e_delta
  f_hh_state_weights <- calc_weights(f_beta, f_delta, s_xmatrix)
  f_hh_weights <- rowSums(f_hh_state_weights)
  f_state_weights <- colSums(f_hh_state_weights)
  f_targets <- (t(f_hh_state_weights) %*% s_xmatrix) / x_scale

  result <- list()
  result$iterations <- i
  result$unscaled_sse <- sseu

  vars <- c("x_scale", "step_scale", "maxdist",
            "f_beta", "f_delta", "f_hh_state_weights",
            "f_hh_weights", "f_targets",
            "hh_weights", "targets", "xmatrix")
  for(var in vars) result[[var]] <- get(var)
  t4 <- proc.time()
  result$total_seconds <- as.numeric((t4 - t1)[3])
  result$iter_seconds <- as.numeric((t3 - t2)[3])

  return(result)
}




poisson_weights <- function(hh_weights, targets, xmatrix, maxiter=200, x_scalefactor=1000, 
                            step_start=nrow(xmatrix), step_stop=1000, step_iter=100){
  # function to solve for state weights for each household, so that:
  #   (a) the weights for each state follow a poisson distribution
  #   (b) the sum of weighted values for each state, for each household characteristic, equals or comes close to
  #       the targets, and
  #   (c) the sum of state weights for each household equals the total national weight for that household
  
  #. definitions ----
  # s prefix means a variable is the scaled version of an input variable (multiplied by a scale factor)
  # i prefix means a variable is an initial value set before iteration
  # e prefix means a variable is an estimate developed during iteration
  # f prefix means a variable is the final value, after iteration
  
  # hh_weights:
  #   vector of national weights (total weight for each household) -- one per household
  
  # targets:
  #   matrix of targets
  #     1 row per state
  #     1 column per characteristic
  
  # xmatrix:
  #   matrix of household characteristics (e.g., wages, interest, dividends)
  #     1 row per household
  #     1 column per characteristic
  
  
  # scale parameters - to scale, we MULTIPLY by these values:
  #   x_scale:      for scaling xmatrix and targets
  #   step_scale:   for scaling the computed step size
  
  
  #. initialization ----
  t1 <- proc.time()
  
  #.  . define indexes ----
  n_states <- nrow(targets) # i.e., s
  n_hh <- nrow(xmatrix) # number of households
  n_characteristics <- ncol(xmatrix)
  
  #.  . compute scaling factors ----
  # to scale, we MULTIPLY by these values
  # if(is.null(x_scale)) x_scale <- 1000 / sum(xmatrix) # seems to be a good default value -- scaled x will sum to 1000
  # create xscale vector so that each target is calibrarted to 1000
  # x_scale <- ifelse(targets==0, 0, 10000 / targets)
  x_scale <- ifelse(colSums(targets)==0, 1, x_scalefactor / colSums(targets))
  # print(x_scale)
  # x_scale <- c(10, 1, 1, 1)
  
  # if(is.null(step_scale)) step_scale <- n_hh * .8
  steps <- rep(step_stop, maxiter)
  isteps <- min(step_iter, maxiter)
  steps[1:isteps] <- seq(step_start, step_stop, length.out = isteps)
  # print(steps)
  
  #.  . scale the problem ----
  s_xmatrix <-sweep(xmatrix, 2, x_scale, FUN = "*") # efficient way to multiply each row of matrix by a vector
  s_targets <- sweep(targets, 2, x_scale, FUN = "*")
  # s_targets <- ifelse(s_targets==0, 1, s_targets)
  # print(x_scale)
  # print(s_targets)
  # targets * x_scale
  
  #.  . compute inverse of xprime-x matrix ----
  #   When we calculate search direction in the loop, we will need the inverse of a derivatives matrix that 
  #   involves the xprime-x matrix. Because the inverse of matrix that has been multiplied by a non-zero scalar equals 
  #   the inverse of the scalar multiplied by inverse of the matrix, we can calculate the inverse of xprime-x once
  #   before entering the loop. Within the loop we will multiply by the inverse of the relevant scalar (computed in the loop).
  xpx <- t(s_xmatrix) %*% s_xmatrix
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible
  
  #.  . define initial values ----
  i_beta <- matrix(0, nrow=n_states, ncol=n_characteristics) # tpc uses 0 as beta starting point
  i_delta <- calc_delta(hh_weights, i_beta, s_xmatrix) # tpc uses initial delta based on initial beta 
  
  i_hh_state_weights <- calc_weights(i_beta, i_delta, s_xmatrix)
  
  i_dist <- s_targets - t(i_hh_state_weights) %*% s_xmatrix # distance from targets using initial values
  
  i_sse <- sum(i_dist^2) # initial sum of squared distances
  
  #.  . set estimated values used in loop to intial values, before entering loop for first time ----
  e_beta <- i_beta
  e_delta <- i_delta
  
  # begin iterations ----
  t2 <- proc.time()
  for(i in 1:maxiter){
    # iterate, searching for best beta and delta coefficients for the poisson model 
    
    # if(i==4){
    #   print("last step:")
    #   print(step)
    #   print(e_beta)
    #   print(e_delta)
    #   print(s_targets)
    #   print(e_targets)
    #   print(dist) 
    #   print(pdist)
    # }
    
    e_hh_state_weights <- calc_weights(e_beta, e_delta, s_xmatrix) # estimated weights for each household for each state
    e_hh_weights <- rowSums(e_hh_state_weights) # estimated total weights for each households
    e_state_weights <- colSums(e_hh_state_weights) # estimated total weights for each state
    
    e_targets <- t(e_hh_state_weights) %*% s_xmatrix
    dist <- e_targets - s_targets # distance from scaled targets
    # print(dist)
    # dist_u <- dist / x_scale
    dist_u <- dist
    pdist <- ifelse(s_targets==0, 1, dist / s_targets)

    sse <- sum(dist^2) # sum of squared errors
    # print(sse)
    # sse <- 1; sseu <- 1
    sseu <- sse
    # sseu <- sse * (x_scale^-2) # unscaled sse
    maxdist <- max(abs(pdist))
    # if(i==5) print(pdist)
    

    # step_scale <- max(n_hh * .75 * (1 - (i / step_scalefactor)/(100)), 100)
    # step_scale <- ifelse(i < step_iter, step_start - (step_start - step_stop) * i / step_iter, step_stop)
    step_scale <- steps[i]
    
    # if(sseu < 1e-6 | (sse < 1e-10 & sseu < 1e-2)) {
    #   print(sprintf("DONE at iteration %i: scaled sse: %.5e, unscaled sse: %.5e", i, sse, sseu))
    #   break
    # }
    
    # if(maxdist < 1e-2) {
    #   print(sprintf("DONE at iteration %i: max error: %.5e", i, maxdist))
    #   break
    # }
    
    if(i <= 10 | (i %% 10)==0) {
      print(sprintf("iteration %i: scaled sse: %.5e, unscaled sse: %.5e max error: %.5e step scale: %.5e", i, sse, sseu, maxdist, step_scale))
    }
    
    # step <- (1 / e_state_weights) * dist %*% invxpx * step_scale # get the step (search direction)
    # step <- (1 / e_state_weights) * dist %*% invxpx * step_scale  # get the step (search direction)
    step <- (1 / e_state_weights) * dist %*% invxpx * 7000 # get the step (search direction)
    # step <- ifelse(is.infinite(step), 0, step)
    # print("e_state_weights"); print(e_state_weights)
    # print("dist"); print(dist)
    # print("step"); print(step)
    
    e_beta <- e_beta + step
    e_delta <- calc_delta(e_hh_weights, e_beta, s_xmatrix)
  }
  t3 <- proc.time()
  # after we exit the loop, e_beta and e_delta have the desired coefficients
  
  # post-loop: calculate weights and return results ----
  f_beta <- e_beta
  f_delta <- e_delta
  f_hh_state_weights <- calc_weights(f_beta, f_delta, s_xmatrix)
  f_hh_weights <- rowSums(f_hh_state_weights)
  f_state_weights <- colSums(f_hh_state_weights)
  f_s_targets <- t(f_hh_state_weights) %*% s_xmatrix
  f_targets <- sweep((t(f_hh_state_weights) %*% s_xmatrix), 2, x_scale, FUN = "/")
  
  result <- list()
  result$iterations <- i
  result$unscaled_sse <- sseu
  
  vars <- c("x_scale", "step_scale", "maxdist",
            "f_beta", "f_delta", "f_hh_state_weights",
            "f_hh_weights", "targets", "f_targets",
            "hh_weights", "xmatrix", "f_s_targets", "s_targets")
  for(var in vars) result[[var]] <- get(var)
  t4 <- proc.time()
  result$total_seconds <- as.numeric((t4 - t1)[3])
  result$iter_seconds <- as.numeric((t3 - t2)[3])
  
  return(result)
}




poisson_weights3 <- function(hh_weights, targets, xmatrix, maxiter=200,
                             x_scale=NULL, 
                             step_start, step_stop, step_iter){
  # function to solve for state weights for each household, so that:
  #   (a) the weights for each state follow a poisson distribution
  #   (b) the sum of weighted values for each state, for each household characteristic, equals or comes close to
  #       the targets, and
  #   (c) the sum of state weights for each household equals the total national weight for that household
  
  #. definitions ----
  # s prefix means a variable is the scaled version of an input variable (multiplied by a scale factor)
  # i prefix means a variable is an initial value set before iteration
  # e prefix means a variable is an estimate developed during iteration
  # f prefix means a variable is the final value, after iteration
  
  # hh_weights:
  #   vector of national weights (total weight for each household) -- one per household
  
  # targets:
  #   matrix of targets
  #     1 row per state
  #     1 column per characteristic
  
  # xmatrix:
  #   matrix of household characteristics (e.g., wages, interest, dividends)
  #     1 row per household
  #     1 column per characteristic
  
  
  # scale parameters - to scale, we MULTIPLY by these values:
  #   x_scale:      for scaling xmatrix and targets
  #   step_scale:   for scaling the computed step size
  
  
  #. initialization ----
  t1 <- proc.time()
  
  #.  . define indexes ----
  n_states <- nrow(targets) # i.e., s
  n_hh <- nrow(xmatrix) # number of households
  n_characteristics <- ncol(xmatrix)
  
  #.  . compute scaling factors ----
  # to scale, we MULTIPLY by these values
  if(is.null(x_scale)) x_scale <- 1000 / sum(xmatrix) # seems to be a good default value -- scaled x will sum to 1000
  # if(is.null(step_scale)) step_scale <- n_hh
  
  #.  . scale the problem ----
  s_xmatrix <- xmatrix * x_scale
  s_targets <- targets * x_scale
  
  #.  . compute inverse of xprime-x matrix ----
  #   When we calculate search direction in the loop, we will need the inverse of a derivatives matrix that
  #   involves the xprime-x matrix. Because the inverse of matrix that has been multiplied by a non-zero scalar equals
  #   the inverse of the scalar multiplied by inverse of the matrix, we can calculate the inverse of xprime-x once
  #   before entering the loop. Within the loop we will multiply by the inverse of the relevant scalar (computed in the loop).
  xpx <- t(s_xmatrix) %*% s_xmatrix
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible
  
  #.  . define initial values ----
  i_beta <- matrix(0, nrow=n_states, ncol=n_characteristics) # tpc uses 0 as beta starting point
  i_delta <- calc_delta(hh_weights, i_beta, s_xmatrix) # tpc uses initial delta based on initial beta
  
  i_hh_state_weights <- calc_weights(i_beta, i_delta, s_xmatrix)
  
  i_dist <- s_targets - t(i_hh_state_weights) %*% s_xmatrix # distance from targets using initial values
  
  i_sse <- sum(i_dist^2) # initial sum of squared distances
  
  steps <- rep(step_stop, maxiter)
  iphasein <- min(step_iter, maxiter) # 
  steps[1:iphasein] <- seq(step_start, step_stop, length.out = iphasein)
  
  #.  . set estimated values used in loop to intial values, before entering loop for first time ----
  e_beta <- i_beta
  e_delta <- i_delta
  
  # begin iterations ----
  t2 <- proc.time()
  for(i in 1:maxiter){
    # iterate, searching for best beta and delta coefficients for the poisson model
    
    e_hh_state_weights <- calc_weights(e_beta, e_delta, s_xmatrix) # estimated weights for each household for each state
    e_hh_weights <- rowSums(e_hh_state_weights) # estimated total weights for each households
    e_state_weights <- colSums(e_hh_state_weights) # estimated total weights for each state
    
    e_targets <- t(e_hh_state_weights) %*% s_xmatrix
    dist <- s_targets - e_targets # distance from scaled targets
    dist_u <- dist / x_scale
    pdist <- ifelse(s_targets==0, 1, dist / s_targets)
    maxdist <- max(abs(pdist)) # maximum absolute distance as a proportion of target
    
    sse <- sum(dist^2) # sum of squared errors
    # print(sse)
    # sse <- 1; sseu <- 1
    sseu <- sse
    sseu <- sse * (x_scale^-2) # unscaled sse
    
    # if(sseu < 1e-6 | (sse < 1e-10 & sseu < 1e-2)) {
    #   print(sprintf("DONE at iteration %i: scaled sse: %.5e, unscaled sse: %.5e", i, sse, sseu))
    #   break
    # }
    
    if(maxdist < .001) {
      print(sprintf("DONE at iteration %i: max error: %.5e", i, maxdist))
      break
    }
    
    if(i <= 10 | (i %% 10)==0) {
      print(sprintf("iteration %i: scaled sse: %.5e, unscaled sse: %.5e max error: %.5e", i, sse, sseu, maxdist))
    }
    
    step_scale <- steps[i]
    step <- (1 / e_state_weights) * dist %*% invxpx * step_scale # get the step (search direction)
    
    e_beta <- e_beta + step
    e_delta <- calc_delta(e_hh_weights, e_beta, s_xmatrix)
  }
  t3 <- proc.time()
  # after we exit the loop, e_beta and e_delta have the desired coefficients
  
  # post-loop: calculate weights and return results ----
  f_beta <- e_beta
  f_delta <- e_delta
  f_hh_state_weights <- calc_weights(f_beta, f_delta, s_xmatrix)
  f_hh_weights <- rowSums(f_hh_state_weights)
  f_state_weights <- colSums(f_hh_state_weights)
  f_s_targets <- t(f_hh_state_weights) %*% s_xmatrix
  f_targets <- (t(f_hh_state_weights) %*% s_xmatrix) / x_scale
  
  result <- list()
  result$iterations <- i
  result$unscaled_sse <- sseu
  
  vars <- c("x_scale", "step_scale", "maxdist",
            "f_beta", "f_delta",
            "targets", "f_targets",
            # "s_targets", "f_s_targets",
            "hh_weights", "f_hh_weights",
            "f_hh_state_weights",
            "xmatrix", "s_xmatrix")
  for(var in vars) result[[var]] <- get(var)
  t4 <- proc.time()
  result$total_seconds <- as.numeric((t4 - t1)[3])
  result$iter_seconds <- as.numeric((t3 - t2)[3])
  
  return(result)
}



poisson_weights4 <- function(hh_weights, targets, xmatrix, maxiter=200,
                             x_scale=NULL, 
                             step_start, step_stop, step_iter){
  # function to solve for state weights for each household, so that:
  #   (a) the weights for each state follow a poisson distribution
  #   (b) the sum of weighted values for each state, for each household characteristic, equals or comes close to
  #       the targets, and
  #   (c) the sum of state weights for each household equals the total national weight for that household
  
  #. definitions ----
  # s prefix means a variable is the scaled version of an input variable (multiplied by a scale factor)
  # i prefix means a variable is an initial value set before iteration
  # e prefix means a variable is an estimate developed during iteration
  # f prefix means a variable is the final value, after iteration
  
  # hh_weights:
  #   vector of national weights (total weight for each household) -- one per household
  
  # targets:
  #   matrix of targets
  #     1 row per state
  #     1 column per characteristic
  
  # xmatrix:
  #   matrix of household characteristics (e.g., wages, interest, dividends)
  #     1 row per household
  #     1 column per characteristic
  
  
  # scale parameters - to scale, we MULTIPLY by these values:
  #   x_scale:      for scaling xmatrix and targets
  #   step_scale:   for scaling the computed step size
  
  
  #. initialization ----
  t1 <- proc.time()
  
  #.  . define indexes ----
  n_states <- nrow(targets) # i.e., s
  n_hh <- nrow(xmatrix) # number of households
  n_characteristics <- ncol(xmatrix)
  
  #.  . compute scaling factors ----
  # to scale, we MULTIPLY by these values
  # if(is.null(x_scale)) x_scale <- 1000 / sum(xmatrix) # seems to be a good default value -- scaled x will sum to 1000
  x_scale <- ifelse(colSums(targets)==0, 1, x_scalefactor / colSums(targets))
  # if(is.null(step_scale)) step_scale <- n_hh
  
  #.  . scale the problem ----
  # s_xmatrix <- xmatrix * x_scale
  # s_targets <- targets * x_scale
  s_xmatrix <-sweep(xmatrix, 2, x_scale, FUN = "*") # efficient way to multiply each row of matrix by a vector
  s_targets <- sweep(targets, 2, x_scale, FUN = "*")
  
  #.  . compute inverse of xprime-x matrix ----
  #   When we calculate search direction in the loop, we will need the inverse of a derivatives matrix that
  #   involves the xprime-x matrix. Because the inverse of matrix that has been multiplied by a non-zero scalar equals
  #   the inverse of the scalar multiplied by inverse of the matrix, we can calculate the inverse of xprime-x once
  #   before entering the loop. Within the loop we will multiply by the inverse of the relevant scalar (computed in the loop).
  xpx <- t(s_xmatrix) %*% s_xmatrix
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible
  
  #.  . define initial values ----
  i_beta <- matrix(0, nrow=n_states, ncol=n_characteristics) # tpc uses 0 as beta starting point
  i_delta <- calc_delta(hh_weights, i_beta, s_xmatrix) # tpc uses initial delta based on initial beta
  
  i_hh_state_weights <- calc_weights(i_beta, i_delta, s_xmatrix)
  
  i_dist <- s_targets - t(i_hh_state_weights) %*% s_xmatrix # distance from targets using initial values
  
  i_sse <- sum(i_dist^2) # initial sum of squared distances
  
  steps <- rep(step_stop, maxiter)
  iphasein <- min(step_iter, maxiter) # 
  steps[1:iphasein] <- seq(step_start, step_stop, length.out = iphasein)
  
  #.  . set estimated values used in loop to intial values, before entering loop for first time ----
  e_beta <- i_beta
  e_delta <- i_delta
  
  # begin iterations ----
  t2 <- proc.time()
  for(i in 1:maxiter){
    # iterate, searching for best beta and delta coefficients for the poisson model
    
    e_hh_state_weights <- calc_weights(e_beta, e_delta, s_xmatrix) # estimated weights for each household for each state
    e_hh_weights <- rowSums(e_hh_state_weights) # estimated total weights for each households
    e_state_weights <- colSums(e_hh_state_weights) # estimated total weights for each state
    
    e_targets <- t(e_hh_state_weights) %*% s_xmatrix
    dist <- s_targets - e_targets # distance from scaled targets
    # dist_u <- dist / x_scale
    pdist <- ifelse(s_targets==0, 1, dist / s_targets)
    maxdist <- max(abs(pdist)) # maximum absolute distance as a proportion of target
    
    sse <- sum(dist^2) # sum of squared errors
    # print(sse)
    # sse <- 1; sseu <- 1
    sseu <- sse
    # sseu <- sse * (x_scale^-2) # unscaled sse
    
    # if(sseu < 1e-6 | (sse < 1e-10 & sseu < 1e-2)) {
    #   print(sprintf("DONE at iteration %i: scaled sse: %.5e, unscaled sse: %.5e", i, sse, sseu))
    #   break
    # }
    
    if(maxdist < .001) {
      print(sprintf("DONE at iteration %i: max error: %.5e", i, maxdist))
      break
    }
    
    if(i <= 10 | (i %% 10)==0) {
      print(sprintf("iteration %i: scaled sse: %.5e, unscaled sse: %.5e max error: %.5e", i, sse, sseu, maxdist))
    }
    
    step_scale <- steps[i]
    step <- (1 / e_state_weights) * dist %*% invxpx * step_scale # get the step (search direction)
    
    e_beta <- e_beta + step
    e_delta <- calc_delta(e_hh_weights, e_beta, s_xmatrix)
  }
  t3 <- proc.time()
  # after we exit the loop, e_beta and e_delta have the desired coefficients
  
  # post-loop: calculate weights and return results ----
  f_beta <- e_beta
  f_delta <- e_delta
  f_hh_state_weights <- calc_weights(f_beta, f_delta, s_xmatrix)
  f_hh_weights <- rowSums(f_hh_state_weights)
  f_state_weights <- colSums(f_hh_state_weights)
  f_s_targets <- t(f_hh_state_weights) %*% s_xmatrix
  f_targets <- sweep(f_s_targets, 2, x_scale, FUN="/")
  result <- list()
  result$iterations <- i
  result$unscaled_sse <- sseu
  
  vars <- c("x_scale", "step_scale", "maxdist",
            "f_beta", "f_delta",
            "targets", "f_targets",
            # "s_targets", "f_s_targets",
            "hh_weights", "f_hh_weights",
            "f_hh_state_weights",
            "xmatrix", "s_xmatrix")
  for(var in vars) result[[var]] <- get(var)
  t4 <- proc.time()
  result$total_seconds <- as.numeric((t4 - t1)[3])
  result$iter_seconds <- as.numeric((t3 - t2)[3])
  
  return(result)
}


