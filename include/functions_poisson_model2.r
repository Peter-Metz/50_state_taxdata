
get_delta <- function(wh, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(wh / colSums(beta_x))
}

get_weights <- function(beta, delta, xmat){
  # get whs: state weights for households, given beta matrix, delta vector, and x matrix
  beta_x <- beta %*% t(xmat)
  # add delta vector to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(mat) mat + delta) 
  exp(beta_xd)
}

scale_problem <- function(problem, scale_goal){
  # problem is a list with at least the following:
  #  targets
  #  x
  # return:
  #   list with the scaled problem, including all of the original elements, plus
  #   scaled versions of x and targets
  #   plus new items scale_goal and scale_factor

  max_targets <- apply(problem$targets, 2, max) # find max target in each row of the target matrix
  scale_factor <- scale_goal / max_targets
  
  scaled_targets <- sweep(problem$targets, 2, scale_factor, "*")
  scaled_x <- sweep(problem$x, 2, scale_factor, "*")
  
  scaled_problem <- problem
  scaled_problem$targets <- scaled_targets
  scaled_problem$x <- scaled_x
  scaled_problem$scale_factor <- scale_factor
  
  scaled_problem
}


jac <- function(ewhs, xmatrix){
  # jacobian of distance vector relative to beta vector, IGNORING delta
  x2 <- xmatrix * xmatrix
  ddiag <- - t(ewhs) %*% x2 # note the minus sign in front
  diag(as.vector(ddiag)) 
}


distances <- function(betavec, wh, xmat, targets){
  # return a distance vector: differences between targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}


solve_poisson <- function(problem, maxiter=100, scale=FALSE, scale_goal=1000, step_scale=1){
  init_step_scale <- step_scale
  
  problem_unscaled <- problem
  if(scale==TRUE) problem <- scale_problem(problem, scale_goal)
  
  # unbundle the problem list and create additional variables needed
  targets <- problem$targets
  wh <- problem$wh
  xmat <- problem$x
  
  xpx <- t(xmat) %*% xmat
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible
  
  beta0 <- matrix(0, nrow=nrow(targets), ncol=ncol(targets)) # tpc uses 0 as beta starting point
  delta0 <- get_delta(wh, beta0, xmat) # tpc uses initial delta based on initial beta 
  
  ebeta <- beta0 # tpc uses 0 as beta starting point
  edelta <- delta0 # tpc uses initial delta based on initial beta 

  sse_vec <- rep(NA_real_, maxiter)
  # iter <- 0
  # step_scale <- 1
  
  for(iter in 1:maxiter){
    # iter <- iter + 1
    ewhs <- get_weights(ebeta, edelta, xmat)
    ews <- colSums(ewhs)
    ewh <- rowSums(ewhs)
    
    etargets <- t(ewhs) %*% xmat
    d <- targets - etargets
    sse <- sum(d^2)
    if(is.na(sse)) break
    sse_vec[iter] <- sse
    best_sse <- min(sse_vec, na.rm=TRUE)
    if(sse==best_sse) best_ebeta <- ebeta
    prior_sse <- sse
    if(iter <=20 | iter %% 20 ==0) print(sprintf("iteration: %i  sse: %.5e ", iter, sse))
    if(sse < 1e-6) {
      # exit if good
      print(sprintf("exit at iteration: %i  sse: %.5e ", iter, sse))
      break
    }
    
    # if sse is > prior sse, adjust step scale downward
    if(sse > best_sse){
      step_scale <- step_scale * .5
      ebeta <- prior_ebeta
    }
    
    prior_ebeta <- ebeta
    
    # ad hoc step
    step <- -(1 / ews) * d %*% invxpx * step_scale
    
    # jval <- jacobian(distances, x=as.vector(ebeta), wh=wh, xmat=xmat, targets=targets, method="simple") # f is differences
    # jval <- jac(ewhs, xmat)
    # step <- solve(jval) %*% as.vector(d) # , tol = 1e-30
    # step <- matrix(step, nrow=nrow(d), byrow=FALSE)
    
    ebeta <- ebeta - step
    edelta <- get_delta(ewh, ebeta, xmat)
  }
  
  best_edelta <- get_delta(ewh, best_ebeta, xmat)
  ewhs <- get_weights(best_ebeta, best_edelta, xmat)
  ewh <- rowSums(ewhs)
  if(scale==TRUE) etargets <- sweep(etargets, 2, problem$scale_factor, "/")
  final_step_scale <- step_scale
  
  keepnames <- c("maxiter", "iter", "sse", "sse_vec", "d", "best_ebeta", "best_edelta", "ewh", "ewhs", "etargets",
                 "problem_unscaled", "scale", "scale_goal", "init_step_scale", "final_step_scale")
  result <- list()
  for(var in keepnames) result[[var]] <- get(var)
  result
}
