
# delta, weights and vtom ----

get_delta <- function(wh, beta, xmat){
  # we cannot let beta %*% xmat get too large!! or exp will be Inf and problem will bomb
  # it will get large when a beta element times an xmat element is large, so either
  # beta or xmat can be the problem
  beta_x <- exp(beta %*% t(xmat))
  log(wh / colSums(beta_x)) # denominator is sum for each person
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
  #  xmat
  # return:
  #   list with the scaled problem, including all of the original elements, plus
  #   scaled versions of x and targets
  #   plus new items scale_goal and scale_factor

  max_targets <- apply(problem$targets, 2, max) # find max target in each row of the target matrix
  scale_factor <- scale_goal / max_targets
  
  scaled_targets <- sweep(problem$targets, 2, scale_factor, "*")
  scaled_x <- sweep(problem$xmat, 2, scale_factor, "*")
  
  scaled_problem <- problem
  scaled_problem$targets <- scaled_targets
  scaled_problem$xmat <- scaled_x
  scaled_problem$scale_factor <- scale_factor
  
  scaled_problem
}


scale_problem_mdn <- function(problem, scale_goal){
  # problem is a list with at least the following:
  #  targets
  #  xmat
  # return:
  #   list with the scaled problem, including all of the original elements, plus
  #   scaled versions of x and targets
  #   plus new items scale_goal and scale_factor
  
  mdn_targets <- apply(problem$targets, 2, median) # find median target in each row of the target matrix
  scale_factor <- scale_goal / mdn_targets
  
  scaled_targets <- sweep(problem$targets, 2, scale_factor, "*")
  scaled_x <- sweep(problem$xmat, 2, scale_factor, "*")
  
  scaled_problem <- problem
  scaled_problem$targets <- scaled_targets
  scaled_problem$xmat <- scaled_x
  scaled_problem$scale_factor <- scale_factor
  
  scaled_problem
}


vtom <- function(vec, nrows){
  # vector to matrix in the same ordering as a beta matrix
  matrix(vec, nrow=nrows, byrow=FALSE)
}


# differences and sse ----

etargs_vec <- function(betavec, wh, xmat, s){
  # return a vector of calculated targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- vtom(betavec, nrows=s)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  as.vector(etargets)
}


diff_vec <- function(betavec, wh, xmat, targets){
  # return a vector of differences between targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}


sse_fn <- function(betavec, wh, xmat, targets){
  # return a single value - sse (sum of squared errors)
  sse <- sum(diff_vec(betavec, wh, xmat, targets)^2)
  sse
}


# step functions ----
step_fd <- function(ebeta, step_inputs){
  # finite differences
  bvec <- as.vector(ebeta)
  gbeta <- numDeriv::grad(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  hbeta <- numDeriv::hessian(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  ihbeta <- solve(hbeta)
  stepvec <- t(ihbeta %*% gbeta)
  step <- vtom(stepvec, nrows=nrow(ebeta))
  step
}

step_fd <- function(ebeta, step_inputs){
  # finite differences -- version to print time
  bvec <- as.vector(ebeta)
  t1 <- proc.time()
  gbeta <- numDeriv::grad(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  t2 <- proc.time()
  print(sprintf("gradient time in seconds: %.1e", (t2-t1)[3]))
  hbeta <- numDeriv::hessian(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  t3 <- proc.time()
  print(sprintf("hessian time in seconds: %.1e", (t3-t2)[3]))
  ihbeta <- solve(hbeta)
  t4 <- proc.time()
  print(sprintf("inverse time in seconds: %.1e", (t4-t3)[3]))
  stepvec <- t(ihbeta %*% gbeta)
  step <- vtom(stepvec, nrows=nrow(ebeta))
  step
}


step_adhoc <- function(ebeta, step_inputs){
  -(1 / step_inputs$ews) * step_inputs$d %*% step_inputs$invxpx * step_inputs$step_scale
}

get_step <- function(step_method, ebeta, step_inputs){
  step <- case_when(step_method=="adhoc" ~ step_adhoc(ebeta, step_inputs),
                    step_method=="finite_diff" ~ step_fd(ebeta, step_inputs),
                    TRUE ~ step_adhoc(ebeta, step_inputs))
  step
}

jac_tpc <- function(ewhs, xmatrix){
  # jacobian of distance vector relative to beta vector, IGNORING delta
  x2 <- xmatrix * xmatrix
  ddiag <- - t(ewhs) %*% x2 # note the minus sign in front
  diag(as.vector(ddiag)) 
}


# solve poisson problem ----
solve_poisson <- function(problem, maxiter=100, scale=FALSE, scale_goal=1000, step_method="adhoc", step_scale=1, tol=1e-3, start=NULL){
  t1 <- proc.time()
  
  if(step_method=="adhoc") step_fn <- step_adhoc else{
    if(step_method=="finite_diff") step_fn <- step_fd
  }
  
  init_step_scale <- step_scale
  
  problem_unscaled <- problem
  if(scale==TRUE) problem <- scale_problem(problem, scale_goal)
  
  # unbundle the problem list and create additional variables needed
  targets <- problem$targets
  wh <- problem$wh
  xmat <- problem$xmat
  
  xpx <- t(xmat) %*% xmat
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible

  if(is.null(start)) beta0 <- matrix(0, nrow=nrow(targets), ncol=ncol(targets)) else # tpc uses 0 as beta starting point
    beta0 <- start
  delta0 <- get_delta(wh, beta0, xmat) # tpc uses initial delta based on initial beta 
  
  ebeta <- beta0 # tpc uses 0 as beta starting point
  edelta <- delta0 # tpc uses initial delta based on initial beta 

  sse_vec <- rep(NA_real_, maxiter)
  
  step_inputs <- list()
  step_inputs$targets <- targets
  step_inputs$step_scale <- step_scale
  step_inputs$xmat <- xmat
  step_inputs$invxpx <- invxpx
  step_inputs$wh <- wh
  
  for(iter in 1:maxiter){
    # iter <- iter + 1
    edelta <- get_delta(wh, ebeta, xmat)
    ewhs <- get_weights(ebeta, edelta, xmat)
    ews <- colSums(ewhs)
    ewh <- rowSums(ewhs)
    step_inputs$ews <- ews
    
    etargets <- t(ewhs) %*% xmat
    d <- targets - etargets
    step_inputs$d <- d
    
    rel_err <- ifelse(targets==0, NA, abs(d / targets))
    max_rel_err <- max(rel_err, na.rm=TRUE)
    sse <- sum(d^2)
    if(is.na(sse)) break # bad result, end it now, we have already saved the prior best result
    
    sse_vec[iter] <- sse
    # sse_vec <- c(seq(200, 100, -1), NA, NA)
    sse_rel_change <- sse_vec / lag(sse_vec) - 1
    # iter <- 5
    # test2 <- ifelse(iter >= 5, !any(sse_rel_change[iter - 0:2] < -.01), FALSE)
    # test2
    # any(sse_rel_change[c(4, 5, 6)] < -.01)
    
    best_sse <- min(sse_vec, na.rm=TRUE)
    if(sse==best_sse) best_ebeta <- ebeta
    prior_sse <- sse
    
    if(iter <=20 | iter %% 20 ==0) print(sprintf("iteration: %i, sse: %.5e, max_rel_err: %.5e", iter, sse, max_rel_err))
    
    #.. stopping criteria ---- iter <- 5
    test1 <- max_rel_err < tol # every distance from target is within our desired error tolerance
    # test2: none the of last 3 iterations had sse improvement of 0.1% or more
    test2 <- ifelse(iter >= 5, !any(sse_rel_change[iter - 0:2] < -.001), FALSE)

    if(test1 | test2) {
      # exit if good
      print(sprintf("exit at iteration: %i, sse: %.5e, max_rel_err: %.5e", iter, sse, max_rel_err))
      break
    }
    
    # if sse is > prior sse, adjust step scale downward
    if(step_method=="adhoc" & (sse > best_sse)){
      step_scale <- step_scale * .5
      ebeta <- best_ebeta # reset and try again
    }
    
    prior_ebeta <- ebeta
    
    # ad hoc step
    # step <- -(1 / ews) * d %*% invxpx * step_scale
    step_inputs$step_scale <- step_scale
    step <- step_fn(ebeta, step_inputs) #  * (1 - iter /maxiter) # * step_scale # * (1 - iter /maxiter)
    # print(step)
    
    ebeta <- ebeta - step
  }
  
  best_edelta <- get_delta(ewh, best_ebeta, xmat)
  ewhs <- get_weights(best_ebeta, best_edelta, xmat)
  ewh <- rowSums(ewhs)
  if(scale==TRUE) etargets <- sweep(etargets, 2, problem$scale_factor, "/")
  final_step_scale <- step_scale
  
  t2 <- proc.time()
  total_seconds <- as.numeric((t2 - t1)[3])
  
  keepnames <- c("total_seconds", "maxiter", "iter", "max_rel_err", "sse", "sse_vec", "d", "best_ebeta", "best_edelta", "ewh", "ewhs", "etargets",
                 "problem_unscaled", "scale", "scale_goal", "init_step_scale", "final_step_scale")
  result <- list()
  for(var in keepnames) result[[var]] <- get(var)
  print("all done")
  result
  # end solve_poisson
}


grad_sse <- function(betavec, wh, xmat, targets){
  # return gradient of the sse function wrt each beta
  
  # get the deltas as we will need them
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  
  # make a data frame for each relevant variable, with h, k, and/or s indexes as needed
  h <- nrow(xmat)
  k <- ncol(xmat)
  s <- nrow(targets)
  
  hstub <- tibble(h=1:h)
  skstub <- expand_grid(s=1:s, k=1:k) %>% arrange(k, s)
  hkstub <- expand_grid(h=1:h, k=1:k) %>% arrange(k, h)
  
  diffs <- diff_vec(betavec, wh, xmat, targets)
  diffsdf <- skstub %>% mutate(diff=diffs)
  
  whdf <- hstub %>% mutate(wh=wh)
  
  xdf <- hkstub %>% mutate(x=as.vector(xmat))
  
  etargets <- etargs_vec(betavec, wh, xmat, s)
  
  targdf <- skstub %>% mutate(target=as.vector(targets))
  etargdf <- skstub %>% mutate(etarget=etargets)
  betadf <- skstub %>% mutate(beta=betavec)
  deltadf <- hstub %>% mutate(delta=delta) 
  
  # now that the data are set up we are ready to calculate the gradient of the sse function
  # break the calculation into pieces using first the chain rule and then the product rule

  # sse = f(beta) = sum over targets [s,k] of (target - g(beta))^2
  #   where g(beta[s,k]) = sum over h(ws[h] * x[h,k]) and ws[h] is the TPC formula
  
  # chain rule for grad, for each beta[s,k] (where gprime is the partial of g wrt beta[s,k]):
  # = 2 * (target - g(beta)) * gprime(beta)
  # = 2 * diffs * gprime(beta[s,k])
  
  # for a single target[s,k]:
  #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
  #     Re-express:
  #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
  #      
  #          = sum[h]  of X * a * b where a=exp(beta *X) and b=exp(delta[h]) and delta is a function of beta
  
  # product rule, still for a single target, gives:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  
  # a = exp(beta * x)
  adf_base <- xdf %>%
    left_join(betadf, by="k") %>%
    mutate(a_exponent=beta * x) %>%
    select(h, s, k, x, beta, a_exponent)
  # adf_base %>% filter(h==1)
  
  adf <- adf_base %>%
    group_by(s, h) %>%
    summarise(a=exp(sum(a_exponent)), .groups="drop") %>% # these are the state weights for each hh BEFORE delta impact
    select(h, s, a)
  # adf %>% filter(h==1)
  
  #    aprime:
  #      deriv of a=exp(beta * x) wrt beta is = x * a
  # this is how much each hh's state weight will change if a beta changes, all else equal, BEFORE delta impact
  aprimedf <- adf %>%
    left_join(xdf, by="h") %>%
    mutate(aprime=x * a) %>%
    select(h, s, k, x, a, aprime) %>%
    arrange(h, s, k)
  # aprimedf %>% filter(h==1)
  
  # b = exp(delta[h])
  bdf <- deltadf %>%
    mutate(b=exp(delta))
  
  # check b: -- good
  # bcheck <- adf %>%
  #   left_join(whdf, by = "h") %>%
  #   group_by(h) %>%
  #   summarise(wh=first(wh), a=sum(a), .groups="drop") %>%
  #   mutate(bcheck=wh / a)
  
  # bprimedf # do this for each hh for each target I think
  # this is how much the delta impact will change if we change a beta - thus we have 1 per h, s, k
  # delta =log(wh/log(sum[s] exp(betas*X))
  # b=exp(delta(h))
  # which is just b = wh / log(sum[s] exp(betas*X))
  # bprime= for each h, for each beta (ie each s-k combination): (according to symbolic differentiation checks):
  #   bprime =  - (wh * xk *exp(bs) / sum[s] exp(bs))^2  where bs is exp(BX) for just that S and just that h
  # note that this bs is the same as a above: the sum, for an s-h combo, of exp(BX)
  
  # for each h, get the sum of their exp(beta * x) as it is a denominator; this is in adf
  # adf %>% filter(h==1)
  asums <- adf %>%
    group_by(h) %>%
    summarise(asum=sum(a), .groups="drop")
  
  bprimedf_base <- adf %>%
    left_join(whdf, by = "h") %>%
    left_join(xdf, by="h") %>%
    left_join(asums, by="h") %>%
    select(h, s, k, wh, x, a, asum)
  # bprimedf_base %>% filter(h==1)
  
  bprimedf <- bprimedf_base %>%
    mutate(bprime= -(wh * x * a / (asum^2)))
  # bprimedf %>% filter(h==1)
  
  # now get gprime:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  #  bprimedf has most of what we need
  gprime_h <- bprimedf %>%
    select(-a) %>% # drop a as it is also in aprime, below
    left_join(bdf, by="h") %>%
    left_join(aprimedf %>% select(-x), by = c("h", "s", "k")) %>% # drop x as it is in bprime
    mutate(gprime_h=x * (a * bprime + b * aprime))
  # gprime_h %>% filter(h==1)
  
  gprime <- gprime_h %>%
    group_by(s, k) %>%
    summarise(gprime=sum(gprime_h), .groups="drop") # sum over the households h
  
  # put it all together to get the gradient by s, k
  # 2 * (target - g(beta)) * gprime(beta)    
  graddf <- diffsdf %>%
    left_join(gprime, by = c("s", "k")) %>%
    mutate(grad=2 * diff * gprime) %>%
    arrange(k, s)
  # graddf <- targdf %>%
  #   left_join(etargdf, by = c("s", "k")) %>%
  #   left_join(gprime, by = c("s", "k")) %>%
  #   mutate(term1=2 * target * gprime,
  #          term2=2 * etarget * gprime)
  
  
  graddf$grad
}


grad_sse_v2 <- function(betavec, wh, xmat, targets){
  # return gradient of the sse function wrt each beta
  
  # get the deltas as we will need them
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  
  # make a data frame for each relevant variable, with h, k, and/or s indexes as needed
  h <- nrow(xmat)
  k <- ncol(xmat)
  s <- nrow(targets)
  
  hstub <- tibble(h=1:h)
  skstub <- expand_grid(s=1:s, k=1:k) %>% arrange(k, s)
  hkstub <- expand_grid(h=1:h, k=1:k) %>% arrange(k, h)
  
  diffs <- diff_vec(betavec, wh, xmat, targets)
  etargets <- etargs_vec(betavec, wh, xmat, s)
  diffsdf <- skstub %>% 
    mutate(target=as.vector(targets),
           etarget=etargets,
           diff=diffs)
  
  whdf <- hstub %>% mutate(wh=wh)
  
  xdf <- hkstub %>% mutate(x=as.vector(xmat))
  
  betadf <- skstub %>% mutate(beta=betavec)
  deltadf <- hstub %>% mutate(delta=delta) 
  
  # now that the data are set up we are ready to calculate the gradient of the sse function
  # break the calculation into pieces using first the chain rule and then the product rule
  
  # sse = f(beta) = sum over targets [s,k] of (target - g(beta))^2
  #   where g(beta[s,k]) = sum over h(ws[h] * x[h,k]) and ws[h] is the TPC formula
  
  # for each target, chain rule for grad, 
  # for each beta[s,k] (where gprime is the partial of g wrt beta[s,k]):
  # = - 2 * (target - g(beta)) * gprime(beta)
  # = - 2 * diffs * gprime(beta[s,k])
  
  # for a single target[s,k]:
  #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
  
  #     Re-express g(beta):
  #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
  #      
  #          = sum[h]  of X * a * b where a=exp(beta *X) and b=exp(delta[h]) and delta is a function of beta
  
  # gprime(beta), still for a single target -- product rule gives:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  
  # a = exp(beta * x)  - this is the same for all
  adf_base <- xdf %>%
    left_join(betadf, by="k") %>%
    mutate(a_exponent=beta * x) %>% # we will raise e to this power
    select(h, s, k, x, beta, a_exponent)
  # adf_base %>% filter(h==1)
  
  adf <- adf_base %>%
    group_by(s, h) %>% # weights are specific to state and household
    # we sum the k elements of the exponent and then raise e to that power
    summarise(a=exp(sum(a_exponent)), .groups="drop") %>% # these are the state weights for each hh BEFORE delta impact
    select(h, s, a)
  # adf %>% filter(h==1)
  
  #    aprime:
  #      deriv of a=exp(beta * x) wrt beta is = x * a
  # this is how much each hh's state weight will change if a beta changes, all else equal, BEFORE delta impact
  # since there is a beta for each s, k combination this will vary with different x[h, k] values
  aprimedf <- adf %>%
    left_join(xdf, by="h") %>%
    mutate(aprime=x * a) %>%
    select(h, s, k, x, a, aprime) %>%
    arrange(h, s, k)
  # aprimedf %>% filter(h==1)
  
  # b = exp(delta[h])
  bdf <- deltadf %>%
    mutate(b=exp(delta))
  
  # bprime - the hardest part -- how much does delta for an h change wrt a change in any beta
  # this is how much the delta impact will change if we change a beta - thus we have 1 per h, s, k
  # delta =log(wh/log(sum[s] exp(betas*X))
  # b=exp(delta(h))
  # which is just b = wh / log(sum-over-s]: exp(beta-for-given-s * x))
  
  # bprime= for each h, for each beta (ie each s-k combination):
  
  # from symbolic differentiation we have:
  # .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2) # which is a for a specific state -- s1 in this case
  # .e9 <- .e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2)
  # so e9 <- exp(b.s1k1 * x1 + b.s1k2 * x2) + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2)
  # which is just asum as calculated below
  #  -(wh * x1 * .e3/(.e9 * log(.e9)^2))
  
  #  .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2)
  #   -(wh * x1 * .e3 / (.e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2))^2)
  
  # which in a, asum notation is:
  #  -(wh * x1 * a / (asum)^2)
  
  # in my a, asum notation below, we have
  # deriv wrt b.s1k1 =
  #    -(wh * x[k=1] * a[s=1] / {asum^2})
  
  # for each h, get the sum of their exp(beta * x) as it is a denominator; this is in adf
  # log(wh / colSums(beta_x)) # denominator is sum for each person
  
  # asum is, for each h, exp(beta-s1k1 +...+ beta-s1kn) + ...+ exp(beta-smk1 + ...+ beta-smkn)
  # each a that we start with is exp(.) for one of the states so this the sum of exp(.) over the states
  asums <- adf %>%
    group_by(h) %>%
    summarise(asum=sum(a), .groups="drop")
  
  bprimedf_base <- adf %>%
    left_join(whdf, by = "h") %>%
    left_join(xdf, by="h") %>%
    left_join(asums, by="h") %>%
    select(h, s, k, wh, x, a, asum)
  # bprimedf_base %>% filter(h==1)
  
  # deriv wrt b.s1k1 =
  #    -(wh * x[k=1] * a[s=1] / {asum^2})
  
  # bprime -- how much does delta for an h change wrt a change in any beta
  bprimedf <- bprimedf_base %>%
    mutate(bprime= -(wh * x * a / {asum^2})) %>%
    arrange(h, k, s)
  # bprimedf %>% filter(h==1)
  
  # now get gprime:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  #  bprimedf has most of what we need
  gprime_h <- bprimedf %>%
    select(-a) %>% # drop a as it is also in aprime, below
    left_join(bdf, by="h") %>%
    left_join(aprimedf %>% select(-x), by = c("h", "s", "k")) %>% # drop x as it is in bprime
    mutate(gprime_h= x * (a * bprime + b * aprime))
  # gprime_h %>% filter(h==1)
  
  gprime <- gprime_h %>%
    group_by(s, k) %>%
    summarise(gprime=sum(gprime_h), .groups="drop") %>% # sum over the households h
    arrange(k, s)

  # put it all together to get the gradient by s, k
  # - 2 * (target - g(beta)) * gprime(beta) FOR EACH TARGET AND ADD THEM UP
  # diffs; gprime now we need to cross each gprime with all distances
  grad_base <- expand_grid(s.d=1:s, s.k=1:k, s.gp=1:s, k.gp=1:k) %>%
    left_join(diffsdf %>% select(s, k, diff) %>% rename(s.d=s, s.k=k), by = c("s.d", "s.k")) %>%
    left_join(gprime %>% select(s, k, gprime) %>% rename(s.gp=s, k.gp=k), by = c("s.gp", "k.gp")) %>%
    mutate(grad=-2 * diff * gprime)
  
  graddf <- grad_base %>%
    group_by(s.gp, k.gp) %>%
    summarise(grad=sum(grad), .groups="drop")
  
  graddf$grad
}


p <- make_problem(h=2, k=2, s=2)
p <- make_problem(h=4, k=2, s=3) # use this
p <- make_problem(h=20, k=4, s=8)

sval <- rep(0, p$s * p$k)
# sval <- runif(p$s*p$k)
betavec <- rep(0, p$s * p$k)
targets <- p$targets
xmat <- p$xmat
wh <- p$wh
beta <- vtom(betavec, p$s)
delta <- get_delta(wh, beta, xmat)
whs <- get_weights(beta, delta, xmat)

# make targets that are all hit, except 1
etargets <- t(whs) %*% xmat
targets <- etargets
row <- 1; col <- 1
targets[row, col] <- etargets[row, col] + 1
diff_vec(betavec, wh, xmat, targets)

grad(sse_fn, x=sval, wh=p$wh, xmat=p$xmat, targets=targets) %>% round(4)
grad_sse_v2(sval, wh=p$wh, xmat=p$xmat, targets=targets) %>% round(4)

jacobian(diff_vec, x=sval, wh=p$wh, xmat=p$xmat, targets=targets)
gprime$gprime
# the diagonal of the jacobian is right...except for the sign



jac <- function(betavec, wh, xmat, targets){
  # return gradient of the sse function wrt each beta
  
  # get the deltas as we will need them
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  
  # make a data frame for each relevant variable, with h, k, and/or s indexes as needed
  h <- nrow(xmat)
  k <- ncol(xmat)
  s <- nrow(targets)
  
  hstub <- tibble(h=1:h)
  skstub <- expand_grid(s=1:s, k=1:k) %>% arrange(k, s)
  hkstub <- expand_grid(h=1:h, k=1:k) %>% arrange(k, h)
  
  # we need functions that get the s or k associated with a given index
  get_s <- function(idx) {skstub$s[idx]}
  get_k <- function(idx) {skstub$k[idx]}
  
  xvec <- as.vector(xmat) # this is in the proper order 
  xdf <- hkstub %>% mutate(x=as.vector(xmat))
  xmat
  diffs <- diff_vec(betavec, wh, xmat, targets)
  
  # create a long form of the jacobian with i indexing rows [differences],
  # and j indexing columns [betas]
  # each element is a partial derivative of d[i] wrt beta[j] and will depend on the h's and k's also
  
  irows_diff <- expand_grid(s.i=1:s, k.i=1:k) %>% 
    arrange(k.i, s.i) %>%
    mutate(i=row_number(), diff=diffs, beta.diff=betavec) %>%
    select(i, s.i, k.i, diff, beta.diff)
  
  jcols_beta <- expand_grid(s.j=1:s, k.j=1:k) %>% 
    arrange(k.j, s.j) %>%
    mutate(j=row_number(), beta.pd=betavec) %>% # we'll need the beta from this iteration
    select(j, s.j, k.j, beta.pd)
  
  # the long version of the jacobian
  jlong_stub <- expand_grid(i=1:(s * k), j=1:(s * k)) %>%
    left_join(irows_diff, by = "i") %>%
    left_join(jcols_beta, by = "j")
  
  
  get_pd <- function(df){
    # get minus gprime(beta) partial deriviative of diff-i (diff) wrt beta-j (beta.pd)
    
    # where
    #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
    
    #     Re-express g(beta):
    #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
    #      
    #          = sum[h]  of X * a * b where a=exp(beta *X) and b=exp(delta[h]) and delta is a function of beta
    
    # gprime(beta), still for a single target -- product rule gives:
    #     gprime(beta)=sum[h] of X[k.d] * (a * bprime + b * aprime)
    # get the x values for all households for this target
    
    
    x_k <- tibble(x=xmat[, df$k.i]) # the x values for this target's k, one value per household
    beta_s <- tibble(k=1:k, beta_s.i=beta[df$s.i, ]) # the beta values for this target's state, one value per k
    # df <- df %>%
    #   mutate(xsum=sum(xmat[, get_k(.$i)])) # we use x[, k] for this target thus we get i, not j
    
    # # a = exp(beta * x)  - this is the same for all
    # adf_base <- xdf %>%
    #   left_join(betadf, by="k") %>%
    #   mutate(a_exponent=beta * x) %>% # we will raise e to this power
    #   select(h, s, k, x, beta, a_exponent)
    # adf_base %>% filter(h==1)
    # 
    # adf <- adf_base %>%
    #   group_by(s, h) %>% # weights are specific to state and household
    #   # we sum the k elements of the exponent and then raise e to that power
    #   summarise(a=exp(sum(a_exponent)), .groups="drop") %>% # these are the state weights for each hh BEFORE delta impact
    #   select(h, s, a)
    # adf %>% filter(h==1)
    amat <- xmat %*% t(beta) # each row has the a_exponent values, 1 per state, for a household
    a <- exp(amat)
    asums <- rowSums(a)
    
    # a <- exp(xmat %*% beta[df$s.i, ]) # one value per h, weight[h, s] calculated without delta
    # we also need asum - the sum of a over all states for this person
    
    # asums <- adf %>%
    #   group_by(h) %>%
    #   summarise(asum=sum(a), .groups="drop")
    # 
    # bprimedf_base <- adf %>%
    #   left_join(whdf, by = "h") %>%
    #   left_join(xdf, by="h") %>%
    #   left_join(asums, by="h") %>%
    #   select(h, s, k, wh, x, a, asum)
    # bprimedf_base %>% filter(h==1)
    
    # deriv wrt b.s1k1 =
    #    -(wh * x[k=1] * a[s=1] / {asum^2})
    
    # bprime -- how much does delta for an h change wrt a change in any beta
    # bprimedf <- bprimedf_base %>%
    #   mutate(bprime= -(wh * x * a / {asum^2})) %>%
    #   arrange(h, k, s)
    
    ab <- x_k %>%
      left_join(df, by = character()) %>%
      mutate(wh=wh,
             a=a[, df$s.i], 
             aprime=x * a,
             b=exp(delta),
             asum=asums,
             bprime = -(wh * x * a / {asum^2}),
             gprime_h = x * (a * bprime + b * aprime))
    # beta_s %>%
    #   left_join(df, by = character())
    ab
  }
  
  df <- jlong_stub %>% filter(i==1, j==1)
  get_pd(df)
  
  
  tmp <- jlong_stub %>%
    filter(i %in% c(1, 4), j %in% 1) %>%
    group_by(i, j) %>%
    group_modify(~get_pd(.x))
  tmp
  tmp %>%
    group_by(i, j) %>%
    summarise(gprime=sum(gprime_h), .groups="drop")
  

  diffs <- diff_vec(betavec, wh, xmat, targets)
  etargets <- etargs_vec(betavec, wh, xmat, s)
  
  # start making data frames
  diffsdf <- skstub %>% 
    mutate(target=as.vector(targets),
           etarget=etargets,
           diff=diffs)
  
  whdf <- hstub %>% mutate(wh=wh)
  
  betadf <- skstub %>% mutate(beta=betavec)
  deltadf <- hstub %>% mutate(delta=delta) 
  
  # now that the data are set up we are ready to calculate the gradient of the sse function
  # break the calculation into pieces using first the chain rule and then the product rule
  
  # sse = f(beta) = sum over targets [s,k] of (target - g(beta))^2
  #   where g(beta[s,k]) = sum over h(ws[h] * x[h,k]) and ws[h] is the TPC formula
  
  # for each target, chain rule for grad, 
  # for each beta[s,k] (where gprime is the partial of g wrt beta[s,k]):
  # = - 2 * (target - g(beta)) * gprime(beta)
  # = - 2 * diffs * gprime(beta[s,k])
  
  # for a single target[s,k]:
  #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
  
  #     Re-express g(beta):
  #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
  #      
  #          = sum[h]  of X * a * b where a=exp(beta *X) and b=exp(delta[h]) and delta is a function of beta
  
  # gprime(beta), still for a single target -- product rule gives:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  
  # a = exp(beta * x)  - this is the same for all
  adf_base <- xdf %>%
    left_join(betadf, by="k") %>%
    mutate(a_exponent=beta * x) %>% # we will raise e to this power
    select(h, s, k, x, beta, a_exponent)
  # adf_base %>% filter(h==1)
  
  adf <- adf_base %>%
    group_by(s, h) %>% # weights are specific to state and household
    # we sum the k elements of the exponent and then raise e to that power
    summarise(a=exp(sum(a_exponent)), .groups="drop") %>% # these are the state weights for each hh BEFORE delta impact
    select(h, s, a)
  # adf %>% filter(h==1)
  
  #    aprime:
  #      deriv of a=exp(beta * x) wrt beta is = x * a
  # this is how much each hh's state weight will change if a beta changes, all else equal, BEFORE delta impact
  # since there is a beta for each s, k combination this will vary with different x[h, k] values
  aprimedf <- adf %>%
    left_join(xdf, by="h") %>%
    mutate(aprime=x * a) %>%
    select(h, s, k, x, a, aprime) %>%
    arrange(h, s, k)
  # aprimedf %>% filter(h==1)
  
  # b = exp(delta[h])
  bdf <- deltadf %>%
    mutate(b=exp(delta))
  
  # bprime - the hardest part -- how much does delta for an h change wrt a change in any beta
  # this is how much the delta impact will change if we change a beta - thus we have 1 per h, s, k
  # delta =log(wh/log(sum[s] exp(betas*X))
  # b=exp(delta(h))
  # which is just b = wh / log(sum-over-s]: exp(beta-for-given-s * x))
  
  # bprime= for each h, for each beta (ie each s-k combination):
  
  # from symbolic differentiation we have:
  # .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2) # which is a for a specific state -- s1 in this case
  # .e9 <- .e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2)
  # so e9 <- exp(b.s1k1 * x1 + b.s1k2 * x2) + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2)
  # which is just asum as calculated below
  #  -(wh * x1 * .e3/(.e9 * log(.e9)^2))
  
  #  .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2)
  #   -(wh * x1 * .e3 / (.e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2))^2)
  
  # which in a, asum notation is:
  #  -(wh * x1 * a / (asum)^2)
  
  # in my a, asum notation below, we have
  # deriv wrt b.s1k1 =
  #    -(wh * x[k=1] * a[s=1] / {asum^2})
  
  # for each h, get the sum of their exp(beta * x) as it is a denominator; this is in adf
  # log(wh / colSums(beta_x)) # denominator is sum for each person
  
  # asum is, for each h, exp(beta-s1k1 +...+ beta-s1kn) + ...+ exp(beta-smk1 + ...+ beta-smkn)
  # each a that we start with is exp(.) for one of the states so this the sum of exp(.) over the states
  asums <- adf %>%
    group_by(h) %>%
    summarise(asum=sum(a), .groups="drop")
  
  bprimedf_base <- adf %>%
    left_join(whdf, by = "h") %>%
    left_join(xdf, by="h") %>%
    left_join(asums, by="h") %>%
    select(h, s, k, wh, x, a, asum)
  # bprimedf_base %>% filter(h==1)
  
  # deriv wrt b.s1k1 =
  #    -(wh * x[k=1] * a[s=1] / {asum^2})
  
  # bprime -- how much does delta for an h change wrt a change in any beta
  bprimedf <- bprimedf_base %>%
    mutate(bprime= -(wh * x * a / {asum^2})) %>%
    arrange(h, k, s)
  # bprimedf %>% filter(h==1)
  
  # now get gprime:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  #  bprimedf has most of what we need
  gprime_h <- bprimedf %>%
    select(-a) %>% # drop a as it is also in aprime, below
    left_join(bdf, by="h") %>%
    left_join(aprimedf %>% select(-x), by = c("h", "s", "k")) %>% # drop x as it is in bprime
    mutate(gprime_h= x * (a * bprime + b * aprime))
  # gprime_h %>% filter(h==1)
  
  gprime <- gprime_h %>%
    group_by(s, k) %>%
    summarise(gprime=sum(gprime_h), .groups="drop") %>% # sum over the households h
    arrange(k, s)
  
  # put it all together to get the gradient by s, k
  # - 2 * (target - g(beta)) * gprime(beta) FOR EACH TARGET AND ADD THEM UP
  # diffs; gprime now we need to cross each gprime with all distances
  grad_base <- expand_grid(s.d=1:s, s.k=1:k, s.gp=1:s, k.gp=1:k) %>%
    left_join(diffsdf %>% select(s, k, diff) %>% rename(s.d=s, s.k=k), by = c("s.d", "s.k")) %>%
    left_join(gprime %>% select(s, k, gprime) %>% rename(s.gp=s, k.gp=k), by = c("s.gp", "k.gp")) %>%
    mutate(grad=-2 * diff * gprime)
  
  graddf <- grad_base %>%
    group_by(s.gp, k.gp) %>%
    summarise(grad=sum(grad), .groups="drop")
  
  graddf$grad
}

