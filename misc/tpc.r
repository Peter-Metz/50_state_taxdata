
source(here::here("include", "libraries.r"))
library(mvtnorm)

<<<<<<< HEAD

# define needed functions ----
get_delta <- function(w, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(w / colSums(beta_x))
}

get_weights <- function(beta, delta, x){
  # get all weights
  beta_x <- beta %*% t(x)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(m) m + delta) 
  exp(beta_xd)
}

=======
>>>>>>> ac916ea2b2f63af4679647d64a1823ff69689c5a
# create the poisson problem ----
# xscale <- 10 and step_scale <- 10 works well with this
h <- 10 # number households
k <- 3 # number characteristics
s <- 5 # number states

# x is an h x k matrix of characteristics
sigma <- matrix(c(1, .8, .6,
                  .8, 1, .4,
                  .6, .4, 1), ncol=3, byrow = TRUE)
sigma

<<<<<<< HEAD
set.seed(1234)
x <- rmvnorm(n=h, mean=c(10, 20, 30), sigma)

=======
>>>>>>> ac916ea2b2f63af4679647d64a1823ff69689c5a
# for the poisson version we need beta and delta

# beta is an s x k matrix of coefficients
beta <- matrix(c(1, 2, 3,
                 4, 5, 6,
                 7, 8, 9,
                 10, 11, 12,
                 13, 14, 15), ncol=3, byrow = TRUE) / 100
beta

# delta is an h-length vector of individual constants in the formula for state weights
delta <- 1:h / 100
<<<<<<< HEAD
=======



# Alterantively, create a large problem ----
# xscale <- 100 and step_scale <- nrow(x) works well with this
h <- 40e3 # number households
k <- 30 # number characteristics
s <- 50 # number states

# alternative way to get x for large or small problem ----
xbar <- seq(10, 100, length.out=k)
xsd <- xbar / 10

set.seed(1234)
x <- matrix(rnorm(n=h * k, mean=xbar, sd=xsd), nrow=h, ncol=k)


# look at x and its correlations ----
x
cor(x)


# define needed functions ----
get_delta <- function(w, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(w / colSums(beta_x))
}
>>>>>>> ac916ea2b2f63af4679647d64a1823ff69689c5a

get_weights(beta, delta, x)


<<<<<<< HEAD
# Alterantively, create a large problem ----
# xscale <- 100 and step_scale <- nrow(x) works well with this
h <- 40e3 # number households
k <- 30 # number characteristics
s <- 50 # number states

# alternative way to get x for large or small problem ----
xbar <- seq(10, 100, length.out=k)
xsd <- xbar / 10

set.seed(1234)
x <- matrix(rnorm(n=h * k, mean=xbar, sd=xsd), nrow=h, ncol=k)


# look at x and its correlations ----
x
cor(x)


=======
>>>>>>> ac916ea2b2f63af4679647d64a1823ff69689c5a
# now define individual by state weights one of 2 ways ----
# either from a poisson model, or randomly (especially for a large problem)
hsweights <- get_weights(beta, delta, x) # poisson model

#.. randomly IF not poisson ----
set.seed(1234)
hsweights <- matrix(runif(n=h * s, min=4, max=8000), nrow=h, ncol=s)
<<<<<<< HEAD
#.. end random ----
=======
>>>>>>> ac916ea2b2f63af4679647d64a1823ff69689c5a

hsweights

# get total 
hweights <- rowSums(hsweights) # sum of weights for each individual
hweights
sweights <- colSums(hsweights) # sum of weights for each state
sweights

# done with data setup ----


# get problem ready to solve ----
xsave <- x # save initial x values in case we scale them


# scale and go ----
xscale <- 100
x <- xsave / xscale

targets <- t(hsweights) %*% x
targets

# now search for the beta coefficients and deltas that are consistent with these data and that meet the targets
# we know the data and we know the aggregate weights, but not the state weights -- e prefix will mean estimated

# define xprime-x and its inverse before entering the loop as it will not change within the loo9p
xpx <- t(x) %*% x
# inverse of matrix that has been multiplied by a non-zero scalar equals inverse of the scalar multiplied by inverse of the matrix,
# so solve xpx once at the start
(ixpx <- solve(xpx)) # this needs to be invertible
# to get the inverse of the jacobian we will multiply by the inverse of the scalar


# define initial values before entering loop ----
ibeta <- matrix(0, nrow=s, ncol=k) # tpc uses 0 as beta starting point
idelta <- get_delta(hweights, ibeta, x) # tpc uses initial delta based on initial beta 

ihsweights <- get_weights(ibeta, idelta, x)
ihsweights

# look at distance from targets using initial values
idist <- targets - t(ihsweights) %*% x
idist

isse <- sum(idist^2)
isse
 
# before start, set values to use when entering loop for first time
ebeta <- ibeta
edelta <- idelta
<<<<<<< HEAD
# step_scale <- 1e4 # note that this has a MAJOR impact on iterations
=======
step_scale <- 1e4 # djb note that this has a MAJOR impact on iterations
>>>>>>> ac916ea2b2f63af4679647d64a1823ff69689c5a
step_scale <- nrow(x)

a <- proc.time()
for(i in 1:500){
  ehsweights <- get_weights(ebeta, edelta, x) # estimated weights for each household for each state
  ehweights <- rowSums(ehsweights) # estimated total weights for each individual
  esweights <- colSums(ehsweights) # estimated total weights for each state
  
  etargets <- t(ehsweights) %*% x
  dist <- targets - etargets
  sse <- sum(dist^2) # a measure of error
  
  if(sse < 1e-8) {
    print(sprintf("DONE at iteration %i: sse is %.5e", i, sse))
    break
  }
  
  if(i<=10 | i %% 5 ==0) {
    print(sprintf("iteration %i: sse is %.5e", i, sse))
  }
  
  # get the step
  step <- matrix(nrow=s, ncol=k)
  for(i in 1:s) step[i, ] <- t((1 /esweights[i]) * ixpx %*% dist[i, ]) * step_scale
  
  ebeta <- ebeta + step
  edelta <- get_delta(ehweights, ebeta, x)
}
b <- proc.time()
b - a

sum(dist^2)
dist
quantile(dist)

fbeta <- ebeta
fdelta <- edelta
fhsweights <- get_weights(fbeta, fdelta, x)

fbeta %>% round(3)
beta %>% round(3)
ibeta %>% round(3)

fdelta %>% round(3)
delta %>% round(3)
idelta %>% round(3)

fhsweights %>% round(3)
hsweights %>% round(3)
ihsweights %>% round(3)
(fhsweights - hsweights) %>% round(6)
fhsweights / hsweights
fhsweights / ihsweights
<<<<<<< HEAD
quantile((fhsweights - hsweights) / hsweights * 100 - 100)
=======
>>>>>>> ac916ea2b2f63af4679647d64a1823ff69689c5a

rowSums(fhsweights)
hweights
(rowSums(fhsweights) / hweights * 100 - 100) %>% round(3)

colSums(fhsweights)
sweights
(colSums(fhsweights) / sweights * 100 - 100) %>% round(3)


t1 <- targets * xscale
t2 <- (t(fhsweights) %*% x) * xscale
t1 %>% round(3)
t2 %>% round(3)
(t2 - t1) %>% round(3)
(t2 / t1 * 100 - 100) %>% round(4)


# regression
# y <- as.vector(hsweights) ...
# x1 <- x
# glm(formula = num_awards ~ prog + math, family = "poisson", data = p)

