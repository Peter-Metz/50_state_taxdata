
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


vtom <- function(vec, nrows){
  matrix(vec, nrow=nrows, byrow=FALSE)
}


# start ----
k <- 2
s <- 3
h <- 4
ks <- k * s

# h x k
xhs <- matrix(c(1, 2, 
               3, 4,
               5, 6,
               7, 8), ncol=k, byrow=TRUE)
xhs

wh <- c(6, 15, 21, 30)

beta <- matrix(c(.1, .2,
                 .3, .4,
                 .5, .6), ncol=k, byrow=TRUE)

delta <- get_delta(wh, beta, xhs)
delta

whs <- get_weights(beta, delta, xhs)
whs # h x s

# 1 x h
rowSums(whs)
# 1 x s
(ws <- colSums(whs))

etargs <- t(whs) %*% xhs
etargs

targs <- matrix(c(4, 8,
                  25, 40,
                  320, 400), ncol=k, byrow=TRUE)
targs


dm <- targs - etargs
dm

dv <- as.vector(dm)
dv

bv <- as.vector(beta)

f <- function(betavec, delta, wh, xmatrix, targets){
  # return the distances from targets as a vector using the OLD delta
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  whs <- calc_weights(beta, delta, xmatrix)
  etargets <- t(whs) %*% xmatrix
  d <- targets - etargets
  as.vector(d)
}

jval <- jacobian(f, x=as.vector(beta), delta=delta, wh=wh, xmatrix=xhs, targets=targs)
jval

# xpx <- t(xhs) %*% xhs
# xxp <- xhs %*% t(xhs)

# this is it ----
(x2 <- xhs * xhs)
(x12 <- xhs[, 1] * xhs[, 2])
x2; x12

(p1 <- (t(whs) %*% x2) %>% as.vector)
(p2 <- t(whs) %*% x12)

jval

# do it efficiently
(p1 <- (t(whs) %*% x2))
(p2 <- t(whs) %*% x12)
t(whs) %>% cbind(p1, p2)

t(whs[, 1]) %*% cbind(x2, x12)
t(whs[, 2]) %*% cbind(x2, x12, x2)

diag(as.vector(t(whs) %*% x2))
p1 <- t(whs) %*% x2
p2 <- t(whs) %*% x12
diag(p1[, 1])
diag(p1[, 2])
diag(p2[, 1])

jac <- bind_rows()

