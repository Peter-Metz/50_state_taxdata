
# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r

# notes ----


# libraries ----
source(here::here("include", "libraries.r"))
# remotes::install_github("tidyverse/dplyr") if needed
library(numDeriv)

devtools::session_info()
(.packages()) %>% sort

# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))

# functions ----
source(here::here("include", "functions_prep_dev.r")) # soon we will replace functions_prep.r with the dev version
source(here::here("include", "functions_poisson_model.r"))


# choose which file to use 
samp1 <- readRDS(here::here("data", fns[1])) %>% 
  select(-nrecs, -pop) # note that we no longer need nrecs; pop ordinarily would not be in the data so drop here and create later
glimpse(samp1)
summary(samp1)
count(samp1, mar)
# djb note that we have nrecs and pop variables -- I want to make sure they are not needed for anything ----
# if you want to target the total number of weighted records we need a variable that is 1 for all records ----

# PREPARE DATA  ----
#.. modify the sample (don't think we need a function for this) ----
# - define income groups
# - create an indicator for each income variable as to whether it is nonzero
# _ expand categoricals into dummies as needed
# if we don't have a variable such as pop where all values are 1, we should create it as it makes it easy to get weighted record counts
samp2 <- samp1 %>%
  mutate(pid=row_number(), # pid -- an id variable for each person in the file
         incgroup=ntile(pincp, 10), # divide the data into 10 income ranges
         pop=1, # it's useful to have a variable that is 1 on every record
         # convert categoricals to dummies if we will base targets upon them
         mar1=ifelse(mar==1, 1, 0), # married
         mar5=ifelse(mar==5, 1, 0), # single
         marx15=ifelse(mar %nin% c(1, 5), 1, 0)
  )
summary(samp2)
ht(samp2)

#.. define the kinds of (weighted) targets we want and prepare the file accordingly ----
# sum:    sum of values
# nnz:    number of nonzero values
# sumneg: sum of negative values
# nneg:   number of zero values
# sumpos: sum of positive value
# npos:   number of positive values

# For the PUF the SOI data provide only the first two kinds of targets, but for the ACS we could have any of them.

# TRY TO AVOID DEPENDENT CONSTRAINTS - redundancy - as they can make the problem very hard to solve.
# For example, suppose there are 3 kinds of income (wages, interest, retirement) plus a total (sum of the 3)
#   -- don't create targets for each of the 3 kinds plus a target for the total -- leave one item out
# Another, less obvious example: don't target the total number of returns plus the number for each marital status - leave one out.

nnz_vars <- c("pop", "mar1", "mar5", "pincp", "wagp") # note that I leave the 3rd marital status out -- marx15
sum_vars <- c("pincp", "wagp", "intp", "pap", "retp", "ssip", "ssp") # DO NOT include total plus all sums - leave one out (otherincp)
sumneg_vars <- "otherincp"

# define a vector of variable names for "possible" targets (a superset) -- we may not target all
possible_target_vars <- make_target_names(
  list(nnz_vars, sum_vars, sumneg_vars),
  c("nnz", "sum", "sumneg"))
possible_target_vars

# prepare data by creating variables with those names:
#   nnz, nneg, and npos will have 1 for rows where the respective variable is nz, neg, or pos, respectively, and 0 otherwise
#   sum will have its respective variable's value
#   sumneg and sumpos will have the variable's value if negative or positive, respectively, and 0 otherwise
samp <- prep_data(samp2, possible_target_vars)
glimpse(samp)


summary_vals <- get_summary_vals(samp, .weight=pwgtp, .sum_vars=possible_target_vars, stabbr, incgroup)
summary_vals

# Create a data frame with all targets for all states and income groups ----
# for the PUF, we will create this using information from Historical Table 2
# for the ACS, we construct the targets from the ACS data itself
all_target_values <- summary_vals


# wrap everything we need for a single income group into a function that returns a list ----

# SINGLE INCOME GROUP ----
#.. define target incgroup, target variable names, and target values for each state in the income group ----
target_incgroup <- 2 # define target income group

possible_target_vars
(target_vars <- possible_target_vars[c(1, 3, 6, 7)])
target_vars <- setdiff(possible_target_vars, c("pap_sum", "ssip_sum", "intp_sum", "otherincp_sumneg"))
target_vars <- setdiff(possible_target_vars, c("pap_sum", "ssip_sum", "intp_sum"))
target_vars <- setdiff(possible_target_vars, c("pap_sum", "ssip_sum"))
target_vars <- setdiff(possible_target_vars, c("pap_sum"))
target_vars <- possible_target_vars

# define target values and states, for this income group
targets_wide <- all_target_values %>%
  filter(incgroup==target_incgroup) %>%
  select(stabbr, incgroup, nrecs, all_of(target_vars)) # a small list of variables to target; we have nrecs because we created it in summary_vals
targets_wide # these are the targets we want to hit
summary(targets_wide)

hweights <- samp %>% filter(incgroup==target_incgroup) %>% .$pwgtp
targets <- targets_wide[, target_vars] %>% as.matrix
xmat <- samp %>% filter(incgroup==target_incgroup) %>% .[, target_vars] %>% as.matrix

targets
xmat

(scale <- 10000 / colSums(targets))
stargets <- targets * scale
sxmat <- sweep(xmat, 2, scale, "*")

scale
xmat[1:5, ]
sxmat[1:5, ]

targets
stargets


get_delta <- function(wh, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(wh / colSums(beta_x))
}

get_weights <- function(beta, delta, x){
  # get all weights
  beta_x <- beta %*% t(x)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(m) m + delta) 
  exp(beta_xd)
}

f <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}

get_result <- function(){
  result <- list()
  result$iter <- i
  result$sse <- sse
  result$d <- d
  result$ebeta <- ebeta
  result$edelta <- edelta
  result$whs <- get_weights(ebeta, edelta, sxmat)
  result
}

stargets <- targets
sxmat <- xmat

xpx <- t(sxmat) %*% sxmat
invxpx <- solve(xpx) # TODO: add error check and exit if not invertible

beta0 <- matrix(0, nrow=nrow(stargets), ncol=ncol(stargets)) # tpc uses 0 as beta starting point
delta0 <- get_delta(hweights, beta0, sxmat) # tpc uses initial delta based on initial beta 

ebeta <- beta0 # tpc uses 0 as beta starting point
edelta <- delta0 # tpc uses initial delta based on initial beta 

maxiter <- 50
for(i in 1:maxiter){
  print(i)
  ewhs <- get_weights(ebeta, edelta, sxmat)
  ews <- colSums(ewhs)
  ewh <- rowSums(ewhs)
  
  etargets <- t(ewhs) %*% sxmat
  d <- stargets - etargets
  sse <- sum(d^2)
  if(sse < 1e-8) {
    # exit if good
    result <- get_result()
    break
  }
  
  jval <- jacobian(f, x=as.vector(ebeta), wh=wh, xmat=sxmat, targets=stargets)
  step <- solve(jval, tol = 1e-30) %*% as.vector(d)
  step <- matrix(step, nrow=nrow(d), byrow=FALSE)
  
  # step <- -(1 / ews) * d %*% invxpx * 10000
  
  ebeta <- ebeta - step
  edelta <- get_delta(ewh, ebeta, sxmat)
  if(i==maxiter) {result <- get_result(); break}
}

# str(result_
result$iter
result$sse
result$d

