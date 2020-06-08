library(Deriv)

delta_h <- expression(log(wh / (exp(b.s1k1*x1 + b.s1k2*x2) + exp(b.s2k1*x1 + b.s2k2*x2))))
Deriv(delta_h, "b.s2k1")


# get derivative of delta_h wrt a b.sk for a given household
# let s=4 and k=3
delta_h <- expression(log(wh / (  exp(b.s1k1*x1 + b.s1k2*x2 + b.s1k3*x3) +
                                  exp(b.s2k1*x1 + b.s2k2*x2 + b.s2k3*x3) +
                                  exp(b.s3k1*x1 + b.s3k2*x2 + b.s3k3*x3) +
                                  exp(b.s4k1*x1 + b.s4k2*x2 + b.s4k3*x3)
                                )))
Deriv(delta_h, "b.s2k1")

b <- expression(wh / (  exp(b.s1k1*x1 + b.s1k2*x2 + b.s1k3*x3) +
                          exp(b.s2k1*x1 + b.s2k2*x2 + b.s2k3*x3) +
                          exp(b.s3k1*x1 + b.s3k2*x2 + b.s3k3*x3) +
                          exp(b.s4k1*x1 + b.s4k2*x2 + b.s4k3*x3)
                        ))
Deriv(b, "b.s2k1")



b <- expression(wh / (  exp(b.s1k1*x1 + b.s1k2*x2) +
                          exp(b.s2k1*x1 + b.s2k2*x2) +
                          exp(b.s3k1*x1 + b.s3k2*x2) +
                          exp(b.s4k1*x1 + b.s4k2*x2)
                        ))
Deriv(b, "b.s2k1")

-(wh * x1 * .e3/(exp(b.s1k1 * x1 + b.s1k2 * x2) +
                   .e3 +
                   exp(b.s3k1 * x1 + b.s3k2 * x2) +
                   exp(b.s4k1 * x1 + b.s4k2 * x2))^2)

Deriv(b, "b.s3k1")
  .e3 <- exp(b.s3k1 * x1 + b.s3k2 * x2)
  -(wh * x1 * .e3/(
    exp(b.s1k1 * x1 + b.s1k2 * x2) +
      exp(b.s2k1 * x1 + b.s2k2 * x2) +
      .e3 + 
      exp(b.s4k1 * x1 + b.s4k2 * x2))^2)

d.s1k1 <- expression((t1 - (x.h1k1*exp(b.s1k1*x.h1k1 + b.s1k2 * x.h1k2 + d.h1))^2))
# d.s1k1 <- expression((t1 - (x.h1k1*exp(b.s1k1*x.h1k1 + b.s1k2 * x.h1k2 + d.h1)

x <- 100
log(x)
exp(log(x))

                            
b <- expression((t1 - x*exp(b1*x))^2 + (t2 - x*exp(b2*x))^2)
Deriv(b, "b1")      

b <- expression((t1 - x^3)^2)
Deriv(b, "x")

b <- expression((t1 - x)^2)
Deriv(b, "x")

                              
                              
#                              x.h2k1*exp(b.s1k1*x.h2k1 + b.s1k2 * x.h2k2 + d.h2))^2))

d.s1k1 <- expression((t1 - (x.h1k1*exp(b.s1k1*x.h1k1 + b.s1k2 * x.h1k2 + 
                                         log(w.h1 / exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2))) +
                              x.h2k1*exp(b.s1k1*x.h2k1 + b.s1k2 * x.h2k2 + d.h2))^2))

d.s1k1 <- expression(
  # difference 1
  (t1 - (
  # household 1
    x.h1k1*exp(b.s1k1*x.h1k1 + b.s1k2 * x.h1k2 +  # h1 weight for s1 before delta1
                 log(w.h1 / (exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2) +
                               exp(b.s2k1*x.h1k1 + b.s2k2*x.h1k2)))) + # delta1
      # household 2
      x.h2k1*exp(b.s1k1*x.h2k1 + b.s1k2 * x.h2k2 +  # h2 weight for s1 before delta2
                   log(w.h2 / (exp(b.s1k1*x.h2k1 + b.s1k2*x.h2k2) +
                                 exp(b.s2k1*x.h2k1 + b.s2k2*x.h2k2)))) # delta2
                              )^2))

Deriv(d.s1k1, "b.s1k1")

# t1 is s1k1 and t2 is s2k1
d.2 <- expression(
  # difference 1 for s1k1
  (t1 - (
    # household 1
    x.h1k1*exp(b.s1k1*x.h1k1 + b.s1k2 * x.h1k2 +  # h1 weight for s1 before delta1
                 log(w.h1 / (exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2) +
                               exp(b.s2k1*x.h1k1 + b.s2k2*x.h1k2)))) + # delta1
      # household 2
      x.h2k1*exp(b.s1k1*x.h2k1 + b.s1k2 * x.h2k2 +  # h2 weight for s1 before delta2
                   log(w.h2 / (exp(b.s1k1*x.h2k1 + b.s1k2*x.h2k2) +
                                 exp(b.s2k1*x.h2k1 + b.s2k2*x.h2k2)))) # delta2
  ))^2 +
    
    # difference 2 for s2k1
    (t2 - (
      # household 1
      x.h1k1*exp(b.s2k1*x.h1k1 + b.s2k2 * x.h1k2 +  # h1 weight for s2 before delta1
                   log(w.h1 / (exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2) +
                                 exp(b.s2k1*x.h1k1 + b.s2k2*x.h1k2)))) + # delta1
        # household 2
        x.h2k1*exp(b.s2k1*x.h2k1 + b.s2k2 * x.h2k2 +  # h2 weight for s2 before delta2
                     log(w.h2 / (exp(b.s1k1*x.h2k1 + b.s1k2*x.h2k2) +
                                   exp(b.s2k1*x.h2k1 + b.s2k2*x.h2k2)))) # delta2
    ))^2)

Deriv(d.2, "b.s1k1")
Deriv(d.2, "b.s2k1")




# 2 households, with 2 states, 2 characteristics = 4 targets
d.4 <- expression(
  # difference 1 for s1k1
  (t1 - (
    # household 1
    x.h1k1*exp(b.s1k1*x.h1k1 + b.s1k2 * x.h1k2 +  # h1 weight for s1 before delta1
                 log(w.h1 / (exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2) +
                               exp(b.s2k1*x.h1k1 + b.s2k2*x.h1k2)))) + # delta1
      # household 2
      x.h2k1*exp(b.s1k1*x.h2k1 + b.s1k2 * x.h2k2 +  # h2 weight for s1 before delta2
                   log(w.h2 / (exp(b.s1k1*x.h2k1 + b.s1k2*x.h2k2) +
                                 exp(b.s2k1*x.h2k1 + b.s2k2*x.h2k2)))) # delta2
  ))^2 +
    
    # difference 2 for s2k1
    (t2 - (
      # household 1
      x.h1k1*exp(b.s2k1*x.h1k1 + b.s2k2 * x.h1k2 +  # h1 weight for s2 before delta1
                   log(w.h1 / (exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2) +
                                 exp(b.s2k1*x.h1k1 + b.s2k2*x.h1k2)))) + # delta1
        # household 2
        x.h2k1*exp(b.s2k1*x.h2k1 + b.s2k2 * x.h2k2 +  # h2 weight for s2 before delta2
                     log(w.h2 / (exp(b.s1k1*x.h2k1 + b.s1k2*x.h2k2) +
                                   exp(b.s2k1*x.h2k1 + b.s2k2*x.h2k2)))) # delta2
    ))^2 +

    # difference 3 for s1k2
    (t3 - (
      # household 1
      x.h1k2*exp(b.s1k1*x.h1k1 + b.s1k2 * x.h1k2 +  # h1 weight for s1 before delta1
                   log(w.h1 / (exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2) +
                                 exp(b.s2k1*x.h1k1 + b.s2k2*x.h1k2)))) + # delta1
        # household 2
        x.h2k2*exp(b.s1k1*x.h2k1 + b.s1k2 * x.h2k2 +  # h2 weight for s1 before delta2
                     log(w.h2 / (exp(b.s1k1*x.h2k1 + b.s1k2*x.h2k2) +
                                   exp(b.s2k1*x.h2k1 + b.s2k2*x.h2k2)))) # delta2
    ))^2  +
    
    # difference 4 for s2k2
    (t4 - (
      # household 1
      x.h1k2*exp(b.s2k1*x.h1k1 + b.s2k2 * x.h1k2 +  # h1 weight for s2 before delta1
                   log(w.h1 / (exp(b.s1k1*x.h1k1 + b.s1k2*x.h1k2) +
                                 exp(b.s2k1*x.h1k1 + b.s2k2*x.h1k2)))) + # delta1
        # household 2
        x.h2k2*exp(b.s2k1*x.h2k1 + b.s2k2 * x.h2k2 +  # h2 weight for s2 before delta2
                     log(w.h2 / (exp(b.s1k1*x.h2k1 + b.s1k2*x.h2k2) +
                                   exp(b.s2k1*x.h2k1 + b.s2k2*x.h2k2)))) # delta2
    ))^2
    )

Deriv(d.4, "b.s1k1")
Deriv(d.4, "b.s2k1")
Deriv(d.4, "b.s1k2")
Deriv(d.4, "b.s2k2")





d <- expression(rowSums(x))
Deriv(d, "x")


f <- function(x){
  (1000 - x^3)^2
}
g <- function(x){
  2 * (1000 - x^3)*(-3*x^2)
}
f(10)
f(11)
g(10)
grad(f, 10)
f(11) - f(10)
f(10.001)


# deriving bprime ----
# delta =log(wh / sum[s] exp(betas*X))
# b=exp(delta(h))
# which is just b = wh / (sum[s] exp(betas*X))

# bprime is deriv of b wrt a beta
# b = wh / (sum-over-s]: exp(beta-for-given-s * x)
b <- expression(wh / (
  exp(b.s1k1 * x1 + b.s1k2 * x2) +
    exp(b.s2k1 * x1 + b.s2k2 * x2) +
    exp(b.s3k1 * x1 + b.s3k2 * x2)
))
Deriv(b, "b.s1k1")
# expression({
#   .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2)
#   -(wh * x1 * .e3/(.e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * 
#                                                                 x1 + b.s3k2 * x2))^2)
# })

#  .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2)
#   -(wh * x1 * .e3 / (.e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2))^2)







