
## Test script

library(EpiModelHIV)

st <- make_nw.hiv()
est <- netest(st$nw,
              formation = st$formation,
              target.stats = st$stats,
              coef.form = -Inf,
              coef.diss = st$coef.diss,
              constraints = ~bd(maxout = 3),
              set.control.ergm = control.ergm(MCMLE.maxit = 500, MPLE.type = "penalized"))

dx <- netdx(est, nsims = 5, nsteps = 250,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
print(dx)
plot(dx)

param <- param.hiv()
init <- init.hiv()
control <- control.hiv()

sim <- netsim(est, param, init, control)
