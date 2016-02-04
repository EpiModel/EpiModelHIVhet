
## Test script

library(EpiModelHIVhet)

st <- make_nw.hiv(part.dur = 2013)
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
init <- init.hiv(i.prev.male = 0.25, i.prev.feml = 0.25)
control <- control.hiv(nsteps = 2600)

sim <- netsim(est, param, init, control)

