
#' @title Deaths Module
#'
#' @description Module for simulating deaths among susceptible and infected
#'              persons within the population.
#'
#' @inheritParams aging.hiv
#'
#' @export
#'
deaths.hiv <- function(dat, at) {

  # Susceptible Deaths ------------------------------------------------------

  ## Variables
  active <- dat$attr$active
  male <- dat$attr$male
  age <- dat$attr$age
  cd4Count <- dat$attr$cd4Count

  di.cd4.aids <- dat$param$di.cd4.aids
  ds.exit.age <- dat$param$ds.exit.age

  ## Eligible are: active uninf, pre-death infected, unhealthy old
  idsEligSus <- which(active == 1 & (is.na(cd4Count) |
                                       cd4Count > di.cd4.aids |
                                       (cd4Count <= di.cd4.aids & age > ds.exit.age)))
  nEligSus <- length(idsEligSus)

  # Set age-sex specific rates
  ds.rates <- dat$param$ds.rates
  if (nEligSus > 0) {
    rates <- ds.rates$mrate[100*male[idsEligSus] + age[idsEligSus]]
  }


  ## Process
  nDeathsSus <- 0; idsDeathsSus <- NULL
  if (nEligSus > 0) {
    vecDeathsSus <- which(rbinom(nEligSus, 1, rates) == 1)
    nDeathsSus <- length(vecDeathsSus)
  }


  ## Update Attributes
  if (nDeathsSus > 0) {
    idsDeathsSus <- idsEligSus[vecDeathsSus]
    dat$attr$active[idsDeathsSus] <- 0
    dat$attr$deathTime[idsDeathsSus] <- at
    dat$attr$deathCause[idsDeathsSus] <- "s"
  }


  # Infected Deaths ---------------------------------------------------------

  ## Variables
  active <- dat$attr$active
  di.cd4.rate <- dat$param$di.cd4.rate

  ## Process
  nDeathsInf <- 0; idsDeathsInf <- NULL

  cd4Count <- dat$attr$cd4Count
  stopifnot(length(active) == length(cd4Count))

  idsEligInf <- which(active == 1 & cd4Count <= di.cd4.aids)
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecDeathsInf <- which(rbinom(nEligInf, 1, di.cd4.rate) == 1)
    if (length(vecDeathsInf) > 0) {
      idsDeathsInf <- idsEligInf[vecDeathsInf]
      nDeathsInf <- length(idsDeathsInf)
    }
  }

  idsDeathsDet <- which(active == 1 & cd4Count <= 0)
  if (length(idsDeathsDet) > 0) {
    idsDeathsInf <- c(idsDeathsInf, idsDeathsDet)
    nDeathsInf <- nDeathsInf + length(idsDeathsDet)
  }


  ## Update Attributes
  if (nDeathsInf > 0) {
    dat$attr$active[idsDeathsInf] <- 0
    dat$attr$deathTime[idsDeathsInf] <- at
    dat$attr$deathCause[idsDeathsInf] <- "i"
  }


  # Update Network ----------------------------------------------------------
  idsDeaths <- c(idsDeathsSus, idsDeathsInf)
  if (length(idsDeaths) > 0) {
    dat$nw <- networkDynamic::deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                                  v = idsDeaths, deactivate.edges = TRUE)
  }


  # Output ------------------------------------------------------------------
  dat$epi$ds.flow[at] <- nDeathsSus
  dat$epi$di.flow[at] <- nDeathsInf

  return(dat)
}
