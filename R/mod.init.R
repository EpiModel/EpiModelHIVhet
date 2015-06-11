
#' @title Initialization Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param.hiv}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.hiv}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.hiv}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#'
initialize.hiv <- function(x, param, init, control, s) {

  # Restarted simulation
  reinit <- ifelse(class(x) == "netsim", TRUE, FALSE)

  if (reinit == FALSE) {
    dat <- list()
    dat$temp <- list()
    if (class(x$fit) == "network") {
      dat$nw <- simulate(x$formation,
                         basis = x$fit,
                         coef = x$coef.form.crude,
                         constraints = x$constraints,
                         control = control.simulate.formula(MCMC.burnin = 1e6))
    } else {
      dat$nw <- simulate(x$fit,
                         control = control.simulate.ergm(MCMC.burnin = 1e6))
    }

    if (class(dat$nw)[1] == "networkDynamic") {
      dat$nw <- network.collapse(dat$nw, at = 1)
    }

    dat$nw <- activate.vertices(dat$nw, onset = 1, terminus = Inf)
    if (control$delete.nodes == TRUE) {
      dat$nw <- network.extract(dat$nw, at = 1)
    }

    ## Network Model Parameters
    dat$nwparam <- list(x[-which(names(x) == "fit")])

    ## Simulation Parameters
    dat$param <- param
    dat$param$modes <- 1

    dat$init <- init
    dat$control <- control

    ## Nodal Attributes
    dat$attr <- list()

    dat$attr$male <- get.vertex.attribute(dat$nw, "male")

    n <- network.size(dat$nw)
    dat$attr$active <- rep(1, n)
    dat$attr$entTime <- rep(1, n)
    dat$attr$deathTime <- rep(NA, n)
    dat$attr$deathCause <- rep(NA, n)

    ## Initialize HIV related attributes
    dat <- initStatus(dat)
    dat <- initAge(dat)
    dat <- initInfTime(dat)
    dat <- initDx(dat)
    dat <- initTx(dat)
    dat <- circ(dat, at = 1)

    ## Stats List
    dat$stats <- list()

    ## Final steps
    dat$epi <- list()
    dat <- prevalence.hiv(dat, at = 1)
    dat <- simnet.hiv(dat, at = 1)

  } else {
    dat <- list()
    dat$nw <- x$network[[s]]
    dat$param <- param
    dat$param$modes <- 1
    dat$control <- control
    dat$nwparam <- x$nwparam
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$stats <- list()
    dat$stats$nwstats <- x$stats$nwstats[[s]]
    dat$temp <- list()
  }

  dat$param$modes <- 1
  class(dat) <- "dat"

  return(dat)
}


initStatus <- function(dat) {

  ## Variables
  status.rand <- dat$init$status.rand
  i.prev.male <- dat$init$i.prev.male
  i.prev.feml <- dat$init$i.prev.feml

  male <- dat$attr$male
  idsMale <- which(male == 1)
  idsFeml <- which(male == 0)
  nMale <- length(idsMale)
  nFeml <- length(idsFeml)
  n <- nMale + nFeml

  ## Process
  if (status.rand == TRUE) {
    status <- rep(NA, n)
    status[idsMale] <- rbinom(nMale, 1, i.prev.male)
    status[idsFeml] <- rbinom(nFeml, 1, i.prev.feml)
    if (sum(status) == 0) {
      status[ssample(1:length(status), 1)] <- 1
    }
    status <- ifelse(status == 1, "i", "s")
  } else {
    status <- rep("s", n)
    status[sample(idsMale, round(i.prev.male*nMale))] <- "i"
    status[sample(idsFeml, round(i.prev.feml*nFeml))] <- "i"
  }
  dat$attr$status <- status

  return(dat)
}


initAge <- function(dat) {

  age <- get.vertex.attribute(dat$nw, "age")
  dat$attr$age <- age


  return(dat)
}


initInfTime <- function(dat) {

  status <- dat$attr$status
  n <- length(status)

  infecteds <- which(status == "i")
  infTime <- rep(NA, n)

  inf.time.dist <- dat$init$inf.time.dist

  if (inf.time.dist == "allacute") {
    max.inf.time <- dat$param$vl.acute.topeak + dat$param$vl.acute.toset
    infTime[infecteds] <- sample(0:(-max.inf.time), length(infecteds), TRUE)
  } else {
    max.inf.time <- dat$init$max.inf.time / dat$param$time.unit
    if (inf.time.dist == "geometric") {
      total.d.rate <- 1/max.inf.time
      infTime[infecteds] <- -rgeom(length(infecteds), total.d.rate)
    }
    if (inf.time.dist == "uniform") {
      infTime[infecteds] <- sample(0:(-max.inf.time), length(infecteds), TRUE)
    }
  }

  ## Enforce that time infected < age
  infTime[infecteds] <- pmax(infTime[infecteds],
                             1 - dat$attr$age[infecteds] * (365 / dat$param$time.unit))

  dat$attr$infTime <- infTime

  timeInf <- 1 - infTime
  dat$attr$ageInf <- pmax(0, dat$attr$age - round(timeInf) * (dat$param$time.unit / 365))

  stopifnot(all(dat$attr$ageInf[infecteds] <= dat$attr$age[infecteds]),
            all(dat$attr$ageInf[infecteds] >= 0))

  return(dat)
}


initDx <- function(dat) {

  n <- sum(dat$attr$active == 1)
  status <- dat$attr$status

  dxStat <- rep(NA, n)
  dxStat[status == "i"] <- 0

  dxTime <- rep(NA, n)

  dat$attr$dxStat <- dxStat
  dat$attr$dxTime <- dxTime

  return(dat)
}


initTx <- function(dat) {

  ## Variables
  status <- dat$attr$status
  n <- sum(dat$attr$active == 1)
  nInf <- sum(status == "i")

  tx.init.cd4.mean <- dat$param$tx.init.cd4.mean
  tx.init.cd4.sd <- dat$param$tx.init.cd4.sd
  tx.elig.cd4 <- dat$param$tx.elig.cd4


  ## Process
  dat$attr$txStat <- rep(NA, n)
  dat$attr$txStartTime <- rep(NA, n)
  dat$attr$txStops <- rep(NA, n)
  dat$attr$txTimeOn <- rep(NA, n)
  dat$attr$txTimeOff <- rep(NA, n)

  txCD4min <- rep(NA, n)
  txCD4min[status == "i"] <- pmin(rnbinom(nInf,
                                          size = nbsdtosize(tx.init.cd4.mean,
                                                            tx.init.cd4.sd),
                                          mu = tx.init.cd4.mean), tx.elig.cd4)
  dat$attr$txCD4min <- txCD4min
  dat$attr$txCD4start <- rep(NA, n)
  dat$attr$txType <- rep(NA, n)

  return(dat)
}
