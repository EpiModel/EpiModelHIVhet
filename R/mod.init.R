
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

  dat <- list()
  dat$temp <- list()
  nw <- simulate(x$fit, control = control.simulate.ergm(MCMC.burnin = 1e6))

  dat$el <- as.edgelist(nw)
  attributes(dat$el)$vnames <- NULL
  p <- tergmLite::ergm_prep(nw, x$formation, x$coef.diss$dissolution, x$coef.form,
                 x$coef.diss$coef.adj, x$constraints, control = tergm::control.simulate.network())
  p$model.form$formula <- NULL
  p$model.diss$formula <- NULL
  dat$p <- p

  ## Network Model Parameters
  dat$nwparam <- list(x[-which(names(x) == "fit")])

  ## Simulation Parameters
  dat$param <- param
  dat$param$modes <- 1

  dat$init <- init
  dat$control <- control

  ## Nodal Attributes
  dat$attr <- list()

  dat$attr$male <- get.vertex.attribute(nw, "male")

  n <- network.size(nw)
  dat$attr$active <- rep(1, n)
  dat$attr$entTime <- rep(1, n)

  dat <- initStatus(dat)

  age <- rep(NA, n)
  age[dat$attr$male == 0] <- sample(init$ages.feml, sum(dat$attr$male == 0), TRUE)
  age[dat$attr$male == 1] <- sample(init$ages.male, sum(dat$attr$male == 1), TRUE)
  dat$attr$age <- age

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

}


#' @title Reinitialization Module
#'
#' @description This function reinitializes the master \code{dat} object on which
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
reinit.hiv <- function(x, param, init, control, s) {
  dat <- list()
  dat$el <- x$el[[s]]
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

  dat$param$modes <- 1
  class(dat) <- "dat"

  return(dat)
}


initStatus <- function(dat) {

  ## Variables
  i.prev.male <- dat$init$i.prev.male
  i.prev.feml <- dat$init$i.prev.feml

  male <- dat$attr$male
  idsMale <- which(male == 1)
  idsFeml <- which(male == 0)
  nMale <- length(idsMale)
  nFeml <- length(idsFeml)
  n <- nMale + nFeml

  ## Process
  status <- rep(0, n)
  status[sample(idsMale, round(i.prev.male * nMale))] <- 1
  status[sample(idsFeml, round(i.prev.feml * nFeml))] <- 1

  dat$attr$status <- status

  return(dat)
}


initInfTime <- function(dat) {

  status <- dat$attr$status
  n <- length(status)

  infecteds <- which(status == 1)
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
  dxStat[status == 1] <- 0

  dxTime <- rep(NA, n)

  dat$attr$dxStat <- dxStat
  dat$attr$dxTime <- dxTime

  return(dat)
}


initTx <- function(dat) {

  ## Variables
  status <- dat$attr$status
  n <- sum(dat$attr$active == 1)
  nInf <- sum(status == 1)

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
  txCD4min[status == 1] <- pmin(rnbinom(nInf,
                                        size = nbsdtosize(tx.init.cd4.mean,
                                                          tx.init.cd4.sd),
                                        mu = tx.init.cd4.mean), tx.elig.cd4)
  dat$attr$txCD4min <- txCD4min
  dat$attr$txCD4start <- rep(NA, n)
  dat$attr$txType <- rep(NA, n)

  return(dat)
}
