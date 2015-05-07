
#' @title Aging Module
#'
#' @description This module ages all active nodes in the population by one time
#'              unit at each time step.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @export
#'
aging.hiv <- function(dat, at) {

  time.unit <- dat$param$time.unit

  age <- dat$attr$age
  active <- dat$attr$active

  age[active == 1] <- age[active == 1] + time.unit/365

  dat$attr$age <- age

  dat$nw <- network::set.vertex.attribute(dat$nw, "age", age)

  return(dat)
}
