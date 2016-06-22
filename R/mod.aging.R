
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
aging_het <- function(dat, at) {

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Attributes
  age <- dat$attr$age
  active <- dat$attr$active

  ## Updates
  age[active == 1] <- age[active == 1] + time.unit/365

  ## Save out
  dat$attr$age <- age

  return(dat)
}
