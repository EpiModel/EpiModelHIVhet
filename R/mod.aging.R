
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

  ## Parameters
  time.unit <- dat$param$time.unit
  agecat.cutoff <- dat$param$agecat.cutoff

  ## Attributes
  age <- dat$attr$age
  male <- dat$attr$male
  active <- dat$attr$active
  agecat <- dat$attr$agecat

  ## Updates
  age[active == 1] <- age[active == 1] + time.unit/365

  agecat[active == 1 & male == 0 & age >= agecat.cutoff] <- 1
  agecat[active == 1 & male == 1 & age >= agecat.cutoff] <- 3

  ## Save out
  dat$attr$age <- age
  dat$attr$agecat <- agecat

  return(dat)
}
