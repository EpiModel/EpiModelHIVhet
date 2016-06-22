
#' @title Circumcision Module
#'
#' @description Module for simulating circumcision status at birth/entry into
#'              the population.
#'
#' @inheritParams aging_het
#'
#' @export
#'
circ <- function(dat, at) {

  circ.prob.birth <- dat$param$circ.prob.birth
  male <- dat$attr$male

  ### Initialize circ status at T1
  if (at == 1) {

    n <- length(male)
    nMales <- sum(male == 1)
    age <- dat$attr$age

    circStat <- circTime <- rep(NA, n)

    circStat[male == 1] <- rbinom(nMales, 1, circ.prob.birth)

    isCirc <- which(circStat == 1)
    circTime[isCirc] <- round(-age[isCirc] * (365 / dat$param$time.unit))

    dat$attr$circStat <- circStat
    dat$attr$circTime <- circTime
  }


  ### Set circ status for incoming nodes
  if (at > 1) {
    circStat <- dat$attr$circStat
    entTime <- dat$attr$entTime

    idsNewMale <- which(male == 1 & entTime == at)

    if (length(idsNewMale) > 0) {
      age <- dat$attr$age[idsNewMale]
      newCirc <- rbinom(length(idsNewMale), 1, circ.prob.birth)
      isCirc <- which(newCirc == 1)

      newCircTime <- rep(NA, length(idsNewMale))
      newCircTime[isCirc] <- round(-age[isCirc] * (365 / dat$param$time.unit))

      dat$attr$circStat[idsNewMale] <- newCirc
      dat$attr$circTime[idsNewMale] <- newCircTime
    }
  }

  return(dat)
}
