
#' @title Prevalence Module
#'
#' @description Module function to calculate and store summary statistics for
#'              disease prevalence, demographics, and other epidemiological
#'              outcomes.
#'
#' @inheritParams aging.hiv
#'
#' @export
#'
prevalence.hiv <- function(dat, at) {

  status <- dat$attr$status
  active <- dat$attr$active
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {

    # Prev vectors
    dat$epi$s.num <-  dat$epi$i.num <- rNA
    dat$epi$num <- dat$epi$cumlNum <- rNA
    dat$epi$cumlInc <- dat$epi$incr <- rNA

    dat$epi$s.num.male <- dat$epi$s.num.feml <- rNA
    dat$epi$i.num.male <- dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- dat$epi$i.prev.feml <- rNA
    dat$epi$incr.male <- dat$epi$incr.feml <- rNA

    dat$epi$num.male <- dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA

  }

  ### Prevalence ###

  ## Overall prevalence/incidence
  dat$epi$s.num[at] <- sum(active == 1 & status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(active == 1 & status == 1, na.rm = TRUE)
  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$cumlNum[at] <- dat$epi$num[1] + sum(dat$epi$b.flow, na.rm = TRUE)
  dat$epi$cumlInc[at] <- sum(dat$epi$si.flow, na.rm = TRUE)
  dat$epi$incr[at] <- (dat$epi$si.flow[at] / dat$epi$s.num[at])*5200

  ## Sex-specific prevalence
  dat$epi$s.num.male[at] <- sum(active == 1 & status == 0 & male == 1, na.rm = TRUE)
  dat$epi$s.num.feml[at] <- sum(active == 1 & status == 0 & male == 0, na.rm = TRUE)
  dat$epi$i.num.male[at] <- sum(active == 1 & status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(active == 1 & status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(active == 1 & status == 1 & male == 1, na.rm = TRUE) /
                             sum(active == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(active == 1 & status == 1 & male == 0, na.rm = TRUE) /
                             sum(active == 1 & male == 0, na.rm = TRUE)
  dat$epi$incr.male[at] <- (dat$epi$si.flow.male[at] / dat$epi$s.num.male[at])*5200
  dat$epi$incr.feml[at] <- (dat$epi$si.flow.feml[at] / dat$epi$s.num.feml[at])*5200

  ### Demographics ###
  dat$epi$num.male[at] <- sum(active == 1 & male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(active == 1 & male == 0, na.rm = TRUE)

  ## Age
  dat$epi$meanAge[at] <- mean(age[active == 1], na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male[active == 1], na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {

  which(attr$active == 1 &
          attr$status == 1 &
          attr$vlLevel <= log10(50) &
          (attr$age - attr$ageInf) * (365 / param$time.unit) >
          (param$vl.acute.topeak + param$vl.acute.toset))

}
