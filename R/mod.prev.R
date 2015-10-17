
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
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(active == 1 & status == 1, na.rm = TRUE)
  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)

  dat$epi$i.num.male[at] <- sum(active == 1 & status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(active == 1 & status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(active == 1 & status == 1 & male == 1, na.rm = TRUE) /
                             sum(active == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(active == 1 & status == 1 & male == 0, na.rm = TRUE) /
                             sum(active == 1 & male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(active == 1 & male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(active == 1 & male == 0, na.rm = TRUE)
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
