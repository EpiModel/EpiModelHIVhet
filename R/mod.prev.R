
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
  dxTime <- dat$attr$dxTime
  vlLevel <- dat$attr$vlLevel
  txTimeOn <- dat$attr$txTimeOn
  txTimeOff <- dat$attr$txTimeOff
  txType <- dat$attr$txType

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
    dat$epi$meanAge <- dat$epi$meanAge.sus <- dat$epi$meanAge.inf <- rNA
    dat$epi$meanAge.male <- dat$epi$meanAge.feml <- rNA
    dat$epi$propMale <- rNA

    dat$epi$newDx <- rNA
    dat$epi$meanTxTimeOn <- dat$epi$meanTxTimeOn.t0 <-
      dat$epi$meanTxTimeOn.t1 <- rNA
    dat$epi$meanVl <- dat$epi$vlSupp <- dat$epi$vlSupp.tx <-
      dat$epi$vlSupp.t0 <- dat$epi$vlSupp.t1 <- rNA
    if (dat$control$clin.array == TRUE) {
      dat$clin <- list()
      dat$clin$txStat <- matrix(NA, length(active), nsteps)
      dat$clin$vlLevel <- matrix(NA, length(active), nsteps)
      dat$clin$cd4Count <- matrix(NA, length(active), nsteps)
    }

    # Incidence vectors
    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA

    dat$epi$txCov <- rNA
    dat$epi$txStart <- rNA
    dat$epi$txStop <- rNA
    dat$epi$txRest <- rNA

  }

  ### Prevalence ###

  ## Overall prevalence/incidence
  dat$epi$s.num[at] <- sum(active == 1 & status == "s")
  dat$epi$i.num[at] <- sum(active == 1 & status == "i")
  dat$epi$num[at] <- sum(active == 1)
  dat$epi$cumlNum[at] <- dat$epi$num[1] + sum(dat$epi$b.flow, na.rm = TRUE)
  dat$epi$cumlInc[at] <- sum(dat$epi$si.flow, na.rm = TRUE)
  dat$epi$incr[at] <- (dat$epi$si.flow[at] / dat$epi$s.num[at])*5200

  ## Sex-specific prevalence
  dat$epi$s.num.male[at] <- sum(active == 1 & status == "s" & male == 1)
  dat$epi$s.num.feml[at] <- sum(active == 1 & status == "s" & male == 0)
  dat$epi$i.num.male[at] <- sum(active == 1 & status == "i" & male == 1)
  dat$epi$i.num.feml[at] <- sum(active == 1 & status == "i" & male == 0)
  dat$epi$i.prev.male[at] <- sum(active == 1 & status == "i" & male == 1) /
                             sum(active == 1 & male == 1)
  dat$epi$i.prev.feml[at] <- sum(active == 1 & status == "i" & male == 0) /
                             sum(active == 1 & male == 0)
  dat$epi$incr.male[at] <- (dat$epi$si.flow.male[at] / dat$epi$s.num.male[at])*5200
  dat$epi$incr.feml[at] <- (dat$epi$si.flow.feml[at] / dat$epi$s.num.feml[at])*5200

  ### Demographics ###
  dat$epi$num.male[at] <- sum(active == 1 & male == 1)
  dat$epi$num.feml[at] <- sum(active == 1 & male == 0)

  ## Age
  dat$epi$meanAge[at] <- mean(age[active == 1])
  dat$epi$meanAge.sus[at] <- mean(age[active == 1 & status == "s"])
  dat$epi$meanAge.inf[at] <- mean(age[active == 1 & status == "i"])
  dat$epi$meanAge.male[at] <- mean(age[active == 1 & male == 1])
  dat$epi$meanAge.feml[at] <- mean(age[active == 1 & male == 0])
  dat$epi$propMale[at] <- mean(male[active == 1])


  ### Diagnosis and treatment ###
  dat$epi$newDx[at] <- sum(dxTime == at)

  ## Treatment adherence
  dat$epi$meanTxTimeOn[at] <- mean(txTimeOn/(txTimeOn + txTimeOff), na.rm = TRUE)
  if (is.nan(dat$epi$meanTxTimeOn[at])) {
    dat$epi$meanTxTimeOn[at] <- NA
  }
  dat$epi$meanTxTimeOn.t0[at] <- mean(txTimeOn[txType == 0]/
                                        (txTimeOn[txType == 0] +
                                           txTimeOff[txType == 0]),
                                      na.rm = TRUE)
  if (is.nan(dat$epi$meanTxTimeOn.t0[at])) {
    dat$epi$meanTxTimeOn.t0[at] <- NA
  }
  dat$epi$meanTxTimeOn.t1[at] <- mean(txTimeOn[txType == 1]/
                                        (txTimeOn[txType == 1] +
                                           txTimeOff[txType == 1]),
                                      na.rm = TRUE)
  if (is.nan(dat$epi$meanTxTimeOn.t1[at])) {
    dat$epi$meanTxTimeOn.t1[at] <- NA
  }

  ## Viral load and suppression
  if (at >= 2) {
    dat$epi$meanVl[at] <- mean(vlLevel[active == 1 & status == "i"])
    dat$epi$vlSupp[at] <- length(whichVlSupp(dat$attr, dat$param)) /
      sum(active == 1 & status == "i", na.rm = TRUE)
    if (is.nan(dat$epi$vlSupp[at])) {
      dat$epi$vlSupp[at] <- NA
    }
    dat$epi$vlSupp.tx[at] <- length(intersect(which(txType %in% 0:1),
                                              whichVlSupp(dat$attr, dat$param))) /
      sum(active == 1 & status == "i" &
            txType %in% 0:1, na.rm = TRUE)
    if (is.nan(dat$epi$vlSupp.tx[at])) {
      dat$epi$vlSupp.tx[at] <- NA
    }
    dat$epi$vlSupp.t0[at] <- length(intersect(which(txType == 0),
                                              whichVlSupp(dat$attr, dat$param))) /
      sum(active == 1 & status == "i" &
            txType == 0, na.rm = TRUE)
    if (is.nan(dat$epi$vlSupp.t0[at])) {
      dat$epi$vlSupp.t0[at] <- NA
    }
    dat$epi$vlSupp.t1[at] <- length(intersect(which(txType == 1),
                                              whichVlSupp(dat$attr, dat$param))) /
      sum(active == 1 & status == "i" &
            txType == 1, na.rm = TRUE)
    if (is.nan(dat$epi$vlSupp.t1[at])) {
      dat$epi$vlSupp.t1[at] <- NA
    }

    ### Clinical array ###
    if (dat$control$clin.array == TRUE) {
      dat$clin$txStat[,at] <- dat$attr$txStat
      dat$clin$cd4Count[,at] <- dat$attr$cd4Count
      dat$clin$vlLevel[,at] <- dat$attr$vlLevel
    }
  }

  ### Supplemental prevalence submod
  if (!is.null(dat$control$getprev.suppl)) {
    dat <- do.call(dat$control[["getprev.suppl"]], list(dat, at))
  }

  return(dat)
}


whichVlSupp <- function(attr, param) {

  which(attr$active == 1 &
          attr$status == "i" &
          attr$vlLevel <= log10(50) &
          (attr$age - attr$ageInf) * (365 / param$time.unit) >
          (param$vl.acute.topeak + param$vl.acute.toset))

}
