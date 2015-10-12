
#' @title Viral Load Module
#'
#' @description Module function for simulating progression of HIV viral load in
#'              natural disease dynamics and in the presence of ART.
#'
#' @inheritParams aging.hiv
#'
#' @export
#'
vl.hiv <- function(dat, at) {

  ## Common variables
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime


  # Assign base VL ----------------------------------------------------------
  if (is.null(dat$attr$vlLevel)) {
    dat$attr$vlLevel <- rep(NA, sum(active == 1))
    dat$attr$vlSlope <- rep(NA, sum(active == 1))
  }
  vlLevel <- dat$attr$vlLevel

  idsEligAsn <- which(active == 1 & status == 1 & is.na(vlLevel))
  if (length(idsEligAsn) > 0) {
    vlLevel[idsEligAsn] <- expected_vl(male = dat$attr$male[idsEligAsn],
                                       age = dat$attr$age[idsEligAsn],
                                       ageInf = dat$attr$ageInf[idsEligAsn],
                                       param = dat$param)
  }


  # Update natural VL -------------------------------------------------------
  txStartTime <- dat$attr$txStartTime
  idsEligUpd <- which(active == 1 & status == 1 &
                        infTime < at & is.na(txStartTime))

  if (length(idsEligUpd) > 0) {
    vlLevel[idsEligUpd] <- expected_vl(male = dat$attr$male[idsEligUpd],
                                       age = dat$attr$age[idsEligUpd],
                                       ageInf = dat$attr$ageInf[idsEligUpd],
                                       param = dat$param)
  }

  # VL decline with ART -----------------------------------------------------
  txStat <- dat$attr$txStat
  idsEligTx <- which(active == 1 & status == 1 & infTime < at & txStat == 1)
  if (length(idsEligTx) > 0) {
    tx.vlsupp.time <- dat$param$tx.vlsupp.time
    tx.vlsupp.level <- dat$param$tx.vlsupp.level

    vlSlope <- dat$attr$vlSlope
    needSlope <- intersect(idsEligTx, which(is.na(vlSlope)))

    vl.slope <- vlSlope
    if (length(needSlope) > 0) {
      vl.diff <- pmin(tx.vlsupp.level - vlLevel[needSlope], 0)
      vl.slope[needSlope] <- vl.diff / tx.vlsupp.time
      dat$attr$vlSlope[needSlope] <- vl.slope[needSlope]
    }

    vlLevel[idsEligTx] <- pmax(vlLevel[idsEligTx] + vl.slope[idsEligTx], tx.vlsupp.level)
  }


  # VL rebound post ART -----------------------------------------------------
  idsEligNoTx <- which(active == 1 & status == 1 &
                         txStat == 0 & !is.na(txStartTime))
  if (length(idsEligNoTx) > 0) {
    tx.vlsupp.time <- dat$param$tx.vlsupp.time

    expVl <- expected_vl(male = dat$attr$male[idsEligNoTx],
                         age = dat$attr$age[idsEligNoTx],
                         ageInf = dat$attr$ageInf[idsEligNoTx],
                         param = dat$param)

    vl.slope <- dat$attr$vlSlope

    vlLevel[idsEligNoTx] <- pmin(vlLevel[idsEligNoTx] - vl.slope[idsEligNoTx], expVl)
  }

  dat$attr$vlLevel <- vlLevel

  return(dat)
}


expected_vl <- function(male, age, ageInf, param) {

  timeInf <- (age - ageInf) * (365 / param$time.unit)

  slope1 <- param$vl.acute.peak / param$vl.acute.topeak
  slope2 <- (param$vl.setpoint - param$vl.acute.peak) /
            (param$vl.acute.toset - param$vl.acute.topeak)

  sl3denom <- expected_cd4(method = "timeto",
                           cd4Count1 = 200, cd4Count2 = 25,
                           male = male, age = age, ageInf = ageInf,
                           time.unit = param$time.unit)
  slope3 <- (param$vl.aidsmax - param$vl.setpoint) / sl3denom

  setptTime <- param$vl.acute.topeak + param$vl.acute.toset
  aidsTime <- expected_cd4(method = "timeto", cd4Count1 = 200,
                           male = male, age = age, ageInf = ageInf,
                           time.unit = param$time.unit)

  gp <- 1 * (timeInf <= param$vl.acute.topeak) +
        2 * (timeInf > param$vl.acute.topeak & timeInf <= setptTime) +
        3 * (timeInf > setptTime & timeInf <= aidsTime) +
        4 * (timeInf > aidsTime)

  vlLevel <- rep(NA, length(timeInf))
  vlLevel[gp == 1] <- timeInf[gp == 1] * slope1
  vlLevel[gp == 2] <- pmax(param$vl.setpoint,
                           param$vl.acute.peak +
                             (timeInf[gp == 2] - param$vl.acute.topeak) * slope2)
  vlLevel[gp == 3] <- param$vl.setpoint
  vlLevel[gp == 4] <- pmin(param$vl.aidsmax,
                           param$vl.setpoint +
                             (timeInf[gp == 4] - aidsTime[gp == 4]) * slope3[gp == 4])

  return(vlLevel)
}
