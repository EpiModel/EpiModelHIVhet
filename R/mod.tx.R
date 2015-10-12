
#' @title HIV Anti-Retroviral Treatment Module
#'
#' @description Module function for simulating HIV therapy after diagnosis,
#'              including adherence and non-adherence to ART.
#'
#' @inheritParams aging.hiv
#'
#' @export
#'
tx.hiv <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  dxStat <- dat$attr$dxStat
  txStat <- dat$attr$txStat
  txStartTime <- dat$attr$txStartTime
  txStops <- dat$attr$txStops
  txTimeOn <- dat$attr$txTimeOn
  txTimeOff <- dat$attr$txTimeOff
  txCD4start <- dat$attr$txCD4start

  cd4Count <- dat$attr$cd4Count
  tx.elig.cd4 <- dat$param$tx.elig.cd4
  tx.coverage <- dat$param$tx.coverage

  txType <- dat$attr$txType
  tx.adhere.full <- dat$param$tx.adhere.full
  tx.adhere.part <- dat$param$tx.adhere.part


  # Start tx for tx naive ---------------------------------------------------

  ## Calculate tx coverage
  allElig <- which(active == 1 & (cd4Count < tx.elig.cd4 | !is.na(txStartTime)))
  txCov <- sum(!is.na(txStartTime[allElig]))/length(allElig)
  if (is.nan(txCov)) {
    txCov <- 0
  }

  idsElig <- which(active == 1 & dxStat == 1 & txStat == 0 &
                     is.na(txStartTime) & cd4Count < tx.elig.cd4)
  nElig <- length(idsElig)
  idsTx <- NULL


  ## Treatment coverage
  nStart <- max(0, min(nElig, round((tx.coverage - txCov) * length(allElig))))
  if (nStart > 0) {
    idsTx <- ssample(idsElig, nStart)
  }


  ## Treatment type assignment
  if (length(idsTx) > 0) {
    needtxType <- which(is.na(txType[idsTx]))
    if (length(needtxType) > 0) {
      txType[idsTx[needtxType]] <- rbinom(length(needtxType), 1, tx.adhere.full)
    }
    if (tx.adhere.part == 0) {
      idsTx <- intersect(idsTx, which(txType == 1))
    }
  }

  if (length(idsTx) > 0) {
    txStat[idsTx] <- 1
    txStartTime[idsTx] <- at
    txStops[idsTx] <- 0
    txTimeOn[idsTx] <- 0
    txTimeOff[idsTx] <- 0
    txCD4start[idsTx] <- cd4Count[idsTx]
  }


  # Stop tx -----------------------------------------------------------------
  idsStop <- NULL
  idsEligStop <- which(active == 1 & dat$attr$txStat == 1 & txType == 0)
  nEligStop <- length(idsEligStop)
  if (nEligStop > 0) {
    vecStop <- which(rbinom(nEligStop, 1, (1 - tx.adhere.part)) == 1)
    if (length(vecStop) > 0) {
      idsStop <- idsEligStop[vecStop]
      txStat[idsStop] <- 0
      txStops[idsStop] <- txStops[idsStop] + 1
    }
  }


  # Restart tx --------------------------------------------------------------
  idsRest <- NULL
  idsEligRest <- which(active == 1 & dat$attr$txStat == 0 & txStops > 0)
  nEligRest <- length(idsEligRest)
  if (nEligRest > 0) {
    vecRes <- which(rbinom(nEligRest, 1, tx.adhere.part) == 1)
    if (length(vecRes) > 0) {
      idsRest <- idsEligRest[vecRes]
      txStat[idsRest] <- 1
      dat$attr$vlSlope[idsRest] <- NA
    }
  }


  # Output ------------------------------------------------------------------
  idsOnTx <- which(txStat == 1)
  idsOffTx <- which(txStat == 0 & !is.na(txStartTime))
  txTimeOn[idsOnTx] <- txTimeOn[idsOnTx] + 1
  txTimeOff[idsOffTx] <- txTimeOff[idsOffTx] + 1

  dat$attr$txStat <- txStat
  dat$attr$txStartTime <- txStartTime
  dat$attr$txStops <- txStops
  dat$attr$txTimeOn <- txTimeOn
  dat$attr$txTimeOff <- txTimeOff
  dat$attr$txType <- txType
  dat$attr$txCD4start <- txCD4start

  return(dat)
}
