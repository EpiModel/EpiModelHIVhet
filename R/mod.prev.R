
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

  if (at == 1) {

    ### Prevalence ###

    ## Overall prevalence/incidence
    dat$epi$s.num <- sum(active == 1 & status == "s")
    dat$epi$i.num <- sum(active == 1 & status == "i")
    dat$epi$num <- sum(active == 1)
    dat$epi$cumNum <- dat$epi$num[1] + sum(dat$epi$b.flow)
    dat$epi$cumlInc <- 0
    dat$epi$incr <- 0

    ## Sex-specific prevalence
    dat$epi$s.num.male <- sum(active == 1 & status == "s" & male == 1)
    dat$epi$s.num.feml <- sum(active == 1 & status == "s" & male == 0)
    dat$epi$i.num.male <- sum(active == 1 & status == "i" & male == 1)
    dat$epi$i.num.feml <- sum(active == 1 & status == "i" & male == 0)
    dat$epi$i.prev.male <- sum(active == 1 & status == "i" & male == 1) /
                           sum(active == 1 & male == 1)
    dat$epi$i.prev.feml <- sum(active == 1 & status == "i" & male == 0) /
                           sum(active == 1 & male == 0)
    dat$epi$incr.male <- 0
    dat$epi$incr.feml <- 0


    ## Age-sex specific prevalence
    if (dat$control$calc.asprev == TRUE) {
      dat$epi$i.prev.1824 <- sum(active == 1 & status == "i" & age >= 18 & age < 25)/
        sum(active == 1 & age >= 18 & age < 25)
      dat$epi$i.prev.2529 <- sum(active == 1 & status == "i" & age >= 25 & age < 30)/
        sum(active == 1 & age >= 25 & age < 30)
      dat$epi$i.prev.3034 <- sum(active == 1 & status == "i" & age >= 30 & age < 35)/
        sum(active == 1 & age >= 30 & age < 35)
      dat$epi$i.prev.3539 <- sum(active == 1 & status == "i" & age >= 35 & age < 40)/
        sum(active == 1 & age >= 35 & age < 40)
      dat$epi$i.prev.4044 <- sum(active == 1 & status == "i" & age >= 40 & age < 45)/
        sum(active == 1 & age >= 40 & age < 45)
      dat$epi$i.prev.45pl <- sum(active == 1 & status == "i" & age >= 45)/
        sum(active == 1 & age >= 45)
      dat$epi$i.prev.1829 <- sum(active == 1 & status == "i" & age < 30)/
        sum(active == 1 & age < 30)
      dat$epi$i.prev.30pl <- sum(active == 1 & status == "i" & age >= 30)/
        sum(active == 1 & age >= 30)

      dat$epi$i.prev.feml.1824 <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 18 & age < 25)/
        sum(active == 1 & age >= 18 & age < 25 & male == 0)
      dat$epi$i.prev.feml.2529 <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 25 & age < 30)/
        sum(active == 1 & age >= 25 & age < 30 & male == 0)
      dat$epi$i.prev.feml.3034 <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 30 & age < 35)/
        sum(active == 1 & age >= 30 & age < 35 & male == 0)
      dat$epi$i.prev.feml.3539 <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 35 & age < 40)/
        sum(active == 1 & age >= 35 & age < 40 & male == 0)
      dat$epi$i.prev.feml.4044 <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 40 & age < 45)/
        sum(active == 1 & age >= 40 & age < 45 & male == 0)
      dat$epi$i.prev.feml.45pl <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 45)/
        sum(active == 1 & age >= 45 & male == 0)
      dat$epi$i.prev.feml.1829 <- sum(active == 1 & status == "i" & male == 0 & age < 30)/
        sum(active == 1 & male == 0 & age < 30)
      dat$epi$i.prev.feml.30pl <- sum(active == 1 & status == "i" & male == 0 & age >= 30)/
        sum(active == 1 & male == 0 & age >= 30)

      dat$epi$i.prev.male.1824 <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 18 & age < 25)/
        sum(active == 1 & age >= 18 & age < 25 & male == 0)
      dat$epi$i.prev.male.2529 <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 25 & age < 30)/
        sum(active == 1 & age >= 25 & age < 30 & male == 0)
      dat$epi$i.prev.male.3034 <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 30 & age < 35)/
        sum(active == 1 & age >= 30 & age < 35 & male == 0)
      dat$epi$i.prev.male.3539 <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 35 & age < 40)/
        sum(active == 1 & age >= 35 & age < 40 & male == 0)
      dat$epi$i.prev.male.4044 <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 40 & age < 45)/
        sum(active == 1 & age >= 40 & age < 45 & male == 0)
      dat$epi$i.prev.male.45pl <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 45)/
        sum(active == 1 & age >= 45 & male == 0)
      dat$epi$i.prev.male.1829 <- sum(active == 1 & status == "i" & male == 0 & age < 30)/
        sum(active == 1 & male == 1 & age < 30)
      dat$epi$i.prev.male.30pl <- sum(active == 1 & status == "i" & male == 0 & age >= 30)/
        sum(active == 1 & male == 1 & age >= 30)
    }


    ### Demographics ###
    dat$epi$num.male <- sum(active == 1 & male == 1)
    dat$epi$num.feml <- sum(active == 1 & male == 0)
    dat$epi$meanAge <- mean(age[active == 1])
    dat$epi$meanAge.sus <- mean(age[active == 1 & status == "s"])
    dat$epi$meanAge.inf <- mean(age[active == 1 & status == "i"])
    dat$epi$meanAge.male <- mean(age[active == 1 & male == 1])
    dat$epi$meanAge.feml <- mean(age[active == 1 & male == 0])
    dat$epi$propMale <- mean(male[active == 1])


    ### Diagnosis and Treatment ###
    dat$epi$newDx <- sum(dxTime == at, na.rm = TRUE)

    ## Tx adherence
    dat$epi$meanTxTimeOn <- NA
    dat$epi$meanTxTimeOn.t0 <- NA
    dat$epi$meanTxTimeOn.t1 <- NA

    ## Viral load
    dat$epi$meanVl <- NA
    dat$epi$vlSupp <- NA
    dat$epi$vlSupp.tx <- NA
    dat$epi$vlSupp.t0 <- NA
    dat$epi$vlSupp.t1 <- NA


    ### Age Mixing ###
    relage <- get_male_relage(dat$nw, at)
    dat$epi$male.relage.mean <- mean(relage)
    dat$epi$male.relage.mode <- dens_mode(relage)


    ### Clinical array for treatment ###
    if (dat$control$clin.array == TRUE) {
      dat$clin <- list()
      dat$clin$txStat <- matrix(NA, length(active), dat$control$nsteps)
      dat$clin$vlLevel <- matrix(NA, length(active), dat$control$nsteps)
      dat$clin$cd4Count <- matrix(NA, length(active), dat$control$nsteps)
    }

  } else {

    ### Prevalence ###

    ## Overall prevalence/incidence
    dat$epi$s.num[at] <- sum(active == 1 & status == "s")
    dat$epi$i.num[at] <- sum(active == 1 & status == "i")
    dat$epi$num[at] <- sum(active == 1)
    dat$epi$cumNum[at] <- dat$epi$num[1] + sum(dat$epi$b.flow)
    dat$epi$cumlInc[at] <- sum(dat$epi$si.flow)
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


    ## Age-sex specific prevalence
    if (dat$control$calc.asprev == TRUE) {
      dat$epi$i.prev.1824[at] <- sum(active == 1 & status == "i" & age >= 18 & age < 25)/
        sum(active == 1 & age >= 18 & age < 25)
      dat$epi$i.prev.2529[at] <- sum(active == 1 & status == "i" & age >= 25 & age < 30)/
        sum(active == 1 & age >= 25 & age < 30)
      dat$epi$i.prev.3034[at] <- sum(active == 1 & status == "i" & age >= 30 & age < 35)/
        sum(active == 1 & age >= 30 & age < 35)
      dat$epi$i.prev.3539[at] <- sum(active == 1 & status == "i" & age >= 35 & age < 40)/
        sum(active == 1 & age >= 35 & age < 40)
      dat$epi$i.prev.4044[at] <- sum(active == 1 & status == "i" & age >= 40 & age < 45)/
        sum(active == 1 & age >= 40 & age < 45)
      dat$epi$i.prev.45pl[at] <- sum(active == 1 & status == "i" & age >= 45)/
        sum(active == 1 & age >= 45)
      dat$epi$i.prev.1829[at] <- sum(active == 1 & status == "i" & age < 30)/
        sum(active == 1 & age < 30)
      dat$epi$i.prev.30pl[at] <- sum(active == 1 & status == "i" & age >= 30)/
        sum(active == 1 & age >= 30)

      dat$epi$i.prev.feml.1824[at] <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 18 & age < 25)/
        sum(active == 1 & age >= 18 & age < 25 & male == 0)
      dat$epi$i.prev.feml.2529[at] <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 25 & age < 30)/
        sum(active == 1 & age >= 25 & age < 30 & male == 0)
      dat$epi$i.prev.feml.3034[at] <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 30 & age < 35)/
        sum(active == 1 & age >= 30 & age < 35 & male == 0)
      dat$epi$i.prev.feml.3539[at] <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 35 & age < 40)/
        sum(active == 1 & age >= 35 & age < 40 & male == 0)
      dat$epi$i.prev.feml.4044[at] <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 40 & age < 45)/
        sum(active == 1 & age >= 40 & age < 45 & male == 0)
      dat$epi$i.prev.feml.45pl[at] <- sum(active == 1 & status == "i" &
                                        male == 0 & age >= 45)/
        sum(active == 1 & age >= 45 & male == 0)
      dat$epi$i.prev.feml.1829[at] <- sum(active == 1 & status == "i" & male == 0 & age < 30)/
        sum(active == 1 & male == 0 & age < 30)
      dat$epi$i.prev.feml.30pl[at] <- sum(active == 1 & status == "i" & male == 0 & age >= 30)/
        sum(active == 1 & male == 0 & age >= 30)

      dat$epi$i.prev.male.1824[at] <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 18 & age < 25)/
        sum(active == 1 & age >= 18 & age < 25 & male == 0)
      dat$epi$i.prev.male.2529[at] <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 25 & age < 30)/
        sum(active == 1 & age >= 25 & age < 30 & male == 0)
      dat$epi$i.prev.male.3034[at] <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 30 & age < 35)/
        sum(active == 1 & age >= 30 & age < 35 & male == 0)
      dat$epi$i.prev.male.3539[at] <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 35 & age < 40)/
        sum(active == 1 & age >= 35 & age < 40 & male == 0)
      dat$epi$i.prev.male.4044[at] <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 40 & age < 45)/
        sum(active == 1 & age >= 40 & age < 45 & male == 0)
      dat$epi$i.prev.male.45pl[at] <- sum(active == 1 & status == "i" &
                                        male == 1 & age >= 45)/
        sum(active == 1 & age >= 45 & male == 0)
      dat$epi$i.prev.male.1829[at] <- sum(active == 1 & status == "i" & male == 0 & age < 30)/
        sum(active == 1 & male == 1 & age < 30)
      dat$epi$i.prev.male.30pl[at] <- sum(active == 1 & status == "i" & male == 0 & age >= 30)/
        sum(active == 1 & male == 1 & age >= 30)
    }


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


    ### Age Mixing ###
    relage <- get_male_relage(dat$nw, at)
    dat$epi$male.relage.mean[at] <- mean(relage)
    dat$epi$male.relage.mode[at] <- dens_mode(relage)

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

dens_mode <- function(data) {
  dens <- density(data)
  dens$x[which.max(dens$y)]
}


get_male_relage <- function(nD, at) {

  el <- get.dyads.active(nD, at = at)

  ages <- get.vertex.attribute(nD, "age")
  ageel <- matrix(ages[el], ncol = 2)

  male <- get.vertex.attribute(nD, "male")
  maleel <- matrix(male[el], ncol = 2)

  male.in.col1 <- which(maleel[, 1] == 1)
  ageel[male.in.col1,] <- ageel[male.in.col1, 2:1]

  out <- apply(ageel, 1, function(x) x[2] - x[1])
  return(out)
}
