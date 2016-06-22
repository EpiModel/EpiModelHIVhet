
#' @title Infection Module
#'
#' @description Module function to simulate transmission over an active discordant
#'              edgelist.
#'
#' @inheritParams aging_het
#'
#' @export
#'
infect_het <- function(dat, at) {

  ## Discordant Edgelist
  del <- discord_edgelist_het(dat, at)

  nInf <- 0
  idsInf <- idsTrans <- NULL

  if (!is.null(del)) {

    ## Acts
    del <- acts(dat, del, at)

    ## Transmission
    dat <- trans(dat, del, at)
    del <- dat$del
    dat$del <- NULL

    ## Update Nodal Attr
    idsInf <- unique(del$sus)
    idsTrans <- unique(del$inf)
    nInf <- length(idsInf)

    if (nInf > 0) {
      dat$attr$status[idsInf] <- 1
      dat$attr$infTime[idsInf] <- at
      dat$attr$ageInf[idsInf] <- dat$attr$age[idsInf]
      dat$attr$dxStat[idsInf] <- 0
      dat$attr$vlLevel[idsInf] <- 0
      dat$attr$txCD4min[idsInf] <-
        pmin(rnbinom(nInf,
                     size = nbsdtosize(dat$param$tx.init.cd4.mean,
                                       dat$param$tx.init.cd4.sd),
                     mu = dat$param$tx.init.cd4.mean),
             dat$param$tx.elig.cd4)
    }

    ## Transmission data frame
    if (dat$control$save.transmat == TRUE) {
      if (nInf > 0) {
        if (at == 2) {
          dat$stats$transmat <- as.data.frame(del)
        } else {
          dat$stats$transmat <- rbind(dat$stats$transmat, as.data.frame(del))
        }
      }
    }

  }

  ## Incidence vector
  dat$epi$si.flow[at] <- nInf
  dat$epi$si.flow.male[at] <- sum(dat$attr$male[idsInf] == 1, na.rm = TRUE)
  dat$epi$si.flow.feml[at] <- sum(dat$attr$male[idsInf] == 0, na.rm = TRUE)

  return(dat)
}


discord_edgelist_het <- function(dat, at) {

  status <- dat$attr$status
  active <- dat$attr$active

  idsInft <- which(active == 1 & status == 1)
  nInft <- length(idsInft)

  del <- NULL

  if (nInft > 0) {

    if (is.null(dat$el)) {
      el <- get.dyads.active(dat$nw, at = at)
    } else {
      el <- dat$el
    }

    if (nrow(el) > 0) {
      el <- el[sample(1:nrow(el)), , drop = FALSE]

      disc <- which(abs(status[el[, 1]] - status[el[, 2]]) == 1)
      if (length(disc) > 0) {
        tmp.del <- el[disc, ]
        tmp.del[status[tmp.del[, 2]] == 1, ] <- tmp.del[status[tmp.del[, 2]] == 1, 2:1]

        del <- list()
        del$sus <- tmp.del[, 2]
        del$inf <- tmp.del[, 1]
      }
    }

  }

  return(del)
}


acts <- function(dat, del, at) {

  # Variables
  nedges <- length(del[[1]])

  act.rate.early <- dat$param$act.rate.early
  act.rate.late <- dat$param$act.rate.late
  act.rate.cd4 <- dat$param$act.rate.cd4

  cd4Count <- dat$attr$cd4Count[del$inf]

  isLate <- which(cd4Count < act.rate.cd4)

  rates <- rep(act.rate.early, nedges)
  rates[isLate] <- act.rate.late


  # Process
  act.rand <- dat$param$acts.rand
  if (act.rand == TRUE) {
    numActs <- rpois(nedges, rates)
  } else {
    numActs <- rates
  }

  cond.prob <- dat$param$cond.prob
  cond.prob <- rep(cond.prob, nedges)

  del$numActs <- numActs

  if (act.rand == TRUE) {
    del$protActs <- rbinom(nedges, rpois(nedges, numActs), cond.prob)
  } else {
    del$protActs <- numActs * cond.prob
  }

  del$protActs <- pmin(numActs, del$protActs)
  del$unprotActs <- numActs - del$protActs

  stopifnot(all(del$unprotActs >= 0))

  return(del)
}


trans <- function(dat, del, at) {

  nedges <- length(del[[1]])
  if (nedges == 0) {
    return(del)
  }

  # Base transmission probability
  vlLevel <- dat$attr$vlLevel[del$inf]
  males <- dat$attr$male[del$sus]
  ages <- dat$attr$age[del$sus]
  circs <- dat$attr$circStat[del$sus]
  prop.male <- dat$epi$propMale[at - 1]
  base.tprob <- hughes_tp(vlLevel, males, ages, circs, prop.male)

  # Acute and aids stage multipliers
  acute.stage.mult <- dat$param$acute.stage.mult
  aids.stage.mult <- dat$param$aids.stage.mult

  isAcute <- which(at - dat$attr$infTime[del$inf] <
                     (dat$param$vl.acute.topeak + dat$param$vl.acute.toset))
  isAIDS <- which(dat$attr$cd4Count[del$inf] < 200)

  base.tprob[isAcute] <- base.tprob[isAcute] * acute.stage.mult
  base.tprob[isAIDS] <- base.tprob[isAIDS] * aids.stage.mult


  # Condoms
  # Probability as a mixture function of protected and unprotected acts
  cond.eff <- dat$param$cond.eff
  prob.stasis.protacts <- (1 - base.tprob*(1 - cond.eff)) ^ del$protActs
  prob.stasis.unptacts <- (1 - base.tprob) ^ del$unprotActs
  prob.stasis <- prob.stasis.protacts * prob.stasis.unptacts
  finl.tprob <- 1 - prob.stasis

  # Transmission
  del$base.tprob <- base.tprob
  del$finl.tprob <- finl.tprob

  stopifnot(length(unique(sapply(del, length))) == 1)

  # Random transmission given final trans prob
  idsTrans <- which(rbinom(nedges, 1, del$finl.tprob) == 1)

  # Subset discord edgelist to transmissions
  del <- keep.attr(del, idsTrans)

  dat$del <- del

  return(dat)
}


hughes_tp <- function(vls, susmales, susages, suscircs, prop.male, fmat = FALSE) {

  suscircs[is.na(suscircs)] <- 0

  sus.hsv2 <- 0.59*prop.male + 0.86*(1 - prop.male)
  sus.gud <- 0.039*prop.male + 0.053*(1 - prop.male)
  sus.tvagin <- 0.068*prop.male + 0.12*(1 - prop.male)
  sus.cerv <- 0.066*(1 - prop.male)

  interc <- -8.3067
  coef.vl <- 1.062566
  coef.male <- 0.6430989
  coef.age <- -0.0403451
  coef.hsv2 <- 0.7625081
  coef.circ <- -0.6377294
  coef.gud <- 0.9749536
  coef.vagin <- 0.9435334
  coef.cerv <- 1.288279

  tp.full <- exp(interc + coef.vl*(vls - 4) +
                   coef.male*susmales + coef.age*(susages - 35) +
                   coef.hsv2*sus.hsv2 + coef.circ*susmales*suscircs +
                   coef.gud*sus.gud + coef.vagin*sus.tvagin +
                   coef.cerv*sus.cerv)

  if (fmat == TRUE) {
    tp.full <- data.frame(tp.full, vls, susmales, susages, suscircs)
  }

  return(tp.full)
}
