
#' @title Births Module
#'
#' @description Module for simulating births/entries into the population, including
#'              initialization of attributes for incoming nodes.
#'
#' @inheritParams aging.hiv
#'
#' @export
#'
births.hiv <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  b.rate.method <- dat$param$b.rate.method
  b.rate <- dat$param$b.rate
  active <- dat$attr$active
  currNwSize <- network::network.size(dat$nw)


  # Process -----------------------------------------------------------------
  nBirths <- 0
  if (b.rate.method == "stgrowth") {
    exptPopSize <- dat$epi$num[1] * (1 + b.rate*at)
    numNeeded <- exptPopSize - sum(active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    }
  }
  if (b.rate.method == "totpop") {
    nElig <- dat$epi$num[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }
  if (b.rate.method == "fpop") {
    nElig <- dat$epi$num.feml[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }


  # Update Attr -------------------------------------------------------------
  if (nBirths > 0) {
    dat <- setBirthAttr(dat, at, nBirths)
  }


  # Update Network ----------------------------------------------------------
  if (nBirths > 0) {
    newIds <- (currNwSize + 1):(currNwSize + nBirths)

    stopifnot(unique(sapply(dat$attr, length)) == (currNwSize + nBirths))

    dat$nw <- networkDynamic::add.vertices.active(x = dat$nw, nv = nBirths,
                                                  onset = at, terminus = Inf)

    dat$nw <- network::set.vertex.attribute(x = dat$nw,
                                   attrname = c("male", "age"),
                                   value = list(male = dat$attr$male[newIds],
                                                age = dat$attr$age[newIds]),
                                   v = newIds)

  }


  # Output ------------------------------------------------------------------
  dat$epi$b.flow[at] <- nBirths

  return(dat)
}



#' @title Assign Vertex Attributes at Network Entry
#'
#' @description Assigns vertex attributes to incoming nodes at birth/entry into
#'              the network.
#'
#' @inheritParams births.hiv
#' @param nBirths Number of new births as determined by \code{\link{births.hiv}}.
#'
#' @export
#'
setBirthAttr <- function(dat, at, nBirths) {

  # Set attributes for new births to NA
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nBirths))
  })
  newIds <- which(is.na(dat$attr$active))


  # Network Status ----------------------------------------------------------
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$entTime[newIds] <- rep(at, nBirths)


  # Demography --------------------------------------------------------------
  prop.male <- ifelse(is.null(dat$param$b.propmale),
                      dat$epi$propMale[1],
                      dat$param$b.propmale)
  dat$attr$male[newIds] <- rbinom(nBirths, 1, prop.male)

  if (dat$param$b.age.const == TRUE) {
    idsMale <- intersect(newIds, which(dat$attr$male == 1))
    nMale <- length(idsMale)
    if (nMale > 0) {
      adm <- dat$temp$age.dens.male
      dat$attr$age[idsMale] <- sample(adm$x, nMale, TRUE, adm$y)
    }

    idsFeml <- intersect(newIds, which(dat$attr$male == 0))
    nFeml <- length(idsFeml)
    if (nFeml > 0) {
      adf <- dat$temp$age.dens.feml
      dat$attr$age[idsFeml] <- sample(adf$x, nFeml, TRUE, adf$y)
    }

  } else {
    dat$attr$age[newIds] <- rep(18, nBirths)
  }


  # Circumcision ------------------------------------------------------------
  dat <- circ(dat, at)


  # Epi/Clinical ------------------------------------------------------------
  dat$attr$status[newIds] <- rep("s", nBirths)

  if (length(unique(sapply(dat$attr, length))) != 1) {
    sapply(dat$attr, length)
    stop("Attribute dimensions not unique")
  }

  if (dat$control$clin.array == TRUE) {
    for (i in 1:3) {
      dat$clin[[i]] <- rbind(dat$clin[[i]], matrix(NA, nBirths, dat$control$nsteps))
    }
  }


  return(dat)
}
