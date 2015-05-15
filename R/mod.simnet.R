
#' @title Network Resimulation Module
#'
#' @description Module function to resimulate the dynamic network forward one
#'              time step conditional on current network structure and vertex
#'              attributes.
#'
#' @inheritParams aging.hiv
#'
#' @export
#'
simnet.hiv <- function(dat, at) {

  resim.int <- dat$control$resim.int
  if (at > 1 & at %% resim.int > 0) {
    return(dat)
  }

  # Delete Nodes ------------------------------------------------------------
  if (at > 1 & dat$control$delete.nodes == TRUE) {
    dat$nw <- network.extract(dat$nw, at = at)
    inactive <- which(dat$attr$active == 0)
    dat$attr <- deleteAttr(dat$attr, inactive)
  }

  # Resimulation ------------------------------------------------------------
  nwparam <- get_nwparam(dat)
  if (at == 1) {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.crude)
  } else {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.adj)
  }

  suppressWarnings(
    dat$nw <- simulate(
      dat$nw,
      formation = nwparam$formation,
      dissolution = nwparam$dissolution,
      coef.form = as.numeric(nwparam$coef.form),
      coef.diss = coef.diss,
      constraints = nwparam$constraints,
      time.start = at,
      time.slices = 1 * resim.int,
      time.offset = 0,
      output = "networkDynamic",
      monitor = dat$control$nwstats.formula))

  if (at == 1) {
    dat$stats$nwstats <- as.data.frame(attributes(dat$nw)$stats)
  } else {
    dat$stats$nwstats <- rbind(dat$stats$nwstats,
                               tail(attributes(dat$nw)$stats, 1 * resim.int))
  }


  return(dat)
}
