
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

  nwparam <- get_nwparam(dat)
  if (at == 1) {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.crude)
  } else {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.adj)
  }

  dat$el <- tergmLite::simulate_network(p = dat$p,
                             el = dat$el,
                             coef.form = nwparam$coef.form,
                             coef.diss = coef.diss,
                             time.start = at)

  # if (at == 1) {
  #   dat$stats$nwstats <- matrix(NA, ncol = 5, nrow = dat$control$nsteps)
  #   colnames(dat$stats$nwstats) <- c("edges", "meandeg", "deg0", "deg1", "concurrent")
  # }
  # n <- attributes(dat$el)$n
  # tab <- table(dat$el)
  # dat$stats$nwstats[at, 1] <- nrow(dat$el)
  # dat$stats$nwstats[at, 2] <- nrow(dat$el)/n
  # dat$stats$nwstats[at, 4] <- sum(tab == 1)/n
  # dat$stats$nwstats[at, 5] <- sum(tab > 1)/n
  # dat$stats$nwstats[at, 3] <- (n - sum(tab == 1) - sum(tab > 1))/n

  return(dat)
}
