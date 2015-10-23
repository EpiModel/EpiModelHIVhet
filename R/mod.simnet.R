
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

  dat <- update_nwp(dat)
  nwparam <- get_nwparam(dat)

  dat$el <- tergmLite::simulate_network(p = dat$p,
                                        el = dat$el,
                                        coef.form = nwparam$coef.form,
                                        coef.diss = nwparam$coef.diss$coef.adj)

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

update_nwp <- function(dat) {

  mf <- dat$p$model.form
  md <- dat$p$model.diss
  mhf <- dat$p$MHproposal.form
  mhd <- dat$p$MHproposal.diss

  n <- attributes(dat$el)$n
  maxdyads <- choose(n, 2)

  ## 1. Update model.form ##

  # edges
  # inputs <- c(0, 1, 0) # not changed
  mf$terms[[1]]$maxval <- maxdyads

  # nodematch
  nodecov <- dat$attr$male
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov, u, nomatch = length(u) + 1)
  inputs <- nodecov
  mf$terms[[2]]$inputs <- c(0, 1, length(inputs), inputs)

  ## Update combined maxval here
  mf$maxval <- c(maxdyads, Inf)


  ## 2. Update model.diss ##
  md$terms[[1]]$maxval <- maxdyads
  md$maxval <- maxdyads


  ## 3. Update MHproposal.form ##
  mhf$arguments$constraints$bd$attribs <-
    matrix(rep(mhf$arguments$constraints$bd$attribs[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxout <-
    matrix(rep(mhf$arguments$constraints$bd$maxout[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxin <- matrix(rep(n, n), ncol = 1)
  mhf$arguments$constraints$bd$minout <-
    mhf$arguments$constraints$bd$minin <- matrix(rep(0, n), ncol = 1)


  ## 4. Update MHproposal.diss ##
  mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd


  ## 5. Output ##
  p <- list(model.form = mf, model.diss = md,
            MHproposal.form = mhf, MHproposal.diss = mhd)

  dat$p <- p
  return(dat)
}
