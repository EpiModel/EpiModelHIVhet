
#' @title Calculate Network Statistics
#'
#' @description This function calculates the target statistics for the formation
#'              and dissolution models estimated in \code{netest}.
#'
#' @param n Population size.
#' @param meandeg Mean degree.
#' @param prop.male Percent of the population that is male.
#' @param ages.male Vector of ages for \code{init_nw}.
#' @param ages.feml Vector of ages for \code{init_nw}.
#' @param start.prev Starting HIV prevalence in the population.
#' @param part.dur Mean duration of partnerships.
#' @param time.unit Time unit used, relative to days.
#'
#' @export
#'
make_nw_het <- function(n = 10000,
                        meandeg = 0.8,
                        prop.male = 0.5,
                        start.prev = 0.05,
                        part.dur = 1000,
                        time.unit = 7) {

  nMale <- round(n * prop.male)

  male <- rep(0, n)
  male[sample(1:n, nMale)] <- 1

  # Set vertex attributes
  nw <- network.initialize(n = n, directed = FALSE)
  nw <- set.vertex.attribute(nw, attrname = "male", value = male)

  # Formation Model
  formation <- ~edges + offset(nodematch("male"))

  # Target stats
  edges.ts <- meandeg * (n/2)

  stats <- edges.ts

  # Dissolution model
  dissolution <- ~offset(edges)
  dur <- part.dur/time.unit
  d.rate <- time.unit * (((1 - start.prev) * 1/(55 - 18)/365) + (start.prev * 1/12/365))
  coef.diss <- dissolution_coefs(dissolution, duration = dur, d.rate = d.rate)

  out <- list()
  out$nw <- nw
  out$time.unit <- time.unit
  out$formation <- formation
  out$stats <- stats
  out$coef.diss <- coef.diss

  return(out)
}
