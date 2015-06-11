
#' @title Calculate Network Statistics
#'
#' @description This function calculates the target statistics for the formation
#'              and dissolution models estimated in \code{netest}.
#'
#' @param n Population size.
#' @param prop.male Percent of the population that is male.
#' @param ages.male Vector of ages for \code{init_nw}.
#' @param ages.feml Vector of ages for \code{init_nw}.
#' @param start.prev Starting HIV prevalence in the population.
#' @param meandeg.male Mean degree of men.
#' @param meandeg.feml Mean degree of females.
#' @param absdiff.offst Mean relative age difference of men over women.
#' @param absdiff.remain Average remainer of absdiffby term after offset.
#' @param prop.conc.male Proportion of men with concurrency.
#' @param prop.conc.feml Proportion of women with concurrency.
#' @param agecat.cutoff Age in years for calculating the \code{agecat} attribute
#'        used in the network model.
#' @param part.dur Mean duration of partnerships.
#' @param time.unit Time unit used, relative to days.
#' @param use.constraints If \code{TRUE}, use the \code{~bd(maxout = 3)} constraint.
#'
#' @export
#'
make_nw.hiv <- function(n = 10000,
                        prop.male = 0.451,
                        ages.male = seq(18, 50, 7/365),
                        ages.feml = seq(18, 50, 7/365),
                        start.prev = 0.05,
                        meandeg.male = 1.005,
                        meandeg.feml = 0.826,
                        absdiff.offst = 5.38,
                        absdiff.remain = 4.11,
                        prop.conc.male = 0.178,
                        prop.conc.feml = 0.028,
                        agecat.cutoff = 35,
                        part.dur = 1125,
                        time.unit = 7,
                        use.constraints = TRUE) {

  # Core demographics

  nMale <- round(n * prop.male)
  nFeml <- n - nMale

  male <- rep(0, n)
  male[sample(1:n, nMale)] <- 1

  age.male <- sample(ages.male, nMale, TRUE)
  age.feml <- sample(ages.feml, nFeml, TRUE)

  age <- rep(NA, n)
  age[male == 1] <- age.male
  age[male == 0] <- age.feml

  agecat <- rep(NA, length(age))
  agecat[male == 0 & age < 30] <- 0
  agecat[male == 0 & age >= 30] <- 1
  agecat[male == 1 & age < 30] <- 2
  agecat[male == 1 & age >= 30] <- 3

  # Set vertex attributes on NW
  nw <- set.vertex.attribute(nw,
                             attrname = c("male", "age", "agecat"),
                             value = list(male = male,
                                          age = age,
                                          agecat = agecat))

  # Formation Model
  formation <- as.formula(paste("~ edges +
                                concurrent(by = 'male') +
                                absdiffby('age', 'male',", absdiff.offst,") +
                                nodefactor('agecat', base = c(1, 3)) +
                                offset(nodematch('male'))"))

  ### Target stats
  edges <- (meandeg.male * nMale/2) + (meandeg.feml * nFeml/2)
  absdiff <- absdiff.remain * edges
  conc <- c(nFeml * prop.conc.feml, nMale * prop.conc.male)
  nf <- unname(table(agecat)) * c(meandeg.feml, meandeg.feml,
                                  meandeg.male, meandeg.male)

  stats <- unname(c(edges = edges, conc = conc, nf = nf, absdiff = absdiff))


  # Constraints
  if (use.constraints == TRUE) {
    constraints <- ~ bd(maxout = 3)
  } else {
    constraints <- ~.
  }

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
  out$constraints <- constraints
  out$coef.diss <- coef.diss

  class(out) <- "basenet"
  return(out)
}
