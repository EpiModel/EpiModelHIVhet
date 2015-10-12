
#' @title Parameters for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description Sets the simulation parameters for the stochastic
#'              network model of HIV-1 Infection among Heterosexuals in
#'              Sub-Saharan Africa for the \code{EpiModelHIV} package.
#'
#' @param time.unit Unit of time relative to one day.
#'
#' @param acute.stage.mult Acute stage multiplier for increased infectiousness
#'        above impact of heightened viral load.
#' @param aids.stage.mult AIDS stage multiplier for increased infectiousness in
#'        AIDS above impact of heightened viral load.
#'
#' @param vl.acute.topeak Time in days to peak viremia during acute infection.
#' @param vl.acute.toset Time in days to viral set point following peak viremia.
#' @param vl.acute.peak Log 10 viral load at acute peak.
#' @param vl.setpoint Log 10 viral load at set point.
#' @param vl.aidsmax Maximum log 10 viral load during AIDS.
#'
#' @param cond.prob Probability of condoms per act with partners.
#' @param cond.eff Efficacy of condoms per act in HIV prevention.
#'
#' @param act.rate.early Daily per-partnership act rate in early disease.
#' @param act.rate.late Daily per-partnership act rate in late disease.
#' @param act.rate.cd4 CD4 count at which the \code{act.rate.late} applies.
#' @param acts.rand If \code{TRUE}, will draw number of total and unprotected
#'        acts from a binomial distribution parameterized by the \code{act.rate}.
#'
#' @param circ.prob.birth Proportion of men circumcised at birth.
#' @param circ.eff Efficacy of circumcision per act in HIV prevention.
#'
#' @param tx.elig.cd4 CD4 count at which a person becomes eligible for treatment.
#' @param tx.init.cd4.mean Mean CD4 count at which person presents for care.
#' @param tx.init.cd4.sd SD of CD4 count at which person presents for care.
#' @param tx.adhere.full Proportion of people who start treatment who are fully
#'        adherent.
#' @param tx.adhere.part Of the not fully adherent proportion, the percent of time
#'        they are on medication.
#' @param tx.vlsupp.time Time in days from treatment initiation to viral suppression.
#' @param tx.vlsupp.level Log 10 viral load level at suppression.
#' @param tx.cd4.recrat.feml Rate of CD4 recovery under treatment for males.
#' @param tx.cd4.recrat.male Rate of CD4 recovery under treatment for females.
#' @param tx.cd4.decrat.feml Rate of CD4 decline under periods of non-adherence
#'        for females.
#' @param tx.cd4.decrat.male Rate of CD4 decline under periods of non-adherence
#'        for males.
#' @param tx.coverage Proportion of treatment-eligible persons who have initiated
#'        treatment.
#' @param tx.prev.eff Proportional amount by which treatment reduces infectivity
#'        of infected partner.
#'
#' @param b.rate General entry rate per day for males and females specified.
#' @param b.rate.method Method for assigning birth rates, with options of "totpop"
#'        for births as a function of the total population size, "fpop" for births
#'        as a function of the female population size, and "stgrowth" for a constant
#'        stable growth rate.
#' @param b.propmale Proportion of entries assigned as male. If NULL, then set
#'        adaptively based on the proportion at time 1.
#' @param agecat.cutoff Age in years for calculating the \code{agecat} attribute
#'        used in the network model.
#'
#' @param ds.exit.age Age at which the age-specific ds.rate is set to 1, with NA
#'        value indicating no censoring.
#' @param ds.rate.mult Simple multiplier for background death rates.
#' @param di.cd4.aids CD4 count at which late-stage AIDS occurs and the risk of
#'        mortality is governed by \code{di.cd4.rate}.
#' @param di.cd4.rate Mortality in late-stage AIDS after hitting a nadir CD4 of
#'        \code{di.cd4.aids}.
#'
#' @details This function sets the parameters for the models.
#'
#' @export
#'
param.hiv <- function(time.unit = 7,

                      acute.stage.mult = 5,
                      aids.stage.mult = 1,

                      vl.acute.topeak = 14,
                      vl.acute.toset = 107,
                      vl.acute.peak = 6.7,
                      vl.setpoint = 4.5,
                      vl.aidsmax = 7,

                      cond.prob = 0.09,
                      cond.eff = 0.78,

                      act.rate.early = 0.362,
                      act.rate.late = 0.197,
                      act.rate.cd4 = 50,
                      acts.rand = TRUE,

                      circ.prob.birth = 0.9,
                      circ.eff = 0.53,

                      tx.elig.cd4 = 350,
                      tx.init.cd4.mean = 120,
                      tx.init.cd4.sd = 40,
                      tx.adhere.full = 0.76,
                      tx.adhere.part = 0.50,
                      tx.vlsupp.time = 365/3,
                      tx.vlsupp.level = 1.5,
                      tx.cd4.recrat.feml = 11.6/30,
                      tx.cd4.recrat.male = 9.75/30,
                      tx.cd4.decrat.feml = 11.6/30,
                      tx.cd4.decrat.male = 9.75/30,
                      tx.coverage = 0.3,
                      tx.prev.eff = 0.96,

                      b.rate = 0.03/365,
                      b.rate.method = "totpop",
                      b.propmale = NULL,
                      agecat.cutoff = 35,

                      ds.exit.age = 55,
                      ds.rate.mult = 1,
                      di.cd4.aids = 50,
                      di.cd4.rate = 2/365) {

  ## Process parameters
  p <- list()
  formal.args <- formals(sys.function())
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }


  ## trans.rate multiplier
  p$trans.rate <- p$trans.rate * p$trans.rate.mult


  ## Death rate transformations
  ltGhana <- EpiModelHIV::ltGhana
  ds.rates <- ltGhana[ltGhana$year == 2011, ]
  ds.rates$mrate <- ds.rates$mrate / 365
  if (is.numeric(ds.exit.age)) {
    ds.rates$mrate[ds.rates$agStart >= ds.exit.age] <- 1
  }
  ds.rates$reps <- ds.rates$agEnd - ds.rates$agStart + 1
  ds.rates$reps[ds.rates$agStart == 100] <- 1
  male <- rep(ds.rates$male, ds.rates$reps)
  mrate <- rep(ds.rates$mrate, ds.rates$reps)
  mrate <- pmin(1, mrate * ds.rate.mult)
  age <- rep(0:100, 2)
  ds.rates <- data.frame(male = male, age, mrate = mrate)
  ds.rates <- ds.rates[ds.rates$age != 0, ]
  p$ds.rates <- ds.rates

  ## Time unit scaling
  if (time.unit > 1) {

    ## Rates multiplied by time unit
    p$act.rate.early <- act.rate.early * time.unit
    p$act.rate.late <- act.rate.late * time.unit
    p$b.rate <- b.rate * time.unit
    p$ds.rates$mrate <- ifelse(p$ds.rates$mrate < 1,
                               p$ds.rates$mrate * time.unit,
                               p$ds.rates$mrate)

    p$dx.prob.feml <- p$dx.prob.feml * time.unit
    p$dx.prob.male <- p$dx.prob.male * time.unit
    p$tx.cd4.recrat.feml <- tx.cd4.recrat.feml * time.unit
    p$tx.cd4.recrat.male <- tx.cd4.recrat.male * time.unit
    p$tx.cd4.decrat.feml <- tx.cd4.decrat.feml * time.unit
    p$tx.cd4.decrat.male <- tx.cd4.decrat.male * time.unit
    p$di.cd4.rate <- di.cd4.rate * time.unit

    ## Intervals divided by time unit
    p$vl.acute.topeak <- vl.acute.topeak / time.unit
    p$vl.acute.toset <- vl.acute.toset / time.unit

    p$tx.vlsupp.time <- tx.vlsupp.time / time.unit

  }

  p$model <- "a2"

  class(p) <- "param.net"
  return(p)
}


#' @title Initial Conditions for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description This function sets the initial conditions for the stochastic
#'              network models in the \code{epimethods} package.
#'
#' @param i.prev.male Prevalence of initially infected males.
#' @param i.prev.feml Prevalence of initially infected females.
#' @param status.rand If \code{TRUE}, set the number of initial infected from a
#'        series of binomial draws, else set the number deterministically.
#' @param inf.time.dist Probability distribution for setting time of infection
#'        for nodes infected at T1, with options of \code{"geometric"} for randomly
#'        distributed on a geometric distribution with a probability of the
#'        reciprocal of the average length of infection, \code{"uniform"} for a
#'        uniformly distributed time over that same interval, or \code{"allacute"} for
#'        placing all infections in the acute stage at the start.
#' @param max.inf.time Maximum infection time in days for infection at initialization,
#'        used when \code{inf.time.dist} is \code{"geometric"} or \code{"uniform"}.
#' @param reinit.age.bystatus If \code{TRUE}, will reinitialize ages by disease
#'        status (infected tend to be older than susceptible) in initialization
#'        module.
#'
#' @details This function sets the initial conditions for the models.
#'
#' @export
#'
init.hiv <- function(i.prev.male = 0.05,
                     i.prev.feml = 0.05,
                     status.rand = FALSE,
                     inf.time.dist = "geometric",
                     max.inf.time = 5 * 365,
                     reinit.age.bystatus = FALSE) {

  ## Process parameters
  p <- list()
  formal.args <- formals(sys.function())
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }


  ## Parameter checks
  if (!(inf.time.dist %in% c("uniform", "geometric", "allacute"))) {
    stop("inf.time.dist must be \"uniform\" or \"geometric\" or \"allacute\" ")
  }

  class(p) <- "init.net"
  return(p)
}


#' @title Control Settings for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description This function sets the control settings for the stochastic
#'              network models in the \code{epimethods} package.
#'
#' @param simno Simulation ID number.
#' @param nsteps Number of time steps to simulate the model over in whatever unit
#'        implied by \code{time.unit}.
#' @param start Starting time step for simulation
#' @param nsims Number of simulations.
#' @param ncores Number of parallel cores to use for simulation jobs, if using
#'        the \code{EpiModel.hpc} package.
#' @param par.type Parallelization type, either of \code{"single"} for multi-core
#'        or \code{"mpi"} for multi-node MPI threads.
#' @param initialize.FUN Module to initialize the model at time 1.
#' @param aging.FUN Module to age active nodes.
#' @param cd4.FUN CD4 progression module.
#' @param vl.FUN HIV viral load progression module.
#' @param dx.FUN HIV diagnosis module.
#' @param tx.FUN HIV treatment module.
#' @param deaths.FUN Module to simulate death or exit.
#' @param births.FUN Module to simulate births or entries.
#' @param edges_correct.FUN Module to adjust the edges coefficient in response
#'        to changes to the population size, with the default function of
#'        \code{\link{edges_correct.hiv}} that preserves mean degree.
#' @param updatenwp.FUN In progress
#' @param resim_nets.FUN Module to resimulate the network at each time step.
#' @param infection.FUN Module to simulate disease infection.
#' @param get_prev.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{prevalence.hiv}}.
#' @param verbose.FUN Module to print simulation progress to screen, with the
#'        default function of \code{\link{verbose.hiv}}.
#' @param resim.int Interval unit for resimulation of network, relative to
#'        \code{time.unit} in \code{\link{param.hiv}}.
#' @param calc.asprev If \code{TRUE}, add age and sex-specific binned prevalence
#'        calculations to epi summary output.
#' @param clin.array If \code{TRUE}, save an array of individual-level attribute
#'        history for treatment status, viral load, and CD4 level.
#' @param delete.nodes If \code{TRUE}, delete vertices upon exit/death, otherwise
#'        deactivate them.
#' @param save.int Interval of time steps at which dat object should be saved to
#'        disk.
#' @param keep.cpdata Keep checkpoint interval data.
#' @param save.network Save out network objects (in \code{simOut} function).
#' @param save.transmat Save out transmission matrix.
#' @param save.nwstats Save out network statistics.
#' @param save.other Other list elements of dat to save out.
#' @param verbose If \code{TRUE}, print progress to console.
#' @param verbose.int Interval for printing progress to console.
#' @param nwstats.formula Formula for network statistics on the network.
#' @param skip.check If \code{TRUE}, skips the error check for parameter values,
#'        initial conditions, and control settings before running the models.
#'        This is suggested only if encountering unnecessary errors when running
#'        new models.
#' @param ... Additional arguments passed to the function.
#'
#' @details This function sets the parameters for the models.
#'
#' @export
#'
control.hiv <- function(simno = 1,
                        nsteps = 100,
                        start = 1,
                        nsims = 1,
                        ncores = 1,
                        par.type = "single",
                        initialize.FUN = initialize.hiv,
                        aging.FUN = aging.hiv,
                        cd4.FUN = cd4.hiv,
                        vl.FUN = vl.hiv,
                        dx.FUN = dx.hiv,
                        tx.FUN = tx.hiv,
                        deaths.FUN = deaths.hiv,
                        births.FUN = births.hiv,
                        edges_correct.FUN = edges_correct.hiv,
                        updatenwp.FUN = NULL,
                        resim_nets.FUN = simnet.hiv,
                        infection.FUN = infect.hiv,
                        get_prev.FUN = prevalence.hiv,
                        verbose.FUN = verbose.hiv,
                        resim.int = 1,
                        calc.asprev = FALSE,
                        clin.array = FALSE,
                        delete.nodes = TRUE,
                        save.int = NULL,
                        keep.cpdata = TRUE,
                        save.transmat = FALSE,
                        save.network = TRUE,
                        save.nwstats = TRUE,
                        save.other = "attr",
                        verbose = TRUE,
                        verbose.int = 100,
                        nwstats.formula = "formation",
                        skip.check = TRUE,
                        ...) {

  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names.dot.args, value = TRUE)

  if (clin.array == TRUE) {
    p$save.other <- c(p$save.other, "clin")
  }

  class(p) <- "control.net"
  return(p)
}

