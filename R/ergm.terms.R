
#' @export
InitErgmTerm.absdiffby <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=FALSE,
                      varnames = c("attrname", "by", "assym"),
                      vartypes = c("character", "character", "numeric"),
                      required = c(TRUE, TRUE, TRUE),
                      defaultvalues = list(NULL, NULL, NULL))

  nodecov <- get.node.attr(nw, a$attrname)
  nodeby <- get.node.attr(nw, a$by)
  coef.names <- paste("absdiffby", a$attrname, a$by, sep=".")

  list(name = "absdiffby",
       coef.names = coef.names,
       pkgname = "EpiModel.hiv",
       inputs = c(a$assym, nodecov, nodeby),
       dependence = FALSE,
       emptynwstats = 0
  )
}
