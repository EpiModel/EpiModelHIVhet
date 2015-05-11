
nbsdtosize <- function(mu, sd) {
  mu ^ 2 / (sd ^ 2 - mu)
}

get_attr <- function(x, sim = 1) {
  if (is.null(x$attr)) {
    stop("No attr on x")
  } else {
    x$attr[[1]]
  }
}

cut_age <- function(age, breaks = c(0, 29, 39, Inf)) {
  cut(age, breaks = breaks, labels = FALSE)
}

keep.attr <- function(attrList, keep) {
  lapply(attrList, function(x) x[keep])
}
