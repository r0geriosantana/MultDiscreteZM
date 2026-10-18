## Objective and analytic gradient for the MZMCMP_mu regression, with a cache
## so that optim(fn, gr) evaluates the likelihood and the gradient in a single
## C++ pass at each parameter value (internal helper).
zmreg_fg <- function(xmat, y1, y2, y3) {
  last_par <- NULL; last <- NULL
  eval_at <- function(par) {
    if (is.null(last_par) || !identical(par, last_par)) {
      last <<- logzmRegIu_vg(par = as.numeric(par), xmat = xmat, y1 = y1, y2 = y2, y3 = y3)
      last_par <<- par
    }
    last
  }
  list(fn = function(par) eval_at(par)[1],
       gr = function(par) eval_at(par)[-1])
}
