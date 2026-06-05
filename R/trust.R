#' @keywords internal
#' @useDynLib blotIt, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats pchisq rnorm
#' @importFrom utils modifyList
"_PACKAGE"

# Column / aesthetic names used via non-standard evaluation in data.table and
# ggplot2 calls.
utils::globalVariables(c("name", "lower", "upper", "ymin", "ymax"))


#' Robust trust-region optimizer
#'
#' A drop-in replacement for \code{trust::trust} backed by a compiled C++
#' kernel (see \code{src/trust_kernel.cpp}). It implements the same trust-region
#' method but does not abort when the Hessian at the current point is
#' ill-conditioned or non-finite: on any failure of the symbolic eigen-solve it
#' falls back to a steepest-descent step within the current trust region and
#' shrinks the region on rejected steps. This avoids the
#' \code{"error code 1 from Lapack routine 'dsyevr'"} crashes that
#' \code{trust::trust} throws for the near-singular fits blotIt produces when a
#' target has many scaling parameters.
#'
#' @param objfun objective function returning a list with components
#'   \code{value}, \code{gradient} and \code{hessian} (as required by
#'   \code{trust::trust}). It may return additional components (e.g.
#'   \code{prediction}, \code{sigma}); these are preserved in the result, taken
#'   from a final evaluation at the solution.
#' @param parinit numeric vector of starting parameters (names are kept).
#' @param rinit,rmax initial and maximum trust-region radius.
#' @param parscale optional positive numeric vector for rescaling parameters.
#'   \code{NULL} or \code{FALSE} disables rescaling.
#' @param iterlim maximum number of iterations.
#' @param fterm,mterm convergence tolerances on the change in objective value
#'   and on the predicted reduction.
#' @param minimize logical; minimize (default) or maximize \code{objfun$value}.
#' @param blather accepted for compatibility with \code{trust::trust}; the
#'   iteration path is not recorded.
#' @param ... further arguments passed on to \code{objfun} at every evaluation.
#'
#' @return A list with \code{argument} (the solution), \code{value},
#'   \code{converged}, \code{iterations}, plus all components returned by a final
#'   call to \code{objfun} at the solution (\code{gradient}, \code{hessian},
#'   and any extras such as \code{prediction}/\code{sigma}).
#'
#' @export
trust <- function(objfun, parinit, rinit, rmax,
                  parscale = NULL,
                  iterlim = 100L,
                  fterm = sqrt(.Machine$double.eps),
                  mterm = sqrt(.Machine$double.eps),
                  minimize = TRUE,
                  blather = FALSE,
                  ...) {

  if (!is.numeric(parinit) || !all(is.finite(parinit))) {
    stop("parinit must be a finite numeric vector")
  }
  if (isFALSE(parscale)) parscale <- NULL
  if (!is.null(parscale)) {
    if (length(parscale) != length(parinit) || any(!is.finite(parscale)) ||
        any(parscale <= 0)) {
      stop("parscale must be a positive finite vector the length of parinit")
    }
  }

  dots <- list(...)
  fn <- if (length(dots) > 0L) {
    function(x) do.call(objfun, c(list(x), dots))
  } else {
    objfun
  }

  res <- trust_kernel(
    fn        = fn,
    parinit   = parinit,
    rinit     = rinit,
    rmax      = rmax,
    iterlim   = as.integer(iterlim),
    fterm     = fterm,
    mterm     = mterm,
    minimize  = minimize,
    parscale_ = parscale
  )

  ## Final evaluation at the solution so the returned list carries every
  ## component objfun produces (gradient, hessian, prediction, sigma, ...).
  out <- objfun(res$argument, ...)
  out$argument   <- res$argument
  out$converged  <- res$converged
  out$iterations <- res$iterations
  out
}
