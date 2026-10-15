#' @keywords internal
#' @useDynLib blotIt, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats pchisq rnorm
#' @importFrom utils modifyList
"_PACKAGE"

# Column names used via non-standard evaluation in data.table and ggplot2
utils::globalVariables(c("name", "lower", "upper", "ymin", "ymax"))


#' Trust-region optimizer
#'
#' Replacement for [trust::trust()] that does not abort on ill-conditioned or
#' non-finite Hessians. See `vignette("blotIt")` for the method.
#'
#' @param objfun objective function returning a list with `value`, `gradient`
#'   and `hessian`. Further components (e.g. `prediction`, `sigma`) are
#'   returned from a final evaluation at the solution.
#' @param parinit numeric vector of starting parameters; names are kept.
#' @param rinit,rmax initial and maximum trust-region radius.
#' @param parscale optional positive numeric vector for rescaling parameters.
#'   `NULL` or `FALSE` disables rescaling.
#' @param iterlim maximum number of iterations.
#' @param fterm,mterm convergence tolerances on the change in objective value
#'   and on the predicted reduction.
#' @param minimize logical; minimize (default) or maximize `objfun$value`.
#' @param blather accepted for compatibility with [trust::trust()]; the
#'   iteration path is not recorded.
#' @param structured logical; if `TRUE`, use the sparse optimizer for large
#'   problems. `objfun` must then return a sparse `hessian` (a `Matrix` object)
#'   and the rank-1 term `constrVec`, the full Hessian being
#'   \eqn{H + 2vv^\top}. `rinit`, `rmax` and `parscale` are ignored.
#' @param ... further arguments passed on to `objfun`.
#'
#' @return list with `argument` (the solution), `value`, `converged`,
#'   `iterations` and all components of `objfun` at the solution.
#'
#' @example inst/examples/exampleTrust.R
#'
#' @export
trust <- function(objfun, parinit, rinit, rmax,
                  parscale = NULL,
                  iterlim = 100L,
                  fterm = sqrt(.Machine$double.eps),
                  mterm = sqrt(.Machine$double.eps),
                  minimize = TRUE,
                  blather = FALSE,
                  structured = FALSE,
                  ...) {

  if (!is.numeric(parinit) || !all(is.finite(parinit))) {
    stop("parinit must be a finite numeric vector")
  }
  if (isFALSE(parscale)) parscale <- NULL
  if (!is.null(parscale) &&
      (length(parscale) != length(parinit) || any(!is.finite(parscale)) ||
       any(parscale <= 0))) {
    stop("parscale must be a positive finite vector the length of parinit")
  }

  if (structured) {
    res <- trustStructured(
      objfun   = objfun,
      parinit  = parinit,
      iterlim  = as.integer(iterlim),
      fterm    = fterm,
      mterm    = mterm,
      minimize = minimize,
      ...
    )
  } else {
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
  }

  ## Final evaluation returns every component of objfun at the solution
  out <- objfun(res$argument, ...)
  out$argument   <- res$argument
  out$converged  <- res$converged
  out$iterations <- res$iterations
  out
}


#' Sparse Levenberg-Marquardt kernel for large problems
#'
#' @param objfun objective returning `value`, `gradient`, a sparse `hessian`
#'   and `constrVec`.
#' @param parinit numeric starting vector; names are kept.
#' @param iterlim,fterm,mterm,minimize as in [trust()].
#' @param ... further arguments passed to `objfun`.
#'
#' @return list with `argument`, `value`, `converged`, `iterations`.
#' @keywords internal
#' @noRd
trustStructured <- function(objfun, parinit, iterlim = 100L,
                            fterm = sqrt(.Machine$double.eps),
                            mterm = sqrt(.Machine$double.eps),
                            minimize = TRUE, ...) {
  sgn <- if (minimize) 1 else -1
  n <- length(parinit)
  In <- Matrix::Diagonal(n)

  ## Solves (Hgn + lambda*I + 2 v v^T) x = b: sparse Cholesky plus
  ## Sherman-Morrison for the rank-1 term. NULL if the factorization fails.
  solveDamped <- function(Hgn, v, lambda, b) {
    A <- Hgn + lambda * In
    ch <- tryCatch(Matrix::Cholesky(A, LDL = FALSE, perm = TRUE),
                   error = function(e) NULL)
    if (is.null(ch)) return(NULL)
    AinvB <- as.numeric(Matrix::solve(ch, b))
    if (is.null(v) || all(v == 0)) return(AinvB)
    AinvV <- as.numeric(Matrix::solve(ch, v))
    denom <- 1 + 2 * sum(v * AinvV)
    AinvB - (2 * sum(v * AinvB) / denom) * AinvV
  }

  ## theta must keep its names: the objective classifies parameters by name
  theta <- parinit
  cur <- objfun(theta, ...)
  f <- sgn * cur$value
  if (!is.finite(f)) stop("trust: objective not finite at the initial parameters")

  ## Nielsen damping, initial lambda relative to the curvature scale
  d0 <- Matrix::diag(cur$hessian)
  lambda <- 1e-3 * max(mean(d0[is.finite(d0)]), 1e-8)
  nu <- 2

  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(iterlim)) {
    g <- sgn * cur$gradient
    if (!all(is.finite(g))) break
    Hgn <- cur$hessian
    v <- cur$constrVec

    ## Optimality test on the curvature-scaled gradient g_i / sqrt(H_ii)
    gScale <- Matrix::diag(Hgn)
    if (!is.null(v)) gScale <- gScale + 2 * v^2
    gScale <- sqrt(pmax(as.numeric(gScale), 1e-12))
    if (max(abs(g) / gScale) < sqrt(.Machine$double.eps)) {
      converged <- TRUE
      break
    }

    stepTaken <- FALSE
    for (inner in seq_len(80L)) {
      p <- solveDamped(Hgn, v, lambda, -g)
      if (is.null(p) || any(!is.finite(p))) {
        lambda <- lambda * nu
        nu <- nu * 2
        next
      }

      ## Reduction predicted by the undamped quadratic model with the full Hessian
      quad <- as.numeric(Matrix::crossprod(p, Hgn %*% p))
      if (!is.null(v)) quad <- quad + 2 * sum(v * p)^2
      predRed <- -(sum(g * p) + 0.5 * quad)

      thetaTry <- theta + p
      trial <- objfun(thetaTry, ...)
      fTry <- sgn * trial$value
      aRed <- f - fTry
      rho <- if (is.finite(fTry) && predRed > 0) aRed / predRed else -Inf

      if (rho <= 0) {
        lambda <- lambda * nu
        nu <- nu * 2
        if (lambda > 1e16) break
        next
      }
      theta <- thetaTry
      cur <- trial
      f <- fTry
      lambda <- lambda * max(1 / 3, 1 - (2 * rho - 1)^3)
      nu <- 2
      stepTaken <- TRUE
      break
    }

    ## No accepted step at any damping: stationary point
    if (!stepTaken) {
      converged <- TRUE
      break
    }
    if (converged) break
  }

  arg <- theta
  names(arg) <- names(parinit)
  list(argument = arg, value = sgn * f, converged = converged,
       iterations = iter)
}
