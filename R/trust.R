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
#' @param structured logical; when \code{TRUE} the optimization is run by the
#'   sparse, structure-exploiting kernel \code{trustStructured()} instead of the
#'   dense C++ kernel. This is required when the number of parameters is large
#'   enough that a dense \eqn{K \times K} Hessian / eigensolve would not fit in
#'   memory. The objective must then return a sparse \code{hessian} plus the
#'   rank-1 \code{constrVec} (see \code{objFunction(..., sparseHessian = TRUE)}).
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
  if (!is.null(parscale)) {
    if (length(parscale) != length(parinit) || any(!is.finite(parscale)) ||
        any(parscale <= 0)) {
      stop("parscale must be a positive finite vector the length of parinit")
    }
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

  ## Final evaluation at the solution so the returned list carries every
  ## component objfun produces (gradient, hessian, prediction, sigma, ...).
  out <- objfun(res$argument, ...)
  out$argument   <- res$argument
  out$converged  <- res$converged
  out$iterations <- res$iterations
  out
}


#' Sparse, structure-exploiting trust-region kernel (Levenberg-Marquardt)
#'
#' Alternative to the dense C++ \code{trust_kernel} for parameter counts at which
#' a dense Hessian or its eigendecomposition is infeasible (dense time-course
#' alignment with \eqn{10^4}-\eqn{10^6} biological parameters). The objective
#' returns a sparse Gauss-Newton \code{hessian} (a symmetric \code{Matrix}
#' object) and the rank-1 normalization coupling \code{constrVec}; the full
#' Hessian is
#' \eqn{H = \mathrm{hessian} + 2\,vv^{\top}} with \eqn{v = }\code{constrVec}.
#'
#' Each step solves the Levenberg-damped system
#' \eqn{(H + \lambda I)\,p = -g} with a sparse Cholesky factorization of the
#' sparse part (the Gauss-Newton Hessian is positive semidefinite, so
#' \eqn{H_{\mathrm{GN}} + \lambda I} is positive definite for \eqn{\lambda > 0}
#' and the factorization never fails) and a Sherman-Morrison rank-1 update for
#' the constraint. \eqn{\lambda} is decreased on accepted steps and increased on
#' rejected ones. For the arrow-shaped Hessians blotIt produces (a diagonal
#' biological block bordered by a handful of scaling/error parameters) the
#' Cholesky fill-in is negligible, so each iteration is \eqn{O(\mathrm{nPars})}.
#'
#' @param objfun objective returning \code{value}, \code{gradient},
#'   \code{hessian} (sparse) and \code{constrVec}.
#' @param parinit numeric starting vector (names preserved).
#' @param iterlim,fterm,mterm,minimize as in \code{trust}.
#' @param ... further arguments passed to \code{objfun}.
#'
#' @return list with \code{argument}, \code{value}, \code{converged},
#'   \code{iterations}.
#' @keywords internal
#' @noRd
trustStructured <- function(objfun, parinit, iterlim = 100L,
                            fterm = sqrt(.Machine$double.eps),
                            mterm = sqrt(.Machine$double.eps),
                            minimize = TRUE, ...) {
  sgn <- if (minimize) 1 else -1
  n <- length(parinit)
  In <- Matrix::Diagonal(n)

  ## Solve (Hgn + lambda*I + 2 v v^T) x = b via sparse Cholesky of the sparse
  ## part plus a Sherman-Morrison rank-1 correction for the constraint.
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

  ## Keep parameter names throughout: objFunction()/resolveFunction() classify
  ## parameters by their names (e.g. "yi_1", "sj_1"), so unnamed vectors break
  ## the fit. Arithmetic theta + p preserves theta's names.
  theta <- parinit
  cur <- objfun(theta, ...)
  f <- sgn * cur$value
  if (!is.finite(f)) stop("trust: objective not finite at the initial parameters")

  ## Nielsen-style adaptive damping: lambda is scaled by the gain ratio between
  ## the actual and the model-predicted reduction. Initial lambda relative to
  ## the curvature scale.
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

    ## First-order optimality on the curvature-scaled gradient g_i / sqrt(H_ii),
    ## which is dimensionless and robust to the very different scales of the
    ## biological, scaling and (constraint-coupled) parameters.
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
      if (is.null(p) || any(!is.finite(p))) { lambda <- lambda * nu; nu <- nu * 2; next }

      ## Predicted reduction from the quadratic model (full Hessian, undamped):
      ## predRed = -(g'p + 1/2 p'H p), with p'H p = p'Hgn p + 2 (v'p)^2.
      quad <- as.numeric(Matrix::crossprod(p, Hgn %*% p))
      if (!is.null(v)) quad <- quad + 2 * sum(v * p)^2
      predRed <- -(sum(g * p) + 0.5 * quad)

      thetaTry <- theta + p
      trial <- objfun(thetaTry, ...)
      fTry <- sgn * trial$value
      aRed <- f - fTry
      rho <- if (is.finite(fTry) && predRed > 0) aRed / predRed else -Inf

      if (rho > 0) {
        theta <- thetaTry
        cur <- trial
        f <- fTry
        lambda <- lambda * max(1 / 3, 1 - (2 * rho - 1)^3)
        nu <- 2
        stepTaken <- TRUE
        break
      } else {
        lambda <- lambda * nu
        nu <- nu * 2
        if (lambda > 1e16) break
      }
    }

    if (!stepTaken) {        # no admissible step found -> stationary
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
