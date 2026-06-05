# ==========================================================================
# Symbolic model inversion via SymPy (Python), reticulate integration
#
# - .onLoad() declares the Python dependency via py_require()
# - A lazy cache holds the imported Python module
# - An accessor function initialises it on first use (delay_load)
# - R functions call into the cached module
# ==========================================================================


# ---- Package initialization -----------------------------------------------

#' @keywords internal
#' @importFrom reticulate py_require
#' @noRd
.onLoad <- function(libname, pkgname) {
  reticulate::py_require("sympy")
}


# ---- Lazy Python module cache ---------------------------------------------

#' @keywords internal
#' @noRd
.blotit_py_cache <- new.env(parent = emptyenv())

#' Get the cached symbolic_inverse Python module
#'
#' On first call the module is imported from \code{inst/python/} with
#' \code{delay_load = TRUE}; subsequent calls return the cached reference.
#'
#' @return A reticulate Python module reference.
#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_symbolic_inverse_py <- function() {
  if (!exists("symbolic_inverse", envir = .blotit_py_cache,
              inherits = FALSE)) {
    .blotit_py_cache$symbolic_inverse <-
      reticulate::import_from_path(
        "symbolic_inverse",
        path = system.file("python", package = "blotIt"),
        delay_load = TRUE
      )
  }
  .blotit_py_cache$symbolic_inverse
}


# ---- symbolicInverse() ----------------------------------------------------

#' Analytically invert a model equation using SymPy via reticulate
#'
#' Given a model string like \code{"yi / sj"} and the name of the
#' biological parameter (e.g. \code{"yi"}), this function calls the
#' Python module \code{symbolic_inverse.invert_model()} to solve the
#' equation \code{model == value} for \code{yi}.  It returns the
#' analytical inverse together with its derivative w.r.t. \code{value}
#' (needed for Gaussian error propagation when scaling sigma).
#'
#' If SymPy cannot find an analytical solution (e.g. transcendental
#' equations), the function returns \code{NULL} and the caller should
#' fall back to numerical inversion.
#'
#' @param model Character string of the model, e.g. \code{"yi / sj"}.
#'   This is the \emph{original} model before log-scale or index
#'   substitutions.
#' @param bioPar Character, the biological parameter to solve for.
#' @param allSymbols Character vector of all symbols in the model
#'   (excluding data columns).
#' @param verbose Logical, print diagnostic information.
#'
#' @return If successful, a list with:
#' \describe{
#'   \item{\code{inverseExpr}}{Character string: R expression for the
#'     inverse.  The variable \code{value} represents the observed
#'     measurement.}
#'   \item{\code{inverseDerivExpr}}{Character string: R expression for
#'     \code{d(inverse)/d(value)}.}
#' }
#' Returns \code{NULL} if analytical inversion fails or reticulate /
#' SymPy are unavailable.
#'
#' @noRd
symbolicInverse <- function(model,
                            bioPar,
                            allSymbols,
                            verbose = FALSE) {

  ## Obtain the Python module (lazy-loaded)
  pymod <- tryCatch(
    get_symbolic_inverse_py(),
    error = function(e) NULL
  )
  if (is.null(pymod)) {
    if (verbose) {
      cat("Python module 'symbolic_inverse' not available,",
          "skipping analytic inversion.\n")
    }
    return(NULL)
  }

  ## Call Python: invert_model(model_str, bio_par, all_symbols)
  result <- tryCatch(
    pymod$invert_model(
      model_str = model,
      bio_par = bioPar,
      all_symbols = as.list(unname(allSymbols))
    ),
    error = function(e) {
      if (verbose) {
        cat("Python invert_model() failed: ",
            conditionMessage(e), "\n")
      }
      NULL
    }
  )

  ## Python returns None (mapped to NULL) when no solution exists
  if (is.null(result)) {
    if (verbose) cat("No analytical inverse found for '", bioPar, "'.\n")
    return(NULL)
  }

  ## Convert reticulate dict to R list
  result <- as.list(result)

  inverseStr <- result$inverse
  inverseDerivStr <- result$inverse_deriv

  if (verbose) {
    cat("Symbolic inverse:       ", bioPar, "=", inverseStr, "\n")
    cat("Symbolic inverse deriv: d/dvalue =", inverseDerivStr, "\n")
  }

  return(list(
    inverseExpr = inverseStr,
    inverseDerivExpr = inverseDerivStr
  ))
}


# ---- prepareAnalyticInverse() ---------------------------------------------

#' Prepare the analytical inverse expressions for use in scaleTarget
#'
#' Takes the raw output of \code{symbolicInverse} and applies the same
#' transformations that \code{alignReplicates} applies to the forward
#' model: (1) log-scale substitution if \code{fitLogscale = TRUE},
#' (2) index-bracket substitution for effect categories.
#'
#' @param symbolicResult Output of \code{symbolicInverse}, or \code{NULL}.
#' @param parameters Named character vector (names = effect category,
#'   values = parameter names), as constructed in \code{alignReplicates}.
#' @param fitLogscale Logical.
#'
#' @return A list with parsed R expressions \code{inverseExpr} and
#'   \code{inverseDerivExpr} ready for \code{eval()}, or \code{NULL}.
#'
#' @noRd
prepareAnalyticInverse <- function(symbolicResult,
                                   parameters,
                                   fitLogscale) {

  if (is.null(symbolicResult)) return(NULL)

  inverseStr <- symbolicResult$inverseExpr
  inverseDerivStr <- symbolicResult$inverseDerivExpr

  ## (1) Apply log-scale transformation: replace each parameter p by exp(p)
  if (fitLogscale) {
    inverseStr <- replaceSymbols(
      parameters,
      paste0("exp(", parameters, ")"),
      inverseStr
    )
    inverseDerivStr <- replaceSymbols(
      parameters,
      paste0("exp(", parameters, ")"),
      inverseDerivStr
    )
  }

  ## (2) Apply index-bracket substitution: p -> p[effect_category]
  for (n in seq_along(parameters)) {
    bracket <- paste0(parameters[n], "[", names(parameters)[n], "]")
    inverseStr <- replaceSymbols(parameters[n], bracket, inverseStr)
    inverseDerivStr <- replaceSymbols(parameters[n], bracket, inverseDerivStr)
  }

  ## Parse to R expressions
  inverseExpr <- tryCatch(
    parse(text = inverseStr),
    error = function(e) NULL
  )
  inverseDerivExpr <- tryCatch(
    parse(text = inverseDerivStr),
    error = function(e) NULL
  )

  if (is.null(inverseExpr) || is.null(inverseDerivExpr)) return(NULL)

  return(list(
    inverseExpr = inverseExpr,
    inverseDerivExpr = inverseDerivExpr,
    inverseStr = inverseStr,
    inverseDerivStr = inverseDerivStr
  ))
}


# ---- evaluateAnalyticInverse() -------------------------------------------

#' Evaluate the analytical inverse to compute scaled values
#'
#' Given fitted parameters and the data, evaluates the analytical inverse
#' expression to obtain the biological parameter values on the common
#' scale, without any numerical root-finding.
#'
#' @param parList Named list of fitted parameter vectors (from
#'   \code{residualsFit[-1]}).
#' @param passParList The standard passParList from scaleTarget.
#' @param passParList2 The standard passParList2 from scaleTarget.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{values}}{Numeric vector of scaled values.}
#'   \item{\code{jacobian}}{Numeric vector of
#'     \code{|d(inverse)/d(value)|}.}
#' }
#' Returns \code{NULL} if analytical inverse expressions are not in
#' \code{passParList} or evaluation fails.
#'
#' @noRd
evaluateAnalyticInverse <- function(parList,
                                    passParList,
                                    passParList2) {

  inverseExpr <- passParList$analyticInverseExpr
  inverseDerivExpr <- passParList$analyticInverseDerivExpr

  if (is.null(inverseExpr) || is.null(inverseDerivExpr)) return(NULL)

  dataFit <- passParList2$dataFit

  ## Evaluate in the same environment as the forward model
  envList <- c(as.list(dataFit), parList)

  values <- tryCatch(
    with(envList, eval(inverseExpr)),
    error = function(e) NULL
  )
  jacobian <- tryCatch(
    with(envList, eval(inverseDerivExpr)),
    error = function(e) NULL
  )

  if (is.null(values) || is.null(jacobian)) return(NULL)
  if (any(!is.finite(values)) || any(!is.finite(jacobian))) return(NULL)

  return(list(
    values = as.numeric(values),
    jacobian = abs(as.numeric(jacobian))
  ))
}
