# Symbolic model inversion via SymPy (Python, through reticulate).


#' @keywords internal
#' @importFrom reticulate py_require
#' @noRd
.onLoad <- function(libname, pkgname) {
  reticulate::py_require("sympy")
}

#' @keywords internal
#' @noRd
.blotit_py_cache <- new.env(parent = emptyenv())

#' Cached `symbolic_inverse` Python module from `inst/python/`
#'
#' @return reticulate module reference, imported lazily on first call.
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


# symbolicInverse() --------------------------------------------------------

#' Solve `model == value` for the biological parameter with SymPy
#'
#' @param model model string before log-scale or index substitutions.
#' @param bioPar biological parameter to solve for.
#' @param allSymbols all model symbols except data columns.
#' @param verbose logical, print diagnostics.
#'
#' @return list with strings `inverseExpr` (inverse in terms of `value`) and
#'   `inverseDerivExpr` (its derivative with respect to `value`), or `NULL` if
#'   no closed form exists or SymPy is unavailable.
#'
#' @noRd
symbolicInverse <- function(model,
                            bioPar,
                            allSymbols,
                            verbose = FALSE) {
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

  if (is.null(result)) {
    if (verbose) cat("No analytical inverse found for '", bioPar, "'.\n")
    return(NULL)
  }

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


# prepareAnalyticInverse() -------------------------------------------------

#' Apply the forward-model substitutions to the analytic inverse
#'
#' @param symbolicResult output of `symbolicInverse()`, or `NULL`.
#' @param parameters character vector of parameters named by effect class.
#' @param fitLogscale logical.
#'
#' @return list with parsed `inverseExpr`, `inverseDerivExpr` and their
#'   strings, or `NULL`.
#'
#' @noRd
prepareAnalyticInverse <- function(symbolicResult,
                                   parameters,
                                   fitLogscale) {
  if (is.null(symbolicResult)) return(NULL)

  inverseStr <- symbolicResult$inverseExpr
  inverseDerivStr <- symbolicResult$inverseDerivExpr

  ## Same substitutions as the forward model: p -> exp(p), then p -> p[class]
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

  for (n in seq_along(parameters)) {
    bracket <- paste0(parameters[n], "[", names(parameters)[n], "]")
    inverseStr <- replaceSymbols(parameters[n], bracket, inverseStr)
    inverseDerivStr <- replaceSymbols(parameters[n], bracket, inverseDerivStr)
  }

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


# evaluateAnalyticInverse() ------------------------------------------------

#' Evaluate the analytic inverse at the fitted parameters
#'
#' @param parList named list of fitted parameter vectors.
#' @param passParList,passParList2 as in `scaleTarget()`.
#'
#' @return list with `values` on the common scale and `jacobian`
#'   \eqn{|d\,\mathrm{inverse}/d\,\mathrm{value}|}, or `NULL` if no analytic
#'   inverse is available or its evaluation is not finite.
#'
#' @noRd
evaluateAnalyticInverse <- function(parList,
                                    passParList,
                                    passParList2) {
  inverseExpr <- passParList$analyticInverseExpr
  inverseDerivExpr <- passParList$analyticInverseDerivExpr

  if (is.null(inverseExpr) || is.null(inverseDerivExpr)) return(NULL)

  dataFit <- passParList2$dataFit

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
