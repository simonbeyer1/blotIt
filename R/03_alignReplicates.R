#' Align replicate data on a common scale
#'
#' Fits an alignment model with biological, scaling and error parameters per
#' target by maximum likelihood and returns the aligned values with their
#' uncertainties. See `vignette("blotIt")` for the method.
#'
#' @param data data frame in long format, usually the output of [readWide()].
#'   Required columns are `name` and `value`, plus the columns named in
#'   `biological`, `scaling` and `error`.
#' @param model character, the alignment model, e.g. `"yi / sj"`.
#' @param errorModel character, model of the standard deviation, e.g.
#'   `"sigma0 + value * sigmaR"`. It may use parameters, data columns and
#'   `value`, which stands for the model prediction.
#' @param biological two-sided formula `par1 + par2 ~ col1 + col2` assigning
#'   the biological parameters of `model` one value per level of the data
#'   columns on the right-hand side.
#' @param scaling two-sided formula as `biological` for the scaling parameters.
#' @param error two-sided formula as `biological` for the error parameters, or
#'   a list of such formulas to give each error parameter its own levels, e.g.
#'   `list(sigmaA ~ name + modality, sigmaB ~ name)`. Use `~ name` for one value
#'   per target.
#' @param fitLogscale logical, fit the logarithm of the parameters.
#' @param normalize logical, constrain the mean of the first biological
#'   parameter to 1.
#' @param averageTechRep logical, average data points that share all
#'   biological and scaling levels.
#' @param verbose logical, print information about each fit.
#' @param namesFactored logical, return `name` columns as factors.
#' @param normalizeInput logical, normalize the data of each target before the
#'   fit: positive data are divided by their mean, other data are centered.
#' @param keepOriginalScale character, `"none"` (default), `"lin"` to rescale
#'   outputs to the mean of the input, or `"log"` to shift them to it.
#' @param outputScale character, scale of the returned aligned values:
#'   `"linear"` (default) or a transformation understood by [plotIt()].
#' @param iterlim integer, iteration limit. `optControls$iterlim` takes
#'   precedence.
#' @param optControls named list of optimizer settings: `rinit` (default 0.1),
#'   `rmax` (10), `parscale` (`FALSE`), `iterlim` (100), `fterm` and `mterm`
#'   (both `sqrt(.Machine$double.eps)`), see [trust()], and
#'   `structuredThreshold` (5000), the parameter count above which the sparse
#'   optimizer is used.
#' @param multistart integer, number of fits per target. Start 1 is the initial
#'   point, the others perturb it randomly. A converged start replaces start 1
#'   if its objective value is lower.
#' @param rootControls list of arguments passed to [nleqslv::nleqslv()] when
#'   the model has no analytic inverse, e.g. `method`, `global`, `control`.
#' @param parinit `NULL` (all parameters start at 1) or a named list of start
#'   values on the scale of the parameter table. Each entry is a number for all
#'   levels or a vector named by level (as in the `level` column of the
#'   parameter table). Parameters and levels not given start at 1.
#' @param multistartSd standard deviation of the normal perturbations of the
#'   start values (on the fit scale) for starts 2 to `multistart`.
#' @param cores integer, number of processes for the multistart fits. Values
#'   \eqn{> 1} fork processes and are not available on Windows. Results do not
#'   depend on `cores`.
#' @param likelihood `"normal"` (default) or `"t"`. With `"t"`, each data point
#'   follows a Student-t distribution with location `model` and scale
#'   `errorModel`, which reduces the influence of outliers. The `sigma` of
#'   `prediction` and `scaled` is then this scale.
#' @param df degrees of freedom of the t distribution (fixed), used only for
#'   `likelihood = "t"`.
#'
#' @return list with
#' \describe{
#' \item{aligned}{one row per biological level with the estimated `value`,
#'     its `sigma` from the Fisher information, `lower` and `upper`, on the
#'     common scale.}
#' \item{scaled}{the measurements mapped to the common scale by the inverse
#'     model, with propagated errors.}
#' \item{prediction}{the data with `value` replaced by the model prediction
#'     and `sigma` by the error model, on the original scale.}
#' \item{original}{the data that were fitted.}
#' \item{originalWithParameters}{the data with one column per fitted
#'     parameter.}
#' \item{parameter}{parameter table with columns `name`, `level`,
#'     `parameter`, `value`, `sigma`, `nll` (\eqn{-2 \log L}), `noPars`,
#'     `noData`, `lower`, `upper`, and attributes `value` (total
#'     \eqn{-2 \log L}) and `df` (degrees of freedom). Parameters without
#'     influence on the likelihood are not fitted and reported as `NA`.}
#' \item{multistart}{final objective value, convergence and iterations of each
#'     start per target, see [plotWaterfall()].}
#' \item{biological, scaling, error}{names of the data columns of each effect.}
#' \item{outputScale}{the `outputScale` argument.}
#' }
#'
#' @example inst/examples/exampleAlignReplicates.R
#' @seealso [readWide()], [plotIt()], [llrTest()]
#'
#' @importFrom stats D
#'
#' @export
alignReplicates <- function(data,
                            model = NULL,
                            errorModel = NULL,
                            biological = NULL,
                            scaling = NULL,
                            error = NULL,
                            fitLogscale = TRUE,
                            normalize = TRUE,
                            averageTechRep = FALSE,
                            verbose = FALSE,
                            namesFactored = TRUE,
                            normalizeInput = TRUE,
                            keepOriginalScale = "none",
                            outputScale = "linear",
                            iterlim = 100,
                            optControls = list(),
                            multistart = 1,
                            rootControls = list(),
                            parinit = NULL,
                            multistartSd = 3,
                            cores = 1,
                            likelihood = c("normal", "t"),
                            df = 5
) {

  ## BLAS runs single-threaded during the fit, restored on exit
  prevBlasThreads <- setBlasThreads(1L)
  if (!is.na(prevBlasThreads)) {
    on.exit(setBlasThreads(prevBlasThreads), add = TRUE)
  }

  ## optControls$iterlim takes precedence over the iterlim argument
  optDefaults <- list(
    rinit    = 0.1,
    rmax     = 10,
    parscale = FALSE,
    iterlim  = iterlim,
    fterm    = sqrt(.Machine$double.eps),
    mterm    = sqrt(.Machine$double.eps)
  )
  optControls <- modifyList(optDefaults, optControls)

  multistart <- as.integer(max(multistart, 1L))
  likelihood <- match.arg(likelihood)
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
    stop("'df' must be a positive number.")
  }
  cores <- as.integer(max(cores, 1L))
  if (cores > 1L && .Platform$OS.type == "windows") {
    warning("cores > 1 is not supported on Windows, using cores = 1.")
    cores <- 1L
  }

  ## A previous result is refitted from its data, starting at its estimates
  parameterData <- NULL
  if (is.list(data) && !is.data.frame(data) &&
      all(c("original", "parameter") %in% names(data))) {
    parameterData <- data$parameter
    data <- as.data.frame(data$original)
    data$sigma <- NULL
  }

  inputCheckReport <- inputCheck(
    data = data,
    model = model,
    errorModel = errorModel,
    biological = biological,
    scaling = scaling,
    error = error,
    fitLogscale = fitLogscale,
    outputScale = outputScale
  )

  if (verbose) {
    cat(inputCheckReport, "\n")
  }

  effects <- identifyEffects(
    biological = biological,
    scaling = scaling,
    error = error
  )

  effectsValues <- effects$effectsValues
  effectsPars <- effects$effectsPars
  errorGroups <- effects$errorGroups

  data <- as.data.frame(data)

  toBeScaled <- splitForScaling(
    data,
    effectsValues,
    normalizeInput
  )

  ## Target names are made unique if a target splits into several scaling blocks
  targets <- make.unique(
    vapply(toBeScaled, function(d) as.character(d$name)[1], ""),
    sep = "_"
  )

  for (n in seq_along(targets)) {
    toBeScaled[[n]]$name <- targets[n]
  }

  parameters <- getSymbols(c(model, errorModel), exclude = colnames(data))

  if (normalize) {
    constraint <- paste("1e3 * (mean(", effectsPars[[1]][1], ") - 1)")
    constStrength <- 1000
  } else {
    constraint <- "0"
    constStrength <- 0
  }

  ## Covariates are data columns used in the (error) model
  covariates <- union(
    getSymbols(
      model,
      exclude = c(effectsPars[[1]], effectsPars[[2]], effectsPars[[3]])
    ),
    getSymbols(
      errorModel,
      exclude = c(effectsPars[[1]], effectsPars[[2]], effectsPars[[3]])
    )
  )
  cat("Covariates:", paste(covariates, sep = ", "), "\n")

  if (
    length(
      setdiff(
        c(
          effectsPars[[1]],
          effectsPars[[2]],
          effectsPars[[3]]
        ),
        parameters
      )
    ) > 0
  ) {
    stop("Not all paramters are defined in either arguments
         'scaling', 'biological' or 'error'")
  }

  ## Names are effect classes; they index the parameter vectors in the model
  names(parameters)[parameters %in% effectsPars[[1]]] <- "biological"
  names(parameters)[parameters %in% effectsPars[[2]]] <- "scaling"
  names(parameters)[parameters %in% effectsPars[[3]]] <- "error"
  for (p in names(errorGroups)) {
    names(parameters)[parameters == p] <- paste0("error_", p)
  }

  if (!is.null(parinit)) {
    if (!is.list(parinit) || is.null(names(parinit)) ||
        !all(names(parinit) %in% parameters)) {
      stop("'parinit' must be a list named by model parameters (",
           paste(parameters, collapse = ", "), ").")
    }
    for (p in names(parinit)) {
      v <- parinit[[p]]
      if (!is.numeric(v) || anyNA(v) || (length(v) != 1L && is.null(names(v)))) {
        stop("parinit$", p, " must be a number or a numeric vector named by level.")
      }
      if (fitLogscale && any(v <= 0)) {
        stop("parinit$", p, " must be positive with fitLogscale = TRUE.")
      }
    }
  }

  ## Must run on the model before the substitutions below
  symbolicResult <- symbolicInverse(
    model = model,
    bioPar = effectsPars[[1]][1],
    allSymbols = parameters,
    verbose = verbose
  )

  ## `value` in the error model is the model prediction
  errorModel <- replaceSymbols(
    "value",
    paste0("(", model, ")"),
    errorModel
  )

  if (fitLogscale == TRUE) {
    model <- replaceSymbols(parameters, paste0(
      "exp(", parameters,
      ")"
    ), model)
    errorModel <- replaceSymbols(parameters, paste0(
      "exp(",
      parameters, ")"
    ), errorModel)
    constraint <- replaceSymbols(parameters, paste0(
      "exp(",
      parameters, ")"
    ), constraint)
    if (verbose) {
      cat("model, errormodel and constraint are scaled: x <- exp(x)\n")
    }
  } else if (verbose == TRUE) {
    cat("model, errormodel and constraint remain linear scaled.\n")
  }

  modelDeriv <- deparse(
    stats::D(parse(text = model), name = effectsPars[[1]][1])
  )

  modelJacobian <- lapply(
    parameters,
    function(p) {
      deparse(stats::D(parse(text = model), name = p))
    }
  )
  errorModelJacobian <- lapply(
    parameters,
    function(p) {
      deparse(stats::D(parse(text = errorModel), name = p))
    }
  )

  constrExpr <- parse(text = constraint)

  ## p -> p[<effect class>]
  for (n in seq_along(parameters)) {
    model <- replaceSymbols(parameters[n], paste0(
      parameters[n],
      "[", names(parameters)[n], "]"
    ), model)
    errorModel <- replaceSymbols(parameters[n], paste0(
      parameters[n],
      "[", names(parameters)[n], "]"
    ), errorModel)
    modelDeriv <- replaceSymbols(parameters[n], paste0(
      parameters[n],
      "[", names(parameters)[n], "]"
    ), modelDeriv)
    modelJacobian <- lapply(modelJacobian, function(myjac) {
      replaceSymbols(
        parameters[n],
        paste0(
          parameters[n], "[", names(parameters)[n],
          "]"
        ), myjac
      )
    })
    errorModelJacobian <- lapply(errorModelJacobian, function(myjac) {
      replaceSymbols(
        parameters[n],
        paste0(
          parameters[n], "[", names(parameters)[n],
          "]"
        ), myjac
      )
    })
  }

  cat("Model:        ", model, "\n", sep = "")
  cat("Error Model:  ", errorModel, "\n", sep = "")

  modelExpr <- parse(text = model)
  modelDerivExpr <- parse(text = modelDeriv)
  errorModelExpr <- parse(text = errorModel)

  modelJacobianExpr <- lapply(
    modelJacobian,
    function(myjac) parse(text = myjac)
  )
  errorModelJacobianExpr <- lapply(
    errorModelJacobian,
    function(myjac) parse(text = myjac)
  )

  analyticInverse <- prepareAnalyticInverse(
    symbolicResult = symbolicResult,
    parameters = parameters,
    fitLogscale = fitLogscale
  )

  if (!is.null(analyticInverse)) {
    cat("Analytic inverse:       ", analyticInverse$inverseStr, "\n")
    cat("Analytic inverse deriv: ", analyticInverse$inverseDerivStr, "\n")
  } else {
    cat("Analytic inverse: not available, using numerical fallback.\n")
  }

  passParList <- list(
    effectsValues = effectsValues,
    parameterData = parameterData,
    averageTechRep = averageTechRep,
    verbose = verbose,
    covariates = covariates,
    parameters = parameters,
    fitLogscale = fitLogscale,
    effectsPars = effectsPars,
    modelExpr = modelExpr,
    errorModelExpr = errorModelExpr,
    constrExpr = constrExpr,
    modelJacobianExpr = modelJacobianExpr,
    errorModelJacobianExpr = errorModelJacobianExpr,
    constStrength = constStrength,
    normalize = normalize,
    modelDerivExpr = modelDerivExpr,
    outputScale = outputScale,
    optControls = optControls,
    multistart = multistart,
    multistartSd = multistartSd,
    likelihood = likelihood,
    df = df,
    cores = cores,
    parinit = parinit,
    errorGroups = errorGroups,
    rootControls = rootControls,
    analyticInverseExpr = if (!is.null(analyticInverse)) analyticInverse$inverseExpr else NULL,
    analyticInverseDerivExpr = if (!is.null(analyticInverse)) analyticInverse$inverseDerivExpr else NULL
  )

  out <- lapply(seq_along(toBeScaled), function(i) {
    try({
      cat("Target ", i, "/", length(toBeScaled), ":", targets[i], "\n", sep = "")
      scaleTarget(currentData = toBeScaled[[i]], passParList = passParList)
    }, silent = FALSE)
  })

  ## Failed targets are dropped from all outputs
  ok <- Filter(function(o) !inherits(o, "try-error"), out)
  parTable <- do.call(rbind, lapply(ok, function(o) o[[6]]))
  attr(parTable, "value") <- do.call(sum, lapply(ok, function(o) attr(o[[6]], "value")))
  attr(parTable, "df") <- do.call(sum, lapply(ok, function(o) attr(o[[6]], "df")))

  outCombined <- lapply(seq_len(5), function(i) do.call(rbind, lapply(ok, function(o) o[[i]])))
  names(outCombined) <- c(
    "prediction",
    "scaled",
    "aligned",
    "original",
    "originalWithParameters"
  )

  if (keepOriginalScale != "none" && normalizeInput) {
    warning("To keep the original scale, use normalizeInput = FALSE")
  }
  if (keepOriginalScale %in% c("lin", "log")) {
    cols <- c("value", "sigma", "lower", "upper")
    if (keepOriginalScale == "log") cols <- c("value", "lower", "upper")
    for (target in unique(outCombined$original$name)) {
      meanOriginal <- mean(subset(outCombined$original, name == target)$value)
      for (tab in c("scaled", "aligned", "prediction")) {
        rows <- which(outCombined[[tab]]$name == target)
        meanTab <- mean(subset(outCombined[[tab]], name == target)$value)
        for (col in cols) {
          x <- outCombined[[tab]][[col]][rows]
          outCombined[[tab]][[col]][rows] <- if (keepOriginalScale == "lin") {
            x * meanOriginal / meanTab
          } else {
            x + meanOriginal - meanTab
          }
        }
      }
    }
  }

  returnList <- list(
    aligned = outCombined$aligned,
    scaled = outCombined$scaled,
    prediction = outCombined$prediction,
    original = outCombined$original,
    originalWithParameters = outCombined$originalWithParameters,
    parameter = parTable,
    multistart = do.call(rbind, lapply(ok, function(o) o$multistart)),
    biological = effectsValues[[1]],
    scaling = effectsValues[[2]],
    error = effectsValues[[3]],
    outputScale = outputScale,
    likelihood = likelihood,
    likelihoodDf = if (likelihood == "t") df else NULL
  )

  if (namesFactored == TRUE) {
    for (tab in c("aligned", "scaled", "prediction", "original", "originalWithParameters")) {
      returnList[[tab]]$name <- factor(
        returnList[[tab]]$name,
        levels = unique(returnList[[tab]]$name)
      )
    }
  }

  return(returnList)
}
