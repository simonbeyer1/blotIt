# getSymbols() -----------------------------------------------------------

#' Symbols in expression strings
#'
#' @param char character vector of expressions; entries `"0"` are ignored.
#' @param exclude character vector of symbols to drop.
#'
#' @return character vector of unique symbol names.
#'
#' @importFrom utils getParseData
#'
#' @noRd
getSymbols <- function(char, exclude = NULL) {
  char <- char[char != "0"]
  out <- parse(text = char, keep.source = TRUE)
  out <- utils::getParseData(out)
  names <- unique(out$text[out$token == "SYMBOL"])
  if (!is.null(exclude)) {
    names <- names[!names %in% exclude]
  }
  return(names)
}


# identifyEffects() ------------------------------------------------------

#' Parameter and column names of the effect formulas
#'
#' @param biological,scaling two-sided formulas.
#' @param error a two-sided formula or a list of them.
#'
#' @return list with `effectsValues`, `effectsPars` and `errorGroups` (named
#'   list error parameter -> data columns, `NULL` for a single error formula)
#'
#' @noRd
identifyEffects <- function(biological = NULL, scaling = NULL, error = NULL) {

  biologicaValues <- getSymbols(as.character(biological)[3])
  scalingValues <- getSymbols(as.character(scaling)[3])
  biologicalPars <- getSymbols(as.character(biological)[2])
  scalingPars <- getSymbols(as.character(scaling)[2])

  errorGroups <- NULL
  if (is.list(error)) {
    errorGroups <- do.call(c, lapply(error, function(f) {
      lhs <- getSymbols(as.character(f)[2])
      stats::setNames(rep(list(getSymbols(as.character(f)[3])), length(lhs)), lhs)
    }))
    if (anyDuplicated(names(errorGroups))) {
      stop("Each error parameter may appear in only one 'error' formula.")
    }
    errorValues <- unique(unlist(errorGroups, use.names = FALSE))
    errorPars <- names(errorGroups)
  } else {
    errorValues <- getSymbols(as.character(error)[3])
    errorPars <- getSymbols(as.character(error)[2])
  }

  effectsValues <- list(
    biologicaValues = biologicaValues,
    scalingValues = scalingValues,
    errorValues = errorValues
  )
  effectsPars <- list(
    biologicalPars = biologicalPars,
    scalingPars = scalingPars,
    errorPars = errorPars
  )

  return(list(
    effectsValues = effectsValues,
    effectsPars = effectsPars,
    errorGroups = errorGroups
  ))
}


# replaceSymbols ----------------------------------------------------------

#' Replace symbols in expression strings
#'
#' @param what character vector of symbols to replace.
#' @param by replacement strings, recycled to the length of `what`.
#' @param x character vector of expressions; entries `"0"` stay unchanged.
#'
#' @return `x` with every symbol in `what` replaced.
#'
#' @importFrom utils getParseData
#'
#' @noRd
replaceSymbols <- function(what, by, x) {
  xOrig <- x
  notZero <- which(x != "0")
  x <- x[notZero]
  useNames <- names(x)
  xParsed <- parse(text = x, keep.source = TRUE)
  data <- utils::getParseData(xParsed)
  by <- rep(by, length.out = length(what))
  names(by) <- what
  data$text[data$text %in% what] <- by[data$text[data$text %in% what]]
  data <- data[data$token != "expr", ]
  breaks <- c(0, which(diff(data$line1) == 1), length(data$line1))
  out <- lapply(
    seq_len((length(breaks) - 1)),
    function(i) {
      paste(
        data$text[seq((breaks[i] + 1), (breaks[i + 1]))],
        collapse = ""
      )
    }
  )
  names(out) <- useNames
  out <- unlist(out)
  xOrig[notZero] <- out
  return(xOrig)
}


# paste_() ----------------------------------------------------------------

#' `paste()` with separator `"_"`
#'
#' @noRd
paste_ <- function(...) paste(..., sep = "_")


# analyzeBlocks() --------------------------------------------------------

#' Groups of biological levels linked by shared scaling levels
#'
#' @param blockMatrix 0/1 matrix, rows are biological levels, columns are
#'   scaling levels, 1 where a level combination was measured.
#'
#' @return list with one integer vector of row indices per connected block,
#'   rows ascending, blocks ordered by their smallest row.
#'
#' @noRd
analyzeBlocks <- function(blockMatrix) {
  out <- which(rowSums(blockMatrix) == 0)
  if (length(out) > 0) {
    blockMatrix <- blockMatrix[-out, , drop = FALSE]
    cat("matrix contains zero rows which have been eliminated\n")
  }

  nR <- nrow(blockMatrix)
  nC <- ncol(blockMatrix)

  ## Union-find over the bipartite graph: nodes 1..nR are rows, nR + 1..nR + nC
  ## are columns, each non-zero entry is an edge.
  edges <- which(blockMatrix == 1, arr.ind = TRUE)
  parent <- seq_len(nR + nC)

  findRoot <- function(x) {
    root <- x
    while (parent[root] != root) root <- parent[root]
    while (parent[x] != root) {
      nxt <- parent[x]
      parent[x] <<- root
      x <- nxt
    }
    root
  }

  for (e in seq_len(nrow(edges))) {
    ra <- findRoot(edges[e, 1L])
    rb <- findRoot(edges[e, 2L] + nR)
    if (ra != rb) parent[ra] <- rb
  }

  rowRoots <- vapply(seq_len(nR), findRoot, integer(1))
  rowOrder <- order(rowRoots, seq_len(nR))
  splitRoots <- rowRoots[rowOrder]
  rComponents <- split(rowOrder, splitRoots)
  rComponents <- rComponents[order(vapply(rComponents, min, integer(1)))]
  names(rComponents) <- NULL

  return(rComponents)
}


# splitForScaling()  ----------------------------------------------------

#' Split data into independently fitted blocks
#'
#' @param data long-format data frame.
#' @param effectsValues list of the biological, scaling and error column names.
#' @param normalizeInput logical; if `TRUE`, each block is divided by its mean
#'   (positive data) or centered by its mean (otherwise).
#'
#' @return list of data frames, one per block.
#'
#' @noRd
splitForScaling <- function(data,
                            effectsValues,
                            normalizeInput) {
  if (!"1" %in% colnames(data)) {
    data["1"] <- 1
    intercept <- FALSE
  } else {
    intercept <- TRUE
  }

  scalingStrings <- Reduce(paste_, data[effectsValues[[2]]])
  biologicalStrings <- Reduce(paste_, data[effectsValues[[1]]])

  scalingStringsUnique <- unique(scalingStrings)
  biologicalStringsUnique <- unique(biologicalStrings)

  blockMatrix <- matrix(
    0,
    ncol = length(scalingStringsUnique),
    nrow = length(biologicalStringsUnique)
  )
  rowIdx <- match(biologicalStrings, biologicalStringsUnique)
  colIdx <- match(scalingStrings, scalingStringsUnique)
  blockMatrix[cbind(rowIdx, colIdx)] <- 1

  scalingsList <- analyzeBlocks(blockMatrix)

  if (!intercept) {
    data <- data[, -which(colnames(data) == "1")]
  }

  outList <- lapply(
    scalingsList,
    function(l) {
      data[biologicalStrings %in% biologicalStringsUnique[l], ]
    }
  )

  if (!normalizeInput) {
    return(outList)
  }
  for (i in seq_len(length(outList))) {
    v <- outList[[i]]$value
    ## Blocks with non-positive values are centered, since their mean can be ~0
    if (all(v > 0)) {
      outList[[i]]$value <- v / mean(v)
    } else {
      outList[[i]]$value <- v - mean(v)
    }
  }

  return(outList)
}


# inputCheck() -----------------------------------------------------------

#' Check the arguments of `alignReplicates()` for structural errors
#'
#' @param data,model,errorModel,biological,scaling,error,fitLogscale,outputScale
#'   as in `alignReplicates()`.
#'
#' @return `"All input checks passed."`, or an error.
#'
#' @noRd
inputCheck <- function(data = NULL,
                       model = NULL,
                       errorModel = NULL,
                       biological = NULL,
                       scaling = NULL,
                       error = NULL,
                       fitLogscale = NULL,
                       outputScale = NULL) {
  if (is.null(model) | is.null(errorModel) | is.null(biological) |
      is.null(scaling) | is.null(error)) {
    stop(
      "All of model, errorModel, biological, scaling, error ",
      "must be set."
    )
  }

  errorList <- if (is.list(error)) error else list(error)
  if (length(errorList) == 0) {
    stop("'error' must be a formula or a non-empty list of formulas.")
  }

  columnNames <- names(data)
  expectedNames <- unique(c(
    getSymbols(as.character(biological)[3]),
    getSymbols(as.character(scaling)[3]),
    unlist(lapply(errorList, function(f) getSymbols(as.character(f)[3])))
  ))

  if (!all(expectedNames %in% columnNames)) {
    stop(
      "Not all column names set in 'biological', 'scaling' and ",
      "'error' are present in 'data'."
    )
  }

  if (!is.logical(fitLogscale)) {
    stop(
      "'fitLogscale' must be logical"
    )
  }

  if (length(as.character(biological)) == 1 |
      length(as.character(scaling)) == 1 |
      any(vapply(errorList, function(f) length(as.character(f)) == 1, TRUE))) {
    stop("Do not pass biological, scaling or error as string.")
  }
  if (length(as.character(biological)) < 3) {
    stop("Left and right-hand side of formula 'biological' is needed")
  }
  if (length(as.character(scaling)) < 3) {
    stop("Left and right-hand side of formula 'scaling' is needed")
  }
  if (any(vapply(errorList, function(f) length(as.character(f)) < 3, TRUE))) {
    stop("Left and right-hand side of formula 'error' is needed")
  }

  returnMsg <- "All input checks passed."
  return(returnMsg)
}


# structuredDiagInverse() ------------------------------------------------

#' Diagonal of the inverse of an arrow plus rank-1 Hessian
#'
#' Exact \eqn{\mathrm{diag}(H^{-1})} for \eqn{H = H_{GN} + 2 v v^T} without
#' forming a dense matrix.
#'
#' @param Hgn sparse symmetric Gauss-Newton Hessian.
#' @param constrVec rank-1 constraint vector \eqn{v}, or `NULL`.
#' @param bioIdx,borderIdx indices of the diagonal biological block and of the
#'   remaining (border) parameters.
#'
#' @return numeric vector, or `NULL` if the biological block is not diagonal.
#'
#' @noRd
structuredDiagInverse <- function(Hgn, constrVec, bioIdx, borderIdx) {
  n <- nrow(Hgn)
  Hbb <- Hgn[bioIdx, bioIdx, drop = FALSE]
  if (!Matrix::isDiagonal(Hbb)) {
    return(NULL)
  }

  Db <- as.numeric(Matrix::diag(Hgn))[bioIdx]
  ## Floor zero or negative curvature at a small positive value
  posScale <- max(Db[is.finite(Db)], na.rm = TRUE)
  floorVal <- sqrt(.Machine$double.eps) * max(posScale, 1)
  Db <- pmax(Db, floorVal)
  DbInv <- 1 / Db

  vb <- if (is.null(constrVec)) numeric(length(bioIdx)) else constrVec[bioIdx]

  ## Sherman-Morrison for P = diag(Db) + 2 vb vb'
  z <- DbInv * vb
  denom <- 1 + 2 * sum(vb * z)
  diagPinv <- DbInv - 2 * z^2 / denom

  m <- length(borderIdx)
  if (m == 0) {
    out <- numeric(n)
    out[bioIdx] <- diagPinv
    return(out)
  }

  C <- Hgn[bioIdx, borderIdx, drop = FALSE]                 # n_b x m (sparse)
  Mb <- as.matrix(Hgn[borderIdx, borderIdx, drop = FALSE])  # m x m (dense)

  ## P^{-1} C = diag(DbInv) C - (2/denom) z (z' C)
  zC <- as.numeric(Matrix::crossprod(z, C))                 # length m
  W <- as.matrix(Matrix::Diagonal(x = DbInv) %*% C) -
    (2 / denom) * outer(z, zC)                              # n_b x m (dense)

  ## Schur complement Sc = M - C' P^{-1} C, ridge-regularized if singular
  Sc <- Mb - as.matrix(Matrix::crossprod(C, W))
  ScInv <- tryCatch(solve(Sc), error = function(e) {
    ridge <- sqrt(.Machine$double.eps) * max(abs(diag(Sc)), 1)
    solve(Sc + diag(ridge, m))
  })

  WS <- W %*% ScInv
  diagBio <- diagPinv + rowSums(WS * W)

  out <- numeric(n)
  out[bioIdx] <- diagBio
  out[borderIdx] <- diag(ScInv)
  out
}


# scaleTarget() ----------------------------------------------------------

#' Fit one block of data and build the output tables
#'
#' @param currentData data frame of one block, see `splitForScaling()`.
#' @param passParList list of settings and parsed expressions assembled by
#'   `alignReplicates()`.
#'
#' @return list with `outPrediction`, `outScaled`, `outAligned`, `original`,
#'   `outOrigParameters`, `parTable` and `multistart`, see the return value of
#'   `alignReplicates()`.
#'
#' @importFrom nleqslv nleqslv
#'
#' @noRd
scaleTarget <- function(currentData,
                        passParList) {
  effectsValues <- passParList$effectsValues
  averageTechRep <- passParList$averageTechRep
  verbose <- passParList$verbose
  covariates <- passParList$covariates
  parameters <- passParList$parameters
  fitLogscale <- passParList$fitLogscale
  effectsPars <- passParList$effectsPars
  normalize <- passParList$normalize
  outputScale <-  passParList$outputScale
  optControls <- passParList$optControls
  multistart <- passParList$multistart
  multistartSd <- passParList$multistartSd
  cores <- passParList$cores
  rootControls <- passParList$rootControls
  errorGroups <- passParList$errorGroups

  currentName <- unique(currentData$name)
  userHasOne <- "1" %in% colnames(currentData)
  currentData$sigma <- NaN
  currentData[["1"]] <- "1"

  dataFitBiological <- do.call(
    paste_, currentData[, effectsValues[[1]], drop = FALSE]
  )
  dataFitScaling <- do.call(
    paste_, currentData[, effectsValues[[2]], drop = FALSE]
  )

  if (averageTechRep) {
    if (verbose) {
      cat("Analyzing technical replicates ... ")
    }
    groups <- interaction(dataFitBiological, dataFitScaling)
    if (any(duplicated(groups))) {
      currentData <- do.call(rbind, lapply(
        unique(groups),
        function(g) {
          subdata <- currentData[groups == g, ]
          outdata <- subdata[1, ]
          outdata$value <- mean(subdata$value)
          return(outdata)
        }
      ))
      cat(
        "data points that could not be distinguished by either ",
        "biological\nor scaling variables have been averaged.\n"
      )
    } else if (verbose) {
      cat("none found.\n")
    }
  }

  dataFit <- data.frame(
    currentData[
      ,
      union(c("name", "time", "value", "sigma"), covariates)
    ],
    biological = do.call(
      paste_,
      currentData[, effectsValues[[1]], drop = FALSE]
    ),
    scaling = do.call(
      paste_,
      currentData[, effectsValues[[2]], drop = FALSE]
    ),
    error = do.call(
      paste_,
      currentData[, effectsValues[[3]], drop = FALSE]
    ),
    stringsAsFactors = FALSE
  )
  for (p in names(errorGroups)) {
    dataFit[[paste0("error_", p)]] <- do.call(
      paste_, currentData[, errorGroups[[p]], drop = FALSE]
    )
  }
  tGroupIdx <- NULL
  if (!is.null(passParList$tGroup)) {
    tGroupLevels <- do.call(paste_, currentData[, passParList$tGroup, drop = FALSE])
    tGroupIdx <- match(tGroupLevels, unique(tGroupLevels))
  }

  ## Levels per effect class, named like the class columns of dataFit
  levelsList <- lapply(
    stats::setNames(nm = unique(c("biological", "scaling", "error", names(parameters)))),
    function(e) unique(as.character(dataFit[[e]]))
  )

  initPars <- genInitPars(
    parameters,
    fitLogscale,
    levelsList
  )
  parLevels <- unlist(lapply(names(parameters), function(e) levelsList[[e]]),
                      use.names = FALSE)
  parSymbols <- rep(unname(parameters), lengths(levelsList[names(parameters)]))
  initPars <- applyPreviousFit(
    initPars, passParList$parameterData, currentName, parSymbols, parLevels,
    fitLogscale
  )
  initPars <- applyParinit(
    initPars, passParList$parinit, parSymbols, parLevels, fitLogscale
  )

  maskIdx <- buildMaskIndex(
    parameters,
    levelsList,
    dataFit
  )

  if (verbose) {
    cat("Starting fit\n")
  }

  passParList2 <- list(
    dataFit = dataFit,
    levelsList = levelsList,
    effectsPars = effectsPars,
    maskIdx = maskIdx,
    initPars = initPars,
    tGroupIdx = tGroupIdx,
    tGroupSize = if (!is.null(tGroupIdx)) tabulate(tGroupIdx)
  )

  # * call of trust() function ----------------------------------------------

  ## Above the threshold a dense Hessian is too large: use the sparse optimizer
  nPars <- length(initPars)
  structuredThreshold <- if (!is.null(optControls$structuredThreshold)) {
    optControls$structuredThreshold
  } else {
    5000L
  }
  useStructured <- nPars > structuredThreshold
  passParList$sparseHessian <- useStructured
  if (useStructured && verbose) {
    cat("Using sparse structured optimizer (", nPars, " parameters).\n", sep = "")
  }

  ## Parameters without influence on the likelihood stay at their start value
  ## and are reported as NA.
  inactive <- findInactivePars(initPars, passParList, passParList2)
  freeIdx <- seq_along(initPars)
  fitInit <- initPars
  if (length(inactive)) {
    freeIdx <- setdiff(freeIdx, inactive)
    fitInit <- initPars[freeIdx]
    passParList2$fixedPars <- initPars[inactive]
    passParList2$freeIdx <- freeIdx
    if (verbose) {
      cat("Fixed without influence on the likelihood:",
          names(initPars)[inactive], "\n")
    }
  }

  trustArgs <- list(
    objfun    = objFunction,
    parinit   = initPars,
    rinit     = optControls$rinit,
    rmax      = optControls$rmax,
    iterlim   = optControls$iterlim,
    fterm     = optControls$fterm,
    mterm     = optControls$mterm,
    blather   = verbose,
    structured = useStructured,
    passParList  = passParList,
    passParList2 = passParList2
  )
  if (!isFALSE(optControls$parscale) && !is.null(optControls$parscale)) {
    trustArgs$parscale <- optControls$parscale
  }

  ## Start 1 is the initial point, the others perturb it. Perturbations are
  ## drawn here, before any fit, so results do not depend on `cores`.
  starts <- c(
    list(fitInit),
    lapply(seq_len(multistart - 1L), function(i) {
      fitInit + rnorm(length(fitInit), mean = 0, sd = multistartSd)
    })
  )
  runStart <- function(parinit) {
    args <- trustArgs
    args$parinit <- parinit
    tryCatch(do.call(trust, args), error = function(e) e)
  }
  fits <- if (cores > 1L && multistart > 1L) {
    parallel::mclapply(starts, runStart, mc.cores = cores)
  } else {
    lapply(starts, runStart)
  }
  if (inherits(fits[[1]], "error")) stop(fits[[1]])
  isFit <- vapply(fits, function(f) is.list(f) && !is.null(f$argument), TRUE)

  fitResult <- fits[[1]]
  for (f in fits[-1][isFit[-1]]) {
    if (f$converged && f$value < fitResult$value) fitResult <- f
  }

  multistartTable <- data.frame(
    name = currentName,
    start = seq_along(fits),
    value = vapply(seq_along(fits), function(i) {
      if (isFit[i]) fits[[i]]$value else NA_real_
    }, 0),
    converged = vapply(seq_along(fits), function(i) {
      isFit[i] && isTRUE(fits[[i]]$converged)
    }, TRUE),
    iterations = vapply(seq_along(fits), function(i) {
      if (isFit[i]) as.integer(fits[[i]]$iterations) else NA_integer_
    }, 0L),
    stringsAsFactors = FALSE
  )

  if (!fitResult$converged) {
    warning(paste("Non-converged fit for target", currentName))
  } else if (verbose == TRUE) {
    cat("Fit converged.\n")
  }

  residualsFit <- resolveFunction(
    currentPars = fitResult$argument,
    passParList = passParList,
    passParList2 = passParList2,
    calcDeriv = FALSE
  )

  if (identical(passParList$likelihood, "t")) {
    info <- objFunction(
      fitResult$argument, passParList,
      c(passParList2, list(expectedInfo = TRUE))
    )
    fitResult$hessian <- info$hessian
    fitResult$constrVec <- info$constrVec
  }

  nFree <- length(fitResult$argument)
  bessel <- sqrt(
    nrow(dataFit) / (nrow(dataFit) - nFree +
                       normalize)
  )

  ## Dense path only, an SVD is infeasible at structured-path sizes
  if (!isTRUE(passParList$sparseHessian)) {
    singleValues <- svd(fitResult[["hessian"]])[["d"]]
    tol <- sqrt(.Machine$double.eps)
    nonIdentifiable <- which(singleValues < tol * singleValues[1])
    if (length(nonIdentifiable) > 0 && verbose) {
      cat("Note: ", length(nonIdentifiable),
          " eigenvalue(s) of Hessian below tolerance.\n")
    }
  }

  ## Diagonal of twice the inverse Fisher information, for all parameter sigmas
  if (isTRUE(passParList$sparseHessian)) {
    effects <- names(parameters)
    nlevPar <- vapply(effects, function(e) length(levelsList[[e]]), integer(1))
    offs <- cumsum(c(0L, nlevPar))[seq_along(parameters)]
    bioIdx <- unlist(lapply(
      which(effects == "biological"),
      function(k) as.integer(offs[k]) + seq_len(nlevPar[k])
    ))
    bioIdx <- stats::na.omit(match(bioIdx, freeIdx))
    borderIdx <- setdiff(seq_along(freeIdx), bioIdx)
    diagInv <- structuredDiagInverse(
      Hgn = fitResult$hessian,
      constrVec = fitResult$constrVec,
      bioIdx = bioIdx,
      borderIdx = borderIdx
    )
    if (is.null(diagInv)) {
      stop("Structured diagonal-of-inverse failed: the biological Hessian ",
           "block is not diagonal (model has multiple biological parameters?). ",
           "The large-parameter path supports a single biological parameter.")
    }
    diag2Inv <- 2 * diagInv
  } else {
    diag2Inv <- diag(2 * MASS::ginv(fitResult$hessian))
  }
  diag2InvFull <- rep(NA_real_, length(initPars))
  diag2InvFull[freeIdx] <- diag2Inv
  parValues <- rep(NA_real_, length(initPars))
  parValues[freeIdx] <- fitResult$argument
  names(parValues) <- names(initPars)

  # * parameter table -------------------------------------------------------

  parTable <- data.frame(
    name = currentName,
    level = parLevels,
    parameter = names(parValues),
    value = parValues,
    nll = fitResult$value,
    noPars = nFree - normalize,
    noData = nrow(dataFit)
  )

  parTable$sigma <- sqrt(diag2InvFull) * bessel

  if (fitLogscale == TRUE) {
    parTable$value <- exp(parTable$value)
    parTable$sigma <- parTable$value * parTable$sigma
  }
  parTable$lower <- parTable$value - parTable$sigma
  parTable$upper <- parTable$value + parTable$sigma

  if (verbose) {
    cat("Estimated parameters on non-log scale:\n")
    print(parTable)
    cat(
      "converged:", fitResult$converged, ", iterations:",
      fitResult$iterations, "\n"
    )
    cat("-2*LL: ", fitResult$value, "on", nrow(currentData) +
          normalize - nFree, "degrees of freedom\n")
  }

  attr(parTable, "value") <- fitResult$value
  attr(parTable, "df") <- nrow(currentData) + normalize - nFree

  # * outPrediction --------------------------------------------------------

  outPrediction <- currentData
  outPrediction$value <- fitResult$prediction
  outPrediction$sigma <- fitResult$sigma * bessel
  outPrediction$lower <- outPrediction$value - outPrediction$sigma
  outPrediction$upper <- outPrediction$value + outPrediction$sigma

  # * outScaled ------------------------------------------------------------

  ## Analytic inverse if available, otherwise numeric inversion with nleqslv
  initialValuesScaled <- rep(0, nrow(currentData))
  if (verbose) {
    cat("Inverting model ... ")
  }

  analyticResult <- evaluateAnalyticInverse(
    parList = residualsFit[-1],
    passParList = passParList,
    passParList2 = passParList2
  )

  if (!is.null(analyticResult)) {
    valuesScaled <- analyticResult$values
    myDeriv <- analyticResult$jacobian
    if (verbose) cat("done (analytic)\n")

  } else {
    if (verbose && !is.null(passParList$analyticInverseExpr)) {
      cat("analytic eval failed, falling back to nleqslv ... ")
    }

    n_data <- length(initialValuesScaled)

    ## nleqslv expects the Jacobian as a matrix
    jacWrapper <- function(x, ...) {
      dv <- as.numeric(evaluateModelJacobian(x, ...))
      n <- length(dv)
      J <- matrix(0, n, n)
      diag(J) <- dv
      J
    }

    reservedNames <- c("x", "fn", "jac", "parList", "passParList",
                       "passParList2")
    safeRootControls <- rootControls[
      !names(rootControls) %in% reservedNames
    ]
    if (is.null(safeRootControls$control)) {
      safeRootControls$control <- list(allowSingular = TRUE)
    } else if (is.null(safeRootControls$control$allowSingular)) {
      safeRootControls$control$allowSingular <- TRUE
    }

    rootArgs <- c(
      list(
        x = initialValuesScaled,
        fn = evaluateModel,
        jac = jacWrapper,
        parList = residualsFit[-1],
        passParList = passParList,
        passParList2 = passParList2
      ),
      safeRootControls
    )

    nleqslvResult <- try(
      do.call(nleqslv::nleqslv, rootArgs),
      silent = TRUE
    )

    nleqslvConverged <- !inherits(nleqslvResult, "try-error") &&
      nleqslvResult$termcd %in% c(1, 2)

    if (nleqslvConverged) {
      valuesScaled <- nleqslvResult$x
      myDeriv <- abs(
        evaluateModelJacobian(
          valuesScaled,
          residualsFit[-1],
          passParList = passParList,
          passParList2 = passParList2
        )
      )
    } else {
      if (!inherits(nleqslvResult, "try-error")) {
        warning("nleqslv did not converge for target ",
                currentName, " (termcd: ", nleqslvResult$termcd,
                " - ", nleqslvResult$message, ")")
      }
      valuesScaled <- try(stop("inversion failed"), silent = TRUE)
    }

    if (verbose) cat("done (numeric)\n")
  }

  if (inherits(valuesScaled, "try-error")) {
    outScaled <- NULL
    warning("Rescaling to common scale not possible for target ",
            currentName, ". Equations not invertible.")
  } else {
    sigmasScaled <- fitResult$sigma * bessel / myDeriv
    if (fitLogscale == TRUE) {
      valuesScaled <- exp(valuesScaled)
      sigmasScaled <- valuesScaled * sigmasScaled
    }
    upperScaled <- valuesScaled + sigmasScaled
    lowerScaled <- valuesScaled - sigmasScaled

    outScaled <- currentData
    outScaled$value <- valuesScaled
    outScaled$sigma <- sigmasScaled
    outScaled$upper <- upperScaled
    outScaled$lower <- lowerScaled
  }

  # * outAligned -----------------------------------------------------------

  ## One row per biological level
  outAligned <- currentData[
    !duplicated(dataFit$biological),
    intersect(effectsValues[[1]], colnames(currentData))
  ]

  outAligned$value <- residualsFit[[effectsPars[[1]][1]]]
  outAligned$parameter <- c(
    do.call(
      rbind,
      lapply(
        seq_len(nrow(outAligned)),
        function (i) {
          out <- parTable[
            which(parTable$level == do.call(
              paste_,
              c(outAligned[i,effectsValues[[1]]])
            )
            ),
          ]$parameter
        }
      )
    )
  )

  outAligned$sigma <- sqrt(diag2InvFull)[parSymbols == effectsPars[[1]][1]] *
    bessel
  if (fitLogscale == TRUE) {
    outAligned$value <- exp(outAligned$value)
    outAligned$sigma <- outAligned$value * outAligned$sigma
  }
  outAligned$lower <- outAligned$value - outAligned$sigma
  outAligned$upper <- outAligned$value + outAligned$sigma

  # * outOrigParameters -------------------------------------------------

  outOrigParameters <- currentData
  for (k in seq_along(parameters)) {
    effect <- names(parameters)[k]
    useLevels <- as.character(dataFit[[effect]])
    index0 <- which(parSymbols == parameters[k])
    index1 <- match(useLevels, as.character(parTable$level[index0]))
    index <- index0[index1]
    outOrigParameters[[parameters[k]]] <- parTable$value[index]
  }
  out <- list(
    outPrediction = outPrediction,
    outScaled = outScaled,
    outAligned = outAligned,
    original = currentData,
    outOrigParameters = outOrigParameters
  )
  if (!userHasOne && !"1" %in% unlist(effectsValues)) {
    out <- lapply(out, function(tab) {
      if (!is.null(tab)) tab[["1"]] <- NULL
      tab
    })
  }
  out$parTable <- parTable
  out$multistart <- multistartTable

  return(out)
}


# genInitPars() -------------------------------------------------

#' Default initial parameters
#'
#' @param parameters named character vector of model symbols, names are the
#'   effect classes.
#' @param fitLogscale logical; start value is 0 if `TRUE`, else 1.
#' @param levelsList named list of levels per effect class.
#'
#' @return numeric vector with one entry per symbol and level, named
#'   `<symbol>_<index>` with a running index.
#'
#' @noRd
genInitPars <- function(parameters,
                        fitLogscale,
                        levelsList) {
  initPars <- do.call(
    "c",
    lapply(
      seq_along(parameters),
      function(n) {
        v <- if (fitLogscale == TRUE) 0 else 1
        l <- length(levelsList[[names(parameters)[n]]])
        p <- as.character(parameters[n])
        structure(rep(v, l), names = rep(p, l))
      }
    )
  )

  names(initPars) <- paste(
    names(initPars),
    seq_along(initPars),
    sep = "_"
  )

  return(initPars)
}


# applyParinit() ---------------------------------------------------------

#' Overwrite default initial values with user-supplied ones
#'
#' @param initPars output of `genInitPars()`.
#' @param parinit `NULL` or named list, see `parinit` in [alignReplicates()].
#' @param parSymbols,parLevels model symbol and level of each entry of
#'   `initPars`.
#' @param fitLogscale logical; values are log-transformed if `TRUE`.
#'
#' @return `initPars` with the supplied entries replaced.
#'
#' @noRd
applyParinit <- function(initPars, parinit, parSymbols, parLevels, fitLogscale) {
  for (p in names(parinit)) {
    idx <- which(parSymbols == p)
    v <- parinit[[p]]
    if (!is.null(names(v))) {
      idx <- idx[parLevels[idx] %in% names(v)]
      v <- v[parLevels[idx]]
    }
    initPars[idx] <- if (fitLogscale) log(v) else v
  }
  initPars
}


# applyPreviousFit() -----------------------------------------------------

#' Start values from the parameter table of a previous fit
#'
#' @param prev `NULL` or the `parameter` table of an [alignReplicates()] result.
#' @inheritParams applyParinit
#' @param currentName name of the current target.
#'
#' @return `initPars` with the entries found in `prev` replaced.
#'
#' @noRd
applyPreviousFit <- function(initPars, prev, currentName, parSymbols,
                             parLevels, fitLogscale) {
  if (is.null(prev)) return(initPars)
  prev <- prev[as.character(prev$name) == currentName & is.finite(prev$value), ]
  if (fitLogscale) prev <- prev[prev$value > 0, ]
  key <- paste(sub("_[0-9]+$", "", prev$parameter), prev$level, sep = "\r")
  idx <- match(paste(parSymbols, parLevels, sep = "\r"), key)
  hit <- !is.na(idx)
  v <- prev$value[idx[hit]]
  initPars[hit] <- if (fitLogscale) log(v) else v
  initPars
}


# findInactivePars() -----------------------------------------------------

#' Indices of parameters without influence on the likelihood
#'
#' A parameter is inactive if its gradient entry and Hessian row vanish at
#' `initPars` and at a shifted point.
#'
#' @return integer vector of indices into `initPars`.
#'
#' @noRd
findInactivePars <- function(initPars, passParList, passParList2) {
  inactive <- rep(TRUE, length(initPars))
  for (x in list(initPars, initPars + 0.1)) {
    o <- objFunction(x, passParList, passParList2)
    if (is.null(o$gradient)) return(integer(0))
    hessRowZero <- Matrix::colSums(abs(o$hessian)) == 0
    if (!is.null(o$constrVec)) hessRowZero <- hessRowZero & o$constrVec == 0
    inactive <- inactive & o$gradient == 0 & hessRowZero
  }
  which(inactive)
}


# genMask() ---------------------------------------------------------

#' Dense 0/1 activation masks, one per parameter
#'
#' Superseded by `buildMaskIndex()`, which holds the same information as index
#' vectors.
#'
#' @param initPars output of `genInitPars()`.
#' @param parameters named character vector of model symbols.
#' @param allLevels level of each entry of `initPars`.
#' @param dataFit data frame with the effect class columns.
#'
#' @return list of numeric 0/1 vectors of length `nrow(dataFit)`.
#'
#' @noRd
genMask <- function(initPars,
                    parameters,
                    allLevels,
                    dataFit) {
  mask <- lapply(
    seq_along(initPars),
    function(k) {
      effect <- getParClass(
        parameters[
          match(
            getParClass(initPars)[k],
            parameters
          )
        ]
      )
      mask_vector <- as.numeric(
        dataFit[[effect]] == allLevels[k]
      )
      return(mask_vector)
    }
  )

  return(mask)
}


# buildMaskIndex() --------------------------------------------------------

#' Index-based parameter activation map
#'
#' Each data point activates exactly one parameter per model symbol. Parameters
#' are ordered as in `genInitPars()`.
#'
#' @param parameters named character vector of model symbols, names are the
#'   effect classes.
#' @param levelsList named list of levels per effect class.
#' @param dataFit data frame with the effect class columns.
#'
#' @return list with `gidx` (per symbol, the global parameter index of each data
#'   point), `nPars`, `constrIdx` (indices of the first biological symbol, which
#'   the normalization constraint acts on) and `constrDenom` (their number).
#'
#' @noRd
buildMaskIndex <- function(parameters, levelsList, dataFit) {
  S <- length(parameters)
  effects <- names(parameters)
  nlev <- vapply(effects, function(e) length(levelsList[[e]]), integer(1))
  offsets <- cumsum(c(0L, nlev))[seq_len(S)]
  nPars <- sum(nlev)

  gidx <- lapply(
    seq_len(S),
    function(k) {
      as.integer(offsets[k]) +
        match(as.character(dataFit[[effects[k]]]), levelsList[[effects[k]]])
    }
  )

  firstBio <- which(effects == "biological")[1]
  constrIdx <- integer(0)
  constrDenom <- 1L
  if (!is.na(firstBio)) {
    constrIdx <- as.integer(offsets[firstBio]) + seq_len(nlev[firstBio])
    constrDenom <- nlev[firstBio]
  }

  list(
    gidx = gidx,
    nPars = nPars,
    constrIdx = constrIdx,
    constrDenom = constrDenom
  )
}


# resolveFunction() -----------------------------------------------------

#' Split the parameter vector by symbol and evaluate the residuals
#'
#' @param currentPars named parameter vector, without fixed parameters.
#' @param passParList,passParList2 settings, see `scaleTarget()`.
#' @param calcDeriv logical, passed to `rssModel()`.
#'
#' @return list with `residuals` from `rssModel()` and one named vector per
#'   model symbol, named by level.
#'
#' @noRd
resolveFunction <- function(currentPars,
                            passParList,
                            passParList2,
                            calcDeriv) {
  parameters <- passParList$parameters
  levelsList <- passParList2$levelsList

  ## Fixed parameters are merged back in global index order
  parsAll <- currentPars
  fixedPars <- passParList2$fixedPars
  if (!is.null(fixedPars)) {
    parsAll <- c(currentPars, fixedPars)
    parsAll <- parsAll[order(as.integer(sub(".*_", "", names(parsAll))))]
  }
  parClass <- getParClass(parsAll)
  parList <- lapply(
    seq_along(parameters),
    function(k) {
      subpar <- parsAll[parClass == parameters[k]]
      names(subpar) <- levelsList[[names(parameters)[k]]]
      subpar
    }
  )

  names(parList) <- parameters

  c(
    list(
      residuals = rssModel(
        parList,
        passParList,
        passParList2,
        calcDeriv = calcDeriv
      )
    ),
    parList
  )
}


# rssModel() -------------------------------------------------------------

#' Weighted residuals, variances and constraint
#'
#' @param parList named list with one vector per model symbol, named by level.
#' @param passParList,passParList2 settings, see `scaleTarget()`.
#' @param calcDeriv logical; if `TRUE`, the derivatives with respect to each
#'   symbol are attached as attributes `residualDeriv` and `varianceDeriv`.
#'
#' @return numeric vector: residuals divided by sigma, variances, constraint.
#'
#' @noRd
rssModel <- function(parList,
                     passParList,
                     passParList2,
                     calcDeriv = TRUE) {
  modelExpr <- passParList$modelExpr
  errorModelExpr <- passParList$errorModelExpr
  constrExpr <- passParList$constrExpr
  modelJacobianExpr <- passParList$modelJacobianExpr
  errorModelJacobianExpr <- passParList$errorModelJacobianExpr

  dataFit <- passParList2$dataFit

  with(
    c(as.list(dataFit), parList), {
      ## Binding for R CMD check, never executed
      if (FALSE) value <- NULL

      ## Assignment into a vector recycles scalar model results to all points
      prediction <- var <- rep(1, length(value))
      prediction[seq_along(prediction)] <- eval(modelExpr)
      res <- prediction - value
      var[seq_along(var)] <- eval(errorModelExpr)^2
      constr <- eval(constrExpr)
      residuals <- c(res / sqrt(var), var, constr)

      residualDeriv <- varianceDeriv <- NULL
      if (calcDeriv) {
        jac <- lapply(
          seq_len(length(parList)),
          function(k) {
            vMod <- vErr <- rep(0, length(value))
            vMod[seq_len(length(value))] <- eval(
              modelJacobianExpr[[k]]
            )
            vErr[seq_len(length(value))] <- eval(
              errorModelJacobianExpr[[k]]
            )

            residualDeriv <- vMod / sqrt(var) - vErr * res / var
            varianceDeriv <- vErr * 2 * sqrt(var)
            list(residualDeriv, varianceDeriv)
          }
        )

        residualDeriv <- lapply(jac, function(j) j[[1]])
        varianceDeriv <- lapply(jac, function(j) j[[2]])
      }
      attr(residuals, "residualDeriv") <- residualDeriv
      attr(residuals, "varianceDeriv") <- varianceDeriv

      return(residuals)
    }
  )
}


# objFunction() ----------------------------------------------------

#' Objective function for `trust()`
#'
#' @param currentPars named parameter vector, without fixed parameters.
#' @param passParList,passParList2 settings, see `scaleTarget()`.
#' @param calcDeriv logical; if `TRUE`, gradient and Gauss-Newton Hessian are
#'   returned.
#'
#' @return list with `value` (\eqn{-2 \log L}), `gradient`, `hessian`,
#'   `constrVec` (rank-1 constraint term, sparse path only), `prediction` and
#'   `sigma`. Non-finite results are returned as `value = Inf` without
#'   derivatives.
#'
#' @noRd
objFunction <- function(currentPars,
                        passParList,
                        passParList2,
                        calcDeriv = TRUE) {
  parameters <- passParList$parameters
  fitLogscale <- passParList$fitLogscale
  constStrength <- passParList$constStrength

  dataFit <- passParList2$dataFit
  levelsList <- passParList2$levelsList

  noData <- nrow(dataFit)

  calculatedResiduals <- resolveFunction(
    currentPars = currentPars,
    passParList = passParList,
    passParList2 = passParList2,
    calcDeriv = calcDeriv
  )$residuals

  residuals <- calculatedResiduals[1:noData]
  variances <- calculatedResiduals[seq((noData + 1), (2 * noData))]
  constraint <- calculatedResiduals[2 * noData + 1]
  residualDeriv <- attr(calculatedResiduals, "residualDeriv")
  varianceDeriv <- attr(calculatedResiduals, "varianceDeriv")

  bessel <- 1
  isT <- identical(passParList$likelihood, "t")
  tIdx <- passParList2$tGroupIdx
  if (isT && is.null(tIdx)) {
    nu <- passParList$df
    tConst <- log(nu * pi) - 2 * (lgamma((nu + 1) / 2) - lgamma(nu / 2))
    tWeight <- (nu + 1) / (nu + residuals^2)
    expWeight <- (nu + 1) / (nu + 3)
    value <- (nu + 1) * sum(log1p(residuals^2 / nu)) + sum(log(variances)) +
      noData * tConst + constraint^2
  } else if (isT) {
    ## Multivariate t per group: one shared scale weight for all its points
    nu <- passParList$df
    nG <- passParList2$tGroupSize
    Q <- as.vector(rowsum(residuals^2, tIdx, reorder = TRUE))
    tConst <- nG * log(nu * pi) - 2 * (lgamma((nu + nG) / 2) - lgamma(nu / 2))
    tWeight <- ((nu + nG) / (nu + Q))[tIdx]
    expWeight <- ((nu + nG) / (nu + nG + 2))[tIdx]
    value <- sum((nu + nG) * log1p(Q / nu)) + sum(log(variances)) +
      sum(tConst) + constraint^2
  } else {
    value <- sum(residuals^2) + bessel * sum(log(pi * variances)) + constraint^2
  }
  gradient <- NULL
  hessian <- NULL
  constrVec <- NULL
  sparseHessian <- isTRUE(passParList$sparseHessian)

  ## Derivatives assemble from the per-symbol index vectors of buildMaskIndex(),
  ## without a dense Jacobian
  if (calcDeriv) {
    maskIdx <- passParList2$maskIdx
    if (is.null(maskIdx)) {
      maskIdx <- buildMaskIndex(parameters, levelsList, dataFit)
    }

    nPars <- maskIdx$nPars
    S <- length(parameters)
    gidx <- maskIdx$gidx
    freeIdx <- passParList2$freeIdx
    fullPars <- currentPars
    if (!is.null(freeIdx)) {
      fullPars <- numeric(nPars)
      fullPars[freeIdx] <- currentPars
      fullPars[-freeIdx] <- passParList2$fixedPars
    }

    ## Gradient: per symbol, segment sums over the data points of each level
    gradient <- numeric(nPars)
    for (k in seq_len(S)) {
      rTerm <- if (isT) tWeight * residuals else residuals
      contrib <- 2 * rTerm * residualDeriv[[k]] +
        (bessel / variances) * varianceDeriv[[k]]
      acc <- rowsum(contrib, gidx[[k]], reorder = FALSE)
      idx <- as.integer(rownames(acc))
      gradient[idx] <- gradient[idx] + acc[, 1L]
    }

    ## The constraint acts on the first biological symbol only
    cj <- numeric(nPars)
    ci <- maskIdx$constrIdx
    if (length(ci)) {
      cjVal <- constStrength / maskIdx$constrDenom
      if (fitLogscale == TRUE) cjVal <- cjVal * exp(fullPars[ci])
      cj[ci] <- cjVal
      gradient <- gradient + 2 * constraint * cj
    }

    ## t: IRLS weights while fitting, expected information for uncertainties
    jacEntries <- residualDeriv[seq_len(S)]
    if (isT) {
      hessWeight <- if (isTRUE(passParList2$expectedInfo)) {
        expWeight
      } else {
        tWeight
      }
      jacEntries <- lapply(jacEntries, function(d) sqrt(hessWeight) * d)
    }
    residualJacobian <- Matrix::sparseMatrix(
      i = rep.int(seq_len(noData), S),
      j = unlist(gidx, use.names = FALSE),
      x = unlist(jacEntries, use.names = FALSE),
      dims = c(noData, nPars)
    )
    gnHessian <- 2 * Matrix::crossprod(residualJacobian)

    if (sparseHessian) {
      ## Full Hessian is gnHessian + 2 * tcrossprod(constrVec)
      hessian <- gnHessian
      constrVec <- cj
    } else {
      ## The rank-1 constraint term is non-zero on the constrained block only
      hessian <- as.matrix(gnHessian)
      if (length(ci)) {
        hessian[ci, ci] <- hessian[ci, ci] + 2 * outer(cj[ci], cj[ci])
      }
    }
    if (!is.null(freeIdx)) {
      gradient <- gradient[freeIdx]
      hessian <- hessian[freeIdx, freeIdx, drop = FALSE]
      if (!is.null(constrVec)) constrVec <- constrVec[freeIdx]
    }
  }

  prediction <- residuals * sqrt(variances) + dataFit$value
  sigma <- sqrt(variances)

  ## Non-finite results become value = Inf, which trust() rejects as a step
  hessianFinite <- if (is.null(hessian)) {
    TRUE
  } else if (sparseHessian) {
    all(is.finite(hessian@x)) && all(is.finite(constrVec))
  } else {
    all(is.finite(hessian))
  }
  badDeriv <- calcDeriv && (any(!is.finite(gradient)) || !hessianFinite)
  if (!is.finite(value) || badDeriv) {
    value <- Inf
    gradient <- NULL
    hessian <- NULL
    constrVec <- NULL
  }

  list(
    value = value, gradient = gradient, hessian = hessian,
    constrVec = constrVec, prediction = prediction, sigma = sigma
  )
}


# evaluateModel() --------------------------------------------------------

#' Model prediction minus data for given biological parameter values
#'
#' @param initPars values of the first biological symbol, one per data point.
#' @param parList named list with one vector per model symbol, named by level.
#' @param passParList,passParList2 settings, see `scaleTarget()`.
#'
#' @return numeric vector, one entry per data point.
#'
#' @noRd
evaluateModel <- function(initPars,
                          parList,
                          passParList,
                          passParList2) {
  effectsPars <- passParList$effectsPars
  modelExpr <- passParList$modelExpr
  dataFit <- passParList2$dataFit

  ## The biological class index is the data point, so entry i applies to point i
  biological <- seq_along(initPars)
  useList <- c(
    list(initPars, biological = biological),
    as.list(dataFit),
    parList
  )
  names(useList)[1] <- effectsPars[[1]][1]

  values <- with(useList, eval(modelExpr) - dataFit$value)

  return(values)
}


# evaluateModelJacobian() -----------------------------------------------

#' Derivative of the model with respect to the first biological symbol
#'
#' @inheritParams evaluateModel
#'
#' @return numeric vector, one entry per data point.
#'
#' @noRd
evaluateModelJacobian <- function(initPars,
                                  parList,
                                  passParList,
                                  passParList2) {
  effectsPars <- passParList$effectsPars
  modelDerivExpr <- passParList$modelDerivExpr
  dataFit <- passParList2$dataFit

  biological <- seq_along(initPars)
  useList <- c(
    list(initPars, biological = biological), as.list(dataFit),
    parList
  )
  names(useList)[1] <- effectsPars[[1]][1]
  derivative_values <- with(useList, eval(modelDerivExpr))

  return(derivative_values)
}


# getMaxY() --------------------------------------------------------------

#' Rounded upper y limit for plot axes
#'
#' @param x numeric scalar, the largest y value.
#'
#' @return numeric scalar \eqn{\ge x}.
#'
#' @noRd
getMaxY <- function(x) {
  rx <- round(x, 1)
  lower <- floor(x)
  upper <- ceiling(x)
  lowDiff <- rx - lower
  uppDiff <- upper - rx

  if (x > 5) {
    return(if (upper %% 2 == 0) upper else lower + 1)
  }
  if (x >= 2) {
    return(if (lowDiff < uppDiff) lower + 0.5 else upper)
  }
  if (rx > x) {
    return(if (rx %% 0.2 == 0) rx else rx + 0.1)
  }
  if (rx %% 0.2 != 0) {
    return(rx + 0.1)
  }
  rx + 0.2
}


# scaleValues() ----------------------------------------------------------

#' Transform values and sigma from the fit scale to the linear scale
#'
#' @param fitLogscale logical; if `TRUE`, values are exponentiated and sigma is
#'   propagated linearly.
#' @param outputScale unused.
#' @param value,sigma numeric vectors on the fit scale.
#' @param upper,lower unused, recomputed as `value +- sigma`.
#'
#' @return list with `value`, `upper`, `lower` and `sigma`.
#'
#' @noRd
scaleValues <- function(fitLogscale,
                        outputScale,
                        value,
                        upper,
                        lower,
                        sigma) {
  if (fitLogscale == TRUE) {
    value <- exp(value)
    sigma <- value * sigma
  }
  upper <- value + sigma
  lower <- value - sigma

  return(list(value = value, upper = upper, lower = lower, sigma = sigma))
}


# getParClass() -------------------------------------------------------

#' Model symbol of indexed parameter names
#'
#' @param paramlist named vector with names of the form `<symbol>_<index>`.
#'
#' @return character vector of symbols.
#'
#' @noRd
getParClass <- function(paramlist) {
  sub("_[0-9]+$", "", names(paramlist))
}
