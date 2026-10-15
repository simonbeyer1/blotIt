# readWide() -------------------------------------------------------------

#' Import measurement data
#'
#' Reads a wide csv file with one column per target and returns it in long
#' format.
#'
#' @param file character, path of the data file.
#' @param description numeric or character vector, the columns describing the
#'   measurements. All other columns are targets.
#' @param time numeric or character of length 1, the time column. It is
#'   renamed to `time`.
#' @param header logical, whether the first line holds the column names.
#' @param ... further arguments passed to [utils::read.csv()].
#'
#' @return data frame with the description columns, `name` (the target) and
#'   `value`. Rows with missing values are dropped.
#'
#' @example inst/examples/exampleReadWide.R
#'
#' @export
#' @importFrom utils read.csv
readWide <- function(file, description = NULL, time = 1, header = TRUE, ...) {
  useData <- read.csv(file, header = header, ...)
  allNames <- colnames(useData)

  if (is.null(description)) {
    stop("Specify columns containing descriptions.")
  }
  if (is.character(description)) {
    noDescription <- which(colnames(useData) %in% description)
  } else {
    noDescription <- description
  }
  if (is.character(time)) {
    noTime <- which(colnames(useData) == time)
  } else {
    noTime <- time
  }

  if (length(noDescription) < length(description)) {
    warning(
      "Not all columns proposed by argument 'description' are available",
      " in file.\nTaking the available ones."
    )
  }

  if (length(noTime) == 0) {
    stop(
      "File did not contain a time column as proposed by 'time' argument."
    )
  }

  descrEntries <- useData[, noDescription]
  restLong <- unlist(useData[, -noDescription])

  newData <- data.frame(
    descrEntries,
    name = rep(allNames[-noDescription], each = dim(useData)[1]),
    value = restLong
  )

  newData <- newData[!is.nan(newData$value), ]
  newData <- newData[!is.na(newData$value), ]

  colnames(newData)[noTime] <- "time"

  return(newData)
}


# plotIt() ---------------------------------------------------------

#' Plot alignment results
#'
#' Plots one table of an [alignReplicates()] result as points with error bars
#' and one as a line with an error band.
#'
#' @param inputList output of [alignReplicates()].
#' @param plotPoints character, the table plotted as points with error bars.
#'   One of `"original"`, `"scaled"`, `"prediction"`, `"aligned"`.
#' @param plotLine character, the table plotted as a line with error band.
#'   Same choices as `plotPoints`.
#' @param spline logical, if `TRUE` the line is a cubic polynomial fit instead
#'   of straight segments between points.
#' @param scales character, passed as `scales` to [ggplot2::facet_wrap()].
#'   `"biological"` is the same as `"free_x"`.
#' @param alignZeros logical, if `TRUE` all facets share the lower y limit.
#' @param plotCaption logical, if `TRUE` a caption names the plotted tables and
#'   error sources.
#' @param ncol numeric, passed as `ncol` to [ggplot2::facet_wrap()].
#' @param useColors character vector of colors, passed as `values` to
#'   [ggplot2::scale_color_manual()]. `NULL` uses the ggplot2 default.
#' @param useOrder character vector of all target names, the facet order.
#' @param plotScaleX,plotScaleY character, transformation of the x and y axis:
#'   `"log"`, `"log2"`, `"log10"`, and for x also `"pseudoLog10"`. `NULL` uses
#'   the `outputScale` of `inputList`.
#' @param doseResponse logical, if `TRUE` the x axis is the dose column
#'   `doseValueName` instead of `time`.
#' @param doseValueName character, the dose column used if `doseResponse` is
#'   `TRUE`.
#' @param labelX,labelY character, axis labels. `NULL` gives `"Time"` or
#'   `"Dose"`, and `"Signal"`.
#' @param ... logical expression for subsetting the plotted tables, e.g.
#'   `name == "pERK1" & time < 60`.
#' @param FacetLabels character, `"simple"` joins the facet values with
#'   underscores, `"full"` prints them as `column: value` lines.
#'
#' @return ggplot object.
#'
#' @example inst/examples/examplePlotIt.R
#'
#' @import ggplot2 data.table
#' @export
#' @author Severin Bang and Svenja Kemmer
plotIt <- function(
    inputList,
    ...,
    plotPoints = "aligned",
    plotLine = "aligned",
    spline = FALSE,
    scales = "free",
    alignZeros = TRUE,
    plotCaption = TRUE,
    ncol = NULL,
    useColors = NULL,
    useOrder = NULL,
    plotScaleY = NULL,
    plotScaleX = NULL,
    doseResponse = FALSE,
    doseValueName = "dose",
    labelX = NULL,
    labelY = NULL,
    FacetLabels = c("simple", "full")[1]
) {
  if (!plotPoints %in% c("original", "scaled", "prediction", "aligned") |
      !plotLine %in% c("original", "scaled", "prediction", "aligned")) {
    stop(
      "\n\t'plotPoints' and 'plotLine' must each be one of
            c('original', 'scaled', 'prediction', 'aligned')\n"
    )
  }

  if (!(FacetLabels == "simple" | FacetLabels == "full")) {
    stop("'FacetLabels' must be either'simple' or 'full'")
  }

  if (!is.null(useOrder)) {
    if (length(setdiff(levels(inputList[[1]]$name), useOrder)) != 0) {
      stop("useOrder doesn't contain all protein names.")
    }
    for (tab in c("aligned", "scaled", "prediction", "original")) {
      inputList[[tab]]$name <- factor(inputList[[tab]]$name, levels = useOrder)
    }
  }

  biological <- inputList$biological
  scaling <- inputList$scaling

  plotList <- inputList

  if (doseResponse == TRUE) {
    xValue <- doseValueName
  } else {
    xValue <- "time"
  }
  if (FacetLabels == "simple") {
    plotList$aligned$biological <- do.call(
      paste,
      c(
        plotList$aligned[
          ,
          biological[!(biological %in% c("name", xValue))],
          drop = FALSE
        ],
        list(sep="_")
      )
    )
    plotList$aligned$scaling <- NA

    ## scaled
    plotList$scaled$biological <- do.call(
      paste,
      c(
        inputList$scaled[
          ,
          biological[!(biological %in% c("name", xValue))],
          drop = FALSE
        ],
        list(sep="_")
      )
    )
    plotList$scaled$scaling <- do.call(
      paste,
      c(
        inputList$scaled[, scaling[scaling != "name"], drop = FALSE],
        list(sep="_")
      )
    )

    ## prediction
    plotList$prediction$biological <- do.call(
      paste,
      c(
        inputList$prediction[
          ,
          biological[!(biological %in% c("name", xValue))],
          drop = FALSE
        ],
        list(sep="_")
      )
    )
    plotList$prediction$scaling <- do.call(
      paste,
      c(
        inputList$prediction[, scaling[scaling != "name"], drop = FALSE],
        list(sep="_")
      )
    )

    ## original
    plotList$original$biological <- do.call(
      paste,
      c(
        inputList$original[
          ,
          biological[!(biological %in% c("name", xValue))],
          drop = FALSE
        ],
        list(sep="_")
      )
    )
    plotList$original$scaling <- do.call(
      paste,
      c(
        inputList$original[, scaling[scaling != "name"], drop = FALSE],
        list(sep="_")
      )
    )
  } else if (FacetLabels == "full") {
    plotList$aligned$biological <- do.call(
      c,
      lapply(
        seq(nrow(plotList$aligned)),
        function(i) {
          paste0(
            lapply(
              biological[!(biological %in% c("name", xValue))],
              function(j) paste0(j, ": ", plotList$aligned[i,j])
            ),
            collapse = "\n"
          )
        }
      )
    )
    plotList$aligned$scaling <- NA

    ## scaled
    plotList$scaled$biological <- do.call(
      c,
      lapply(
        seq(nrow(plotList$scaled)),
        function(i) {
          paste0(
            lapply(
              biological[!(biological %in% c("name", xValue))],
              function(j) paste0(j, ": ", plotList$scaled[i,j])
            ),
            collapse = "\n"
          )
        }
      )
    )
    plotList$scaled$scaling <- do.call(
      c,
      lapply(
        seq(nrow(plotList$scaled)),
        function(i) {
          paste0(
            lapply(
              scaling[!(scaling %in% c("name"))],
              function(j) paste0(j, ": ", plotList$scaled[i,j])
            ),
            collapse = "\n"
          )
        }
      )
    )

    ## prediction
    plotList$prediction$biological <- do.call(
      c,
      lapply(
        seq(nrow(plotList$prediction)),
        function(i) {
          paste0(
            lapply(
              biological[!(biological %in% c("name", xValue))],
              function(j) paste0(j, ": ", plotList$prediction[i,j])
            ),
            collapse = "\n"
          )
        }
      )
    )
    plotList$prediction$scaling <- do.call(
      c,
      lapply(
        seq(nrow(plotList$prediction)),
        function(i) {
          paste0(
            lapply(
              scaling[!(scaling %in% c("name"))],
              function(j) paste0(j, ": ", plotList$prediction[i,j])
            ),
            collapse = "\n"
          )
        }
      )
    )

    ## original
    plotList$original$biological <- do.call(
      c,
      lapply(
        seq(nrow(plotList$original)),
        function(i) {
          paste0(
            lapply(
              biological[!(biological %in% c("name", xValue))],
              function(j) paste0(j, ": ", plotList$original[i,j])
            ),
            collapse = "\n"
          )
        }
      )
    )
    plotList$original$scaling <- do.call(
      c,
      lapply(
        seq(nrow(plotList$original)),
        function(i) {
          paste0(
            lapply(
              scaling[!(scaling %in% c("name"))],
              function(j) paste0(j, ": ", plotList$original[i,j])
            ),
            collapse = "\n"
          )
        }
      )
    )
  }

  plotListPoints <- plotList[[plotPoints]]
  plotListLine <- plotList[[plotLine]]

  plotListPoints <- subset(plotListPoints, ...)
  plotListLine <- subset(plotListLine, ...)

  usedErrors <- list(
    aligned = "Fisher Information",
    scaled = "Propagated error model to common scale",
    prediction = "Error model",
    original = "None"
  )

  usedData <- list(
    aligned = "Estimated true values",
    scaled = "Original data scaled to common scale",
    prediction = "Predictions from model evaluation on original scale",
    original = "Original data"
  )

  captionText <- paste0(
    "Datapoints: ", usedData[[plotPoints]], "\n",
    "Errorbars: ", usedErrors[[plotPoints]], "\n",
    "Line: ", usedData[[plotLine]], "\n",
    if (plotPoints != plotLine) {
      paste0("Errorband: ", usedErrors[[plotLine]], "\n")
    },
    "\n",
    "Date: ", Sys.Date()
  )

  if (scales == "biological") {
    scales <- "free_x"
  }

  if (doseResponse == TRUE) {
    Xlabel <-  "Dose"
    xVariable <- doseValueName
  } else {
    Xlabel <-  "Time"
    xVariable <- "time"
  }

  Ylabel <- "Signal"

  if (!is.null(labelX)){
    Xlabel <- labelX
  }

  if (!is.null(labelY)){
    Ylabel <- labelY
  }

  if (!is.null(plotScaleX)) {
    errWidth <- 0
  } else {
    errWidth <- max(as.numeric(plotListPoints[[xVariable]]))/50
  }

  plotListPoints[[xVariable]] <- as.numeric(plotListPoints[[xVariable]])
  plotListLine[[xVariable]] <- as.numeric(plotListLine[[xVariable]])

  if (plotPoints == "aligned" & plotLine == "aligned") {
    g <- ggplot(
      data = plotListPoints,
      aes(
        x = .data[[xVariable]],
        y = .data$value,
        group = .data$biological,
        color = .data$biological,
        fill = .data$biological
      )
    )
    g <- g + facet_wrap(~name, scales = scales, ncol = ncol, labeller = labeller(.default = label_both))
    if (!is.null(useColors)) {
      useColors <- c(useColors, rep("gray", 100))
      g <- g + scale_color_manual("Biological", values = useColors) +
        scale_fill_manual("Biological", values = useColors)
    }
  } else {
    g <- ggplot(
      data = plotListPoints,
      aes(
        x = .data[[xVariable]],
        y = .data$value,
        group = .data$scaling,
        color = .data$scaling,
        fill = .data$scaling
      )
    )
    g <- g + facet_wrap(
      ~ name * biological,
      scales = scales,
      ncol = ncol,
      labeller = labeller(.default = label_both)
    )
    if (!is.null(useColors)) {
      useColors <- c(useColors, rep("gray", 100))
      g <- g + scale_color_manual("Scaled", values = useColors) +
        scale_fill_manual("Scaled", values = useColors)
    }
  }


  g <- g + geom_point(data = plotListPoints, size = 2.5)
  g <- g + geom_line(data = plotListLine, linewidth = 1)

  ## Error bars on points and bands on lines only where the table has sigma
  pointBars <- plotPoints != "original" && plotLine != plotPoints
  lineBand <- switch(plotPoints,
    original = plotLine != "original",
    aligned = plotLine == "aligned",
    TRUE
  )
  if (pointBars) {
    g <- g + geom_errorbar(
      data = plotListPoints,
      aes(ymin = lower, ymax = upper),
      linewidth = 0.5,
      width = errWidth,
      alpha = 0.5
    )
  }
  if (lineBand && spline) {
    g <- g + geom_smooth(
      data = plotListLine,
      se = FALSE,
      method = "lm",
      formula = y ~ poly(x, 3)
    )
  }
  if (lineBand && !spline && plotLine != "original") {
    g <- g + geom_ribbon(
      data = plotListLine,
      aes(ymin = lower, ymax = upper),
      alpha = 0.3,
      lty = 0
    )
  }

  g <- g + theme_bw(base_size = 20) +
    theme(
      legend.position = "top",
      legend.key = element_blank(),
      strip.background = element_rect(color = NA, fill = NA),
      axis.line.x = element_line(linewidth = 0.3, colour = "black"),
      axis.line.y = element_line(linewidth = 0.3, colour = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    )
  g <- g + xlab(paste0("\n",Xlabel)) + ylab(paste0(Ylabel,"\n"))

  if (alignZeros && plotPoints != "original") {
    ## Common lower y limit for all facets, upper limit per facet
    plotListPoints <- as.data.table(plotListPoints)
    blankData <- plotListPoints[
      ,
      list(ymax = max(upper), ymin = min(lower)),
      by = c("name", "biological", "scaling")
    ]
    blankData[, ":="(ymin = min(ymin))]
    blankData[
      ,
      ":="(ymax = getMaxY(ymax)),
      by = c("name", "biological", "scaling")
    ]
    blankData <- melt(
      blankData,
      id.vars = c("name", "biological", "scaling"),
      measure.vars = c("ymax", "ymin"),
      value.name = "value"
    )
    blankData[, ":="(xVariable = 0, variable = NULL)]
    setnames(blankData, "xVariable", xVariable)
    g <- g + geom_blank(
      data = as.data.frame(blankData),
      aes(x = .data[[xVariable]], y = .data$value)
    )
  }

  if (plotCaption) {
    g <- g + labs(caption = captionText)
  }

  if (is.null(plotScaleY)) {
    if (inputList$outputScale != "linear") {
      g <- g + coord_trans(y = inputList$outputScale)
    }
  } else if (plotScaleY %in% c("log", "log2", "log10")) {
    g <- g + scale_y_continuous(trans = plotScaleY)
  }

  if (is.null(plotScaleX)) {
    if (inputList$outputScale != "linear") {
      g <- g + coord_trans(x = inputList$outputScale)
    }
  } else if (plotScaleX == "pseudoLog10") {
    g <- g + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))
  } else if (plotScaleX %in% c("log", "log2", "log10")) {
    g <- g + scale_x_continuous(trans = plotScaleX)
  }

  return(g)
}


# llrTest() --------------------------------------------------------------

#' Likelihood ratio test
#'
#' Tests two nested [alignReplicates()] results against each other, e.g.
#' whether a buffer effect is negligible or a measurement is an outlier.
#'
#' @param H0 output of [alignReplicates()] under the null hypothesis, a special
#'   case of `H1`.
#' @param H1 output of [alignReplicates()] under the alternative.
#' @param check logical, if `TRUE` stop unless the effects of `H0` are a subset
#'   of those of `H1`.
#'
#' @return list with `llr` (difference of \eqn{-2 \log L}), `statistic` (a
#'   description of the reference distribution) and `p.value` from the
#'   \eqn{\chi^2} distribution with the difference of the degrees of freedom.
#'
#' @example inst/examples/exampleLlrTest.R
#'
#' @export
llrTest <- function(H0, H1, check = TRUE) {
  if (!identical(H0$likelihood, H1$likelihood) ||
      !identical(H0$likelihoodDf, H1$likelihoodDf) ||
      !identical(H0$tGroup, H1$tGroup)) {
    stop("H0 and H1 must use the same likelihood.")
  }

  biological0 <- union(
    H0$biological[!(H0$biological %in% c("name", "time"))],
    "1"
  )
  biological1 <- union(
    H1$biological[!(H1$biological %in% c("name", "time"))],
    "1"
  )

  scaling0 <- union(H0$scaling[H0$scaling != "name"], "1")
  scaling1 <- union(H1$scaling[H1$scaling != "name"], "1")

  error0 <- union(H0$error[H0$error != "name"], "1")
  error1 <- union(H1$error[H1$error != "name"], "1")

  if (check) {
    if (
      !all(biological0 %in% biological1) |
      !all(scaling0 %in% scaling1) | !all(error0 %in% error1)
    ) {
      stop("H0 is not a special case of H1.")
    }
  }

  value0 <- attr(H0$parameter, "value")
  value1 <- attr(H1$parameter, "value")
  df0 <- attr(H0$parameter, "df")
  df1 <- attr(H1$parameter, "df")

  list(
    llr = value0 - value1,
    statistic = paste0(
      "chisquare with ", df0 - df1, " degrees of freedom."
    ),
    p.value = pchisq(value0 - value1, df = df0 - df1, lower.tail = FALSE)
  )
}


# plotWaterfall() ---------------------------------------------------------

#' Waterfall plot of multistart fits
#'
#' Shows the final objective value of every start per target, sorted
#' increasingly. Starts with equal values reached the same optimum.
#'
#' @param out output of [alignReplicates()].
#' @param tol numeric, starts with an objective value within `tol` of the best
#'   value are colored as best fits.
#'
#' @return ggplot object.
#'
#' @example inst/examples/examplePlotWaterfall.R
#'
#' @export
plotWaterfall <- function(out, tol = 1) {
  ms <- out$multistart
  if (is.null(ms) || nrow(ms) == 0) stop("'out' has no multistart table.")
  ms <- ms[!is.na(ms$value), ]
  ms <- ms[order(ms$name, ms$value), ]
  ms$index <- stats::ave(ms$value, ms$name, FUN = seq_along)
  best <- stats::ave(ms$value, ms$name, FUN = min)
  ms$fit <- ifelse(!ms$converged, "not converged",
                   ifelse(ms$value - best <= tol, "best", "other"))
  ggplot(ms, aes(x = .data$index, y = .data$value, color = .data$fit)) +
    geom_point() +
    facet_wrap(~name, scales = "free") +
    scale_color_manual(values = c(
      best = "#D55E00", other = "black", `not converged` = "grey60"
    )) +
    labs(x = "start (sorted)", y = "-2 log L", color = NULL) +
    theme_bw()
}
