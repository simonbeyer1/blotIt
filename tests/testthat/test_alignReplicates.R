# alignReplicates() - model parsing ----------------------------------------------

test_that("alignReplicates() - model parsing", {
    simDataWideFile <- system.file(
        "extdata", "simDataWide.csv",
        package = "blotIt"
    )

    input_data <- readWide(simDataWideFile, description = seq_len(3))

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = "one_parameter",
            scaling = sj ~ name + experiment,
            error = sigmaR ~ name + 1
        ),
        "Do not pass biological, scaling or error as string."
    )

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel <- "value * sigmaR",
            biological <- yi ~ name + time + condition,
            error = sigmaR ~ name + 1
        ),
        "All of model, errorModel, biological, scaling, error must be set"
    )

    expect_error(
        alignReplicates(
            data = input_data,
            biological = NULL,
            error = left ~ right2
        ),
        "All of model, errorModel, biological, scaling, error must be set."
    )

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = yi ~ condition,
            scaling = wrong ~ experiment,
            error = sigmaR ~ name + 1,
            fitLogscale = TRUE
        ),
        "Not all paramters are defined in either arguments
         'scaling', 'biological' or 'error'"
    )

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = yi ~ condition,
            scaling = sj ~ experiment,
            error = sigmaR ~ name + 1,
            fitLogscale = "string"
        ),
        "'fitLogscale' must be logical"
    )
})


# alignReplicates() - unchanged default behavior ---------------------------------

fitSimLog2 <- function(data = NULL, error = sigmaA ~ name, ...) {
    if (is.null(data)) {
        data <- readWide(
            system.file("extdata", "simDataWide.csv", package = "blotIt"),
            description = seq_len(3)
        )
        data$value <- log2(data$value)
    }
    suppressWarnings(alignReplicates(
        data, model = "yi + sj", errorModel = "exp(sigmaA)",
        biological = yi ~ name + time + condition,
        scaling = sj ~ name + experiment, error = error,
        normalizeInput = FALSE, fitLogscale = FALSE, ...
    ))
}

## The fixtures contain the internal column "1", which results no longer have
readFixture <- function(file) {
    lapply(readRDS(test_path("fixtures", file)), function(x) {
        if (is.data.frame(x)) x[["1"]] <- NULL
        x
    })
}

test_that("alignReplicates() - default call reproduces v0.2.0", {
    ref <- readFixture("default_v0.2.0.rds")
    out <- fitSimLog2()
    refMs <- readFixture("multistart_v0.2.0.rds")
    set.seed(1)
    outMs <- fitSimLog2(multistart = 4, likelihood = "normal")
    for (n in names(ref)) {
        expect_equal(out[[n]], ref[[n]], tolerance = 1e-8)
        expect_equal(outMs[[n]], refMs[[n]], tolerance = 1e-8)
    }
    ## Bitwise identity holds only on the platform that produced the fixtures
    skip_on_ci()
    skip_on_cran()
    for (n in names(ref)) {
        expect_identical(out[[n]], ref[[n]])
        expect_identical(outMs[[n]], refMs[[n]])
    }
})

test_that("alignReplicates() - multistart result does not depend on cores", {
    skip_on_os("windows")
    set.seed(5)
    seq1 <- fitSimLog2(multistart = 4)
    set.seed(5)
    par2 <- fitSimLog2(multistart = 4, cores = 2)
    expect_identical(seq1, par2)
    expect_equal(nrow(par2$multistart), 4 * length(unique(par2$original$name)))
    expect_s3_class(plotWaterfall(par2), "ggplot")
})


# alignReplicates() - parinit ---------------------------------------------------

simShifted <- function(seed, offset) {
    set.seed(seed)
    d <- expand.grid(
        time = seq_len(8), gel = paste0("G", seq_len(5)),
        stringsAsFactors = FALSE
    )
    d$name <- "pX"
    gelEffect <- rnorm(5, 0, 2)
    d$value <- offset + 2 * sin(d$time) + gelEffect[match(d$gel, unique(d$gel))] +
        rnorm(nrow(d), 0, 0.2)
    d
}

test_that("alignReplicates() - parinit on offset data matches centered fit", {
    fit <- function(d, ...) {
        suppressWarnings(alignReplicates(
            d, model = "yi + sj",
            errorModel = "sqrt(exp(2*sigmaA) + exp(2*sigmaB) * 2^(-2*value))",
            biological = yi ~ name + time, scaling = sj ~ name + gel,
            error = sigmaA + sigmaB ~ name,
            normalizeInput = FALSE, fitLogscale = FALSE, ...
        ))
    }
    centered <- fit(simShifted(1, 0))
    offset <- fit(
        simShifted(1, 10),
        parinit = list(sj = 9, sigmaA = log(0.2), sigmaB = log(0.2))
    )
    ## Equal up to the normalization of the biological parameters
    ctr <- function(x) x - mean(x)
    expect_equal(
        ctr(offset$aligned$value), ctr(centered$aligned$value),
        tolerance = 1e-4
    )
})

test_that("alignReplicates() - parinit input validation and level names", {
    expect_error(fitSimLog2(parinit = list(nope = 1)), "parinit")
    expect_error(fitSimLog2(parinit = list(sj = c(1, 2, 3))), "named by level")
    out <- fitSimLog2(
        parinit = list(sigmaA = c(pAKT = 0.1, pERK = 0.2)), iterlim = 0
    )
    sig <- out$parameter[grepl("^sigmaA", out$parameter$parameter), ]
    expect_equal(sig$value, c(0.1, 0.2))
})


# alignReplicates() - error formula list and inactive parameters ------------------

simTwoModalities <- function() {
    set.seed(3)
    d <- expand.grid(
        time = 1:6, gel = c("W1", "W2", "W3", "P1", "P2", "P3"),
        stringsAsFactors = FALSE
    )
    d$name <- "pX"
    d$modality <- ifelse(substr(d$gel, 1, 1) == "W", "WB", "PCR")
    d$isWB <- as.numeric(d$modality == "WB")
    d$value <- sin(d$time) + rnorm(6)[match(d$gel, unique(d$gel))] +
        rnorm(nrow(d), 0, ifelse(d$isWB == 1, 0.3, 0.1))
    d
}

test_that("alignReplicates() - list of error formulas gives own levels", {
    out <- alignReplicates(
        simTwoModalities(), model = "yi + sj",
        errorModel = "sqrt(exp(2*sigmaA)^2 + exp(2*sigmaB))",
        biological = yi ~ name + time, scaling = sj ~ name + gel,
        error = list(sigmaA ~ name + modality, sigmaB ~ name),
        normalizeInput = FALSE, fitLogscale = FALSE
    )
    par <- out$parameter
    expect_setequal(par$level[grepl("^sigmaA", par$parameter)], c("pX_WB", "pX_PCR"))
    expect_equal(par$level[grepl("^sigmaB", par$parameter)], "pX")
    expect_equal(out$error, c("name", "modality"))
})

test_that("alignReplicates() - parameter without influence is NA", {
    out <- expect_no_error(alignReplicates(
        simTwoModalities(), model = "yi + sj",
        errorModel = "sqrt(exp(2*sigmaA) + isWB * exp(2*sigmaB))",
        biological = yi ~ name + time, scaling = sj ~ name + gel,
        error = list(sigmaA ~ name + modality, sigmaB ~ name + modality),
        normalizeInput = FALSE, fitLogscale = FALSE
    ))
    par <- out$parameter
    off <- par[par$level == "pX_PCR" & grepl("^sigmaB", par$parameter), ]
    on <- par[par$level == "pX_WB" & grepl("^sigmaB", par$parameter), ]
    expect_true(is.na(off$value) && is.na(off$sigma))
    expect_true(is.finite(on$value) && is.finite(on$sigma))
    expect_true(all(is.finite(out$prediction$sigma)))
    expect_equal(attr(par, "df"), nrow(out$original) + 1 - sum(!is.na(par$value)))
})

test_that("alignReplicates() - duplicate error parameter is rejected", {
    expect_error(
        alignReplicates(
            simTwoModalities(), model = "yi + sj", errorModel = "exp(sigmaA)",
            biological = yi ~ name + time, scaling = sj ~ name + gel,
            error = list(sigmaA ~ name, sigmaA ~ modality)
        ),
        "only one"
    )
})


# alignReplicates() - Student-t likelihood ----------------------------------------

simOutlier <- function(outlier) {
    set.seed(11)
    d <- expand.grid(
        time = seq_len(8), gel = paste0("G", seq_len(4)),
        stringsAsFactors = FALSE
    )
    d$name <- "pX"
    d$value <- 2 * sin(d$time) + rnorm(4)[match(d$gel, unique(d$gel))] +
        rnorm(nrow(d), 0, 0.1)
    d$value[3] <- d$value[3] + outlier
    d
}

fitOutlier <- function(d, ...) {
    out <- suppressWarnings(alignReplicates(
        d, model = "yi + sj", errorModel = "exp(sigmaA)",
        biological = yi ~ name + time, scaling = sj ~ name + gel,
        error = sigmaA ~ name, normalizeInput = FALSE, fitLogscale = FALSE, ...
    ))
    out$aligned$value
}

test_that("alignReplicates() - t likelihood limits the shift by an outlier", {
    clean <- fitOutlier(simOutlier(0))
    shiftNormal <- max(abs(fitOutlier(simOutlier(5)) - clean))
    shiftT <- max(abs(fitOutlier(simOutlier(5), likelihood = "t", df = 3) - clean))
    expect_lt(shiftT, shiftNormal / 5)
})

test_that("alignReplicates() - t likelihood output and validation", {
    out <- suppressWarnings(alignReplicates(
        simOutlier(5), model = "yi + sj", errorModel = "exp(sigmaA)",
        biological = yi ~ name + time, scaling = sj ~ name + gel,
        error = sigmaA ~ name, normalizeInput = FALSE, fitLogscale = FALSE,
        likelihood = "t", df = 3
    ))
    expect_identical(out$likelihood, "t")
    expect_identical(out$likelihoodDf, 3)
    expect_true(all(is.finite(out$aligned$sigma) & out$aligned$sigma > 0))
    expect_true(all(out$multistart$converged))
    expect_error(fitSimLog2(likelihood = "t", df = 0), "'df'")
    expect_error(fitSimLog2(likelihood = "cauchy"), "should be one of")
    expect_error(llrTest(fitSimLog2(), out), "same likelihood")
})


# alignReplicates() - previous result as data -------------------------------------

test_that("alignReplicates() - previous result starts at its estimates", {
    first <- fitSimLog2()
    second <- fitSimLog2(data = first, iterlim = 0)
    expect_equal(second$parameter$value, first$parameter$value, tolerance = 1e-12)
    expect_false("1" %in% names(first$prediction))
})

test_that("llrTest() - checks nesting of the error effects", {
    h1 <- fitSimLog2()
    h0 <- fitSimLog2(error = sigmaA ~ name + experiment)
    expect_error(llrTest(h0, h1), "not a special case")
})


# alignReplicates() - grouped t likelihood ----------------------------------------

simNoisyGel <- function() {
    set.seed(1)
    d <- expand.grid(
        time = seq_len(8), Exp_ = paste0("E", 1:5), stringsAsFactors = FALSE
    )
    d$name <- "pX"
    d$value <- 2 * sin(d$time) + rnorm(5)[match(d$Exp_, unique(d$Exp_))] +
        rnorm(nrow(d), 0, ifelse(d$Exp_ == "E5", 1.5, 0.1))
    d
}

fitNoisyGel <- function(d, ...) {
    out <- suppressWarnings(alignReplicates(
        d, model = "yi + sj", errorModel = "exp(sigmaA)",
        biological = yi ~ name + time, scaling = sj ~ name + Exp_,
        error = sigmaA ~ name, normalizeInput = FALSE, fitLogscale = FALSE, ...
    ))
    out$aligned$value
}

test_that("alignReplicates() - tGroup down-weights a noisy gel as a whole", {
    d <- simNoisyGel()
    ref <- fitNoisyGel(d[d$Exp_ != "E5", ])
    dist <- function(y) sqrt(mean((y - ref)^2))
    distNormal <- dist(fitNoisyGel(d))
    distT <- dist(fitNoisyGel(d, likelihood = "t"))
    distGroup <- dist(fitNoisyGel(d, likelihood = "t", tGroup = "Exp_"))
    expect_lt(distGroup, distT)
    expect_lt(distGroup, distNormal)
})

test_that("alignReplicates() - tGroup of single points equals pointwise t", {
    d <- simNoisyGel()
    d$point <- as.character(seq_len(nrow(d)))
    expect_equal(
        fitNoisyGel(d, likelihood = "t", tGroup = "point"),
        fitNoisyGel(d, likelihood = "t"),
        tolerance = 1e-6
    )
    expect_error(fitNoisyGel(d, tGroup = "Exp_"), "requires likelihood")
    expect_error(fitNoisyGel(d, likelihood = "t", tGroup = "nope"), "columns")
})


# alignReplicates() - analytic inverse on the log scale -----------------------------

test_that("alignReplicates() - analytic and numeric inverse agree with fitLogscale", {
    set.seed(3)
    d <- expand.grid(time = seq_len(6), gel = paste0("G", 1:4), stringsAsFactors = FALSE)
    d$name <- "pX"
    d$value <- exp(sin(d$time) + rnorm(4)[match(d$gel, unique(d$gel))] +
        rnorm(nrow(d), 0, 0.2))
    out <- suppressWarnings(alignReplicates(
        d, model = "yi * sj", errorModel = "sigmaR * value",
        biological = yi ~ name + time, scaling = sj ~ name + gel,
        error = sigmaR ~ name, normalizeInput = FALSE, fitLogscale = TRUE
    ))
    p <- out$parameter
    sj <- p$value[grepl("^sj", p$parameter)][match(d$gel, unique(d$gel))]
    sigmaR <- p$value[grepl("^sigmaR", p$parameter)]
    expect_equal(out$scaled$value, d$value / sj, tolerance = 1e-10)
    expect_equal(out$scaled$sigma / out$prediction$sigma, 1 / sj, tolerance = 1e-10)
    expect_gt(sigmaR, 0)
})
