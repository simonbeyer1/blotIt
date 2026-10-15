exampleData <- readWide(
  system.file("extdata", "INFdata_WB_wide.csv", package = "blotIt"),
  description = seq_len(3)
)

out <- alignReplicates(
  data = exampleData,
  model = "yi / sj",
  errorModel = "value * sigmaR",
  biological = yi ~ name + time + condition,
  scaling = sj ~ name + experiment,
  error = sigmaR ~ name
)

## Original data with the model prediction
plotIt(out, plotPoints = "original", plotLine = "prediction")

## Data scaled to the common scale with the aligned values
plotIt(out, plotPoints = "scaled", plotLine = "aligned")

## Aligned values only
plotIt(out, plotPoints = "aligned", plotLine = "aligned")
