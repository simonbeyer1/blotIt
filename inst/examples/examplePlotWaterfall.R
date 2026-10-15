exampleData <- readWide(
  system.file("extdata", "INFdata_WB_wide.csv", package = "blotIt"),
  description = seq_len(3)
)
exampleData <- exampleData[exampleData$name == exampleData$name[1], ]

out <- alignReplicates(
  data = exampleData,
  model = "yi / sj",
  errorModel = "value * sigmaR",
  biological = yi ~ name + time + condition,
  scaling = sj ~ name + experiment,
  error = sigmaR ~ name,
  multistart = 10
)
plotWaterfall(out)
