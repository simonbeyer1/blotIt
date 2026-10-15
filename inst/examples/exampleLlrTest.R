llrData <- readWide(
  system.file("extdata", "exampleLlrTest.csv", package = "blotIt"),
  description = seq_len(4)
)

## H0: the buffer column is not a biological effect
H0 <- alignReplicates(
  data = llrData,
  model = "yi / sj",
  errorModel = "value * sigmaR",
  biological = yi ~ name + time + stimmulus,
  scaling = sj ~ name + ID,
  error = sigmaR ~ name + 1,
  fitLogscale = FALSE
)

## H1: the buffer column is a biological effect
H1 <- alignReplicates(
  data = llrData,
  model = "yi / sj",
  errorModel = "value * sigmaR",
  biological = yi ~ name + time + stimmulus + buffer,
  scaling = sj ~ name + ID,
  error = sigmaR ~ name + 1,
  fitLogscale = FALSE
)

llrTest(H0, H1)
