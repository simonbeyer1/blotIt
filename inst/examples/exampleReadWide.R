simDataWideFile <- system.file("extdata", "simDataWide.csv", package = "blotIt")
readWide(simDataWideFile, description = seq_len(3))
