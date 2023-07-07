## Set seed.
set.seed(1234)

## load the data.
data(drops_GE)
data(drops_GnE)
data(drops_K)

## Remove identifiers that we don't need.
drops_GE_GnE <- rbind(drops_GE[, -c(2, 4)], drops_GnE[, -c(2, 4)])

## Restrict to 10 genotypes.
testDat <- drops_GE_GnE[drops_GE_GnE$Variety_ID %in%
                          levels(drops_GE_GnE$Variety_ID)[1:10], ]
testDat <- droplevels(testDat)

testK <- drops_K[levels(testDat$Variety_ID), levels(testDat$Variety_ID)]

## Restrict to 11 environments.
# At least 10 are needed for training to prevent removal of all genotypes.
# 1 used for testing.
testDat <- testDat[testDat$Experiment %in% levels(testDat$Experiment)[1:11], ]
testDat <- droplevels(testDat)

## Select indices and create indicesDat for testing.
indices <- c("Tnight.Early", "Tnight.Flo")

indicesDat <- testDat[testDat$Variety_ID == "A3",
                      c("Tnight.Early", "Tnight.Flo")]
rownames(indicesDat) <- substring(rownames(indicesDat), first = 1, last = 6)

## Create partionDat.
partitionDat <- data.frame(E = levels(testDat$Experiment), partition = 1:11)


### Check that input checks work as expected.
expect_error(GnE(dat = 1, Y = 1, G = 1, E = 1),
             "dat should be a data.frame")
expect_error(GnE(dat = testDat, Y = "a", G = 1, E = 1),
             "a should be a column in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "a", E = 1),
             "a should be a column in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID", E = "a"),
             "a should be a column in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment"),
             "Either indices or indicesData should be provided")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment"),
             "Either indices or indicesData should be provided")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = 1, indicesData = 1),
             "Either indices or indicesData should be provided")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = 1),
             "indices should be a vector of length > 1 of columns in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices[1]),
             "indices should be a vector of length > 1 of columns in dat")

## Check that columns can be refered to by names and by numbers.
modBase <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
               E = "Experiment", indices = indices)
modBaseNum <- GnE(dat = testDat, Y = 3, G = 2, E = 1, indices = indices)

expect_equal(modBase, modBaseNum)

## Check output structure.
expect_inherits(modBase, "list")
expect_equal(names(modBase),
             c("predTrain", "predTest", "resTrain", "resTest", "mu",
               "envInfoTrain", "envInfoTest", "parGeno", "trainAccuracyEnv",
               "testAccuracyEnv", "trainAccuracyGeno", "testAccuracyGeno",
               "lambda", "lambdaSequence", "RMSEtrain", "RMSEtest", "Y", "G",
               "E", "indices", "quadratic"))

## Check full output object.
expect_equal_to_reference(modBase, "modBaseGnE")

## Check that test environments can be added.

expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, testEnv = 1),
             "testEnv should be a vector of environments present in dat")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, testEnv = "bla"),
             "testEnv should be a vector of environments present in dat")


modTest <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
               E = "Experiment", indices = indices, testEnv = "Cam12R")

## Check full output object.
expect_equal_to_reference(modTest, "modTestGnE")


## Check that indicesData can be used.
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indicesData = 1),
             "indicesData should be a data.frame with all environments")
expect_error(GnE(dat = testDat[testDat$Experiment != "Cam12R", ],
                 Y = "grain.yield", G = "Variety_ID", E = "Experiment",
                 indicesData = indicesDat),
             "All environments in indicesData should be in dat")

expect_warning(modIndDat <- GnE(dat = testDat, Y = "grain.yield",
                                G = "Variety_ID", E = "Experiment",
                                indicesData = indicesDat),
               "The following columns in indicesDat are already in dat")

expect_equal(mean(modIndDat$trainAccuracyEnv$r), 0.674521689495499)


## Check that redundant genotypes are removed.
expect_warning(GnE(dat = testDat[!testDat$Experiment %in% c("Gai12W", "Gai13R"), ],
                   Y = "grain.yield", G = "Variety_ID",
                   E = "Experiment", indices = indices),
               "the following genotypes have < 10 observations, and are removed")
expect_error(GnE(dat = testDat[!testDat$Experiment %in% c("Gai12W", "Gai13R"), ],
                 Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices),
             "No data left after removing genotypes with < 10 observations")


## Check that specifying lambda works correctly.
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = "a"),
             "lambda should be a numeric vector")

modLambda <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda)
expect_equal(modLambda$lambda, modBase$lambda)
expect_equal(modLambda$lambdaSequence, modBase$lambda)


## Check that option quadratic works correctly.
modQuad <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
               E = "Experiment", indices = indices, lambda = modBase$lambda,
               quadratic = TRUE)
expect_equal(modQuad$indices, c("Tnight.Early", "Tnight.Flo",
                                "Tnight.Early_quad", "Tnight.Flo_quad"))
expect_equal(mean(modQuad$trainAccuracyEnv$r), 0.898854022315224)


## Check that option K works correctly.
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 K = 1),
             "K should be a matrix with all genotypes in dat in its row")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 K = testK[1:3, 1:3]),
             "K should be a matrix with all genotypes in dat in its row")

modK <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
            E = "Experiment", indices = indices, lambda = modBase$lambda,
            K = testK)
expect_equal(mean(modK$trainAccuracyEnv$r), 0.75100036563115)


## Check that option scale works correctly.
modAll <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
              E = "Experiment", indices = indices, lambda = modBase$lambda,
              scaling = "all")
expect_equal(mean(modAll$trainAccuracyEnv$r), 0.791322119880585)


## Check that option partition works correctly.
partitionDatNonnum <- partitionDat
partitionDatNonnum$partition <- as.character(partitionDatNonnum$partition)
partitionDatFewFolds <- partitionDat
partitionDatFewFolds$partition <- 1

expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 partition = 1),
             "partition should be a data.frame")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 partition = partitionDat[, 1, drop = FALSE]),
             "partition should have columns E and partition")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 partition = partitionDatNonnum),
             "Column partition in partition should be a numeric column with")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 partition = partitionDatFewFolds),
             "Column partition in partition should be a numeric column with")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 partition = partitionDat[-1, ]),
             "All training environments should be in partition")

modPart <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
               E = "Experiment", indices = indices, lambda = modBase$lambda,
               partition = partitionDat)

expect_equal(mean(modPart$trainAccuracyEnv$r), 0.791322119880585)

expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 partition = NULL, nfolds = "a"),
             "nfolds should be a numeric value of 4 or more")
expect_error(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                 E = "Experiment", indices = indices, lambda = modBase$lambda,
                 partition = NULL, nfolds = 3),
             "nfolds should be a numeric value of 4 or more")

modPartNULL <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
                   E = "Experiment", indices = indices, lambda = modBase$lambda,
                   partition = NULL, nfolds = 4)

expect_equal(mean(modPartNULL$trainAccuracyEnv$r), 0.791322119880585)


## Check that output can be written to file.

tmpFile <- tempfile()

modOut <- GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
              E = "Experiment", indices = indices, lambda = modBase$lambda,
              testEnv = "Cam12R", outputFile = tmpFile)

outTrain <- read.csv(paste0(tmpFile, "_perEnv_Train.csv"))
outTest <- read.csv(paste0(tmpFile, "_perEnv_Test.csv"))

expect_equal(modOut$trainAccuracyEnv, outTrain)
expect_equal(modOut$testAccuracyEnv, outTest)


## Check that output is written to console with verbose.

expect_stdout(GnE(dat = testDat, Y = "grain.yield", G = "Variety_ID",
               E = "Experiment", indices = indices, lambda = modBase$lambda,
               testEnv = "Cam12R", verbose = TRUE),
               "Training environments ( grain.yield )", fixed = TRUE)





