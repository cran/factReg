#' @title
#' Genomic prediction using glmnet, with a genotype-specific penalized
#' regression model.
#'
#' @description
#' \loadmathjax
#' .... These models can be fitted either for the original
#' data, or on the residuals of a model with only main effects.
#'
#' @inheritParams GnE
#'
#' @param useRes Indicates whether the genotype-specific regressions are to be
#' fitted on the residuals of a model with main effects. If \code{TRUE}
#' residuals of a model with environmental main effects are used, if
#' \code{FALSE} the regressions are fitted on the original data.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{predTrain}{Vector with predictions for the training set (to do: Add
#'   the factors genotype and environment; make a data.frame)}
#'   \item{predTest}{Vector with predictions for the test set (to do: Add the
#'   factors genotype and environment; make a data.frame). To do: add estimated
#'   environmental main effects, not only predicted environmental main effects}
#'   \item{mu}{the estimated overall (grand) mean}
#'   \item{envInfoTrain}{The estimated environmental main effects, and the
#'   predicted effects, obtained when the former are regressed on the averaged
#'   indices, using penalized regression.}
#'   \item{envInfoTest}{The predicted environmental main effects for the test
#'   environments, obtained from penalized regression using the estimated
#'   main effects for the training environments and the averaged indices.}
#'   \item{parGeno}{data.frame containing the estimated genotypic
#'   main effects (first column) and sensitivities (subsequent columns)}
#'   \item{testAccuracyEnv}{a \code{data.frame} with the accuracy (r) for each
#'   test environment}
#'   \item{trainAccuracyEnv}{a \code{data.frame} with the accuracy (r) for each
#'   training environment}
#'   \item{trainAccuracyGeno}{a \code{data.frame} with the accuracy (r) for
#'   each genotype, averaged over the training environments}
#'   \item{testAccuracyGeno}{a \code{data.frame} with the accuracy (r) for
#'   each genotype, averaged over the test environments}
#'   \item{RMSEtrain}{The root mean squared error on the training environments}
#'   \item{RMSEtest}{The root mean squared error on the test environments}
#'   \item{Y}{The name of the trait that was predicted, i.e. the column name
#'   in dat that was used}
#'   \item{G}{The genotype label that was used, i.e. the argument G that was
#'   used}
#'   \item{E}{The environment label that was used, i.e. the argument E that was
#'   used}
#'   \item{indices}{The indices that were used, i.e. the argument indices that
#'   was used}
#'   \item{lambdaOpt}{}
#'   \item{pargeno}{}
#'   \item{quadratic}{The quadratic option that was used}
#' }
#'
#' @export
perGeno <- function(dat,
                    Y,
                    G,
                    E,
                    indices = NULL,
                    indicesData = NULL,
                    testEnv = NULL,
                    weight = NULL,
                    useRes = TRUE,
                    outputFile = NULL,
                    corType = c("pearson", "spearman"),
                    partition = data.frame(),
                    nfolds = 10,
                    alpha = 1,
                    scaling = c("train", "all", "no"),
                    quadratic = FALSE,
                    verbose = FALSE) {
  ## Input checks.
  if (!inherits(dat, "data.frame")) {
    stop("dat should be a data.frame.\n")
  }
  ## Get column names.
  traitName <- if (is.numeric(Y)) names(dat)[Y] else Y
  genoName <- if (is.numeric(G)) names(dat)[G] else G
  envName <- if (is.numeric(E)) names(dat)[E] else E
  ## Rename data columns for Y, G and E.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  scaling <- match.arg(scaling)
  corType <- match.arg(corType)
  ## Either indices or indicesData should be provided.
  if ((is.null(indices) && is.null(indicesData)) ||
      (!is.null(indices) && !is.null(indicesData))) {
    stop("Either indices or indicesData should be provided.\n")
  }
  if (!is.null(indices) && (!is.character(indices) ||
                            length(indices) <= 1 ||
                            !all(hasName(x = dat, name = indices)))) {
    stop("indices should be a vector of length > 1 of columns in dat.\n")
  }
  if (!is.null(indicesData)) {
    if (!inherits(indicesData, "data.frame") ||
        !all(levels(dat$E) %in% rownames(indicesData))) {
      stop("indicesData should be a data.frame with all environments in its ",
           "rownames.\n")
    }
    if (!all(rownames(indicesData) %in% levels(dat$E))) {
      stop("All environments in indicesData should be in dat.\n")
    }
    presCols <- colnames(indicesData)[colnames(indicesData) %in%
                                        colnames(dat)]
    if (length(presCols) > 0) {
      warning("The following columns in indicesDat are already in dat. Values ",
              "in dat will be overwritten:\n",
              paste(presCols, collapse = ", "), ".\n")
    }
  }
  ## Check testEnv.
  if (!is.null(testEnv) && (!is.character(testEnv) || length(testEnv) < 1 ||
                            !all(testEnv %in% levels(dat$E)))) {
    stop("testEnv should be a vector of environments present in dat.\n")
  }
  ## Get training envs.
  trainEnv <- setdiff(levels(dat$E), testEnv)
  ## Remove missing values from training envs.
  dat <- dat[!(is.na(dat$Y) & dat$E %in% trainEnv), ]
  redundantGeno <- names(which(table(dat$G[!is.na(dat$Y) &
                                             dat$E %in% trainEnv]) < 10))
  if (length(redundantGeno) > 0) {
    warning("the following genotypes have < 10 observations, and are removed:",
            "\n\n", paste(redundantGeno, collapse = ", "), "\n")
    dat <- dat[!(dat$G %in% redundantGeno), ]
    if (nrow(dat) == 0) {
      stop("No data left after removing genotypes with < 10 observations.\n")
    }
  }
  if (!is.null(indicesData)) {
    indices <- colnames(indicesData)
    ## Remove columns from dat that are also in indices and then merge indices.
    dat <- dat[, setdiff(names(dat), indices)]
    dat <- merge(dat, indicesData, by.x = "E", by.y = "row.names")
  }
  dat <- droplevels(dat)
  ## Scale environmental variables.
  if (scaling == "train") {
    muTr <- colMeans(dat[dat$E %in% trainEnv, indices])
    sdTr <- sapply(X = dat[dat$E %in% trainEnv, indices], sd)
    dat[, indices] <- scale(dat[, indices], center = muTr, scale = sdTr)
  } else if (scaling == "all") {
    dat[, indices] <- scale(dat[, indices])
  }
  if (quadratic) {
    ## Add quadratic environmental columns to data.
    datIndQuad <- dat[, indices]^2
    colnames(datIndQuad) <- paste0(indices, "_quad")
    dat <- cbind(dat, datIndQuad)
    ## Add the quadratic columns to the indices.
    indices <- c(indices, paste0(indices, "_quad"))
  }
  if (is.null(weight)) {
    dat$W <- 1
  } else {
    stopifnot(length(weight) == nrow(dat))
    dat$W <- weight
  }
  ## Split dat into training and test set.
  dTrain <- dat[dat$E %in% trainEnv, ]
  #dTrainNA <- .....
  #dTrain <- dat[dat$E %in% trainEnv & !is.na(dat$Y), ]
  dTrain <- droplevels(dTrain)
  if (!is.null(testEnv)) {
    dTest <- dat[dat$E %in% testEnv, ]
    dTest <- droplevels(dTest)
    nGenoTest <- nlevels(dTest$G)
  }
  nEnvTest <- length(testEnv)
  nEnvTrain <- nlevels(dat$E) - nEnvTest
  nGenoTrain <- nlevels(dTrain$G)
  ## When partition == data.frame() (the default),
  ## do leave one environment out cross-validation.
  if (!is.null(partition)) {
    if (!inherits(partition, "data.frame")) {
      stop("partition should be a data.frame.\n")
    }
    if (ncol(partition) == 0) {
      partition <- unique(data.frame(E = dTrain$E,
                                     partition = as.numeric(dTrain$E)))
    } else {
      if (!setequal(colnames(partition), c("E", "partition"))) {
        stop("partition should have columns E and partition.\n")
      }
      if (!is.numeric(partition$partition) ||
          length(unique(partition$partition)) < 4) {
        stop("Column partition in partition should be a numeric column with",
             "at least 4 different values.\n")
      }
      if (!all(levels(dTrain$E) %in% partition$E)) {
        stop("All training environments should be in partition.\n")
      }
      partition <- partition[partition$E %in% trainEnv, ]
    }
    ## Construct foldid from partition
    foldDat <- merge(dTrain, partition)
    foldid <- foldDat$partition
    nfolds <- length(unique(foldid))
  } else {
    foldid <- NULL
    if (!is.numeric(nfolds) || length(nfolds) > 1 || nfolds < 4) {
      stop("nfolds should be a numeric value of 4 or more.\n")
    }
  }
  mm <- Matrix::sparse.model.matrix(Y ~ E + G, data = dTrain)
  modelMain <- glmnet::glmnet(y = dTrain$Y, x = mm, thresh = 1e-18, lambda = 0)
  cf <- c(modelMain$a0, modelMain$beta[-1])
  names(cf) <- c("(Intercept)", paste0("E", levels(dTrain$E)[-1]),
                 paste0("G", levels(dTrain$G)[-1]))
  ## Even if useRes is FALSE, genoMain will still be used for comparison with a
  ## main effects only model
  genoMain <- setNames(c(0, cf[(nEnvTrain + 1):(nEnvTrain + nGenoTrain - 1)]),
                       levels(dTrain$G))
  ## Extract from the output the estimated environmental main effects
  envMain <- matrix(c(0, cf[2:nEnvTrain]), ncol = 1,
                    dimnames = list(levels(dTrain$E)))
  dTrain$yRes <- dTrain$Y -
    as.numeric(predict(modelMain, newx = mm, thresh = 1e-18, lambda = 0))
  ## Define the design matrix for the factorial regression model.
  geFormulaTrain <- as.formula(paste0("yRes ~ -1 +",
                                      paste(paste0(indices, ":G"),
                                            collapse = " + ")))
  ## Construct design matrix for training set.
  mTrain <- Matrix::sparse.model.matrix(geFormulaTrain, data = dTrain)
  if (!is.null(testEnv)) {
    ## Construct design matrix for test set.
    geFormulaTest <- update(geFormulaTrain, new = NULL ~ .)
    mTest <- Matrix::sparse.model.matrix(geFormulaTest, data = dTest)
  }
  predTrain <- rep(NA, nrow(dTrain))
  names(predTrain) <- rownames(dTrain)
  if (!is.null(testEnv)) {
    predTest <- rep(NA, nrow(dTest))
    names(predTest) <- rownames(dTest)
  } else {
    predTest <- NULL
  }
  nIndices <- length(indices)
  ## Create a data.frame that will contain the estimated genotypic
  ## main effects (first column), and the estimated environmental
  ## sensitivities (subsequent) columns.
  parGeno <- matrix(nrow = nGenoTrain, ncol = nIndices + 1,
                    dimnames = list(levels(dTrain$G), c("main", indices)))
  parGeno <- as.data.frame(parGeno)
  lambdaOpt <- setNames(rep(NA, nGenoTrain), levels(dTrain$G))
  for (gg in levels(dTrain$G)) {
    ggInd <- which(levels(dTrain$G) == gg)
    cn <- ggInd + (0:(nIndices - 1)) * nGenoTrain
    scn <- colnames(mTrain)[cn]
    gnInd <- which(dTrain$G == gg)
    mgg <- mTrain[gnInd, scn, drop = FALSE]
    mTrainGeno <- scale(mgg)
    if (useRes) {
      yTemp <- dTrain$yRes[gnInd]
    } else {
      yTemp <- dTrain$Y[gnInd]
    }
    grouped <- if (!is.null(foldid)) max(table(foldid[gnInd])) > 2 else
      length(yTemp) / nfolds >= 3
    glmnetOut <-
      glmnet::cv.glmnet(x = mTrainGeno,
                        y = yTemp,
                        weights = dTrain$W[gnInd],
                        foldid = if (!is.null(foldid)) seq_along(foldid[gnInd]),
                        nfolds = nfolds,
                        standardize = TRUE,
                        intercept = TRUE,
                        alpha = alpha,
                        grouped = grouped)
    lambda <- glmnetOut$lambda
    lambdaIndex <- which(lambda == glmnetOut$lambda.min)
    cfe <- as.numeric(glmnetOut$glmnet.fit$beta[, lambdaIndex])
    mu <- as.numeric(glmnetOut$glmnet.fit$a0[lambdaIndex])
    parGeno[ggInd, ] <- c(mu, cfe)
    lambdaOpt[ggInd] <- glmnetOut$lambda.min
    predTrain[gnInd] <- as.numeric(predict(object = glmnetOut,
                                           newx = mTrainGeno,
                                           s = "lambda.min"))
    if (useRes) {
      predTrain[gnInd] <- predTrain[gnInd] + genoMain[gg]
    }
    ## Predictions for test set.
    if (!is.null(testEnv) && gg %in% levels(dTest$G)) {
      mugg <- Matrix::colMeans(mgg)
      sdgg <- apply(X = mgg, MARGIN = 2, FUN = sd)
      mggTest <- mTest[dTest$G == gg, scn, drop = FALSE]
      mTest[dTest$G == gg, scn] <- scale(mggTest, center = mugg, scale = sdgg)
      glmnetPred <- as.numeric(predict(object = glmnetOut,
                                       newx = mTest[dTest$G == gg,
                                                    scn, drop = FALSE],
                                       s = "lambda.min"))
      predTest[dTest$G == gg] <- glmnetPred
      if (useRes) {
        predTest[dTest$G == gg] <- predTest[dTest$G == gg] + genoMain[gg]
      }
    }
  }
  ## Compute the mean of each environmental index, in each environment
  indFrame <- aggregate(dat[, indices], by = list(E = dat$E), FUN = mean)
  indFrameTrain <- indFrame[indFrame$E %in% trainEnv, ]
  indFrameTrain <- merge(indFrameTrain, envMain, by.x = "E", by.y = "row.names")
  colnames(indFrameTrain)[ncol(indFrameTrain)] <- "envMainFitted"
  if (!is.null(partition)) {
    indFrameTrain <- merge(indFrameTrain, partition)
  }
  rownames(indFrameTrain) <- indFrameTrain$E
  glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                 y = indFrameTrain$envMainFitted,
                                 alpha = alpha,
                                 foldid = indFrameTrain$partition,
                                 grouped = nrow(indFrameTrain) / nfolds >= 3)
  parEnvTrain <- predict(object = glmnetOut,
                         newx = as.matrix(indFrameTrain[, indices]),
                         s = "lambda.min")
  if (!is.null(testEnv)) {
    indFrameTest <- indFrame[indFrame$E %in% testEnv, ]
    rownames(indFrameTest) <- indFrameTest$E
    parEnvTest  <- predict(object = glmnetOut,
                           newx = as.matrix(indFrameTest[, indices]),
                           s = "lambda.min")
  }
  ## add the predicted environmental main effects:
  predTrain <- predTrain + as.numeric(parEnvTrain[as.character(dTrain$E), ])
  indicesTrain <- data.frame(indFrameTrain,
                             envMainPred = as.numeric(parEnvTrain))
  if (!is.null(testEnv)) {
    ## add the predicted environmental main effects:
    predTest <- predTest + as.numeric(parEnvTest[as.character(dTest$E), ])
    indicesTest <- data.frame(indFrameTest,
                              envMainPred = as.numeric(parEnvTest))
  }
  ## Compute statistics for training data.
  predMain <- genoMain[as.character(dTrain$G)]
  trainAccuracyEnv <- getAccuracyEnv(datNew = dTrain[, "Y"],
                                     datPred = predTrain,
                                     datPredMain = predMain,
                                     datE = dTrain[, "E"],
                                     corType = corType,
                                     rank = TRUE)
  ## Compute accuracies for genotypes.
  trainAccuracyGeno <- getAccuracyGeno(datNew = dTrain[, "Y"],
                                       datPred = predTrain,
                                       datG = dTrain[, "G"],
                                       corType = corType)
  if (!is.null(testEnv)) {
    ## Compute statistics for test data.
    predMainTest <- genoMain[as.character(dTest$G)]
    testAccuracyEnv <- getAccuracyEnv(datNew = dTest[, "Y"],
                                      datPred = predTest,
                                      datPredMain = predMainTest,
                                      datE = dTest[, "E"],
                                      corType = corType,
                                      rank = TRUE)
    ## Compute accuracies for genotypes.
    testAccuracyGeno <- getAccuracyGeno(datNew = dTest[, "Y"],
                                        datPred = predTest,
                                        datG = dTest[, "G"],
                                        corType = corType)
  } else {
    indicesTest <- NULL
    testAccuracyEnv  <- NULL
    testAccuracyGeno <- NULL
  }
  ## Create RMSE.
  RMSEtrain <- sqrt(mean((dTrain$Y - predTrain) ^ 2, na.rm = TRUE))
  if (!is.null(testEnv)) {
    RMSEtest  <- sqrt(mean((dTest$Y - predTest) ^ 2, na.rm = TRUE))
  } else {
    RMSEtest <- NULL
  }
  if (!is.null(outputFile)) {
    ## Write output to csv.
    resultFilePerEnvTrain <- paste0(outputFile, "_perEnv_Train.csv")
    write.csv(trainAccuracyEnv, file = resultFilePerEnvTrain,
              row.names = FALSE, quote = FALSE)
    if (!is.null(testEnv)) {
      resultFilePerEnvTest <- paste0(outputFile, "_perEnv_Test.csv")
      write.csv(testAccuracyEnv, file = resultFilePerEnvTest,
                row.names = FALSE, quote = FALSE)
    }
  }
  if (verbose) {
    ## Print output to console.
    if (!is.null(testEnv)) {
      cat("\n\n", "Test environments (", traitName, ")", "\n\n")
      print(format(testAccuracyEnv, digits = 2, nsmall = 2))
    }
    cat("\n\n", "Training environments (", traitName,")", "\n\n")
    print(format(trainAccuracyEnv, digits = 2, nsmall = 2))
  }
  ## Create output object.
  out <- list(predTrain = predTrain,
              predTest = predTest,
              envInfoTrain = indicesTrain,
              envInfoTest  = indicesTest,
              testAccuracyEnv = testAccuracyEnv,
              trainAccuracyEnv = trainAccuracyEnv,
              trainAccuracyGeno = trainAccuracyGeno,
              testAccuracyGeno = testAccuracyGeno,
              RMSEtrain = RMSEtrain,
              RMSEtest = RMSEtest,
              Y = traitName,
              G = genoName,
              E = envName,
              indices = indices,
              genotypes = levels(dTrain$G),
              lambdaOpt = lambdaOpt,
              parGeno = parGeno,
              quadratic = quadratic)
  return(out)
}
