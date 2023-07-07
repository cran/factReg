#' @title
#' Penalized factorial regression using glmnet
#'
#' @description
#' \loadmathjax
#' Based on multi-environment field trials, fits the factorial regression model
#' \mjeqn{Y_{ij} = \mu + e_j + g_i + \sum_{k=1}^s \beta_{ik} x_{ij} +
#' \epsilon_{ij},}{ascii} with environmental main effects \mjeqn{e_j}{ascii},
#' genotypic main effects \mjeqn{g_{i}}{ascii} and genotype-specific
#' environmental sensitivities \mjeqn{\beta_{ik}}{ascii}. See e.g. Millet
#' et al 2019 and Bustos-Korts et al 2019. There are \mjeqn{s}{ascii}
#' environmental indices with values \mjeqn{x_{ij}}{ascii}. Optionally,
#' predictions can be made for a set of test environments, for which
#' environmental indices are available. The new environments must contain the
#' same set of genotypes, or a subset.
#'
#' Penalization: the model above is fitted using glmnet, simultaneously
#' penalizing \mjeqn{e_j}{ascii}, \mjeqn{g_i}{ascii} and
#' \mjeqn{\beta_{ik}}{ascii}. If \code{penG = 0} and \code{penE = 0}, the main
#' effects \mjeqn{g_i}{ascii} and \mjeqn{e_j}{ascii} are not penalized. If these
#' parameters are 1, the the main effects are penalized to the same degree as
#' the sensitivities. Any non negative values are allowed. Cross validation is
#' performed with each fold containing a number of environments (details below).
#'
#' After fitting the model, it is possible to replace the estimated
#' genotypic main effects and sensitivities by
#' their predicted genetic values. Specifically, if a kinship matrix \code{K}
#' is assigned the function performs genomic prediction with g-BLUP for the
#' genotypic main effect and each of the sensitivities in turn.
#'
#' Predictions for the test environments are first constructed using the
#' estimated genotypic main effects and sensitivities; next, predicted
#' environmental main effects are added. The latter are obtained by
#' regressing the estimated environmental main effects for the training
#' environments on the average values of the indices in these environments,
#' as in Millet et al. 2019.
#'
#' @param dat A \code{data.frame} with data from multi-environment trials.
#' Each row corresponds to a particular genotype in a particular environment.
#' The data do not need to be balanced, i.e. an environment does not need to
#' contain all genotypes. \code{dat} should contain the training as well as the
#' test environments (see testEnv)
#' @param Y The trait to be analyzed: either of type character, in which case
#' it should be one of the column names in \code{dat}, or numeric, in which
#' case the Yth column of \code{dat} will be analyzed.
#' @param G The column in \code{dat} containing the factor genotype (either
#' character or numeric).
#' @param E The column in \code{dat} containing the factor environment
#' (either character or numeric).
#' @param K A kinship matrix. Used for replacing the estimated genotypic main
#' effect and each of the sensitivities by genomic prediction from a g-BLUP
#' model for each of them. If \code{NULL}, the estimated effects from the model
#' are returned and used for constructing predictions.
#' @param indices The columns in \code{dat} containing the environmental
#' indices (vector of type character). Alternatively, if the indices are always
#' constant within environments (i.e. not genotype dependent), the
#' environmental data can also be provided using the argument \code{indicesData}
#' (see below).
#' @param indicesData An optional \code{data.frame} containing environmental
#' indices (covariates); one value for each environment and index. It should
#' have the environment names as row names (corresponding to the names
#' contained in \code{dat$E}); the column names are the indices. If
#' \code{indices} (see before) is also provided, the latter will be ignored.
#' @param testEnv vector (character). Data from these environments are not used
#' for fitting the model. Accuracy is evaluated for training and test
#' environments separately. The default is \code{NULL}, i.e. no test
#' environments, in which case the whole data set is training. It is also
#' possible that there are test environments, but without any data; in this
#' case, no accuracy is reported for test environments (CHECK correctness).
#' @param weight Numeric vector of length \code{nrow(dat)}, specifying the
#' weight (inverse variance) of each observation, used in glmnet. Default
#' \code{NULL}, giving constant weights.
#' @param outputFile The file name of the output files, without .csv extension
#' which is added by the function. If not \code{NULL} the prediction accuracies
#' for training and test environments are written to separate files. If
#' \code{NULL} the output is not written to file.
#' @param corType type of correlation: Pearson (default) or Spearman rank sum.
#' @param partition \code{data.frame} with columns E and partition. The column
#' E should contain the training environments (type character); partition
#' should be of type integer. Environments in the same fold should have
#' the same integer value. Default is \code{data.frame()}, in which case the
#' function uses a leave-one-environment out cross-validation. If \code{NULL},
#' the (inner) training sets used for cross-validation will be drawn randomly
#' from all observations, ignoring the environment structure. In the latter
#' case, the number of folds (nfolds) can be specified.
#' @param nfolds Default \code{NULL}. If \code{partition == NULL}, this can be
#' used to specify the number of folds to be used in glmnet.
#' @param alpha Type of penalty, as in glmnet (1 = LASSO, 0 = ridge; in between
#'  = elastic net). Default is 1.
#' @param lambda Numeric vector; defines the grid over which the penalty lambda
#' is optimized in cross validation. Default: NULL (defined by glmnet).
#' Important special case: lambda = 0 (no penalty).
#' @param penG numeric; default 0. If 1, genotypic main effects are
#' penalized. If 0, they are not. Any non negative real number is allowed.
#' @param penE numeric; default 0. If 1, environmental main effects are
#' penalized. If 0, they are not. Any non negative real number is allowed.
#' @param scaling determines how the environmental variables are scaled.
#' "train" : all data (test and training environments) are scaled
#' using the mean and and standard deviation in the training environments.
#' "all" : using the mean and standard deviation of all environments.
#' "no" : No scaling.
#' @param quadratic boolean; default \code{FALSE}. If \code{TRUE}, quadratic
#' terms (i.e., squared indices) are added to the model.
#' @param verbose boolean; default \code{FALSE}. If \code{TRUE}, the accuracies
#' per environment are printed on screen.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{predTrain}{A data.frame with predictions for the training set}
#'   \item{predTest}{A data.frame with predictions for the test set}
#'   \item{resTrain}{A data.frame with residuals for the training set}
#'   \item{resTest}{A data.frame with residuals for the test set}
#'   \item{mu}{the estimated overall mean}
#'   \item{envInfoTrain}{The estimated environmental main effects, and the
#'   predicted effects, obtained when the former are regressed on the averaged
#'   indices, using penalized regression}
#'   \item{envInfoTest}{The predicted environmental main effects for the test
#'   environments, obtained from penalized regression using the estimated
#'   main effects for the training environments and the averaged indices}
#'   \item{parGeno}{data.frame containing the estimated genotypic main effects
#'   (first column) and sensitivities (subsequent columns)}
#'   \item{trainAccuracyEnv}{a data.frame with the accuracy (r) for each
#'   training environment, as well as the root mean square error (RMSE), mean
#'   absolute deviation (MAD) and rank (the latter is a proportion: how many
#'   of the best 5 genotypes are in the top 10). To be removed or further
#'   developed. All these quantities are also evaluated for a model with only
#'   genotypic and environmental main effects (columns rMain, RMSEMain and
#'   rankMain)}
#'   \item{testAccuracyEnv}{A data.frame with the accuracy for each test
#'   environment, with the same columns as trainAccuracyEnv}
#'   \item{trainAccuracyGeno}{a data.frame with the accuracy (r) for each
#'   genotype, averaged over the training environments}
#'   \item{testAccuracyGeno}{a data.frame with the accuracy (r) for each
#'   genotype, averaged over the test environments}
#'   \item{lambda}{The value of lambda selected using cross validation}
#'   \item{lambdaSequence}{The values of lambda used in the fits of glmnet. If
#'   \code{lambda} was provided as input, the value of \code{lambda} is
#'   returned}
#'   \item{RMSEtrain}{The root mean squared error on the training environments}
#'   \item{RMSEtest}{The root mean squared error on the test environments}
#'   \item{Y}{The name of the trait that was predicted, i.e. the column name
#'   in \code{dat} that was used}
#'   \item{G}{The genotype label that was used, i.e. the argument G that was
#'   used}
#'   \item{E}{The environment label that was used, i.e. the argument E that
#'   was used}
#'   \item{indices}{The indices that were used, i.e. the argument indices that
#'   was used}
#'   \item{quadratic}{The quadratic option that was used}
#' }
#'
#' @examples
#' ## load the data, which are contained in the package
#' data(drops_GE)
#' data(drops_GnE)
#'
#' ## We remove identifiers that we don't need.
#' drops_GE_GnE <- rbind(drops_GE[, -c(2, 4)], drops_GnE[, -c(2, 4)])
#'
#' ## Define indeces.
#' ind <- colnames(drops_GE)[13:23]
#'
#' ## Define test environments.
#' testenv <- levels(drops_GnE$Experiment)
#'
#' ## Additive model, only main effects (set the penalty parameter to a large value).
#' Additive_model <- GnE(drops_GE_GnE, Y = "grain.yield", lambda = 100000,
#'                       G = "Variety_ID", E = "Experiment", testEnv = testenv,
#'                       indices = ind, penG = FALSE, penE = FALSE,
#'                       alpha = 0.5, scaling = "train")
#' \donttest{
#' ## Full model, no penalization (set the penalty parameter to zero).
#' Full_model <- GnE(drops_GE_GnE, Y = "grain.yield", lambda = 0,
#'                   G = "Variety_ID", E = "Experiment", testEnv = testenv,
#'                   indices = ind, penG = FALSE, penE = FALSE,
#'                   alpha = 0.5, scaling = "train")
#'
#' ## Elastic Net model, set alpha parameter to 0.5.
#' Elnet_model <- GnE(drops_GE_GnE, Y = "grain.yield", lambda = NULL,
#'                    G = "Variety_ID", E = "Experiment", testEnv = testenv,
#'                    indices = ind, penG = FALSE, penE = FALSE,
#'                    alpha = 0.5, scaling = "train")
#'
#' ## Lasso model, set alpha parameter to 1.
#' Lasso_model <- GnE(drops_GE_GnE, Y = "grain.yield", lambda = NULL,
#'                    G = "Variety_ID", E = "Experiment", testEnv = testenv,
#'                    indices = ind, penG = FALSE, penE = FALSE,
#'                    alpha = 1, scaling = "train")
#'
#' ## Ridge model, set alpha parameter to 0.
#' Ridge_model <- GnE(drops_GE_GnE, Y = "grain.yield", lambda = NULL,
#'                    G = "Variety_ID", E = "Experiment", testEnv = testenv,
#'                    indices = ind, penG = FALSE, penE = FALSE,
#'                    alpha = 0, scaling = "train")
#' }
#'
#' @references Millet, E.J., Kruijer, W., Coupel-Ledru, A. et al. Genomic
#' prediction of maize yield across European environmental conditions. Nat
#' Genet 51, 952–956 (2019). \doi{10.1038/s41588-019-0414-y}
#'
#' @export
GnE <- function(dat,
                Y,
                G,
                E,
                K = NULL,
                indices = NULL,
                indicesData = NULL,
                testEnv = NULL,
                weight = NULL,
                outputFile = NULL,
                corType = c("pearson", "spearman"),
                partition = data.frame(),
                nfolds = 10,
                alpha = 1,
                lambda = NULL,
                penG = 0,
                penE = 0,
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
  ## Rename data columns for Y, G and E. - this includes checks.
  dat <- renameYGE(dat = dat, Y = Y, G = G, E = E)
  stopifnot(penG >= 0) # also > 1 is possible!
  stopifnot(penE >= 0)
  scaling <- match.arg(scaling)
  corType <- match.arg(corType)
  if (is.null(lambda)) {
    lambdaProvided <- FALSE
  } else {
    if (!is.numeric(lambda)) {
      stop("lambda should be a numeric vector.\n")
    }
    lambdaProvided <- TRUE
    lambda <- sort(lambda, decreasing = TRUE)
  }
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
  ## Check kinship.
  if (!is.null(K)) {
    if (!is.matrix(K) || !setequal(levels(dat$G), colnames(K)) ||
        !setequal(levels(dat$G), rownames(K))) {
      stop("K should be a matrix with all genotypes in dat in its row and ",
           "column names.\n")
    }
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
    ## Scaling when sd = 0 is not possible. Use a very small value instead.
    sdTr <- pmax(sapply(X = dat[dat$E %in% trainEnv, indices], sd), 1e-10)
    dat[, indices] <- scale(dat[, indices], center = muTr, scale = sdTr)
  } else if (scaling == "all") {
    ## Scaling when sd = 0 is not possible. Use a very small value instead.
    sdFull <- pmax(sapply(X = dat[, indices], sd), 1e-10)
    dat[, indices] <- scale(dat[, indices], scale = sdFull)
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
  dTrain <- droplevels(dTrain)
  if (!is.null(testEnv)) {
    dTest <- dat[dat$E %in% testEnv, ]
    dTest <- droplevels(dTest)
    nGenoTest <- nlevels(dTest$G)
  } else{
    dTest <- NULL
  }
  nEnv <- nlevels(dat$E)
  nEnvTest <- length(testEnv)
  nEnvTrain <- nEnv - nEnvTest
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
  mainOnly <- rep(0, nGenoTrain)
  names(mainOnly) <- levels(dTrain$G)
  mainTemp <- cf[(nEnvTrain + 1):(nEnvTrain + nGenoTrain - 1)]
  names(mainTemp) <- substring(names(mainTemp), first = 2)
  mainOnly[names(mainTemp)] <- mainTemp
  ## Define the design matrix for the factorial regression model.
  geFormula <- as.formula(paste0("Y ~ -1 + E + G +",
                                 paste(paste0(indices, ":G"),
                                       collapse = " + ")))
  ## Construct design matrix for training + test set.
  # important when there are NA's in the test set.
  opts <- options(na.action = "na.pass")
  on.exit(options(opts), add = TRUE)
  ma <- Matrix::sparse.model.matrix(geFormula, data = rbind(dTrain, dTest))
  colnames(ma)[1:(nEnv + nGenoTrain - 1)] <-
    substring(colnames(ma)[1:(nEnv + nGenoTrain - 1)], first = 2)
  # define the vector indicating to which extent parameters are to be penalized.
  penaltyFactorA <- rep(1, ncol(ma))
  penaltyFactorA[1:nEnv] <- penE
  penaltyFactorA[(nEnv + 1):(nEnv + nGenoTrain - 1)] <- penG
  # note: even if unpenalized, the estimated main effects change,
  # depending on lambda!
  # run glmnet, either (if provided) for a single value of lambda
  # (using glmnet), or using cv.glmnet
  thr <- 1e-07
  if (lambdaProvided && length(lambda) == 1) {
    if (lambda == 0) {
      thr <- 1e-11
    }
    glmnetOutA <- glmnet::glmnet(x = ma[1:nrow(dTrain),],
                                 y = dTrain$Y,
                                 lambda = lambda,
                                 weights = dTrain$W,
                                 alpha = alpha, standardize = TRUE,
                                 penalty.factor = penaltyFactorA,
                                 intercept = TRUE,
                                 thres = thr)
    lambdaIndex <- 1
    lambdaSequence <- lambda
    cfe <- glmnetOutA$beta
    mu <- as.numeric(glmnetOutA$a0)
  } else {
    glmnetOutA <- glmnet::cv.glmnet(x = ma[1:nrow(dTrain), ], y = dTrain$Y,
                                    lambda = lambda,
                                    weights = dTrain$W,
                                    foldid = foldid, nfolds = nfolds,
                                    alpha = alpha, standardize = TRUE,
                                    penalty.factor = penaltyFactorA,
                                    grouped = nrow(dTrain) / nfolds >= 3,
                                    intercept = TRUE)
    lambda <- glmnetOutA$lambda
    lambdaSequence <- lambda
    lambdaIndex <- which(lambda == glmnetOutA$lambda.min)
    cfe <- glmnetOutA$glmnet.fit$beta
    mu <- glmnetOutA$glmnet.fit$a0[lambdaIndex]
  }
  ## Extract from the output the estimated environmental main effects
  envMain <- as.matrix(cfe[1:nEnvTrain, lambdaIndex, drop = FALSE])
  envMain2 <- envMain
  envMain2[, 1] <- cf[1:nEnvTrain]
  ## Create a data.frame that will contain the estimated genotypic
  ## main effects (first column), and the estimated environmental
  ## sensitivities (subsequent) columns.
  tempInd2 <- (nEnv + nGenoTrain):ncol(ma)
  ## assign the genotypic main effects
  parGeno <-
    data.frame(main = c(0, cfe[(nEnv + 1): (nGenoTrain + nEnv - 1),
                               lambdaIndex]), row.names = levels(dTrain$G))
  ## assign the genotype specific sensitivities (subsequent columns)
  cfeIndRows <- rownames(cfe[tempInd2, , drop = FALSE])
  cfeGenotypes <- sapply(X = cfeIndRows, FUN = function(cfeIndRow) {
    substring(strsplit(x = cfeIndRow, split = ":")[[1]][1], first = 2)
  })
  cfeIndices <- sapply(X = cfeIndRows, FUN = function(cfeIndRow) {
    strsplit(x = cfeIndRow, split = ":")[[1]][2]
  })
  cfeDf <- data.frame(geno = cfeGenotypes, index = cfeIndices,
                      val = cfe[tempInd2, lambdaIndex],
                      stringsAsFactors = FALSE)
  cfeDf$geno <- factor(cfeDf$geno, levels = rownames(parGeno))
  cfeDf$index <- factor(cfeDf$index, levels = indices)
  parGenoIndices <- reshape(cfeDf, direction = "wide", timevar = "index",
                            idvar = "geno")
  colnames(parGenoIndices)[-1] <- levels(cfeDf$index)
  parGeno <- merge(parGeno, parGenoIndices, by.x = "row.names", by.y = "geno")
  rownames(parGeno) <- parGeno[["Row.names"]]
  parGeno <- parGeno[, c("main", indices)]
  ## If kinship provided replace by kinship predictions.
  if (!is.null(K)) {
    parGenoDat <- parGeno
    parGenoDat$G <- factor(rownames(parGenoDat))
    K <- K[parGenoDat$G, parGenoDat$G]
    for (j in 1:(length(indices) + 1)) {
      kinPred <- rrBLUP::kin.blup(data = parGenoDat, geno = "G",
                                  pheno = names(parGenoDat)[j], K = K)
      parGeno[names((kinPred$pred)), j] <- as.numeric(kinPred$pred)
    }
  }
  ## Make predictions for training set.
  cfePred <- cfe[, lambdaIndex]
  cfePred[(nEnv + nGenoTrain):ncol(ma)] <- as.matrix(parGeno[, -1])
  predTrain <- data.frame(E = dTrain$E, G = dTrain$G,
                          pred = as.numeric(ma[1:nrow(dTrain), ] %*% cfePred + mu))
  resTrain <- data.frame(E = dTrain$E, G = dTrain$G,
                         res = dTrain$Y - predTrain$pred)
  if (!is.null(testEnv)) {
    ## Make predictions for test set.
    predTest <- data.frame(E = dTest$E, G = dTest$G,
               pred = as.numeric(ma[(nrow(dTrain)+ 1):(nrow(dTrain) + nrow(dTest)), ] %*%
                                   cfePred + mu))
    resTest <- data.frame(E = dTest$E, G = dTest$G,
                          res = dTest$Y - predTest$pred)
  } else {
    predTest <- NULL
    resTest <- NULL
  }
  ## Compute the mean of each environmental index, in each environment
  indFrame <- aggregate(dat[, indices], by = list(E = dat$E), FUN = mean)
  rownames(indFrame) <- indFrame$E
  indFrameTrain <- indFrame[trainEnv, ]
  if (!is.null(testEnv)) indFrameTest  <- indFrame[testEnv, ]
  indFrameTrain <- merge(indFrameTrain, envMain, by.x = "E", by.y = "row.names")
  indFrameTrain2 <- merge(indFrameTrain, envMain2, by.x = "E",
                          by.y = "row.names")
  colnames(indFrameTrain)[ncol(indFrameTrain)] <- "envMainFitted"
  colnames(indFrameTrain2)[ncol(indFrameTrain2)] <- "envMainFitted"
  if (!is.null(partition)) {
    indFrameTrain <- merge(indFrameTrain, partition)
    indFrameTrain2 <- merge(indFrameTrain2, partition)
  }
  rownames(indFrameTrain) <- indFrameTrain$E
  rownames(indFrameTrain2) <- indFrameTrain2$E
  ## predict env. main effects.
  ## note : parEnvTrain and parEnvTest will now be matrices; not vectors.
  if (is.null(partition)) {
    glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                   y = indFrameTrain$envMainFitted,
                                   alpha = alpha, nfolds = nfolds,
                                   grouped = nrow(dTrain) / nfolds < 3)
    glmnetOut2 <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain2[, indices]),
                                    y = indFrameTrain2$envMainFitted,
                                    alpha = alpha, nfolds = nfolds,
                                    grouped = nrow(dTrain) / nfolds < 3)
  } else {
    glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                   y = indFrameTrain$envMainFitted,
                                   alpha = alpha,
                                   foldid = indFrameTrain$partition,
                                   grouped = max(table(indFrameTrain$partition)) > 2)
    glmnetOut2 <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain2[, indices]),
                                    y = indFrameTrain2$envMainFitted,
                                    alpha = alpha,
                                    foldid = indFrameTrain2$partition,
                                    grouped = max(table(indFrameTrain2$partition)) > 2)
  }
  parEnvTrain <- predict(object = glmnetOut,
                         newx = as.matrix(indFrameTrain[, indices]),
                         s = "lambda.min")
  if (!is.null(testEnv)) {
    parEnvTest  <- predict(object = glmnetOut,
                           newx = as.matrix(indFrameTest[, indices]),
                           s = "lambda.min")
    parEnvTest2 <- predict(object = glmnetOut2,
                           newx = as.matrix(indFrameTest[, indices]),
                           s = "lambda.min")
  }
  indicesTest <- NULL
  if (!is.null(testEnv)) {
    ## add the predicted environmental main effects:
    predTest$pred <- predTest$pred +
      as.numeric(parEnvTest[as.character(dTest$E), ])
    indicesTest <- data.frame(indFrameTest,
                              envMainPred = as.numeric(parEnvTest))
  }
  indicesTrain <- data.frame(indFrameTrain,
                             envMainPred = as.numeric(parEnvTrain))
  ## Compute statistics for training data.
  predMain <- as.numeric(predict(modelMain, newx = mm))
  trainAccuracyEnv <- getAccuracyEnv(datNew = dTrain[, "Y"],
                                     datPred = predTrain$pred,
                                     datPredMain = predMain,
                                     datE = dTrain[, "E"],
                                     corType = corType,
                                     rank = TRUE)
  ## Compute accuracies for genotypes.
  trainAccuracyGeno <- getAccuracyGeno(datNew = dTrain[, "Y"],
                                       datPred = predTrain$pred,
                                       datG = dTrain[, "G"],
                                       corType = corType)
  if (!is.null(testEnv)) {
    ## Compute statistics for test data.
    predMainTest <- mainOnly[as.character(dTest$G)] +
      as.numeric(parEnvTest2[as.character(dTest$E), ])
    testAccuracyEnv <- getAccuracyEnv(datNew = dTest[, "Y"],
                                      datPred = predTest$pred,
                                      datPredMain = predMainTest,
                                      datE = dTest[, "E"],
                                      corType = corType,
                                      rank = TRUE)
    ##################
    if (is.null(partition)) {
      glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                     y = trainAccuracyEnv$r, alpha = alpha,
                                     foldid = nfolds)
    } else {
      glmnetOut <- glmnet::cv.glmnet(x = as.matrix(indFrameTrain[, indices]),
                                     y = trainAccuracyEnv$r, alpha = alpha,
                                     foldid = indFrameTrain$partition,
                                     grouped = max(table(indFrameTrain$partition)) > 2)
    }
    rTest  <- predict(object = glmnetOut,
                      newx = as.matrix(indFrameTest[, indices]),
                      s = "lambda.min")
    testAccuracyEnv$rEst <- as.numeric(rTest)
    ## Compute accuracies for genotypes.
    testAccuracyGeno <- getAccuracyGeno(datNew = dTest[, "Y"],
                                        datPred = predTest$pred,
                                        datG = dTest[, "G"],
                                        corType = corType)
  } else {
    testAccuracyEnv  <- NULL
    testAccuracyGeno <- NULL
  }
  ## Create RMSE.
  RMSEtrain <- sqrt(mean((dTrain$Y - predTrain$pred) ^ 2, na.rm = TRUE))
  if (!is.null(testEnv)) {
    RMSEtest  <- sqrt(mean((dTest$Y - predTest$pred) ^ 2, na.rm = TRUE))
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
  if (verbose == TRUE) {
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
              resTrain = resTrain,
              resTest = resTest,
              mu = mu,
              envInfoTrain = indicesTrain,
              envInfoTest  = indicesTest,
              parGeno = parGeno,
              trainAccuracyEnv = trainAccuracyEnv,
              testAccuracyEnv = testAccuracyEnv,
              trainAccuracyGeno = trainAccuracyGeno,
              testAccuracyGeno = testAccuracyGeno,
              lambda = lambda[lambdaIndex],
              lambdaSequence = lambdaSequence,
              RMSEtrain = RMSEtrain,
              RMSEtest = RMSEtest,
              Y = traitName,
              G = genoName,
              E = envName,
              indices = indices,
              quadratic = quadratic)
  return(out)
}
