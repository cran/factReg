#' Helper function for renaming Y, G, and E
#'
#' Helper function for renaming columns Y, G and Y in dat
#'
#' @importFrom utils hasName
#' @noRd
#' @keywords internal
renameYGE <- function(dat,
                      Y,
                      G,
                      E) {
  traitName <- ifelse(is.numeric(Y), names(dat)[Y], Y)
  for (var in c("Y", "G", "E")) {
    varVal <- get(var)
    if (is.numeric(varVal)) {
      names(dat)[varVal] <- var
    } else if (is.character(varVal)) {
      if (!hasName(x = dat, varVal)) {
        stop(varVal, " should be a column in dat.\n")
      }
      names(dat)[names(dat) == varVal] <- var
    }
  }
  ## Convert G and E to factor.
  if (!is.factor(dat[["G"]])) {
    dat[["G"]] <- as.factor(dat[["G"]])
  }
  if (!is.factor(dat[["E"]])) {
    dat[["E"]] <- as.factor(dat[["E"]])
  }
  dat <- droplevels(dat)
  return(dat)
}

#' Helper function
#'
#' Helper function for ....
#'
#' @noRd
#' @keywords internal
exRank <- function(y0,
                   yPred,
                   topOnly = TRUE,
                   k = 5,
                   m = 10) {
  b <- sort(y0, decreasing = TRUE)[k]
  topK <- which(y0 >= b)

  a <- sort(yPred, decreasing = TRUE)[m]
  topM <- which(yPred >= a)

  return(length(intersect(topK, topM)) / k)
}

#' Helper function
#'
#' Helper function for ....
#'
#' @noRd
#' @keywords internal
corFun <- function(x,
                   y,
                   corType) {
  cor(x, y, use = "na.or.complete", method = corType)
}


#' Helper function
#'
#' Helper function for ....
#'
#' @noRd
#' @keywords internal
getAccuracyEnv <- function(datNew,
                           datPred,
                           datPredMain,
                           datE,
                           corType,
                           rank = FALSE) {
  ## Split data by environments.
  sNew <- split(datNew, datE)
  sPred <- split(datPred, datE)
  sPredMain <- split(datPredMain, datE)
  ## Define helper functions.
  RMSEFun <- function(x, y) {
    sqrt(mean((x - y) ^ 2))
  }
  MADFun <- function(x, y) {
    mean(abs(x - y))
  }
  ## Compute statistics for test data.
  accuracy <- data.frame(Env = levels(datE),
                         r = mapply(FUN = corFun, sNew, sPred,
                                    MoreArgs = list(corType = corType)),
                         rMain = mapply(FUN = corFun, sNew, sPredMain,
                                        MoreArgs = list(corType = corType)),
                         RMSE = mapply(FUN = RMSEFun, sNew, sPred),
                         RMSEMain = mapply(FUN = RMSEFun, sNew, sPredMain),
                         MAD = mapply(FUN = MADFun, sNew, sPred),
                         row.names = NULL)
  ## Add rank.
  if (rank) {
    rankDat <- data.frame(
      rank = mapply(FUN = exRank, sNew, sPred),
      rankMain = mapply(FUN = exRank, sNew, sPredMain),
      row.names = NULL)
    accuracy <- cbind(accuracy, rankDat)
  }
  return(accuracy)
}

#' Helper function
#'
#' Helper function for ....
#'
#' @noRd
#' @keywords internal
getAccuracyGeno <- function(datNew,
                            datPred,
                            datG,
                            corType) {
  ## Split data by genotypes
  sNew <- split(datNew, datG)
  sPred <- split(datPred, datG)
  accuracy <- data.frame(Geno = levels(datG),
                         r = mapply(FUN = corFun, sNew, sPred,
                                    MoreArgs = list(corType = corType)))
  return(accuracy)
}
