#' @title LDA_activity_single
#'
#' @description calculation of clonogenic activity from data collected by a
#'   limiting dilution assay (LDA) experiment (i.e. numbers of: cells seeded,
#'   wells, positive wells).
#'
#' @param x numeric data.frame or matrix with three columns (cells,
#'   wells, positive wells)
#' @param name optional: experiment name (e.g. name of cell line)
#' @param treat optional: treatment (e.g. irradiation dose in Gy)
#'
#' @return list object with estimated activity, 95%-confidence interval,
#'   84%-confidence interval, estimated parameters, corresponding covariance
#'   matrix, fit-object and p-value for cooperativity-test
#'
#' @examples
#' x <- data.frame("cells" = c(10,50,100,250),
#'                 "wells" = rep(25,4),
#'                 "positive" = c(2,5,10,20))
#' act <- LDA_activity_single(x)
#' data(LDAdata)
#' cell.line <- unique(LDAdata$name)[1]
#' x <- subset.data.frame(
#'       LDAdata,
#'       subset = (name==cell.line) & (Group == 0))
#' LDA_activity_single(x[,4:6])
#' @importFrom stats "glm" "binomial" "predict" "qnorm" "pnorm"
#' @export
#'
LDA_activity_single <- function(x,
                         name = "cell line a",
                         treat = "no") {
  if (!(class(x)[1] %in% c("data.frame","matrix"))){
    stop("error: x must be of class data.frame or matrix")
  }
  x <- as.data.frame(x)
  if (ncol(x) != 3){
    stop("error: number of columns must be 3 (cells, wells, positive)")
  }
  colnames(x) <- c("cells","wells","positive")
  if (!is.numeric(x$cells) | !is.numeric(x$wells) | !is.numeric(x$positive)){
    stop("error: all elements of x must be numeric")
  }
  if (min(x) < 0){
    stop("error: x must be non-negative")
  }
  if (min(x$cells) == 0){
    stop("error: cells must be positive")
  }
  if (min(x$wells - x$positive) < 0){
    stop("error: number of wells must not be smaller than
         number of positive wells")
  }

  single <- data.frame("cells" = rep(NA,sum(x$wells)),"positive" = NA)
  flnr <- 1
  for (i in seq_along(rownames(x))){
    single$cells[flnr:(flnr+x$wells[i]-1)] <- x$cells[i]
    single$positive[flnr:(flnr+x$wells[i]-1)] <-
      c(rep(1,x$positive[i]),rep(0,x$wells[i]-x$positive[i]))
    flnr <- flnr + x$wells[i]
  }

  X.glm <- data.frame("y" = single$positive,
                      "x" = log(single$cells))
  fit.mod <- suppressWarnings(glm(y ~ x,
                 family = binomial(link = "cloglog"),
                 data = X.glm))

  sum.fit <- summary(fit.mod)
  est <- sum.fit$coefficients
  Sig <- sum.fit$cov.unscaled

  new.data <- data.frame("x" = log((10^seq(-4,
                                           max(log10(10*x$cells)),
                                           1/10000))))
  pred <- predict(object = fit.mod,
                  newdata = new.data,
                  type = "response",
                  se.fit = TRUE)

  x.est <- exp(-est[1,1]/est[2,1])
  # approximate 95% confidence interval
  d.c <-  exp(new.data$x)
  calc_act_CI <- function(pred, alpha){
    x.uc <- (1-(pred$fit+qnorm(1-alpha/2)*pred$se.fit))
    x.uc[x.uc<0] <- 0
    x.uc <- log(x.uc)
    x.lc <- log(1-(pred$fit+qnorm(alpha/2)*pred$se.fit))

    x.lb.1 <- d.c[which(x.lc<(-1))[1]]
    x.lb.0 <- d.c[(which(x.lc<(-1))[1])-1]
    x.lb <- (x.lb.0+x.lb.1)/2
    rm(x.lb.0,x.lb.1)

    x.ub.1 <- d.c[which(x.uc<(-1))[1]]
    x.ub.0 <- d.c[(which(x.uc<(-1))[1])-1]
    x.ub <- (x.ub.0+x.ub.1)/2
    rm(x.ub.0,x.ub.1)
    return(c(x.ub,x.lb))
  }
  CI.95 <- calc_act_CI(pred,0.05)
  # approximate 83.5% confidence interval
  CI.84 <- calc_act_CI(pred,0.165)

  # non-linearity
  b <- est[2,1]
  if (b < 1){
    p.b <- 2 * pnorm(q = b,mean = 1,sd = est[2,2])
  } else {
    p.b <- 2 * (1-pnorm(q = b,mean = 1,sd = est[2,2]))
  }

  results <- list("name" = name,
                  "treatment" = treat,
                  "act" = x.est,
                  "CI" = CI.95,
                  "CI.appr.SF" = CI.84,
                  "est" = est,
                  "Sigma" = Sig,
                  "model" = fit.mod,
                  "p.lin.Model" = p.b)
  class(results) <- "LDA_activity_object"
  return(results)
}
