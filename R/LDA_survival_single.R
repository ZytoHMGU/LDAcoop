#' @title LDA_suvival_single
#'
#' @description calculate clonogenic survival fraction from LDA_activity objects
#'
#' @param act.0 reference activity
#' @param act.x activity after treatment
#'
#' @return list object with survival fraction, estimated confidence intervals
#'   (by error propagation through first order Taylor series approximation and
#'   by combination of 84%-uncertainty-intervals of activity estimates)
#'
#' @examples
#' x <- data.frame("dose" = c(10,50,100,250,10,50,100,250),
#'                 "wells" = rep(25,8),
#'                 "positive" = c(2,5,10,20,1,2,6,11),
#'                 "group" = c(rep("A",4),rep("B",4)))
#' act <- LDA_activity(x)
#' sf <- LDA_survival(act)
#' @export
#'
LDA_survival_single <- function(act.0,act.x){
  result <- list("treat" = act.x$treatment,
              "sf" = NA,
              "CI.ep" = c(NA,NA),
              "CI.act" = c(NA,NA))
  est.0 <- act.0$est
  Sigma.0 <- act.0$Sigma
  est.x <- act.x$est
  Sigma.x <- act.x$Sigma

  # g = log(sf)
  g <- est.x[1]/est.x[2] - est.0[1]/est.0[2]
  result$sf <- exp(g)

  sig.g2 <- 1/est.0[2]^2 * ( Sigma.0[1,1] - 2 * est.0[1]/est.0[2] * Sigma.0[1,2] +
                               (est.0[1]/est.0[2])^2 * Sigma.0[2,2] )  +
    1/est.x[2]^2 * ( Sigma.x[1,1] - 2 * est.x[1]/est.x[2] * Sigma.x[1,2] +
                       (est.x[1]/est.x[2])^2 * Sigma.x[2,2] )
  sig.g <- sqrt(sig.g2)

  result$CI.ep <- exp(g+qnorm(p = c(0.025,0.975))*sig.g)

  result$CI.act <- c(act.0$CI.appr.SF[1]/act.x$CI.appr.SF[2],
                     act.0$CI.appr.SF[2]/act.x$CI.appr.SF[1])
  return(result)
}
