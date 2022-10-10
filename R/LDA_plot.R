#' @title LDA_plot
#'
#' @description plot activity from LDA data
#'
#' @param LDA_tab LDA data.frame
#'      ("cells", "wells", "positive", "group", "replicate")
#'
#' @return none
#'
#' @examples
#' x <- data.frame("cells" = rep(c(10,50,100,250),times = 4),
#'                 "wells" = rep(25,16),
#'                 "positive" = c(2,5,10,20,1,2,6,11,3,4,8,22,1,1,7,12),
#'                 "group" = rep(c(rep("A",4),rep("B",4)),times = 2),
#'                 "replicate" = c(rep(1,8),rep(2,8)))
#' LDA_plot(x)
#' data(LDAdata)
#' Z1 <- subset.data.frame(LDAdata,subset = name == unique(LDAdata$name)[1])
#' LDA_plot(Z1[,c("S-value","# Tested","# Clonal growth","Group","replicate")])
#' @importFrom grDevices "colorRampPalette"
#' @importFrom graphics "abline" "par" "plot" "title"
#' @importFrom Hmisc "errbar" "ceil"
#' @importFrom graphics "polygon"
#' @export
#'
LDA_plot <- function(LDA_tab,
                     uncertainty = "ep",
                     xlim = NULL,
                     uncertainty.band = FALSE){
  if (class(LDA_tab)[1] != "data.frame"){
    stop("error: input must be of class data.frame")
  }
  if (!(uncertainty %in% c("ep","act"))){
    stop("error: uncertainty must be either
         'ep' (error propagation) or
         'act' (activity CIs)")
  }
  if (ncol(LDA_tab) == 3){
    LDA_tab$Group <- 0
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(
    mar = c(2.5, 3.5, 0.5, 0.5),
    mgp = c(1.5, 0.5, 0)
  )

  # get data per group and replicate
  grps <- unique(LDA_tab$Group)
  N_treat <- length(grps)

  colhex <- colorRampPalette(c("#43E08700", "#00612A00"))(N_treat) # D23264
  colors <- col2rgb(colhex) / 255
  alpha <- 0.42

  x_sp <- NULL
  for (gi in seq_along(grps)){
    d_g <- subset.data.frame(x = LDA_tab,
                             subset = Group == grps[gi])
    rplcts <- unique(d_g$replicate)
    for (ri in rplcts){
      ddd <- subset.data.frame(x = d_g,
                               subset = replicate == ri)
      d_a <- LDA_activity_single(x = ddd[,1:3])
      d <- d_a$model$data
      p <- aggregate.data.frame(x = d$y,
                                by = list(d$x),
                                FUN = mean)
      n <- aggregate.data.frame(x = !is.na(d$y),
                                by = list(d$x),
                                FUN = sum)
      x <- merge(p,n,by = "Group.1")
      colnames(x) <- c("x","p","n")
      x$ind <- x$p > (1 - 1e-14)
      x$pch <- c(3,6)[(x$ind)+1] # '+' for each replicate
      x$p[x$ind] <- ((x$p*(x$n*length(rplcts))-0.5)/(x$n*length(rplcts)))[x$ind]
      x$col <- colhex[gi]
      x_sp <- rbind(x_sp,x)
    }
    d_a <- LDA_activity_single(x = d_g[,1:3])
    d <- d_a$model$data
    p <- aggregate.data.frame(x = d$y,
                              by = list(d$x),
                              FUN = mean)
    n <- aggregate.data.frame(x = !is.na(d$y),
                              by = list(d$x),
                              FUN = sum)
    x <- merge(p,n,by = "Group.1")
    colnames(x) <- c("x","p","n")
    x$ind <- x$p > (1 - 1e-14)
    x$pch <- c(19,6)[(x$ind)+1] # mean-symbol
    x$p[x$ind] <- ((x$p*x$n-0.5)/x$n)[x$ind]
    x$col <- colhex[gi]
    x_sp <- rbind(x_sp,x)
  }

  if (sum(x_sp$ind==TRUE) > 0){
    x_sp$p[x_sp$ind] <- min(x_sp$p[x_sp$ind])
  }

  # +------------------------+
  # | plot Part I : activity |
  # +------------------------+
  if (is.null(xlim)) {
    xlim <- c(0,max(exp(x_sp$x)))
  }
  if (N_treat > 1){
    par(mfrow = c(2,1))
  }
  plot(x = exp(x_sp$x),
       xlim = xlim,
       y = log(1-x_sp$p), #ylim = c(1.1*min(log(1-x$p)),0),
       ylab = "log fract. nonresp.",
       xlab = "cells seeded",
       pch = x_sp$pch,col = x_sp$col)

  xlim <- c(0, max(x_sp$x))
  new.data <- data.frame(
    "x" = log(seq(1/4000,
                  exp(xlim[2]),
                  (exp(xlim[2]))/4000)))

  for (gi in seq_along(grps)){
    d_g <- subset.data.frame(x = LDA_tab,
                             subset = Group == grps[gi])
    d_a <- LDA_activity_single(x = d_g[,1:3])
    pred <- predict(object = d_a$model,
                    newdata = new.data,
                    type = "response",
                    se.fit = TRUE)
    x.uc <- (1-(pred$fit+qnorm(1-alpha/2)*pred$se.fit))
    x.uc[x.uc<0] <- 0
    x.uc <- log(x.uc)
    x.lc <- log(1-(pred$fit+qnorm(alpha/2)*pred$se.fit))
    lines(exp((new.data$x)),log(1-pred$fit),lwd=2,col = colhex[gi])
    if (uncertainty.band == TRUE){
      z <- qnorm(0.975)
      #z <- qnorm(1-0.165/2)
      polygon_y <- c(1-(pred$fit+z*pred$se.fit),
                     rev(1-(pred$fit-z*pred$se.fit)))
      polygon_y[(polygon_y <= 0)] <- NaN
      polygon_y <- log(polygon_y)
      polygon_y[is.nan(polygon_y)] <- min(polygon_y,na.rm = TRUE)
      polygon(x = c(exp(new.data$x),rev(exp(new.data$x))),
              y = polygon_y,border = NA,
              col = rgb(red = col2rgb(colhex[gi])[1],
                        green = col2rgb(colhex[gi])[2],
                        blue = col2rgb(colhex[gi])[3],alpha = 127,maxColorValue = 255))
    }
    abline(h = -1,lty = 3,col = "grey42")
  }

  # +------------------------+
  # | Part II :      SF      |
  # +------------------------+
  if (N_treat > 1){
    collect_sf <- data.frame(
      "Exp" = 0,
      "treat" = 0,
      "sf" = 1,
      "sf.msd" = 1,
      "sf.psd" = 1
    )
    SF <- LDA_survival(LDA_tab[,1:4])
    for (t in seq_along(SF)) {
      CurSF <- SF[[t]]
      if (uncertainty == "ep"){
        collect_sf <- rbind(collect_sf, c(t, CurSF$treat, CurSF$sf, CurSF[[3]]))
      } else {
        collect_sf <- rbind(collect_sf, c(t, CurSF$treat, CurSF$sf, CurSF[[4]]))
      }
    }
    #  collect_sf <- collect_sf[-1,]

    #par(mar = c(2.5, 3.75, 2.5, 0.5), mgp = c(1.5, 0.5, 0))
    #for (sfi in seq_along(SF)) {
    #if (SF[[sfi]]$name == "no name") {
    #  SF[[sfi]]$name <- paste0("Experiment ", sfi)
    #}
    # PD <-
    #    subset.data.frame(x = collect_sf, subset = (collect_sf$"Exp" == sfi))
    plot(
      x = collect_sf$treat,
      y = log10(collect_sf$sf),
      #main = SF[[sfi]]$name,
      col.main = rgb(
        red = 0, green = 148/255, blue = 64/255, alpha = 1,
        maxColorValue = 1
      ),
      las = 1,
      ylab = "",
      xaxt = "n",
      yaxt = "n",
      xlab = "treatment",
      ylim = c(floor(min(log10(collect_sf$sf.msd))),0),
      col = colhex,
      axes = FALSE,
      pch = 19
    )
    ytick <- floor(min(log10(collect_sf$sf.msd), na.rm = TRUE)):
      ceil(max(log10(collect_sf$sf.psd), na.rm = TRUE))
    axis(
      side = 2,
      at = ytick,
      las = 1,
      labels = paste0(10^ytick * 100, "%")
    )
    title(
      ylab = "clonogenic survival",
      line = 2.5
    )
    axis(
      side = 1,
      at = collect_sf$treat,
      labels = collect_sf$treat
    )
    with(data = collect_sf, errbar(treat, log10(sf), log10(sf.msd), log10(sf.psd),
                                   col = colhex,
                                   add = TRUE, pch = 1, errbar.col = colhex
    ))
  }
}
