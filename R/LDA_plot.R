#' @title LDA_plot
#'
#' @description plot activity from LDA data object
#'
#' @param LDA_act output from LDA_activity()
#'
#' @return none
#'
#' @examples
#' x <- data.frame("dose" = c(10,50,100,250,10,50,100,250),
#'                 "wells" = rep(25,8),
#'                 "positive" = c(2,5,10,20,1,2,6,11),
#'                 "group" = c(rep("A",4),rep("B",4)))
#' LDA_plot(LDA_activity(x))
#'
#' @importFrom grDevices "colorRampPalette"
#' @importFrom graphics "abline" "par" "plot" "title"
#' @importFrom graphics "polygon"
#' @export
#'
LDA_plot <- function(LDA_act){
  if (!(class(LDA_act) %in% c("LDA_activity_object","LDA_activity_list"))){
    stop("error: class of input must be of type LDA_activity, e.g. as returned
         by LDA_activity().")
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  plot_single <- function(LDA,xlim = NA,alpha = 0.05,col = "black"){
    d <- LDA$model$data
    if (is.na(xlim)){
      xlim <- c(0,1*max(d$x))
    }
    p <- aggregate.data.frame(x = d$y,by = list(d$x),mean)
    n <- aggregate.data.frame(x = !is.na(d$y),by = list(d$x),sum)
    x <- merge(p,n,by = "Group.1")
    colnames(x) <- c("x","p","n")
    x$ind <- x$p > (1 - 1e-14)
    x$pch <- c(1,6)[(x$ind)+1]
    x$p[x$ind] <- ((x$p*x$n-0.5)/x$n)[x$ind]

    new.data <- data.frame("x" = log(seq(1/1000,exp(xlim[2]),(exp(xlim[2]))/1000)))
    pred <- predict(object = LDA$model,
                    newdata = new.data,
                    type = "response",
                    se.fit = T)
    x.uc <- (1-(pred$fit+qnorm(1-alpha/2)*pred$se.fit))
    x.uc[x.uc<0] <- 0
    x.uc <- log(x.uc)
    x.lc <- log(1-(pred$fit+qnorm(alpha/2)*pred$se.fit))

    plot(x = exp(x$x), xlim = exp(xlim),
         y = log(1-x$p),ylim = c(1.1*min(log(1-x$p)),0),
         ylab = "log fract. nonresp.",
         xlab = "dose [#cells]",
         pch = x$pch,col = col)
    lines(exp((new.data$x)),log(1-pred$fit),lwd=2)
    polygon_y <- c(1-(pred$fit+1.96*pred$se.fit),rev(1-(pred$fit-1.96*pred$se.fit)))
    polygon_y[(polygon_y <= 0)] <- NaN
    polygon_y <- log(polygon_y)
    polygon_y[is.nan(polygon_y)] <- min(polygon_y,na.rm = T)
    polygon(x = c(exp(new.data$x),rev(exp(new.data$x))),
            y = polygon_y,
            col = rgb(0.2,0.2,0.2,0.2,1))
    abline(h = -log(exp(1)),lty = 1)
    abline(v = LDA$act,lty = 1)
    abline(v = LDA$CI,lty = 3)
  }

  if (identical(class(LDA_act),"LDA_activity_object")){
    plot_single(LDA_act,col = "#43E08700")
    title(paste0(LDA_act$name," ",LDA_act$treatment))
  }
  if (identical(class(LDA_act),"LDA_activity_list")){
    n_sp <- length(LDA_act)
    colhex <- colorRampPalette(c("#43E08700", "#00612A00"))(n_sp) # D23264
    colors <- col2rgb(colhex) / 255
    alpha <- 0.42
    if (n_sp > 3){
      par(mfrow = c(ceiling(n_sp/2),2))
    }
    for (i in 1:n_sp){
      plot_single(LDA_act[[i]],col = colhex[i])
      title(paste0(LDA_act[[i]]$name,", treatment: ",LDA_act[[i]]$treatment))
    }
  }
}
