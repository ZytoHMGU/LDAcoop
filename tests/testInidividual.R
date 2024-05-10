# test cases
library(LDAcoop)
data(LDAdata)
cell.line <- unique(LDAdata$name)[1]
x <- subset.data.frame(
     LDAdata,
     subset = (name==cell.line) & (Group == 0))
summary(x)

x <- x[,4:6]
colnames(x) <- c("cells","wells","positive")
#LDA_activity_single(x[,4:6])

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
summary(fit.mod)

sum.fit <- summary(fit.mod)
est <- sum.fit$coefficients
Sig <- sum.fit$cov.unscaled

x.est <- exp(-est[1,1]/est[2,1])

hl10 <- log10(x.est)
R5.min <- seq(hl10 - 10, hl10 -  4, 1/10  )[-1]
R4.min <- seq(hl10 -  4, hl10 -  3, 1/10  )[-1]
R3.min <- seq(hl10 -  3, hl10 -  2, 1/33  )[-1]
R2.min <- seq(hl10 -  2, hl10 -  1, 1/100 )[-1]
R1.min <- seq(hl10 -  1, hl10 -1/10,1/333 )[-1]
R0     <- seq(hl10 -1/10,hl10 +1/10,1/1000)[-1]
R1.max <- seq(hl10 +1/10,hl10 +  1, 1/333 )[-1]
R2.max <- seq(hl10 +  1, hl10 +  2, 1/100 )[-1]
R3.max <- seq(hl10 +  2, hl10 +  3, 1/33  )[-1]
R4.max <- seq(hl10 +  3, hl10 +  4, 1/10  )[-1]
R5.max <- seq(hl10 +  4, hl10 + 10, 1/10  )[-1]
Svec <- 10^c(R5.min,R4.min,R3.min,R2.min,R1.min,
             R0,
             R1.max,R2.max,R3.max,R4.max,R5.max)
rm(R5.min,R4.min,R3.min,R2.min,R1.min,R0,R1.max,R2.max,R3.max,R4.max,R5.max)

new.data <- data.frame("x" = log(Svec))
pred <- predict(object = fit.mod,
                newdata = new.data,
                type = "response",
                se.fit = TRUE)

d.c <-  exp(new.data$x)
#calc_act_CI <- function(pred, alpha){

alpha <- 0.05

# predicted:
x.uc <- (1-(pred$fit+qnorm(1-alpha/2)*pred$se.fit))
x.lc <- (1-(pred$fit+qnorm(alpha/2)*pred$se.fit))

par(mfrow = c(1,2))
plot(x = new.data$x,y = pred$fit,type = 'l',xlim = c(-10,10),ylim = c(0,1))
par(new = T)
plot(x = new.data$x,y = 1-x.uc,type = 'l',lty = 2,col = "blue",xlim = c(-10,10),ylim = c(0,1))
par(new = T)
plot(x = new.data$x,y = 1-x.lc,type = 'l',lty = 2,xlim = c(-10,10),ylim = c(0,1))
abline(h = c(0,1),col = 'red')

x.lc[x.lc<0] <- 0
x.lc[x.lc>1] <- 1
x.uc[x.uc<0] <- 0
x.uc[x.uc>1] <- 1
x.m <- 1-pred$fit
plot(x = new.data$x,y = 1-x.m,type = 'l',xlim = c(-10,10),ylim = c(0,1))
par(new = T)
plot(x = new.data$x,y = 1-x.uc,type = 'l',lty = 2,col = "blue",xlim = c(-10,10),ylim = c(0,1))
par(new = T)
plot(x = new.data$x,y = 1-x.lc,type = 'l',lty = 2,xlim = c(-10,10),ylim = c(0,1))
abline(h = c(0,1),col = 'red')

# lin plot
x.lc <- log(x.lc)
x.uc <- log(x.uc)
x.m <- log(x.m)

plot(x = new.data$x,y = x.m,type = 'l',xlim = c(-10,10),ylim = c(-10,0))
par(new = T)
plot(x = new.data$x,y = x.uc,type = 'l',lty = 2,col = "blue",xlim = c(-10,10),ylim = c(-10,0))
par(new = T)
plot(x = new.data$x,y = x.lc,type = 'l',lty = 2,xlim = c(-10,10),ylim = c(-10,0))
abline(h = -1,col = 'red')

xl <- c(-0.42,0.42)
yl <- c(-1.2,-0.8)
plot(x = new.data$x,y = x.m,type = 'b',pch = '.',xlim = xl,ylim = yl)
par(new = T)
plot(x = new.data$x,y = x.uc,type = 'b',pch = '.',lty = 2,col = "blue",xlim = xl,ylim = yl)
par(new = T)
plot(x = new.data$x,y = x.lc,type = 'b',pch = '.',lty = 2,xlim = xl,ylim = yl)
abline(h = -1,col = 'red')

# idea and detail
xl <- c(0,1.42)
yl <- c(-1.2,-0.8)
plot(x = exp(new.data$x),y = x.m,type = 'o',pch = '.',lty = 0,xlim = xl,ylim = yl)
abline(h = -1,col = 'red')
abline(v = x.est,col = 'red')
xl <- x.est+c(-0.02,0.02)
yl <- c(-1.02,-0.98)
plot(x = exp(new.data$x),y = x.m,type = 'o',pch = '.',lty = 0,xlim = xl,ylim = yl)
abline(h = -1,col = 'red')
abline(v = x.est,col = 'red')

# lower curve: x.uc 'coming from left'
x.r <- d.c[which(x.uc<(-1))[1]]
x.l <- d.c[which(x.uc<(-1))[1]-1]
y.r <- x.uc[which(x.uc<(-1))[1]]
y.l <- x.uc[which(x.uc<(-1))[1]-1]
x.s <- x.l + (x.r-x.l) * (1+y.l)/(y.l-y.r)

# upper curve: x.lc 'coming from right'
x.lc <- x.lc[length(x.lc):1]
d.c <- d.c[length(d.c):1]
x.r <- d.c[which(x.lc>(-1))[1]]
x.l <- d.c[which(x.lc>(-1))[1]-1]
y.r <- x.lc[which(x.lc>(-1))[1]]
y.l <- x.lc[which(x.lc>(-1))[1]-1]
x.s <- x.l - abs(x.r-x.l) * abs(y.l+1)/abs(y.l-y.r)
