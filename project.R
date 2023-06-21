install.packages("fitdistrplus")
library(fitdistrplus)
install.packages("lmomco")
library(lmomco)
install.packages("evd")
library(evd)
install.packages("FAdist")
library(FAdist)
install.packages("evd")
library(evd)
install.packages("MASS")
library(MASS)
install.packages("optimx")
library(optimx)
install.packages("survival")
library(survival)
install.packages("goftest")
library(goftest)
install.packages("stats")
library(stats)
install.packages("pracma")
library(pracma)
install.packages("ExtDist")
library(ExtDist)
install.packages("readxl")
library("readxl")
rice = read_excel("D:\\Statistic_TU\\project\\raw data (rice).xlsx")
total = read_excel("D:\\Statistic_TU\\project\\rice total.xlsx")
head(rice)
dim(rice)
summary(rice)
x1 <- rice$Q1
x2 <- rice$Q2
x3 <- rice$Q3
x4 <- rice$Q4

x1 <- x1/1000
x2 <- x2/1000
x3 <- x3/1000
x4 <- x4/1000
x1
x2
x3
x4
total <- total$total
total <- total/1000
total
summary(x1)
median(x1)
median(x2)
median(x3)
median(x4)

#weibull
fitw<-fitdist(x4, "weibull")
plot(fitw)
summary(fitw)
gofstat(fitw)
ad.test (x1, "weibull", shape = 3.843397, scale = 2142.618166)

#lognormal
fitln1 <- fitdist(x1, "lnorm")
fitln2 <- fitdist(x2, "lnorm")
fitln3 <- fitdist(x3, "lnorm")
fitln4 <- fitdist(x4, "lnorm")
fitln
plot(fitln)
summary(fitln1)
summary(fitln2)
summary(fitln3)
summary(fitln4)
gofstat(fitln1)
fitln3 <- fitdist(x3, "lnorm")
summary(fitln3)
ad.test (x1, "lnorm", meanlog = fitln1$estimate[1], sdlog = fitln1$estimate[2])
ad.test (x2, "lnorm", meanlog = fitln2$estimate[1], sdlog = fitln2$estimate[1])
ad.test (x3, "lnorm", meanlog = fitln3$estimate[1], sdlog = fitln3$estimate[1])
ad.test (x4, "lnorm", meanlog = fitln4$estimate[1], sdlog = fitln4$estimate[1])
fitln1$estimate[1]
fitln1$estimate[2]

#gamma
fitg <- fitdist(x4, "gamma")
plot(fitg)
summary(fitg)
gofstat(fitg)
ad.test (x4, "gamma", shape=14.434225996, rate=0.006308761)



#gev
dgev <- function(x,a,b,c) 1/b*(1+c*((x-a)/b))^(-1/c-1)*exp(-(1+c*((x-a)/b))^(-1/c))
pgev <- function(q,a,b,c) exp(-(1+c * (q-a)/b)^(-1/c))
fitgev <- fitdist(x1, "gev",start=list(a=1,b=1,c=0.5))
plot(fitgev)
fitgev <- fitdist(x4, "gev",start=list(a=1500,b=100,c=20), optim.method = "L-BFGS-B", lower = c(0.001, 0.001,0.1), upper = c(1000, 400,5))
summary(fitgev)

gofstat(fitgev)
ad.test (x1, "gev", a = 1697.51190277, b = 441.80929727, c = -0.04350634)
ad.test (x1, "gev", a = fitgev$estimate[1], b = fitgev$estimate[2], c = fitgev$estimate[3])

#burr
mu <- mean(x1)
sigma2 <- var(x1)
lambda <- (sqrt(sigma2)/mu)^(-2/2.5)
k <- (2.5/lambda)*(2^(2/lambda)-1)

dBurr <- function(x,c,k) (c*k*x^(c-1))/(1+x^c)^(k+1)
pBurr <- function(q,c,k) 1-(1+q^c)^-k
fitb <- fitdist(x4,"Burr", start = list(c=0.5,k=1), optim.method = "L-BFGS-B", lower = c(0.1,0.001), upper = c(100, 10))
summary(fitb)
ad.test(x1, "Burr", c=fitb$estimate[1], k=fitb$estimate[2])
ad.test(x2, "Burr", c=fitb$estimate[1], k=fitb$estimate[2])
ad.test(x3, "Burr", c=fitb$estimate[1], k=fitb$estimate[2])
ad.test(x4, "Burr", c=fitb$estimate[1], k=fitb$estimate[2])
ad.test(total, "Burr", c=fitb$estimate[1], k=fitb$estimate[2])

#jsb

djsbb <- function(x,a,b,c,d) (b/sqrt(2*pi))*(d/((x-c)*(c+d-x)))*exp((-1/2)*(a+(b*log((x-c)/(c+d-x))))^2 ) 
pjsbb <- function(q,a,b,c,d)  {
  res <- numeric(length(q))
  for (i in 1:length(q)) {
    if (q[i] >=0) {
      res[i] <- pnorm(a+(b*log((q[i]-c)/(d-q[i]+c))))
    } else {
      res [i] <- 0
    }
  }
  return (res)
}

fitJs <- fitdist(x4, "jsbb",start=list(a=0.5,b=1,c=800,d=3000), optim.method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(10, 10,1500,3500))
summary(fitJs)
ad.test (x4, "jsbb", a = fitJs$estimate[1], b = fitJs$estimate[2], c = fitJs$estimate[3], d = fitJs$estimate[4])


#loglogistic 3p
dL3p <- function(x,a,b,y) (a/b)*(x-y)^(a-1)*(1+(x-y)^a/b^a)^(-2)
pL3p <- function(q,a,b,y)  1/(1+(b/y)*(q/c) ^(-a))
fitL3p <- fitdist(x1, "L3p",start=list(a=1,b=1000,y=1),optim.method = "L-BFGS-B",lower = c(0.001, 0.001,0.1), upper = c(2,3,4))

fitL3p <- fitdist(x4, "L3p",start=list(a=0.01,b=100,y=100),optim.method = "L-BFGS-B",lower = c(0,0,0),upper = c(1,1800,500))
summary(fitL3p)
ad.test (x1, "L3p", a = fitL3p$estimate[1], b = fitL3p$estimate[2], y = fitL3p$estimate[3])



#gamma3p

dgammat <- function(x,a,b,c) (((x-c)^(a-1))/(b^a)*gamma(a)) *(exp(-((x-c)/b)))
pgammat <- function(q,a,b,c) {
  res <- numeric(length(q))
  for (i in 1:length(q)) {
    if (q[i] >= 0) {
      lowerincg <- gammainc(b/q[i],a)
      res[i] <- 1-(lowerincg[1]/gamma(a))
    } else {
      res[i] <- 0
    }
  }
  return (res)
}
fitgamma3 <- fitdist(x4, "gammat",start=list(a=3,b=20,c=500), optim.method = "L-BFGS-B", lower = c(0.1,0.1,0.1), upper = c(5, 50,1000))
summary(fitgamma3)

ad.test(x4, "gammat", a = fitgamma3$estimate[1], b =fitgamma3$estimate[2], c = fitgamma3$estimate[3])



# bootstrap
set.seed(123)

#para bootstrap Q1
ParBS.CIQ1 <- function (x, N, samp.size, alpha=0.05){
  # generate the bootstrap samples
  n <- length(x)
  theta.star <- c ()
  for (j in 1:N) {
    x.star <- rlnorm(27, meanlog = fitln1$estimate[1], sdlog = fitln1$estimate[2])
    fits <- fitdist(x.star , "lnorm")
    theta.hat <- exp(fits$estimate[1]+((fits$estimate[2]^2)/2))
    theta.star[j] <- theta.hat
  }
  t <- exp(fitln1$estimate[1]+((fitln1$estimate[2]^2)/2))
  l <- ceiling ((1.2)*alpha*N)
  u <- ceiling ((1-(1/2)*alpha)*N)
  theta.star <- sort(theta.star)
  return(c(2*t-theta.star[u],2*t-theta.star[l]))
}

#para bootstrap Q2
ParBS.CIQ2 <- function (x, N, samp.size, alpha=0.05){
  # generate the bootstrap samples
  n <- length(x)
  theta.star <- c ()
  for (j in 1:N) {
    x.star <- rlnorm(27, meanlog = fitln2$estimate[1], sdlog = fitln$estimate[2])
    fits <- fitdist(x.star , "lnorm")
    theta.hat <- exp(fits$estimate[1]+((fits$estimate[2]^2)/2))
    theta.star[j] <- theta.hat
  }
  t <- exp(fitln2$estimate[1]+((fitln2$estimate[2]^2)/2))
  l <- ceiling ((1.2)*alpha*N)
  u <- ceiling ((1-(1/2)*alpha)*N)
  theta.star <- sort(theta.star)
  return(c(2*t-theta.star[u],2*t-theta.star[l]))
}

#para bootstrap Q3
ParBS.CIQ3 <- function (x, N, samp.size, alpha=0.05){
  # generate the bootstrap samples
  n <- length(x)
  theta.star <- c ()
  for (j in 1:N) {
    x.star <- rlnorm(27, meanlog = fitln3$estimate[1], sdlog = fitln3$estimate[2])
    fits <- fitdist(x.star , "lnorm")
    theta.hat <- exp(fits$estimate[1]+((fits$estimate[2]^2)/2))
    theta.star[j] <- theta.hat
  }
  t <- exp(fitln3$estimate[1]+((fitln3$estimate[2]^2)/2))
  l <- ceiling ((1.2)*alpha*N)
  u <- ceiling ((1-(1/2)*alpha)*N)
  theta.star <- sort(theta.star)
  return(c(2*t-theta.star[u],2*t-theta.star[l]))
}

#para bootstrap Q4
ParBS.CIQ4 <- function (x, N, samp.size, alpha=0.05){
  # generate the bootstrap samples
  n <- length(x)
  theta.star <- c ()
  for (j in 1:N) {
    x.star <- rlnorm(27, meanlog = fitln4$estimate[1], sdlog = fitln4$estimate[2])
    fits <- fitdist(x.star , "lnorm")
    theta.hat <- exp(fits$estimate[1]+((fits$estimate[2]^2)/2))
    theta.star[j] <- theta.hat
  }
  t <- exp(fitln4$estimate[1]+((fitln4$estimate[2]^2)/2))
  l <- ceiling ((1.2)*alpha*N)
  u <- ceiling ((1-(1/2)*alpha)*N)
  theta.star <- sort(theta.star)
  return(c(2*t-theta.star[u],2*t-theta.star[l]))
}

ParBS.CIQ1(x1,10000, samp.size = 100,alpha=0.05)
ParBS.CIQ2(x2,10000, samp.size = 100,alpha=0.05)
ParBS.CIQ3(x3,10000, samp.size = 100,alpha=0.05)
ParBS.CIQ4(x4,10000, samp.size = 100,alpha=0.05)


#plot

# estimate the density function
dens <- density(x4)

max(x4)
# plot a histogram of the density
breaks <- seq(1100, 3500, length.out = 9)
hist(x4, breaks = breaks, col = "lightblue", main = "Quarter 4", xlab = "Data Values (kiloton)", ylab = "f(x)", freq = FALSE)

# add the PDFs of the 8 distributions
x <- seq(min(x4), max(x4), length = 100)
lines(x, dweibull(x, fitw$estimate[1], fitw$estimate[2]), col = "red", lwd = 2)
lines(x, dlnorm(x, fitln3$estimate[1], fitln3$estimate[2]), col = "green", lwd = 2)
lines(x, dgamma(x, fitg$estimate[1], fitg$estimate[2]), col = "purple", lwd = 2)
lines(x, dburr(x, fitb$estimate[1], fitb$estimate[2]), col = "orange", lwd = 2)
lines(x, dgev(x, fitgev$estimate[1], fitgev$estimate[2],fitgev$estimate[3]), col = "brown", lwd = 2)
lines(x, djsbb(x, fitJs$estimate[1], fitJs$estimate[2], fitJs$estimate[3], fitJs$estimate[4]), col = "blue", lwd = 2)
lines(x, dgammat(x, fitgamma3$estimate[1], fitgamma3$estimate[2], fitgamma3$estimate[3]), col = "black", lwd = 2)
lines(x, dL3p(x, fitL3p$estimate[1], fitL3p$estimate[2], fitL3p$estimate[3]), col = "gray", lwd = 2)

# add a legend
legend("topright", legend = c("Weibull", "Lognormal", "Gamma", "Burr", "GEV", "Jsb", "Gam3p", "LLD3"), col = c("red", "green", "purple", "orange", "brown", "blue", "black", "gray"), lwd = 2)

# find probbality 
#Q1
plnorm(1500, meanlog=fitln1$estimate[1], sdlog = fitln1$estimate[2])
a <- plnorm(2000,meanlog=fitln1$estimate[1],sdlog=fitln1$estimate[2])
b <- plnorm(1500,meanlog=fitln1$estimate[1],sdlog=fitln1$estimate[2])
result1 <- a - b
result1
c <- plnorm(2500,meanlog=fitln1$estimate[1],sdlog=fitln1$estimate[2])
d <- plnorm(2000,meanlog=fitln1$estimate[1],sdlog=fitln1$estimate[2])
result2 <- c - d
result2
plnorm(2500, meanlog=fitln1$estimate[1], sdlog = fitln1$estimate[2],lower.tail=FALSE)

#Q2
plnorm(1500, meanlog=fitln2$estimate[1], sdlog = fitln2$estimate[2])
a <- plnorm(2000,meanlog=fitln2$estimate[1],sdlog=fitln2$estimate[2])
b <- plnorm(1500,meanlog=fitln2$estimate[1],sdlog=fitln2$estimate[2])
result1 <- a - b
result1
c <- plnorm(2500,meanlog=fitln2$estimate[1],sdlog=fitln2$estimate[2])
d <- plnorm(2000,meanlog=fitln2$estimate[1],sdlog=fitln2$estimate[2])
result2 <- c - d
result2
plnorm(2500, meanlog=fitln2$estimate[1], sdlog = fitln2$estimate[2],lower.tail=FALSE)

#Q3
plnorm(1500, meanlog=fitln3$estimate[1], sdlog = fitln3$estimate[2])
a <- plnorm(2000,meanlog=fitln3$estimate[1],sdlog=fitln3$estimate[2])
b <- plnorm(1500,meanlog=fitln3$estimate[1],sdlog=fitln3$estimate[2])
result1 <- a - b
result1
c <- plnorm(2500,meanlog=fitln3$estimate[1],sdlog=fitln3$estimate[2])
d <- plnorm(2000,meanlog=fitln3$estimate[1],sdlog=fitln3$estimate[2])
result2 <- c - d
result2
plnorm(2500, meanlog=fitln3$estimate[1], sdlog = fitln3$estimate[2],lower.tail=FALSE)

#Q4
plnorm(1500, meanlog=fitln4$estimate[1], sdlog = fitln4$estimate[2])
a <- plnorm(2000,meanlog=fitln4$estimate[1],sdlog=fitln4$estimate[2])
b <- plnorm(1500,meanlog=fitln4$estimate[1],sdlog=fitln4$estimate[2])
result1 <- a - b
result1
c <- plnorm(2500,meanlog=fitln4$estimate[1],sdlog=fitln4$estimate[2])
d <- plnorm(2000,meanlog=fitln4$estimate[1],sdlog=fitln4$estimate[2])
result2 <- c - d
result2
plnorm(2500, meanlog=fitln4$estimate[1], sdlog = fitln4$estimate[2],lower.tail=FALSE)
