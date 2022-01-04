rm(list = ls(all=TRUE))
setwd("~/Documents/Masters 2020:2021/Bayesian DM/BDM Assignment 2")
library(Rfast)
library(actuar)
set.seed(2020)

####################################
##Question 1
# data
N <- c(4,7,3)
i <- c(1,2,3)
X1j <- c(2.588, 2.364, 0.153, 0.418)
X2j <- c(0.016, 0.243, 0.330, 0.167, 2.507, 1.091, 1.311)
X3j <- c(1.431, 0.444, 0.053)
X <- list(X1j,X2j,X3j)

# Z, Y
YY <- function(i){
  prod(X[[i]]+1)
}

ZZ <- function(i){
  log(YY(i))
}

Y <- c(YY(1), YY(2), YY(3))
Z <- c(ZZ(1), ZZ(2), ZZ(3))

# MLE estimation 
# log Likelihood function
LL <- function(alphs){
  # alphs = vector of 3 alphas
  
  # value of likelihood 
  val <- 1
  for (i in 1:3){
    val <- val*alphs[i]^N[i]/(Y[i])^(alphs[i]+1)
  }
  
  # take log
  # negative since optim minimizes
  return(-log(val))
}

# optim performs minimization
mle <- optim(par = runif(3), fn = LL)
mle$par #1.341256 1.998522 2.295341
N/Z # from differentiating - same

###################
# functions for accrej gibbs 
gamm <- function(a){
  ni <- N[n]
  zi <- Z[n]
  
  ((zi+lam)^(ni-1))*a^ni*exp(-a*(zi+lam))/gamma(ni-1)
}

pareto <- function(a){
  # in form given by actuar
  2*(sum(alphsk))^2/(a+sum(alphsk))^3
}

# prob. integral transform 
randpar <- function(){
  # random uniform no 
  p <- runif(1)
  
  # return inv cdf of p
  sum(alphsk)/sqrt(1-p) - sum(alphsk)
}

c.optim <- function(a){
  # return -log(pi*g/h) to minimize
  val <- log(gamm(a))-log(pareto(a))
  
  return(-val)
}

# starting values
alphs <- mle$par
lam <- 1/mean(alphs)

# how many samples?
ns <- 4000

# storage of alpha series 
A <- matrix(0, nrow = ns+1, ncol = 3)
A[1,] <- alphs

####### Gibbs with Accept-Reject
set.seed(2020)
for (s in 2:(ns+1)){
  for (n in 1:3){
    # set acc indicator
    acc <- FALSE
    # set alpha and alphsk
    alphsk <- alphs[-n]
    alpha <- alphs[n]
    
    # find C to minimize gamma ratio
    alphastar <- nlm(c.optim, 2)$estimate
    C <- gamm(alphastar)/pareto(alphastar)
    
    while(!acc){ # loop until sample is accepted
      # generate candidate from pareto
      x <- randpar()
      #x <- rpareto(n = 1, shape = 2, scale = sum(alphsk))
      
      # gamma ratio
      g <- gamm(x)/(C*pareto(x))
      
      # generate u
      u <- runif(1)
      
      # accept-reject
      if(g>=u){
        A[s,n] <- x         # add to matrix
        alphs[n] <- x       # update alpha vector  
        lam <- 1/mean(alphs)# update lambda
        acc <- TRUE         # terminate if accepted
      }
    }
  }
}

# estimates
# burn-in
burn <- 100
# sample with burn-in removed 
samp <- A[-c(1:(burn+1)),]
# posterior mean
alphahat <- matrix(colmeans(samp), ncol = 3, nrow=1)
ones <- matrix(1, nrow = ns-burn, ncol = 1)
# standard error
means <- ones%*%alphahat
SEs <- colSums((samp-means)^2)/(nrow(samp)-1)

# traceplots
pdf("q1trace.pdf", width = 12, compress = FALSE)
par(mfrow=c(3,1))
for (i in 1:3){
  plot(A[,i], type = "l", ylab = expression(alpha),
       main = paste(expression(i), "=",i, sep = " "), col = "red")
  abline(v=100, col = "darkgrey")
  abline(h=alphahat[i], lwd = 2, col = "blue")
  abline(h = alphahat[i]+SEs[i], lty = 2, col = "grey", lwd = 2)
  abline(h = alphahat[i]-SEs[i], lty = 2, col = "grey", lwd = 2)
}
dev.off()

# autocorrelations
pdf("q1acf.pdf", width = 10, compress = F)
par(mfrow=c(3,1))
for (i in 1:3){
  acf(samp[,i], main = paste("i =", i, sep = " "))  
}
dev.off()

# histograms
pdf("q1hists.pdf", width = 12, compress = FALSE)
par(mfrow=c(1,3))
for(i in 1:3){
  hist(samp[,i], 
       main = paste(expression(i), "=",i,sep = " "), 
       xlab = bquote(alpha),
       breaks = 40, freq = F, col = "lightblue", ylim = c(0,0.8))
  lines(density(samp[,i]), col = "red", lwd = 2)
  abline(v = alphahat[i], lty = 2, col = "blue", lwd = 2)
  abline(v = alphahat[i]+SEs[i], lty = 2, col = "grey", lwd = 2)
  abline(v = alphahat[i]-SEs[i], lty = 2, col = "grey", lwd = 2)
}
dev.off()

#####################################
# Question 2 
library(rstan)

#data 
N <- c(3, 7, 4, 8, 5, 9)
X <- c(10, 33, 3, 39, 5, 50)
lamdahatdata <- X/N # scale by # weeks
mean(lamdahatdata) # sample mean

# data as named list
brkdwns <- list(J=length(N),X=X, N=diag(N))

# diffused prior mu
# uniform prior sigma
sfile1 <- c(
  "data{
  int<lower=0> J;         // number of facilities
  int<lower=0> X[J];      // number of breakdowns at facility J
  matrix[J,J] N;         // number of weeks observed
  }
  parameters{
  real mu;
  real<lower=0> sigma;
  vector[J] lambda;      // lambda_i of poisson 
  }
  transformed parameters {
  vector[J] theta;
  theta = N * lambda;
  }
  model{
  sigma ~ uniform(0.1, 0.5);
  lambda ~ lognormal(mu, sigma);
  X ~ poisson(theta);
  }
  generated quantities{
  real lamhat;
  real xpred;
  lamhat = lognormal_rng(mu, sigma);
  xpred = poisson_rng(6*lamhat);
  }")

set.seed(2020)
fit1 <- stan(
  model_code = sfile1,
  data = brkdwns,
  chains = 6,
  warmup = 2000,
  iter = 4000,
  cores = 2,
  refresh = 0
)

#estimates and CI
res <- summary(fit1, pars = "lambda")$summary # lambda estimates 
#xtable(res)

# joint dist of sigma and mu
pdf("q1sigmu.pdf", width = 10, compress = FALSE)
pairs(fit1, pars = c("sigma", "mu"))
dev.off()

#traceplots
pdf(file = 'q2sigma.pdf', compress = F, width = 10)
traceplot(fit1, pars = c("sigma"), include = TRUE, inc_warmup = TRUE)
dev.off()
#
pdf(file = 'q2lambda.pdf', compress = F, width = 10)
traceplot(fit1, pars = c("lambda"), include = TRUE, inc_warmup = TRUE)
dev.off()
#
pdf(file = 'q2mu.pdf', compress = F, width = 10)
traceplot(fit1, pars = c("mu"), include = TRUE, inc_warmup = TRUE)
dev.off()

# predictive dist
# extract pars
mu <- summary(fit1, pars = c("mu", "sigma"))$summary[1]
sigma <- summary(fit1, pars = c("mu", "sigma"))$summary[2]

# plotting 
nx <-1000000
lams <- rlnorm(n = nx, meanlog = mu, sdlog = sigma)
x7 <- rpois(n = nx, lambda = 6*lams)

pdf(file = 'q2preddist.pdf', compress = F, width = 10)
hist(x7, breaks = 60, freq = F, main = "", 
     xlab = expression("x"), col ="lightblue", 
     ylim = c(0,0.05), xlim = c(0, 70))
lines(density(x7), col = "red", lwd = 2)
dev.off()

# use generated values to verify
xpred <- extract(fit1, "xpred")
#pdf(file = 'q2preddist.pdf', compress = F, width = 10)
hist(xpred$xpred, xlim=c(0,100), breaks=60, 
     freq=FALSE, col = "lightblue", ylim = c(0,0.05))
lines(density(xpred$xpred), col = "red", lwd = 2)
#lines(density(x7), col = "red", lwd = 2)
#dev.off()

################
### sensitivity analysis 
# uninformative prior for sigma, N(1,1) for mu
sfile2<- c(
  "data{
  int<lower=0> J;         // number of facilities
  int<lower=0> X[J];      // number of breakdowns at facility J
  matrix[J,J] N;         // number of weeks observed
  }
  parameters{
  real mu;
  real<lower=0> sigma;
  vector[J] lambda;      // lambda_i of poisson 
  
  }
  transformed parameters {
  vector[J] theta;
  theta = N * lambda;
  }
  model{
  mu ~ normal(1,1);
  sigma ~ uniform(0.1, 1);
  lambda ~ lognormal(mu, sigma);
  X ~ poisson(theta);
  }")

set.seed(2020)
fit1 <- stan(
  #file = sfile.stan, # what goes here?
  model_code = sfile2,
  data = brkdwns,
  chains = 8,
  warmup = 2000,
  iter = 4000,
  cores = 2,
  refresh = 0
)

#estimates and CI
plot(fit1)
print(fit1)

#results 
res <- summary(fit1, pars = "lambda")$summary
#xtable(res)

plot(fit1, pars = c("lambda"))
plot(fit1, pars = c("sigma"))
pairs(fit1, pars = c("sigma", "mu"))

#traceplots to check convergence
pdf("trcase1.pdf", width = 10, compress = F)
traceplot(fit1, pars = c("sigma"))
dev.off()

traceplot(fit1, pars = c("lambda"))
traceplot(fit1, pars = c("mu"))

##### mu N(1,1), sigma gamma
#choosing Gamma pars
#which pars maximise probability between 0.1 and 0.5
max.prob <- function(pars){
  s <- pars[1]
  r <- pars[2]
  
  diff<-pgamma(q = 0.5, shape = s, rate = r) - pgamma(q = 0.1, shape = s, rate = r)
  return(-diff)
}
max <- nlm(f = max.prob, p = c(2,8))
s <- max$estimate[1] #106.1765
r <- max$estimate[2] #347.4632

xser <- seq(0, 0.6, 0.01)
par(mfrow=c(1,1))
plot(x=xser, y=sapply(X = xser, FUN = dgamma, shape = 2, rate =8),
     xlim=c(0,0.6), type = "l")

abline(v =0.1)
abline(v=0.5)
#what prob do we get 
-max.prob(c(2,8))
####

#data 
N <- c(3, 7, 4, 8, 5, 9)
X <- c(10, 33, 3, 39, 5, 50)

# data as named list
brkdwns <- list(J=length(N),X=X, N=diag(N), R=8, S=2)

sfile3 <- c(
  "data{
  int<lower=0> J;         // number of facilities
  int<lower=0> X[J];      // number of breakdowns at facility J
  matrix[J,J] N;         // number of weeks observed
  real<lower=0> R;      // pars for gamma prior of sigma
  real<lower=0> S;
  }
  parameters{
  real mu;
  real<lower=0> sigma;
  vector[J] lambda;      // lambda_i of poisson 
  }
  transformed parameters {
  vector[J] theta;
  theta = N * lambda;
  }
  model{
  mu ~ normal(1,1);
  sigma ~ gamma(S, 1/R); // gamma dist defined differently R vs Stan
  lambda ~ lognormal(mu, sigma);
  X ~ poisson(theta);
  }
  generated quantities{
  real lamhat;
  real xpred;
  lamhat = lognormal_rng(mu, sigma);
  xpred = poisson_rng(6*lamhat);
  }")

set.seed(2020)
fit3 <- stan(
  #file = sfile.stan, # what goes here?
  model_code = sfile3,
  data = brkdwns,
  chains = 6,
  warmup = 2000,
  iter = 4000,
  cores = 2,
  refresh = 0
)

#estimates and CI
res <- summary(fit3, pars = "lambda")$summary # lambda estimates 
#xtable(res)

##############################################
# Question 3
# data
X <- c(0.32, 0.63, 0.73, 0.38, 0.54, 0.95, 0.60, 1.03, 
       0.66, 0.41, 0.48, 0.27, 0.60, 0.21, 0.39, 0.28, 
       1.03, 0.68, 0.10, 0.69)

equipdata <- list(X = X, J = length(X), 
                  anu = -log(0.05)/log(3))
#using generated quantities 
equiprun2 <- c(
  "data{
  int<lower=0> J;         // number of observations
  real<lower=0> X[J];     // 
  real anu;               // alpha_nu constant
  }
  parameters{
  real<lower=0> nu;         // nu of weibull (failure rate, dist X)
  real sigma;               // sigma of weibull (dist X)
  }
  transformed parameters {
  real mu;                  // mean of x
  mu = sigma * tgamma(1 + 1/nu);
  }
  model{
  nu ~ pareto(1, anu);        // (y_min, alpha)
  sigma ~ uniform(0.2/(log(2)) ^ (1/nu), 1.5/(log(2)) ^ (1/nu));
  X ~ weibull(nu, sigma);     // (alpha, sigma)
  }
  generated quantities{
  real<lower=0> p;          // Pr[X <= 0,5]
  p = weibull_cdf(0.5, nu, sigma);
  }")

set.seed(2020)
fitequip2 <- stan(
  model_code = equiprun2,
  data = equipdata,
  chains = 6,
  warmup = 2000,
  iter = 4000,
  cores = 2,
  refresh = 0
)

# results 
#traceplots
#check convergence
traceplot(fitequip2, pars = "nu")
traceplot(fitequip2, pars = "sigma")

#estimates
xtable(summary(fitequip2)$summary)

# tracepplots for mu and p
pdf("q3mutrace.pdf", width = 10, compress = F)
par(mfrow=c(1,2))
traceplot(fitequip2, pars = "mu", include = TRUE, inc_warmup = TRUE)
dev.off()

pdf("qptrace.pdf", width = 10, compress = F)
traceplot(fitequip2, pars = "p", include = TRUE, inc_warmup = TRUE)
dev.off()

# required plots
# expectation and CI 
mures <- summary(fitequip2)$summary[3,]
pres <- summary(fitequip2)$summary[4,]
#xtable(rbind(mures, pres))

#histograms 
#mu and prob(X<0.5)
pdf("q3mu.pdf", width = 10, compress = F)
par(mfrow=c(1,2))
hist(extract(fitequip2)$mu, breaks = 40, col = "lightblue",
     xlab = expression(mu), freq = F, main = "", ylim = c(0,7))
lines(density(extract(fitequip2)$mu), col = "red")
hist(extract(fitequip2)$p, breaks = 40, col = "lightblue",
     xlab = expression(p), freq = F, main = "", ylim = c(0,7), xlim = c(0,1))
lines(density(extract(fitequip2)$p), col = "red")
dev.off()
