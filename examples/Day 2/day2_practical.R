#################################################
### Bayesian causal inference day 2 practical ###
#################################################

set.seed(1)
library(BART)
setwd("~/Downloads/LeuvenBCI-main") # change filepath accordingly
load("data/lalonde.RData")

## define variables
Y <- "re78"
treat <- "treat"
covar <- c("age", "education", "black", "hispanic", "married",
           "nodegree", "re74", "re75", "u74", "u75")
n <- nrow(ldw_psid) # sample size

## benchmark estimate: difference-in-means estimator on experimental data
diff_means <- mean(ldw[ldw$treat == 1,Y]) - mean(ldw[ldw$treat == 0,Y])

## run probit-link BART regression to fit propensity score model
Mskip <- 2000 # number of burn-in samples
M <- 5000 # number of posterior samples
pi_post <- pbart(x.train = ldw_psid[,covar], y.train = ldw_psid[,treat], ndpost = M, nskip = Mskip, printevery = 500L)

## extract posterior samples of pi
pi_draws <- pi_post$prob.train

## prepare test data for fitting BART outcome regression model
test_dat_tr <- cbind(rep(1, nrow(ldw_psid)), ldw_psid[,covar])
test_dat_co <- cbind(rep(0, nrow(ldw_psid)), ldw_psid[,covar])
names(test_dat_tr)[1] <- "treat"
names(test_dat_co)[1] <- "treat"
test_dat <- rbind(test_dat_tr, test_dat_co)

## fit BART outcome regression model
mu_post <- wbart(x.train = ldw_psid[,c(treat,covar)], y.train = ldw_psid[,Y], x.test = test_dat, ndpost = M, nskip = Mskip, printevery = 500L)

## extract posterior samples for mu
mu_draws_tr <- mu_post$yhat.test[,1:n]
mu_draws_co <- mu_post$yhat.test[,(n+1):(2*n)]

###################################
### Posterior inference for ATE ###
###################################

## Question (i) 
# construct CATE posterior sample for BART
# compute the posterior mean and 95% central credible interval



## obtain one-step posterior draws for the ATE
# trim propensity score posterior draws
pi_trim <- pi_draws
pi_trim[pi_trim < 0.05] <- 0.05
pi_trim[pi_trim > 0.95] <- 0.95

chi_tilde <- numeric(M)
for (b in 1:M) {
  ## generate dirichlet weights by normalizing iid exponential weights
  wts <- rexp(n)
  wts <- wts/sum(wts)
  
  ## evaluate corrected parameter for ATE
  chi_tilde[b] <- sum(wts*(ldw_psid[,treat]*(ldw_psid[,Y]-mu_draws_tr[b,])/pi_trim[b,]+mu_draws_tr[b,] - (1-ldw_psid[,treat])*(ldw_psid[,Y] - mu_draws_co[b,])/(1-pi_trim[b,])-mu_draws_co[b,]))
}

## Question (ii) 
# compute the posterior mean and 95% central credible interval for chi_tilde


##################################
### Evaluate degree of overlap ###
##################################

## compute posterior mean of pi from probit BART model
pi_postmean <- colMeans(pi_draws)

ldw_psid_prop <- ldw_psid
ldw_psid_prop$pi_est <- pi_postmean

## define overlap histogram plot function from 
plot_hist <- function(data, var, treat, main = NULL, odds = FALSE,
                      breaks = 40, density = TRUE, xlim = NULL, ylim = NULL,
                      xlab = NULL, text.size = 0.8) {
  ntr <- sum(data[, treat] == 1)
  nco <- sum(data[, treat] == 0)
  if (odds == TRUE) {
    data[, var] <- log(data[, var]/(1-data[, var]))
    if (is.null(xlab) == TRUE) {xlab <- "Log Odds"}
  } else {
    if (is.null(xlab) == TRUE) {xlab <- "Propensity Score"}
  }
  if (is.null(xlim)) {
    if (odds == TRUE) {
      xlim <- range(data[, var])
      cat(xlim)
    } else {
      xlim <- c(0,1)
    }
  }
  intervals <- seq(xlim[1], xlim[2], length.out = breaks + 1)
  h0 <- as.numeric(table(cut(data[data[,treat]==0, var],
                             breaks = intervals, include.lowest = TRUE)))
  h1 <- as.numeric(table(cut(data[data[,treat]==1, var],
                             breaks = intervals, include.lowest = TRUE)))
  if (density == TRUE) {
    h0 <- h0/sum(h0); h1 <- h1/sum(h1)
  }
  s <- cbind.data.frame(h0, h1)
  if (is.null(ylim)) {
    ylim.max <- max(s$h0, s$h1) * 1.2
    ylim <- c(-ylim.max, ylim.max)
  }
  par(mar = c(4, 4, 1, 1))
  barplot(s$h0 * -1, names.arg = sprintf("%.2f", intervals[-1]),
          col = "#AAAAAA80", main = main, cex.lab = 1.3,
          ylim = ylim, xlab = xlab, cex.axis = 1.2, cex.names = 1.2,
          ylab = "Density", border = NA, axes = TRUE)
  barplot(s$h1, col = "#ff000080", add = TRUE,
          border = NA, cex.axis = 1.2)
  abline(h = 0, col = "gray60", lty = 2, lwd = 1.5)
  axis(1, at = seq(1, 60, length.out = breaks/2), labels = FALSE)
  usr <- par("usr")
  user_x <- usr[1] + 0.03 * (usr[2] - usr[1])
  user_y <- usr[3] + 0.92 * (usr[4] - usr[3])
  text(user_x, user_y, paste("Ntr = ", ntr), pos = 4, cex = text.size)
  text(user_x, user_y - 0.05 * (usr[4] - usr[3]), paste("Nco = ", nco),
       pos = 4, cex = text.size)
  box()
}

# overlap plot for propensity score
plot_hist(ldw_psid_prop, "pi_est", treat, odds = FALSE, breaks = 50,
          density = TRUE, main = "", xlim = c(0,1), ylim = NULL)

# overlap plot on log odds scale
plot_hist(ldw_psid_prop, "pi_est", treat, odds = TRUE, breaks = 50,
          density = TRUE, main = "", xlim = NULL, ylim = NULL)

###################################
### Posterior inference for ATT ###
###################################

## Question (iv)
# obtain posterior samples of CATT
# compute the posterior mean and 95% central credible interval

## Question (v)
# obtain one-step posterior draws for the ATT
# trim propensity score posterior draws

## Question (vi)
# Posterior inference for ATT with estimated propensity score
