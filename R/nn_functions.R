#library(doParallel)

#' Computes the recommended sample sizes for a null hypothesis test
#'
#' @param s vector of length 1 (symmetric) or 2 (asymmetric) indicating
#' response standard deviations
#' @param d minimum detectable difference between treatments
#' @param conf 1 - type I error rate
#' @param power 1 - type II error rate
#' @param N finite deployment population, if NULL no finite population correction is used
#'
#' @return recommended sample sizes as a vector
#' @export
#'
#' @examples test_size_nht(s=c(0.5,0.10), d=0.2, conf=0.95, power=0.8, N=NULL)
#' test_size_nht(s=0.5, d=0.2, conf=0.95, power=0.8, N=NULL)
test_size_nht <- function(s, d, conf=0.95, power=0.8, N=NULL) {
  # computes the recommended sample size for a null hypothesis test
  # comparing two treatments with finite population correction
  # s is a vector of length 1 (symmetric) or 2 (asymmetric)
  #   indicating response std deviation(s)
  # d is the minimum detectable difference between treatments
  # conf is 1 - rate of type I error of the null hypothesis test
  # power is 1 - rate of type II error of the null hypothesis test
  # N is the finite population. If N=NULL, then no finite population correction is used.
  if (d <= 0) {warning("The minimum detectable difference between treatments should be positive")}
  z_alpha <- qnorm(1 - (1-conf)/2)
  z_beta <- qnorm(power)
  if(length(s) == 1) { # symmetric response variance
    if (is.null(N)) {
      out <- (z_alpha + z_beta)^2 * (2 * s^2) / d^2 #eq (2)
    } else {
      out <- (z_alpha + z_beta)^2 * (2 * s^2) * N /
        (d^2 * (N-1) + (z_alpha + z_beta)^2 * 4 * s^2)
    }
  } else { # asymmetric response variance
    if(is.null(N)) {
      n1 <- (z_alpha + z_beta)^2 * (s[1]^2 + s[1]*s[2]) / d^2
      n2 <- (z_alpha + z_beta)^2 * (s[1]*s[2] + s[2]^2) / d^2
      out <- c(n1, n2)
    } else {
      n1 <- (z_alpha + z_beta)^2 * N * (s[1]^2 + s[1]*s[2]) /
        (d^2 * (N-1) + (z_alpha + z_beta)^2 * (s[1] + s[2])^2)
      n2 <- (z_alpha + z_beta)^2 * N * (s[2]^2 + s[1]*s[2]) /
        (d^2 * (N-1) + (z_alpha + z_beta)^2 * (s[1] + s[2])^2)
      out <- c(n1, n2)
    }

  }
  out
}

# FUNCTIONS FOR 2-ARM TEST & ROLL (SYMMETRIC AND ASYMMETRIC) =====
#' Computes the per-customer profit for test & roll with 2 arms
#'
#' @param n vector of length 2 containing the sample sizes
#' @param N size of deployment population
#' @param s vector of length 2 containing the known standard deviations of the outcome
#' @param mu vector of length 2 containing the means of the prior on the mean response
#' @param sigma vector of length 2 containing the standard deviations of the prior on the mean response
#' @param log_n whether or not log(n) is an input rather than n (to avoid negative solutions), 'TRUE' or 'FALSE'
#'
#' @return per-customer profit for N customers
#' @export
#'
#' @examples profit_nn(n=100, N=10000, s=.1, mu=c(.7,.5), sigma=c(.2,.2))
#' profit_nn(n=c(100,200), N=10000, s=.1, mu=c(.7,.5), sigma=c(.2,.2))
profit_nn <- function(n, N, s, mu, sigma, log_n=FALSE) {
  # computes the per-customer profit for test & roll with 2 arms
  # where response is normal with (possibly asymmetric) normal priors
  # n is a vector of length 2 of sample sizes
  # N is the size of the deployment population
  # s is a vector of length 2 of the (known) std dev of the outcome
  # mu is a vector of length 2 of the means of the prior on the mean response
  # sigma is a vector of length 2 of the std dev of the prior on the mean response
  # if length(n)=1, equal sample sizes are assumed
  # if length(s)=1 symmetric priors are assumed and only the first elements of mu and sigma are used
  if (length(n) == 1) n <- rep(n, 2)
  if (log_n) n <- exp(n)

  # Catch error where priors indicate treatments are very different
  if (n[1]>N) stop("The algorithm recommends allocating almost all of the sample to treatment 1. This is probably because your priors indicate that there is a very high chance that treatment 1 performs better than treatment 2. If this is true, then you should deploy treatment 1 without an A/B test.")
  if (n[2]>N) stop("The algorithm recommends allocating almost all of the sample to treatment 2. This is probably because your priors indicate that there is a very high chance that treatment 2 performs better than treatment 1. If this is true, then you should deploy treatment 2 without an A/B test.")

  stopifnot(N >= sum(n), n[1]>0, n[2]>0, sum(sigma <= 0) == 0, sum(s <= 0) == 0)
  if (length(s) == 1) { # symmetric
    deploy <- (N - n[1] - n[2]) * (mu[1] + (2 * sigma[1]^2) /
                                     (sqrt(2*pi) * sqrt(s[1]^2 * (n[1] + n[2]) / (n[1]*n[2]) + 2 * sigma[1]^2))) #eq (9)
    test <- mu[1] * (n[1] + n[2])
  } else { # asymmetric
    e <- mu[1] - mu[2] #eq (18)
    v <- sqrt (sigma[1]^4 / (sigma[1]^2 + s[1]^2 / n[1] ) + sigma[2]^4 / (sigma[2]^2 + s[2]^2 / n[2] ) )
    deploy <- (N - n[1] - n[2]) * ( mu[2] + e * pnorm(e/v) + v * dnorm(e/v) )
    test <- mu[1] * n[1] + mu[2] * n[2]
  }
  #c(profit = test + deploy, profit_test = test, profit_deploy = deploy)
  (test + deploy)/N
}

#' Computes the profit-maximizing test size for test & roll with 2 arms
#' (asymmetric case not yet implemented)
#'
#' @param N size of deployment population
#' @param s vector of length 2 containing the standard deviations of the outcome
#' @param mu vector of length 2 containing the means of the prior on the mean response
#' @param sigma
#'
#' @return a vector containing the sample sizes
#' @export
#'
#' @examples test_size_nn(N=10000, s=.1, mu=c(.7,.7), sigma=c(.05,.05))
#' test_size_nn(N=10000, s=c(.1,.2), mu=c(.7,.7), sigma=c(.05,.05))
test_size_nn <- function(N, s, mu, sigma) {
  # computes the profit-maximizing test size for a 2-armed Test & Roll
  # where response is normal with normal priors (possibly asymmetric)
  # N is the size of the deployment population
  # s is a vector of length 2 of the (known) std dev of the outcome
  # mu is a vector of length 2 of the means of the prior on the mean response
  # sigma is a vector of length 2 of the std dev of the prior on the mean response
  # if length(s)=1 symmetric priors are assumed and only the first elements of mu and sigma are used
  if (any(mu <= 0)) {warning("The mean response should be positive")}
  stopifnot(N>2, sum(s <= 0) == 0, sum(sigma <= 0) == 0)
  if (length(s)==1) { # symmetric
    n <- ( -3 * s[1]^2 + sqrt( 9*s[1]^4 + 4*N*s[1]^2*sigma[1]^2 ) ) / (4 * sigma[1]^2) #eq (10)
    n <- rep(n, 2)
  } else { #asymmetric
    n <- optim(par=log(rep(round(N*0.05), 2 ) - 1), fn=profit_nn, control=list(fnscale=-1),
               N=N, s=s, mu=mu, sigma=sigma, log_n=TRUE)$par
    n <- exp(n)
  }
  n
}

#' Computes per-customer profit with perfect information
#' Where response is normal with symmetric normal priors
#'
#' @param mu means of the prior on the mean response
#' @param sigma standard deviations of the prior on the mean response
#'
#' @return per-customer profit with perfect information
#' @export
#'
#' @examples profit_perfect_nn(mu=.7, sigma=.02)
profit_perfect_nn <- function(mu, sigma){
  # computes the per-customer profit with perfect information
  # where response is normal with symmetric normal priors
  # todo: adapt to asymmetric case (closed form may not be possible)
  stopifnot(sigma > 0)
  (mu + sigma/sqrt(pi)) #eq (13)
}

#' Computes the rate of incorrect deployments with symmetric normal priors
#' Where response is normal with symmetric normal priors
#' @param n vector of length 2 containing the sample sizes
#' @param s vector of length 1 (symmetrical) containing the standard deviations of the outcome
#' @param sigma vector of length 1 (symmetrical) containing the standard deviations of the prior on the mean response
#'
#' @return rate of incorrect deployments
#' @export
#'
#' @examples error_rate_nn(n=100, s=.5, sigma=.2)
#' error_rate_nn(n=c(100,200), s=.5, sigma=.2)
error_rate_nn <- function(n, s, sigma) {
  # computes the rate of incorrect deployments
  # where response is normal with symmetric normal priors
  if (length(n)==1) n <- rep(n, 2)
  stopifnot(n[1]>0, n[2]>0, s>0, sigma>0)
  stopifnot(length(s)==1, length(sigma)==1)
  1/2 - 2*atan( sqrt(2)*sigma*sqrt(n[1]*n[2] / (n[1] + n[2]) ) / s ) / (2*pi) #eq (12)
}

#' Plots prior densities against mean response (profit per customer)
#'
#' @param mu means of the prior on the mean response
#' @param sigma standard deviations of the prior on the mean response
#'
#' @return graph plotting prior density vs. profit per customer
#' @export
#'
#' @examples #plot_prior_mean_resp_nn(mu=c(5,10), sigma=c(10,10))
plot_prior_mean_resp_nn <- function(mu, sigma) {
  if (length(mu)==1) {mu <- rep(mu, 2); sigma <- rep(sigma,2)}
  xlim <- c(min(mu[1] - 4*sigma[1], mu[2] - 4*sigma[2]),
            max(mu[1] + 4*sigma[1], mu[2] + 4*sigma[2]))
  ylim <- c(0, 1.05*max(dnorm(mu[1], mean=mu[1], sd=sigma[1]),
                        dnorm(mu[2], mean=mu[2], sd=sigma[2])))
  xs <- seq(xlim[1], xlim[2], length.out=100)
  plot(xs, dnorm(xs, mean=mu[1], sd=sigma[1]), type="l", col="red", xlim=xlim, ylim=ylim,
       main="Prior on Mean Response", xlab="Profit Per Customer", ylab="Prior Density")
  lines(xs, dnorm(xs, mean=mu[2], sd=sigma[2]), type="l", col="blue", lty=3)
  legend("topright", legend=c("Treatment 1", "Treatment 2"), col=c("red", "blue"), lty=c(1,3), bty="n")
}

#' Plots prior densities against response (profit per customer)
#'
#' @param s known standard deviations of the response
#' @param mu means of the prior on the mean response
#' @param sigma standard deviations of the prior on the mean response
#'
#' @return graph plotting prior density vs. profit per customer
#' @export
#'
#' @examples plot_prior_resp_nn(s=c(10,20), mu=c(5,10), sigma=c(10,10))
plot_prior_resp_nn <- function(s, mu, sigma) {
  if (length(mu)==1) {mu <- rep(mu, 2); sigma <- rep(sigma,2); s <- rep(s, 2)}
  sigma_resp1 <- sqrt(sigma[1]^2 + s[1]^2)
  sigma_resp2 <- sqrt(sigma[2]^2 + s[2]^2)
  xlim <- c(min(mu[1] - 4*sigma_resp1, mu[2] - 4*sigma_resp2),
            max(mu[1] + 4*sigma_resp1, mu[2] + 4*sigma_resp2))
  ylim <- c(0, 1.05*max(dnorm(mu[1], mean=mu[1], sd=sigma_resp1),
                        dnorm(mu[2], mean=mu[2], sd=sigma_resp2)))
  xs <- seq(xlim[1], xlim[2], length.out=100)
  plot(xs, dnorm(xs, mean=mu[1], sd=sigma_resp1), type="l", col="red", xlim=xlim, ylim=ylim,
       main="Prior on Response", xlab="Profit Per Customer", ylab="Prior Density")
  lines(xs, dnorm(xs, mean=mu[2], sd=sigma_resp2), type="l", col="blue", lty=3)
  legend("topright", legend=c("Treatment 1", "Treatment 2"), col=c("red", "blue"), lty=c(1,3), bty="n")
}

#' Plots prior densities against treatment effect (difference in profit per customer)
#'
#' @param mu means of the prior on the mean response
#' @param sigma standard deviations of the prior on the mean response
#' @param abs whether or not to take the absolute difference, 'TRUE' or 'FALSE'
#'
#' @return graph plotting prior density vs. difference in profit per customer
#' @export
#'
#' @examples plot_prior_effect_nn(mu=c(5,10), sigma=c(10,10), abs=FALSE)
plot_prior_effect_nn <- function(mu, sigma, abs=FALSE) {
  if (length(mu)==1) {mu <- rep(mu, 2); sigma <- rep(sigma,2)}
  sigma_effect <- sqrt(sigma[1]^2 + sigma[2]^2)
  xlim <- c(mu[2] - mu[1] - 4*sigma_effect, mu[2] - mu[1] + 4*sigma_effect)
  if (!abs) {
    xs <- seq(xlim[1], xlim[2], length.out=100)
    plot(xs, dnorm(xs, mean=mu[2] - mu[1], sd=sigma_effect), type="l", xlim=xlim,
         main="Prior on Treatment Effect (m2 - m1)",
         xlab="Difference in Profit per Customer", ylab="Prior Density")
  } else {
    xs <- seq(0, max(-xlim[1], xlim[2]), length.out=100)
    ys <- dnorm(xs, mean=mu[2] - mu[1], sd=sigma_effect) +
      dnorm(-xs, mean=mu[2] - mu[1], sd=sigma_effect)
    plot(xs, ys, type="l", main="Prior on Treatment Effect |m2 - m1|",
         xlab="Absolute Difference in Profit Per Customer", ylab="Prior Density")
  }
}

# FUNCTIONS FOR K-ARM TEST & ROLL (REQUIRES SIMULATION) =====
#' Utility function used in function 'profit_nn_sim()' to simulate one set of potential outcomes
#' Draws a true mean for each arm and generates N observations from each arm
#' Returns profits and error rates under perfect information, test & roll, and Thomson sampling
#'
#' @param n vector of length K containing the sample sizes
#' @param N the deployment population
#' @param s vector of length K containing the standard deviations of the outcome
#' @param mu vector of length K containing the means of the prior on the mean response
#' @param sigma vector of length K containing the standard deviations of the prior on the mean response
#' @param K the number of arms (treatments)
#' @param TS whether or not to run Thomson sampling, 'TRUE' or 'FALSE'
#'
#' @return The profits and error rates under perfect information, test & roll, and Thomson sampling under one simulation
#' @export
#'
#' @examples one_rep_profit(n=c(100,100), N=1000, s=c(.1,.1), mu=c(.1,.1), sigma=c(.05,.05), K=2, TS=FALSE)
#' one_rep_profit(n=c(100,200,300), N=1000, s=c(.1,.2,.3), mu=c(.1,.2,.3), sigma=c(.01,.03,.05), K=3, TS=FALSE)

one_rep_profit <-function(n, N, s, mu, sigma, K, TS=FALSE) {
  if (any(mu <= 0)) {warning("The mean response should be positive")} # warning: mu <= 0
  stopifnot(N >= sum(n), all(n > 0), sum(s <= 0) == 0, sum(sigma <=0) == 0) # input validation
  # utility function used in profit_nn_sim() to simulate one set of potential outcomes
  m <- rnorm(K, mu, sigma) # draw a true mean for each arm
  y <- matrix(rnorm(N*K, m, s), nrow=N, ncol=K, byrow=TRUE) # N observations from each arm

  # perfect information profit
  perfect_info <- sum(y[,which.max(m)]) # Perfect information: sum observations from arm with highest m

  # test and roll profit with sample sizes n
  postmean <- rep(NA, K)
  for (k in 1:K) {
    postvar <- 1/(1/sigma[k]^2 + n[k]/s[k]^2)
    postmean[k] <- postvar*(mu[k]/sigma[k]^2 + sum(y[1:n[k], k])/s[k]^2) #TODO, input validation
  }
  delta <- which.max(postmean) # pick the arm with the highest posterior mean
  error <- delta != which.max(m)
  deploy_1 <- delta == 1
  test_roll <- 0
  for (k in 1:K) test_roll <- test_roll + sum(y[1:n[k], k]) # profit from first n observations for each arm
  test_roll <- test_roll + sum(y[(sum(n)+1):N, delta]) # profit from remaining observations for selected arm

  # Thompson sampling profit
  thom_samp <- NA
  if (TS==TRUE) {
    n <- rep(0, K) # Initialize each arm with zero observations
    # note mu and sigma are initialized at the priors
    postvar <- t(1/(1/sigma^2 + t(matrix(1:N, nrow=N, ncol=K)) * 1/s^2))
    postmean <- postvar * t(mu/sigma^2 + t(apply(y, 2, cumsum)) * 1/s^2)
    for (i in 1:N) {
      k <- which.max(rnorm(K, mu, sigma)) # draw a sample from the current posterior for each arm and choose arm with largest draw
      n[k] <- n[k] + 1 # increase the number of observations used from sampled arm
      mu[k] <- postmean[n[k], k] # advance mu and sigma
      sigma[k] <- sqrt(postvar[n[k], k])
    }
    thom_samp <- 0
    for (k in 1:K) thom_samp <- thom_samp + sum(y[1:n[k], k])
  }
  #return(c(perfect_info=perfect_info, test_roll=test_roll, thom_samp=thom_samp,
           #error = error, deploy_1=deploy_1))
  return(c(perfect_info=perfect_info, test_roll=test_roll, thom_samp=thom_samp,
           error = error))
}

#' Computes the per-customer profit for test & roll with K arms
#'
#'
#' @param n sample sizes for test & roll (vector length 1 or K)
#' @param N deployment population
#' @param s standard deviations of the outcome (vector of length 1 or K)
#' @param mu means of the priors on the mean response (vector of length 1 or K)
#' @param sigma standard deviations of the priors on the mean response (vector of length 1 or K)
#' @param K number of arms
#' @param TS whether or not to run Thomson sampling, 'TRUE' or 'FALSE'
#' @param R number of simulation replications
#'
#' @return A list containing the profit, regret, and error rates
#' @export
#'
#' @examples profit_nn_sim(n=c(100,100), N=1000, s=c(.1,.1), mu=c(.1,.1), sigma=c(.05,.05))
#' profit_nn_sim(n=c(100,200,300), N=1000, s=c(.1,.2,.3), mu=c(.1,.2,.3), sigma=c(.01,.03,.05), K=3, TS=FALSE, R=100)
#'
profit_nn_sim <- function(n, N, s, mu, sigma, K=2, TS=FALSE, R=1000) {
  if (any(mu <= 0)) {warning("The mean response should be positive")} # warning: mu > 0
  stopifnot(N >= sum(n), all(n > 0), sum(s <= 0) == 0, sum(sigma <=0) == 0) # input validation
  # computes the per-customer profit for test & roll with K arms
  # where response is normal with (asymmetric) normal priors
  # R is the number of simulation replications
  # n is the sample size for test & roll
  # N is the size of the total population
  # s are the (known) std devs of the outcome (vector of length 1 or K)
  # mu are the means of the priors on the mean response (vector of length 1 or K)
  # sigma are the std devs of the priors on the mean response (vector of length 1 or K)
  # TS is a switch for computing profit for Thompson Sampling
  # if s is length 1, then arms are assumed to be symmetric
  if (length(s) == 1) { s <- rep(s, K); mu <- rep(mu, K); sigma <- rep(sigma, K) }
  if (length(n) == 1) n <- rep(n, K)
  stopifnot(length(s) == K, length(mu) == K, length(sigma) == K)
  stopifnot(sum(sigma <= 0) == 0, sum(s <= 0) == 0)
  stopifnot(N > K)
  reps <- foreach(i=1:R) %dopar% one_rep_profit(n, N, s, mu, sigma, K, TS)
  reps <- as.data.frame(do.call(rbind, reps))
  profit <- apply(reps[,1:3], 2, mean) #column 1 (perfect info), column 2 (test&roll), column 3 (thomson sampling)
  profit <- rbind(exp_profit=profit,
                  apply(reps[,1:3], 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE))
  profit <- profit / N
  regret_draws <- 1 - reps[,1:3] / reps$perfect_info
  regret <- apply(regret_draws, 2, mean)
  regret <- rbind(exp_regret=regret,
                  apply(regret_draws, 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE))
  error <- mean(reps[,4])
  #deploy_1 <- mean(reps[,5])
  return(list(profit=profit, regret=regret, error_rate=error, #deploy_1_rate = deploy_1,
              profit_draws=reps, regret_draws=regret_draws))
}

#' Utility function used in test_size_nn_sim() to simulate one set of potential outcomes
#' Returns profits for all possible equal sample sizes
#'
#'
#' @param n_vals potential values for n (1:(floor(N/K)-1))
#' @param N deployment population
#' @param s standard deviations of the outcome
#' @param mu means of the priors on the mean response
#' @param sigma standard deviations of the priors on the mean response
#' @param K number of arms (treatments)
#'
#' @return A 2-column matrix with values of n in the first column and profits in the second column
#' @export
#'
#' @examples one_rep_test_size(1:(floor(10/2)-1), N=10, s=c(10,10), mu= 20, sigma=10, K=2)
one_rep_test_size <- function(n_vals, N, s, mu, sigma, K) {
  if (any(mu <= 0)) {warning("The mean response should be positive")} # warning: mu > 0
  stopifnot(sum(s <= 0) == 0, sum(sigma <=0) == 0) # input validation
  # utility function used in test_size_nn_sim() to simulate one set of potential outcomes
  # and profits for all possible equal sample sizes

  # potential outcomes
  m <- rnorm(K, mu, sigma) # draw a true mean for the arm
  y <- matrix(rnorm(N*K, m, s), nrow=N, ncol=K, byrow=TRUE) # N observations from each arm

  postmean <- matrix(NA, nrow=length(n_vals), ncol=K)
  for (k in 1:K) {
    postvar <- 1/(1/sigma[k]^2 + n_vals/s[k]^2)
    postmean[,k] <- postvar*(mu[k]/sigma[k]^2 + cumsum(y[1:(floor(N/K)-1),k])/s[k]^2)
  }
  delta <- apply(postmean, 1, which.max) # pick the arm with the highest posterior mean
  error <- delta != which.max(m)
  profit <- rep(0, length(n_vals))
  for (i in seq_along(n_vals)) {
    for (k in 1:K) profit[i] <- profit[i] + sum(y[1:n_vals[i], k]) # profit from first n[i] observations for each arm
    profit[i] <- profit[i] + sum(y[(n_vals[i]*K + 1):N, delta[i]]) # profit from remaining observations for selected arm
  }
  return(cbind(n=n_vals, profit))
}

#' Computes the profit-maximizing test size and profits for a multi-armed test & roll
#'
#' @param N deployment population
#' @param s standard deviations of the response (length 1(symmetric) or K)
#' @param mu vector of length K containing means of the priors on the mean response
#' @param sigma vector of length K containing the standard deviations of the priors on the mean response
#' @param K number of arms (treatments)
#' @param R number of simulation replications
#'
#' @return a list with the sample sizes and expected profit per customer
#' @export
#'
#' @examples test_size_nn_sim(N=1000, s=.1, mu=.1, sigma=.05, K=2, R=1000)
#'
test_size_nn_sim <- function(N, s, mu, sigma, K=2, R=1000) {
  if (any(mu <= 0)) {warning("The mean response should be positive")} # warning: mu > 0
  # computes the profit-maximizing test size for a multi-armed test & roll
  # where response is normal with normal priors (possibly asymmetric)
  # N is the size of the deployment population
  # K is the number of arms
  # s is a K-vector of (known) std devs of the response, if length(s)==1, symmetric priors are used
  # mu is a K-vector of length K the means of the priors on the outcome
  # sigma is a K-vector of std devs of the priors on the mean response
  stopifnot(N > 2, sum(s <= 0) == 0, sum(sigma <= 0) == 0)
  if (length(s==1)) { # symmetric arms
    # n is same for all arms; solve by enumeration
    s <- rep(s, K); mu <- rep(mu, K); sigma <- rep(sigma, K)
    n_vals <- 1:(floor(N/K)-1) # potential values for n
    reps <- foreach(i=1:R) %dopar% one_rep_test_size(n_vals, N, s, mu, sigma, K)
    reps <- as.data.frame(do.call(rbind, reps))
    exp_profit <- xtabs(profit ~ n, data=reps) / R
    n <- rep(n_vals[which.max(exp_profit)], K)
  } else { # asymmetric
    stopifnot(length(mu)==K, length(s) == K, length(sigma) == K)
    # todo:finish this
    # best option is to use a numeric optimization, but then noise in profit_nn_sim
    # becomes critical. comments out is one strategy
    # start values based on two-arm symmetric
    #n <- ( -3 * s[1]^2 + sqrt( 9*s[1]^4 + 4*N*s[1]^2*sigma[1]^2 ) ) / (4 * sigma[1]^2)
    #n <- optim(par=log(n), fn=profit_nn_sim, control=list(fnscale=-1),
    #           N=N, s=s, mu=mu, sigma=sigma, K=K, R=1000, log_n=TRUE)$par
    #n <- optim(par=n, fn=profit_nn_sim, control=list(fnscale=-1),  # more precise simulation
    #           N=N, s=s, mu=mu, sigma=sigma, K=K, R=5000, log_n=TRUE)$par
    n <- NA
  }
  return(list(n=n, max(exp_profit)/N))
}

# SUMMARY FUNCTION =====
#' Provides summary of a test & roll plan
#'
#' @param n vector of length 2 containing the sample sizes
#' @param N deployment population
#' @param s known standard deviations of the outcome
#' @param mu means of the priors on the mean response
#' @param sigma standard deviations of the priors on the mean response
#'
#' @return a data frame containing summary statistics such as profit per customer,
#' profits from test phase, error rates, etc.
#' @export
#'
#' @examples test_eval_nn(n=c(100,100), N=1000, s=.1, mu=.1, sigma=.05)
#' test_eval_nn(n=c(100,200), N=1000, s=c(.1,.2), mu=c(.1,.2), sigma=c(.05,.1))
test_eval_nn <- function(n, N, s, mu, sigma) {
  # provides a complete summary of a test & roll plan
  # n is a vector of length 2 of sample sizes
  # N is the size of the deployment population
  # s is a vector of length 2 of the (known) std dev of the outcome
  # mu is a vector of length 2 of the means of the prior on the mean response
  # sigma is a vector of length 2 of the std dev of the prior on the mean response
  # if length(n)=1, equal sample sizes are assumed
  # if length(s)=1 symmetric priors are assumed and only the first elements of mu and sigma are used
  if (any(mu <= 0)) {warning("The mean response should be positive")}
  stopifnot(N >= n[1] + n[2], n[1] > 0, n[2] >0, sum(s <= 0) == 0, sum(sigma <=0) == 0)
  profit <- profit_nn(n, N, s, mu, sigma)*N
  if (length(s)==1) { # symmetric
    test <- mu[1] * (n[1] + n[2])
    deploy <- profit - test
    rand <- mu[1]*N # choose randomly
    perfect <- profit_perfect_nn(mu, sigma)*N
    error_rate <- error_rate_nn(n, s, sigma)
    deploy_1_rate <- 0.5
    gain <- profit - rand
  } else { # asymmetric
    test <- mu[1] * n[1] + mu[2] *n[2]
    deploy <- profit - test
    rand <- ((mu[1] + mu[2])*0.5)*N
    out <- profit_nn_sim(n, N, s, mu, sigma, R=10000)
    perfect <- out$profit["exp_profit", "perfect_info"]*N
    error_rate <- out$error_rate
    gain <- profit - rand
    #deploy_1_rate <- out$deploy_1_rate
  }
  gain <- (profit - rand) / (perfect - rand)
  data.frame(n1=n[1], n2=n[2],
             profit_per_cust = profit/N,
             profit = profit,
             profit_test = test,
             profit_deploy = deploy,
             profit_rand = rand,
             profit_perfect = perfect,
             profit_gain = gain,
             regret = 1-profit/perfect,
             error_rate = error_rate,
             #deploy_1_rate = deploy_1_rate,
             tie_rate = 0
             )
}
