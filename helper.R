# Author: Gord Burtch
# Content: Helper functions for clustered regressions simulations.

## Let's load some libraries.
library('arm')
library('simstudy')
library('mvtnorm')
library('lme4')
library('multiwayvcov')
library('clusterSEs') # for clustered SEs and wild bootstrap clustered SEs.
library('ggplot2')
library('dplyr')
library('haven')
library('RColorBrewer')
library('ggthemes')
library('fishmethods') # for intra-cluster correlations.
library('fwildclusterboot')
library('patchwork')
library('viridis')

### FUNCTION TO SIMULATE A DATA FRAME COMPRISED OF CLUSTERED DATA

# This function simulates clustered data.
# n = observations.
# n_cluster = number of clusters they are associated with.
# rho = intra-cluster correlation.
# param is a set of two values that reflect our regression beta coefficients.

# This function is not perfect, by any means. 
# The shared cluster contribution to variance in x and e leads to a correlation between x and e,
# which can bias estimates?
gen_cluster <-
  function(param = c(.1, .5),
           n = 1000,
           n_cluster = 50,
           rho = .5) {
    # Required package: mvtnorm
    
    while (n %% n_cluster != 0){
      n <- n+1
    }
    
    # Here we are simulating two vectors, from a bivariate normal distribution.
    # These two variables will serve as our variable, X, and the error term.
    # The two vectors (variables) will have means of 0, and will not covary, i.e., we don't have endogeneity here!
    
    # NOTE about how rho works in this setup:
    
    # The simulation will essentially make the variance (sd) of x and e will both depend on cluster membership.
    # However, notice that x and e are simulated to be uncorrelated with one another - we do not have OVB problem.  
    
    # The individual level contribution to the error will have a variance of 1-rho.
    # So, if rho = 1, the individual level draw essentially contributes nothing to X.
    # As you will see below, cluster contribution will drive all the variation in the observed variable in that situation.
    Sigma_i <- matrix(c(1, 0, 0, 1 - rho), ncol = 2)
    values_i <- rmvnorm(n = n, sigma = Sigma_i)
    x_i <- values_i[, 1]
    error_i <- values_i[, 2]
    
    # Now we are going to assign the above datapoints to clusters, evenly.
    cluster_name <- rep(1:n_cluster, each = n / n_cluster)
    
    # For X and error, simulated above, we assign each observation to a cluster.
    # The assigned cluster will then shift each observation's pair X and error by a pair of values drawn from another bivariate normal distribution.
    # The cluster-specific contribution to error will have a variance of rho, again just because?
    # We are now simulating distributional draws that will serve as the cluster-specific contributions to X and error.
    Sigma_cl <- matrix(c(1, 0, 0, rho), ncol = 2)
    values_cl <- rmvnorm(n = n_cluster, sigma = Sigma_cl)
    x_cl <- rep(values_cl[, 1], each = n / n_cluster)
    error_cl <- rep(values_cl[, 2], each = n / n_cluster)
    
    # Now we stick together the individual level values and their cluster-specific contributions.
    x <- x_i + x_cl
    error <- error_i + error_cl
    
    # Finally, we simulate y as a function of our betas, and the X and error we just simulated.
    y <- param[1] + param[2] * x + error
    
    df_tmp <-
      data.frame(x, y, cluster = cluster_name, x_i, x_cl, error, error_i, error_cl)
  return(df_tmp)
}

### SINGLE ITERATION OF THE SIM: CREATE A SAMPLE, ESTIMATE MODEL, STORE RESULTS

# This function nests the data generation function.
# For each simulation iteration, simulate clustered data with the provided params.
# Then estimate the linear regression model, store the beta, and the standard errors + confidence intervals.
# Calculate the SEs naively, with regular clustering, or fast and wild cluster bootstrap (good for few clusters).
cluster_sim <- function(param = c(.1, .5),
                        n = 1000,
                        n_cluster = 50,
                        rho = .5,
                        SEs = "regular", FE = FALSE) {
  # Required packages: multiwayvcov, fwildclusterboot
  
  # Note we have to do some weirdness here because of the way the boottest function is coded.
  # We need df_it to be a global variable, because boottest looks for it implicitly in global memory.
  df_it <<-
    gen_cluster(
      param = param,
      n = n ,
      n_cluster = n_cluster,
      rho = rho
    )
  
  x_err_cor <- cor(df_it$x,df_it$error)
  
  if (FE == TRUE) {
    fit <- lm(data = df_it, y ~ x + factor(cluster))
  } else {
    fit <- lm(data = df_it, y ~ x)
  }
  
  b1 <- coef(fit)[2]
  
  if (SEs == "regular") {
    Sigma <- vcov(fit)
    b1_ci95 <- confint(fit)[2, ]
  } else if (SEs == "clustered") {
    Sigma <- cluster.vcov(fit, ~ cluster)
    se <- sqrt(diag(Sigma)[2])
    t_critical <- qt(.025, df = n - 2, lower.tail = FALSE)
    lower <- b1 - t_critical * se
    upper <- b1 + t_critical * se
    b1_ci95 <- c(lower, upper)
  } else if (SEs == "fast_n_wild") {
    # Per Roodman et al. 2019, employ Webb weights if cluster count <= 12.
    if (n_cluster <= 12) {
      weights = "webb"
    } else {
      weights = "rademacher"
    }
    boot <-
      boottest(
        fit,
        B = 999,
        param = "x",
        clustid = "cluster",
        type = weights
      )
    lower <- boot$conf_int[1]
    upper <- boot$conf_int[2]
    b1_ci95 <- c(lower, upper)
  }
  
  return(c(b1, b1_ci95, x_err_cor))
}

# Lastly, a function to iterate / repeat the simulation.
# A data frame is returned containing all the betas, se's, CIs, etc.
# We are also creating, for each sim, a binary indicator of whether the CIs contain the true parameter for b1.
run_cluster_sim <-
  function(n_sims = 1000,
           param = c(.1, .5),
           n = 1000,
           n_cluster = 50,
           rho = .5,
           SEs = "regular", FE = FALSE) {
    # Required packages: dplyr
    
    df_res <-
      replicate(n_sims,
                cluster_sim(
                  param = param,
                  n = n,
                  rho = rho,
                  n_cluster = n_cluster,
                  SEs = SEs,
                  FE = FE
                ))
    
    df_res <- as.data.frame(t(df_res))
    names(df_res) <- c('b1', 'ci95_lower', 'ci95_upper','x_err_cor')
    df_res <- df_res %>%
      mutate(id = 1:n(),
             param_caught = ci95_lower <= param[2] &
               ci95_upper >= param[2])
    return(df_res)
  }
