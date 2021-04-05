# Author: Gord Burtch
# Date: March 31st, 2021
# Credit to Yuki Yanai for basic helper functions, slightly modified to accommodate fast and wild cluster bootstrap errors:
# http://yukiyanai.github.io/teaching/rm1/contents/R/clustered-data-analysis.html

# Set a seed so you replicate the same random draws / results.
set.seed(1234)

# Load helper functions: one to generate a clustered data set for y~x;
# One to estimate an lm() and recover estimates with requested SEs; 
# One to nest nest draws / estimations, to implemented a set of simulations.
source("helper.R")

# Here's one simulated dataset.
df_test <- gen_cluster(n = 10000, n_cluster = 50, rho = 0.7)

# Here's what 'x' looks like ignoring cluster assignment.
# Nice and normal, right?
ggplot(data = df_test, aes(x = y)) +
  geom_density(fill = "red", alpha = 0.5) +
  xlab("True Error") + ylab("Density") +
  scale_fill_manual(name = "Cluster", labels = as.character(seq(1:50)), values = mypalette) +
  theme_bw() +
  theme(text = element_text(family = "Economica")) +
  NULL

# In reality, what we have is something like mixture of many normals.
mypalette <- colorRampPalette(brewer.pal(11, "Spectral"))(50)
ggplot(data = df_test, aes(x = error, group = cluster)) +
  geom_density(aes(fill = factor(cluster)), alpha = 0.5) +
  xlab("True Error") + ylab("Density") +
  scale_fill_manual(name = "Cluster", labels = as.character(seq(1:50)), values = mypalette) +
  theme_bw() +
  theme(text = element_text(family = "Economica")) +
  NULL

# Here is the intra-cluster correlation coefficient for x and the true error.
clus.rho(df_test$x, df_test$cluster)
clus.rho(df_test$error, df_test$cluster)
clus.rho(df_test$y, df_test$cluster)

### SIMULATIONS

# Let's start with un-clustered data and see how we do.
# We have 50 clusters by default, but in the un-clustered sim that cluster contributes nothing to x / error.
# That is, the intracluster correlation, rho, will be 0. In the 'clustered' case we set intra-cluster correlation to 0.7.
sim_unclustered <- run_cluster_sim(rho = 0, SEs = "regular")
sim_clustered <- run_cluster_sim(rho = 0.7, SEs = "regular")

# Let's look at our distribution of estimates from the unclustered simulations.
ci95_unclustered <- ggplot(
  sim_unclustered,
  aes(
    x = as.numeric(id),
    y = b1,
    ymin = ci95_lower,
    ymax = ci95_upper,
    color = param_caught
  )
) +
  geom_errorbar(alpha = 0.5, color = "gray") +
  geom_point() +
  labs(x = 'Simulation ID', y = 'Estimate (95% CI)', title = 'Simulation Estimates with Regular SEs (Rho = 0)') +
  scale_color_manual(
    name = 'CIs Contain Estimate?',
    labels = c('Yes', 'No'),
    values = c("red", "blue")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Economica"),
    axis.title = element_text(size = 16, margin = margin(
      t = 10,
      b = 0,
      r = 20,
      l = 0
    )),
    axis.text = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0.5,
             linetype = 'dashed',
             color = "chartreuse4") +
  NULL

# Now, How did we do when clustering was present? Badly!
# Let's look at our distribution of estimates from the unclustered simulations. 
ci95_clustered <- ggplot(
  sim_clustered,
  aes(
    x = as.numeric(id),
    y = b1,
    ymin = ci95_lower,
    ymax = ci95_upper,
    color = param_caught
  )
) +
  geom_errorbar(alpha = 0.5, color = "gray") +
  geom_point() +
  labs(x = 'Simulation ID', y = 'Estimate (95% CI)', title = 'Simulation Estimates with Regular SEs (50 Clusters; Rho = 0.7)') +
  scale_color_manual(guide=FALSE,
    name = 'Type I Error',
    labels = c('Yes', 'No'),
    values = c("red", "blue")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Economica"),
    axis.title = element_text(size = 16, margin = margin(
      t = 10,
      b = 0,
      r = 20,
      l = 0
    )),
    axis.text = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0.5,
             linetype = 'dashed',
             color = "chartreuse4") +
  NULL

ci95_unclustered / ci95_clustered

# We see that across all simulations, in ~45% of instances, the 95% CI did *not* contain the true parameter.
mean(sim_unclustered$param_caught)
mean(sim_clustered$param_caught)

### LET'S NOW USE CLUSTERED SE'S TO RESOLVE THIS.

# Here we run sims with the default rho of 0.5 across our default of 50 clusters.
sim_clustered_fix <- run_cluster_sim(rho = 0.7, SEs = "clustered")

# Now, How did we do when clustering was present? Fuckin terribly!
# Let's look at our distribution of estimates from the unclustered simulations. We can take 100 of them.
ci95_clustered_fix <- ggplot(
  sim_clustered_fix,
  aes(
    x = as.numeric(id),
    y = b1,
    ymin = ci95_lower,
    ymax = ci95_upper,
    color = param_caught
  )
) +
  geom_errorbar(alpha = 0.5, color = "gray") +
  geom_point() +
  labs(x = 'Simulation ID', y = 'Estimate (95% CI)', title = 'Simulation Estimates with Clustered SEs (50 Clusters; Rho = 0.7)') +
  scale_color_manual(guide=FALSE,
                     name = 'Type I Error',
                     labels = c('Yes', 'No'),
                     values = c("red", "blue")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Economica"),
    axis.title = element_text(size = 16, margin = margin(
      t = 10,
      b = 0,
      r = 20,
      l = 0
    )),
    axis.text = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0.5,
             linetype = 'dashed',
             color = "chartreuse4") +
  NULL
ci95_clustered_fix

# We see that we do much better; the 95% CI contains the correct parameter about 92-93% of the time.
mean(sim_clustered_fix$param_caught)

ci95_clustered / ci95_clustered_fix

### NOW, WHAT HAPPENS WHEN WE HAVE TOO FEW CLUSTERS?
# Here we run sims with the rho of 0.7 across a small number of (10) clusters.
sim_clustered_few <-
  run_cluster_sim(rho = 0.7,
                  n_cluster = 6,
                  SEs = "regular")
sim_clustered_few_fix <-
  run_cluster_sim(rho = 0.7,
                  n_cluster = 6,
                  SEs = "clustered")

# First, let's confirm that we still botch things if we ignore clustering here.
# Indeed, we do. In fact, we do worse.
ci95_clustered_few <- ggplot(
  sim_clustered_few,
  aes(
    x = as.numeric(id),
    y = b1,
    ymin = ci95_lower,
    ymax = ci95_upper,
    color = param_caught
  )
) +
  geom_errorbar(alpha = 0.5, color = "gray") +
  geom_point() +
  labs(x = 'Simulation ID', y = 'Estimate (95% CI)', title = 'Simulation Estimates with Regular SEs (Clusters = 6; Rho = 0.7)') +
  scale_color_manual(guide=FALSE,
                     name = 'Type I Error',
                     labels = c('Yes', 'No'),
                     values = c("red", "blue")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Economica"),
    axis.title = element_text(size = 16, margin = margin(
      t = 10,
      b = 0,
      r = 20,
      l = 0
    )),
    axis.text = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0.5,
             linetype = 'dashed',
             color = "chartreuse4") +
  NULL
ci95_clustered_few

# The CIs are way too tight; they omit the true estimate 80%+ of the time.
mean(sim_clustered_few$param_caught)

# Now let's see if clustering fixes it with a small number of clusters?
ci95_clustered_few_fix <- ggplot(
  sim_clustered_few_fix,
  aes(
    x = as.numeric(id),
    y = b1,
    ymin = ci95_lower,
    ymax = ci95_upper,
    color = param_caught
  )
) +
  geom_errorbar(alpha = 0.5, color = "gray") +
  geom_point() +
  labs(x = 'Simulation ID', y = 'Estimate (95% CI)', title = 'Simulation Estimates with Clustered SEs (Clusters = 6; Rho = 0.7)') +
  scale_color_manual(guide=FALSE,
                     name = 'Type I Error',
                     labels = c('Yes', 'No'),
                     values = c("red", "blue")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Economica"),
    axis.title = element_text(size = 16, margin = margin(
      t = 10,
      b = 0,
      r = 20,
      l = 0
    )),
    axis.text = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0.5,
             linetype = 'dashed',
             color = "chartreuse4") +
  NULL
ci95_clustered_few / ci95_clustered_few_fix

# We capture the true estimate just 85% of the time using 95% CIs. Not great...
mean(sim_clustered_few_fix$param_caught)

### SO WHAT TO DO? LET'S TRY WILD CLUSTER BOOTSTRAP SEs.

# Fast and wild cluster bootstrap errors will be a bit slow. Let it run.
sim_clustered_few_fwfix <-
  run_cluster_sim(rho = 0.7,
                  n_cluster = 6,
                  SEs = "fast_n_wild")

# Do fast and wild clustered SEs resolve this? 
ci95_clustered_few_fwfix <- ggplot(
  sim_clustered_few_fwfix,
  aes(
    x = as.numeric(id),
    y = b1,
    ymin = ci95_lower,
    ymax = ci95_upper,
    color = param_caught
  )
) +
  geom_errorbar(alpha = 0.5, color = "gray") +
  geom_point() +
  labs(x = 'Simulation ID', y = 'Estimate (95% CI)', title = 'Simulation Estimates with Wild Bootstrap Clustered SEs (Clusters = 6; Rho = 0.7)') +
  scale_color_manual(guide=FALSE,
                     name = 'Type I Error',
                     labels = c('Yes', 'No'),
                     values = c("red", "blue")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Economica"),
    axis.title = element_text(size = 16, margin = margin(
      t = 10,
      b = 0,
      r = 20,
      l = 0
    )),
    axis.text = element_text(size = 14)
  ) +
  geom_hline(yintercept = 0.5,
             linetype = 'dashed',
             color = "chartreuse4") +
  NULL
ci95_clustered_few_fwfix

# You see that it helps somewhat; we capture the true estimate about 93.7% of the time.
# Note that we shift to using 'webb' weights below 12 clusters, rather than rademacher weights.
mean(sim_clustered_few_fwfix$param_caught)

### LET'S NOW EXPLORE MORE SYSTEMATICALLY HOW THIS VARIES OVER CLUSTER COUNTS.

# We are going to look at how we do with regular standard errors, between 2 and 50 clusters. 
mean_CI_coverage_reg <- rep(0,49)
sd_CI_coverage_reg <- rep(0,49)
for (j in 2:50) {
  sim_clustered_fix_tmp <-
    run_cluster_sim(rho = 0.7,
                    n_cluster = j,
                    SEs = "regular")
  mean_CI_coverage_reg[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
  sd_CI_coverage_reg[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
}

regular_SE_by_clustN <- data.frame(coverage = mean_CI_coverage_reg, sd_CI_coverage_reg, clusters = seq(2:50))

# We are going to look at how we do with clustered standard errors, between 2 and 50 clusters. 
mean_CI_coverage <- rep(0,49)
sd_CI_coverage <- rep(0,49)
for (j in 2:50) {
  sim_clustered_fix_tmp <-
    run_cluster_sim(rho = 0.7,
                    n_cluster = j,
                    SEs = "clustered")
  mean_CI_coverage[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
  sd_CI_coverage[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
}

cluster_SE_by_clustN <- data.frame(coverage = mean_CI_coverage, sd_CI_coverage, clusters = seq(2:50))

# Now we are going to repeat this process, but employing fast and wild bootstrap clustered standard errors.
# Note that Fast'n'Wild can break down with errors if you have really few clusters (e.g., 2-4). 
# As such, I'm doing 5 through 50 clusters here. If you have fewer than 5 clusters, then you may be SoL.  
mean_CI_coverage_fw <- rep(0,49)
sd_CI_coverage_fw <- rep(0,49)
for (j in 2:50) {
  tryCatch(
    {
    sim_clustered_fix_tmp <-
      run_cluster_sim(rho = 0.7,
                      n_cluster = j,
                      SEs = "fast_n_wild")
    mean_CI_coverage_fw[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
    sd_CI_coverage_fw[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
  },
  error=function(cond) {
    message(paste("\nFast'n'Wild Failed.\n"))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  })
      
}

cluster_SE_by_clustN_fw <- data.frame(coverage = mean_CI_coverage_fw, sd_CI_coverage_fw, clusters = seq(2:50))

perfPalette_fw <- colorRampPalette(brewer.pal(9, "Greens"))(length(unique(cluster_SE_by_clustN_fw$coverage))) 
unique_coverages_fw <- cluster_SE_by_clustN_fw %>% group_by(coverage) %>% summarize(coverage = first(coverage)) %>% arrange(desc(coverage))
unique_coverages_fw$heat <- perfPalette_fw

cluster_SE_by_clustN_fw <- cluster_SE_by_clustN_fw %>% merge(unique_coverages_fw, by="coverage")

# Now we can plot the SE performance side-by-side. 
ggplot(data = cluster_SE_by_clustN_fw, aes(y=coverage*100,x=(clusters+1),alpha=0.7)) + 
  geom_hline(yintercept=95,color="red",size=1,linetype="dashed")+
  #geom_vline(xintercept=30,color="black",size=1)+
  geom_line(aes(color="green"),size=2) + 
  geom_line(data=cluster_SE_by_clustN, aes(color="blue"),size=2)+
  geom_line(data=regular_SE_by_clustN, aes(color="purple"),size=2)+
  xlab("# of Clusters (1000 Sims Each)") +
  ylab("Avg 95% CI Coverage of True Beta") +
  theme_classic() +
  theme(
    text = element_text(size=18,family = "Economica"),
    axis.title = element_text(size = 16, margin = margin(
      t = 10,
      b = 0,
      r = 20,
      l = 0
    )),
    axis.text = element_text(size = 14)
  ) +
  ylim(-5,97) +
  scale_x_continuous(breaks=c(2,10,20,30,40,50))+
  scale_y_continuous(breaks=seq(5,95,by=10))+
  scale_alpha_continuous(guide=FALSE)+
  ggtitle("Downward Bias in Standard Errors Under (Few) Clusters") +
  scale_color_manual(name="SE Type",values=c("blue","green","orange"),labels=c("Clustered","Wild Clustered","Regular"))
NULL

### PLOT THEM ATOP ONE ANOTHER TO SEE RELATIVE PERFORMANCE...

