library(extraDistr)
library(ggthemes)
library(ggplot2)
library(knitr)
library(plyr)
library(EnvStats)
library(MASS)
library(copula)
library(reshape2)
library(viridis)  
library(randomForest)
library(dplyr)
library(combinat)
library(reshape2)
library(tidyr)
##################################################################
# Function to calculate test statistic

m_whitney_ts = function(sample1, sample2){
  count = 0
  for (i in 1:length(sample1)){
    for (j in 1:length(sample2)){
      if (sample1[i] > sample2[j]){
        count = count + 1
      }
    }
  }
  return(count)
}

# Function to draw sample

draw_sample = function(dist,size,par1,par2,...){
  if (dist == 'normal') {
    mean=par1;sd=par2;sample <- rnorm(size, mean, sd)
  } 
  else if (dist == 'uniform') {
    min_val= par1;max_val=par2;sample=runif(size, min_val, max_val)
  } 
  else if (dist == 'poisson') {
    lambda <- par1;sample <- rpois(size, lambda)
  } 
  else if (dist == 'cauchy') {
    location=par1;scale=par2;sample <- rcauchy(size, location,scale)
  } 
  else if (dist == 'exponential') {
    rate=par1;sample <- rexp(size,rate=1)
  }
  else if (dist == 'weibull') {
    rate=par1;scale= par2;sample <- rweibull(size,rate,scale)
  }
  else if (dist == 'shifed exponential') {
    rate=par1;sample <- rexp(size,rate=1)+par2
  }
  else if (dist == 'gamma') {
    shape=par1;rate=par2;sample <- rcauchy(size,shape,rate)
  } 
  else if (dist == 'beta') {
    shape1=par1;shape2=par2;sample <- rbeta(size,shape1,shape2)
  } 
  else if (dist == 'laplace') {
    mu=par1;sigma=par2;sample <- rlaplace(size,mu,sigma)
  } 
  else if (dist == 'binomial') {
    trials=par1;prob=par2;sample <- rbinom(size,trials,prob)
  } 
  else if (dist == 'negative binomial') {
    s=par1;prob=par2;sample <- rbinom(size,s,prob)
  }
  else if (dist == 'geometric') {
    prob=par1;sample <- rgeom(size,prob)
  } 
  else if (dist == 'logistic') {
    location=par1;scale=par2;sample <- rlogis(size,location,scale)
  } 
  else if (dist == 'lognormal') {
    location=par1;scale=par2;sample <- rlnorm(size,location,scale)
  }
  else if (dist == 'beta I') {
    location=par1;scale=par2;sample <- rbeta(size,location,scale)
  }
  else if (dist == 't') {
    sample <- rt(size,par1)
  } 
  else if (dist == 'chisq') {
    sample <- rchisq(size,par1)
  } 
  else if (dist == 'f') {
    location=par1;scale=par2;sample <- rf(size,location,scale)
  }
  else if (dist == "beta II"){
    sample = numeric(size)
    location = par1;scale = par2;x<- rbeta(size,location,scale)
    for (i in 0:size){
      sample[i] = x[i]/(1-x[i])
    }
  }
  else {
    stop("Unsupported distribution")
  }
  
  return(sample)
}

# Function to test distribution free of test statistic

TS_H = function(iter = 1000,dist,para1=0,para2=1,size1,size2,theta){
  ts_val = numeric(iter)
  for (i in 1:iter){
    x = draw_sample(dist,size1,para1,para2)
    y = draw_sample(dist,size2,para1+theta,para2)
    ts_val[i] = m_whitney_ts(x,y)
  }
  return(ts_val)
}

# To plot Mann Whitney Test statistic with the distributions

plot_H0 <- function(m, n) {
  dist_settings <- list(
    beta        = list(p1 = 1/2, p2 = 1/2, col = "cyan"),
    cauchy      = list(p1 = 1/2, p2 = 1, col = "orange"),
    normal = list(p1 = 0, p2 = 1, col = "purple"),
    chisq = list(p1 = 5, col = "yellow")
  )
  
  # 2. Generate the data using a loop (lapply) instead of manual typing
  plot_data <- do.call(rbind, lapply(names(dist_settings), function(d) {
    set <- dist_settings[[d]]
    data.frame(
      Distribution = d,
      values = TS_H(dist = d, size1 = m, size2 = n, para1 = set$p1, para2 = set$p2, theta = 0)
    )
  }))

  # 3. Create a named vector for the fill colors
  fill_colors <- sapply(dist_settings, `[[`, "col")

  # 4. The Plot
  ggplot(plot_data, aes(x = values, fill = Distribution)) +
    geom_density(alpha = 0.4, color = "black", size = 0.9) + # Set border color here once
    geom_vline(xintercept = m * n / 2, linetype = "dashed", color = "black") +
    scale_fill_manual(values = fill_colors) +
    labs(
      title = "Distribution under H0",
      x = "Mann-Whitney Statistic",
      y = "Density"
    ) +
    theme_stata()
}

# Plotted figures under H_0
plot_H0(20,25)
plot_H0(10,50)
plot_H0(65,55)
plot_H0(200,150)
plot_H0(400,450)

plot_H1 <- function(m, n) {
  # 1. Setup our variables
  thetas <- c(0.1, 0.5, 1, -1)
  dists <- c("normal", "cauchy", "logistic", "laplace")
  
  # Map distributions to their specific colors
  fill_colors <- c("normal" = "cyan", "cauchy" = "red", "logistic" = "blue", "laplace" = "yellow")
  
  # 2. Generate data for all combinations of Distribution and Theta
  # We use expand.grid and lapply to avoid manual typing
  plot_data <- do.call(rbind, lapply(thetas, function(th) {
    do.call(rbind, lapply(dists, function(d) {
      data.frame(
        Distribution = d,
        Theta = th,
        # Standard parameters (0,1) for these three distributions
        values = TS_H(dist = d, size1 = m, size2 = n, para1 = 0, para2 = 1, theta = th)
      )
    }))
  }))

  # 3. The Plot
  ggplot(plot_data, aes(x = values, fill = Distribution)) +
    geom_density(alpha = 0.4, color = "black", size = 0.7) +
    geom_vline(xintercept = m * n / 2, linetype = "dashed", color = "black") +
    # This creates a separate box for each Theta value
    facet_wrap(~Theta, labeller = label_both, scales = "free_y") + 
    scale_fill_manual(values = fill_colors) +
    labs(
      title = "Distribution Comparison across Different Theta Values",
      subtitle = paste("Sample sizes: m =", m, ", n =", n),
      x = "Mann-Whitney Statistic",
      y = "Density"
    ) +
    theme_stata() +
    theme(legend.position = "bottom")
}

# Plotted figures under H_1
plot_H1(5,7)
plot_H1(19,15)
plot_H1(45,57)
plot_H1(50,70)

##################################################################
#Exact null dist
# Function to compute the exact distribution of the Mann-Whitney statistic
exact_mann_whitney_distribution <- function(m, n) {
  # Generate all possible ways to choose m elements from (m + n) ranks
  total_ranks <- 1:(m + n)
  all_combinations <- combn(total_ranks, m, simplify = FALSE)
  
  # Compute the Mann-Whitney statistic for each combination
  mw_stats <- sapply(all_combinations, function(rankings) {
    sum(rankings) - m * (m + 1) / 2  # Compute U statistic
  })
  
  # Compute frequency table
  freq_table <- table(mw_stats) / length(mw_stats)
  
  # Convert to a data frame
  df <- data.frame(U_stat = as.numeric(names(freq_table)), Probability = as.numeric(freq_table))
  
  return(df)
}

# Example: m = 3, n = 4
m <- 7
n <- 5
exact_distribution <- exact_mann_whitney_distribution(m, n)
print(exact_distribution)

##################################################################
# Asymptotic Distributions
assymp_normal_qq <- function(distn, m, n) {
  # 1. Calculate the sample proportion
  total_n <- m + n
  lambda_ <- m / total_n
  
  # 2. Generate and Standardize the Statistic
  # (Standardization formula based on asymptotic variance theory)
  # under H_0
  raw_stat <- TS_H(dist = distn, size1 = m, size2 = n, para1 = 0, para2 = 1, theta = 0)
  centered_stat <- (raw_stat / (m * n)) - 0.5
  
  # Scaling factor (Standard Error)
  std_error <- sqrt((3 - 2 * lambda_) / (12 * lambda_ * (1 - lambda_)))
  z_scores <- (centered_stat * sqrt(total_n)) / std_error
  
  # 3. Create Plot
  ggplot(data.frame(z = z_scores), aes(sample = z)) +
    geom_qq(color = "blue", alpha = 0.6) +
    geom_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(
      title = paste("Asymptotic Normality Check:", distn),
      subtitle = paste("Sample sizes: m =", m, "n =", n),
      x = "Theoretical Normal Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_stata()
}
# Plotting QQplot
assymp_normal_qq("cauchy",20,40)
assymp_normal_qq("normal",50,70)
assymp_normal_qq("laplace",30,40)
assymp_normal_qq("logistic",60,70)
assymp_normal_qq("exponential",80,110)
##################################################################
qq_mwu_TS_H <- function(dists = c("cauchy", "normal", "exponential")) {
  set.seed(123)
  
  get_U <- function(x, y) {
    r <- rank(c(x, y), ties.method = "average")
    Rx <- sum(r[1:length(x)])
    Rx - length(x) * (length(x) + 1) / 2
  }
  
  m <- 3
  i_vals <- 0:9999
  
  out <- do.call(rbind, lapply(dists, function(dist_name) {
    z <- sapply(i_vals, function(i) {
      n <- 7 + 10 * i
      
      x <- draw_sample(dist_name, m, par1 = 0, par2 = 1)
      y <- draw_sample(dist_name, n, par1 = 0, par2 = 1)
      
      mu_U <- m * n / 2
      sd_U <- sqrt(m * n * (m + n + 1) / 12)
      
      (get_U(x, y) - mu_U) / sd_U
    })
    
    data.frame(
      dist = dist_name,
      iteration = i_vals,
      n = 7 + 10 * i_vals,
      z = z
    )
  }))
  
  ggplot(out, aes(sample = z)) +
    stat_qq(size = 1, alpha = 0.5, color = "steelblue") +
    stat_qq_line(color = "red", linewidth = 0.8) +
    facet_wrap(~ dist, scales = "free") +
    labs(
      title = "QQ Plots of Standardized Mann-Whitney U Statistic",
      subtitle = "m = 3, n = 7 + 10*i for i = 0,1,...,9999",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal(base_size = 12)
}

qq_mwu_TS_H <- function(dists = c("cauchy", "normal", "exponential")) {
  set.seed(123)
  
  get_U <- function(x, y) {
    r <- rank(c(x, y), ties.method = "average")
    Rx <- sum(r[1:length(x)])
    Rx - length(x) * (length(x) + 1) / 2
  }
  
  i_vals <- 1:500
  
  out <- do.call(rbind, lapply(dists, function(dist_name) {
    z <- sapply(i_vals, function(i) {
      m <- 10 * ceiling(sqrt(i))
      n <- 10 * i^2
      
      x <- draw_sample(dist_name, m, par1 = 0, par2 = 1)
      y <- draw_sample(dist_name, n, par1 = 0, par2 = 1)
      
      mu_U <- m * n / 2
      sd_U <- sqrt(m * n * (m + n + 1) / 12)
      
      (get_U(x, y) - mu_U) / sd_U
    })
    
    data.frame(
      dist = dist_name,
      iteration = i_vals,
      m = 10 * ceiling(sqrt(i_vals)),
      n = 10 * i_vals,
      z = z
    )
  }))
  
  ggplot(out, aes(sample = z)) +
    stat_qq(size = 1, alpha = 0.5, color = "steelblue") +
    stat_qq_line(color = "red", linewidth = 0.8) +
    facet_wrap(~ dist, scales = "free") +
    labs(
      title = "QQ Plots of Standardized Mann-Whitney U Statistic",
      subtitle = "m = 10*ceiling(sqrt(i)), n = 10*i^2 for i = 1,...,500",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal(base_size = 12)
}

qq_mwu_TS_H <- function(dists = c("cauchy", "normal", "exponential")) {
  set.seed(123)
  
  get_U <- function(x, y) {
    r <- rank(c(x, y), ties.method = "average")
    Rx <- sum(r[1:length(x)])
    Rx - length(x) * (length(x) + 1) / 2
  }
  
  p <- 0.01
  i_vals <- 1:1000
  
  out <- do.call(rbind, lapply(dists, function(dist_name) {
    z <- sapply(i_vals, function(i) {
      x_mod <- i %% 2
      m <- 100 * i * (p^x_mod) * ((1 - p)^(1 - x_mod))
      n <- 100 * i - m
      
      x <- draw_sample(dist_name, m, par1 = 0, par2 = 1)
      y <- draw_sample(dist_name, n, par1 = 0, par2 = 1)
      
      mu_U <- m * n / 2
      sd_U <- sqrt(m * n * (m + n + 1) / 12)
      
      (get_U(x, y) - mu_U) / sd_U
    })
    
    x_mod_vals <- i_vals %% 2
    m_vals <- 10 * i_vals * (p^x_mod_vals) * ((1 - p)^(1 - x_mod_vals))
    n_vals <- 10 * i_vals - m_vals
    
    data.frame(
      dist = dist_name,
      iteration = i_vals,
      m = m_vals,
      n = n_vals,
      z = z
    )
  }))
  
  ggplot(out, aes(sample = z)) +
    stat_qq(size = 1, alpha = 0.5, color = "steelblue") +
    stat_qq_line(color = "red", linewidth = 0.8) +
    facet_wrap(~ dist, scales = "free") +
    labs(
      title = "QQ Plots of Standardized Mann-Whitney U Statistic",
      subtitle = "i odd:m = 10*i*p;i even:m=10*i*(1-p), n = 10*i - m; i = 1,...,1000",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal(base_size = 12)
}
# Example
qq_mwu_TS_H("cauchy")
qq_mwu_TS_H("exponential")
qq_mwu_TS_H("normal")

# histograms under H_0
plot_asymptotic_histograms_H0 <- function(m, n, iter = 1000) {
  # 1. Define the distributions to test
  dists <- c("laplace", "normal", "cauchy", "logistic")
  
  # 2. Generate and Standardize Data
  # We use the Z-score formula: (U/mn - 0.5) / SE
  total_n <- m + n
  lambda_ <- m / total_n
  std_error <- sqrt((3 - 2 * lambda_) / (12 * lambda_ * (1 - lambda_) * total_n))
  
  plot_data <- do.call(rbind, lapply(dists, function(d) {
    # Generate raw statistics
    raw_stats <- TS_H(dist = d, size1 = m, size2 = n, para1 = 0, para2 = 1, theta = 0)
    
    # Standardize to Z-scores
    z_scores <- (raw_stats / (m * n) - 0.5) / std_error
    
    data.frame(values = z_scores, Distribution = d)
  }))
  
  # 3. Create the Grid Plot
  ggplot(plot_data, aes(x = values)) +
    geom_histogram(bins = 15, fill = "aquamarine", color = "blue", alpha = 0.7) +
    facet_wrap(~Distribution, scales = "free", labeller = label_both) +
    labs(
      title = "Asymptotic Histograms of Standardized Mann-Whitney Statistic",
      subtitle = paste("Large Sample Size: m =", m, ", n =", n),
      x = "Standardized Z-values",
      y = "Frequency"
    ) +
    theme_stata() + # Matches the styling in your image
    theme(strip.text = element_text(size = 12, face = "bold"))
}

# histograms under H_1
plot_asymptotic_histograms_H1 <- function(m, n, iter = 1000) {
  # 1. Define the distributions to test
  dists <- c("laplace", "normal", "cauchy", "logistic")
  
  # 2. Generate and Standardize Data
  # We use the Z-score formula: (U/mn - 0.5) / SE
  total_n <- m + n
  lambda_ <- m / total_n
  std_error <- sqrt((3 - 2 * lambda_) / (12 * lambda_ * (1 - lambda_) * total_n))
  
  plot_data <- do.call(rbind, lapply(dists, function(d) {
    # Generate raw statistics
    raw_stats <- TS_H(dist = d, size1 = m, size2 = n, para1 = 0, para2 = 1, theta = 5)
    
    # Standardize to Z-scores
    z_scores <- (raw_stats / (m * n) - 0.5) / std_error
    
    data.frame(values = z_scores, Distribution = d)
  }))
  
  # 3. Create the Grid Plot
  ggplot(plot_data, aes(x = values)) +
    geom_histogram(bins = 15, fill = "aquamarine", color = "blue", alpha = 0.7) +
    facet_wrap(~Distribution, scales = "free", labeller = label_both) +
    labs(
      title = "Asymptotic Histograms of Standardized Mann-Whitney Statistic",
      subtitle = paste("Large Sample Size: m =", m, ", n =", n),
      x = "Standardized Z-values",
      y = "Frequency"
    ) +
    theme_stata() + # Matches the styling in your image
    theme(strip.text = element_text(size = 12, face = "bold"))
}


# --- Execute for large sample size ---
plot_asymptotic_histograms_H0(m = 100, n = 100)
plot_asymptotic_histograms_H0(m = 1000, n = 1000)
plot_asymptotic_histograms_H1(m = 10, n = 20)
plot_asymptotic_histograms_H1(m = 50, n = 30)
plot_asymptotic_histograms_H1(m = 100, n = 150)

# Function to simulate and standardize Mann-Whitney U statistics
################################################################
simulate_mannwhitney_asymptotic <- function(dist, par1, par2, m, n, iter = 5000) {
  mw_stats <- numeric(iter)
  for (i in 1:iter) {
    x <- draw_sample(dist, m, par1, par2)
    y <- draw_sample(dist, n, par1, par2)  # Under H0, same distribution
    
    u_stat <- m_whitney_ts(x, y)  # Custom U-statistic
    expected_u <- m * n / 2
    sd_u <- sqrt(m * n * (m + n + 1) / 12)
    
    # Standardized U-statistic
    mw_stats[i] <- (u_stat - expected_u) / sd_u
  }
  return(mw_stats)
}

######################################################
# Visualizing asymptotic normality for increasing m, n
######################################################
# Prepare data frame for plotting
plot_data <- do.call(rbind, lapply(names(results_list), function(name) {
  data.frame(Standardized_U = results_list[[name]], Sample_Size = name)
}))

# Plot density of standardized Mann-Whitney U-statistic
p <- ggplot(plot_data, aes(x = Standardized_U, fill = Sample_Size, color = Sample_Size)) +
  geom_density(alpha = 0.3) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), 
                color = "black", linetype = "dashed", size = 1.2, 
                aes(linetype = "Standard Normal")) +
  labs(
    title = "Asymptotic Normality of Mann-Whitney U-statistic",
    subtitle = "Comparison of standardized U-statistic to N(0,1) for increasing sample sizes",
    x = "Standardized Mann-Whitney U-statistic",
    y = "Density",
    linetype = NULL
  ) +
  scale_linetype_manual(values = c("Standard Normal" = "dashed")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(p)

#######################################################
# Asymptotic Normality Check for Cauchy and Exponential
#######################################################

# Function to compute standardized Mann-Whitney statistics for given distribution
asymptotic_mw_plot <- function(dist_name, par1, par2 = NULL, m, n, iter = 1000) {
  
  mw_z_scores <- numeric(iter)
  
  for (i in 1:iter) {
    # Generate samples
    sample1 <- draw_sample(dist = dist_name, size = m, par1 = par1, par2 = par2)
    sample2 <- draw_sample(dist = dist_name, size = n, par1 = par1, par2 = par2)
    
    # Mann-Whitney U-statistic
    mw_stat <- m_whitney_ts(sample1, sample2)
    
    # Standardization under H0
    expected_mw <- m * n / 2
    sd_mw <- sqrt(m * n * (m + n + 1) / 12)
    
    # Z-score
    z_mw <- (mw_stat - expected_mw) / sd_mw
    mw_z_scores[i] <- z_mw
  }
  
  # Return dataframe
  data.frame(z_mw = mw_z_scores, Distribution = dist_name, m = m, n = n)
}

####################################
# Run for Cauchy and Exponential
####################################

# Sample sizes
m <- 100
n <- 150
iter <- 3000  # number of simulations

# Cauchy with location = 0, scale = 1
cauchy_res <- asymptotic_mw_plot("cauchy", par1 = 0, par2 = 1, m = m, n = n, iter = iter)

# Exponential with rate = 1
exp_res <- asymptotic_mw_plot("exponential", par1 = 1, m = m, n = n, iter = iter)

# Combine results
combined_res <- rbind(cauchy_res, exp_res)

############################################
# Plotting Asymptotic Normality
############################################

# Density Plot
p1 <- ggplot(cauchy_res, aes(x = z_mw)) +
  geom_density(fill = "brown", alpha = 0.4, color = "black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), 
                color = "blue", linetype = "dashed", size = 1.2) +
  labs(
    title = paste0("Asymptotic Normality of Mann-Whitney U (Cauchy(0,1))"),
    subtitle = paste0("Sample Sizes: m = ", m, ", n = ", n, "; Iter = ", iter),
    x = "Standardized Mann-Whitney U (Z)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

p2 <- ggplot(exp_res, aes(x = z_mw)) +
  geom_density(fill = "cyan", alpha = 0.4, color = "black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), 
                color = "blue", linetype = "dashed", size = 1.2) +
  labs(
    title = paste0("Asymptotic Normality of Mann-Whitney U (Exponential(1))"),
    subtitle = paste0("Sample Sizes: m = ", m, ", n = ", n, "; Iter = ", iter),
    x = "Standardized Mann-Whitney U (Z)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

# Combine plots
grid.arrange(p1, p2, ncol = 2)
##################################################################
# Checking the level maintainance
set.seed(123)

# Sample sizes
m <- 11
n <- 13

# Number of simulations
Nsim <- 10000

# Store p-values
pvals <- numeric(Nsim)

for (i in 1:Nsim) {
  
  # Under H0: same distribution -> Normal
  x <- rnorm(m)
  y <- rnorm(n)
  
  # Mann-Whitney test (one-sided: theta > 0)
  test <- wilcox.test(x, y, alternative = "greater", exact = TRUE)
  
  pvals[i] <- test$p.value
}

# Estimated size at alpha = 0.05
size_05 <- mean(pvals < 0.05)

# Estimated size at alpha = 0.01
size_01 <- mean(pvals < 0.01)

# Print results
cat("Estimated size at alpha = 0.05:", size_05, "\n")
cat("Estimated size at alpha = 0.01:", size_01, "\n")
############################################################
# Cauchy
for (i in 1:Nsim) {
  # Under H0: same distribution
  x <- rcauchy(m)
  y <- rcauchy(n)
  # Mann-Whitney test (one-sided: theta > 0)
  test <- wilcox.test(x, y, alternative = "greater", exact = TRUE)
  
  pvals[i] <- test$p.value
}

# Estimated size at alpha = 0.05
size_05 <- mean(pvals < 0.05)

# Estimated size at alpha = 0.01
size_01 <- mean(pvals < 0.01)

# Print results
cat("Estimated size at alpha = 0.05:", size_05, "\n")
cat("Estimated size at alpha = 0.01:", size_01, "\n")
#########################################################
# Poisson
for (i in 1:Nsim) {
  # Under H0: same distribution
  x <- rpois(m, lambda = 10)
  y <- rpois(n, lambda = 10)
  # Mann-Whitney test (one-sided: theta > 0)
  test <- wilcox.test(x, y, alternative = "greater", exact = TRUE)
  
  pvals[i] <- test$p.value
}

# Estimated size at alpha = 0.05
size_05 <- mean(pvals < 0.05)

# Estimated size at alpha = 0.01
size_01 <- mean(pvals < 0.01)

# Print results
cat("Estimated size at alpha = 0.05:", size_05, "\n")
cat("Estimated size at alpha = 0.01:", size_01, "\n")
##################################################################
# HLE and CI
set.seed(123)

# ------------------------------------------------------------
# One replication: Mann-Whitney/Wilcoxon rank-sum inference
# ------------------------------------------------------------
one_run <- function(n1, n2, dist, theta, conf.level = 0.95, exact = FALSE) {
  
  # Generate baseline samples
  gen_sample <- function(n, dist) {
    switch(dist,
           normal      = rnorm(n, mean = 0, sd = 1),
           poisson     = rpois(n, lambda = 3),
           cauchy      = rcauchy(n, location = 0, scale = 1),
           exponential = rexp(n, rate = 1),
           stop("Unknown distribution")
    )
  }
  
  x <- gen_sample(n1, dist)
  z <- gen_sample(n2, dist)
  y <- z + theta   # location shift model
  
  # Wilcoxon rank-sum / Mann-Whitney test with CI and HL estimator
  wt <- suppressWarnings(
    wilcox.test(x, y,
                alternative = "two.sided",
                conf.int = TRUE,
                conf.level = conf.level,
                exact = exact)
  )
  
  # In R, estimate = location shift for x - y.
  # Since our theta is defined through y = x + theta,
  # true theta corresponds to (y - x), so we flip the sign.
  hl_theta <- -as.numeric(wt$estimate)
  ci_theta <- -rev(as.numeric(wt$conf.int))
  
  reject <- (wt$p.value < 0.05)
  cover  <- (ci_theta[1] <= theta && theta <= ci_theta[2])
  ci_len <- diff(ci_theta)
  
  c(
    hl        = hl_theta,
    ci_low    = ci_theta[1],
    ci_high   = ci_theta[2],
    p_value   = wt$p.value,
    reject    = reject,
    cover     = cover,
    ci_length = ci_len
  )
}

# ------------------------------------------------------------
# Simulation over many replications
# ------------------------------------------------------------
simulate_mw <- function(n1 = 20, n2 = 20,
                        dist = "normal",
                        theta = 0,
                        nsim = 2000,
                        conf.level = 0.95,
                        exact = FALSE) {
  
  out <- replicate(
    nsim,
    one_run(n1, n2, dist, theta, conf.level, exact)
  )
  
  out <- t(out)
  out <- as.data.frame(out)
  
  summary_row <- data.frame(
    distribution   = dist,
    theta_true     = theta,
    n1             = n1,
    n2             = n2,
    nsim           = nsim,
    mean_HL        = mean(out$hl, na.rm = TRUE),
    bias_HL        = mean(out$hl - theta, na.rm = TRUE),
    rmse_HL        = sqrt(mean((out$hl - theta)^2, na.rm = TRUE)),
    rejection_rate = mean(out$reject, na.rm = TRUE),
    coverage       = mean(out$cover, na.rm = TRUE),
    avg_CI_length  = mean(out$ci_length, na.rm = TRUE)
  )
  
  list(summary = summary_row, raw = out)
}

# ------------------------------------------------------------
# Run all cases
# ------------------------------------------------------------
dists  <- c("normal", "poisson", "cauchy", "exponential")
thetas <- c(0, -1, 1)

results_list <- list()
summary_table <- data.frame()

idx <- 1
for (dd in dists) {
  for (th in thetas) {
    cat("Running:", dd, "theta =", th, "\n")
    sim <- simulate_mw(
      n1 = 20, n2 = 20,
      dist = dd,
      theta = th,
      nsim = 2000,
      conf.level = 0.95,
      exact = FALSE
    )
    results_list[[idx]] <- sim
    summary_table <- rbind(summary_table, sim$summary)
    idx <- idx + 1
  }
}

print(summary_table)

##################################################################
# Discrete Distribution
plot_asum1_disc <- function(m, n) {
  dist_list <- list(
    binomial            = list(p1 = 7,    p2 = 0.35, col = "red"),
    poisson             = list(p1 = 4.9,  p2 = 1,    col = "cyan"),
    'negative binomial'   = list(p1 = 8,    p2 = 0.79, col = "brown"),
    geometric           = list(p1 = 0.37, p2 = 1,    col = "yellow")
  )
  plot_data <- do.call(rbind, lapply(names(dist_list), function(d) {
    set <- dist_list[[d]]
    data.frame(
      Distribution = d,
      values = TS_H(dist = d, size1 = m, size2 = n, para1 = set$p1, para2 = set$p2, theta = 0)
    )
}))
  
  # 3. Extract colors
  fill_colors <- sapply(dist_list, `[[`, "col")
  
  # 4. Create the plot
  ggplot(plot_data, aes(x = values, fill = Distribution)) +
    geom_density(alpha = 0.4, color = "black", size = 0.7) +
    geom_vline(xintercept = (m * n) / 2, linetype = "dashed", color = "black", lwd = 0.8) +
    scale_fill_manual(values = fill_colors) +
    labs(
      title = "Distribution for Discrete Data",
      subtitle = paste("Reference Line at Mean (m*n/2) =", (m * n) / 2),
      x = "Mann-Whitney Statistic",
      y = "Density"
    ) +
    theme_stata()
}

# Run the function
plot_asum1_disc(19, 13)
###################################################
# Function to compute Mann-Whitney test statistic
m_whitney_ts <- function(sample1, sample2) {
  test <- wilcox.test(sample1, sample2)$statistic
  return(as.numeric(test))
}

# Generate 4 bivariate samples (single run)
generate_samples <- function(n) {
  mu <- c(0, 0)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  
  # 1. Bivariate Normal
  normal_samples <- mvrnorm(n, mu, Sigma)
  
  # 2. Bivariate Exponential using Copula
  cop <- normalCopula(param = 0.3, dim = 2)
  exp_samples <- rCopula(n, cop)
  exp_samples <- cbind(qexp(exp_samples[,1]), qexp(exp_samples[,2]))
  
  # 3. Bivariate Poisson using Cholesky
  Z <- mvrnorm(n, mu = c(0, 0), Sigma)
  pois_samples <- cbind(
    qpois(pnorm(Z[,1]), lambda = 7),
    qpois(pnorm(Z[,2]), lambda = 11)
  )
  
  # 4. Bivariate Cauchy
  cauchy_samples <- cbind(
    rcauchy(n, location = 0, scale = 1),
    rcauchy(n, location = 1, scale = 2.5)
  )
  
  return(list(normal_samples, exp_samples, pois_samples, cauchy_samples))
}

# Single-run Mann-Whitney statistics
n <- 75
samples <- generate_samples(n)

mann_whitney_results <- data.frame(
  values = c(
    m_whitney_ts(samples[[1]][,1], samples[[1]][,2]),
    m_whitney_ts(samples[[2]][,1], samples[[2]][,2]),
    m_whitney_ts(samples[[3]][,1], samples[[3]][,2]),
    m_whitney_ts(samples[[4]][,1], samples[[4]][,2])
  ),
  Distribution = rep(
    c("Bivariate Normal", "Bivariate Exponential",
      "Bivariate Poisson", "Bivariate Cauchy"),
    each = 1
  )
)

# Replicate for density plot
mann_whitney_results <- mann_whitney_results[rep(1:4, each = 1000), ]

ggplot(mann_whitney_results,
       aes(x = values, col = Distribution, fill = Distribution)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = n * n / 2,
             linetype = "dashed", col = "black", lwd = 0.9) +
  labs(
    title = "Mann-Whitney Test Statistic Density for Different Bivariate Distributions",
    x = "Mann-Whitney Statistic",
    y = "Density"
  ) +
  scale_fill_manual(values = c(
    "Bivariate Normal"      = "cyan",
    "Bivariate Exponential" = "red",
    "Bivariate Poisson"     = "blue",
    "Bivariate Cauchy"      = "brown"
  )) +
  theme_minimal()

############################################################
# Power Curves
# Consistency
# right tailed test
MW_right_tail<-function(distn,size1,size2,alpha=0.025,iter=1000,rge=4,...){
  effect_sizes <- seq(0,rge, by = 0.1);power_values=numeric(length(effect_sizes))
  for(i in 1:length(effect_sizes)){
    reject<-0
    for(j in 1:iter){
      sample1=draw_sample(dist=distn,par1=0,par2=1,size=size1);sample2=draw_sample(dist=distn,par1=(0+effect_sizes[i]),par2=1,size=size1)
      if(wilcox.test(sample1,sample2)$p.value <alpha) reject=reject+1
    }
    power_values[i]=reject/iter
  }
  return(list(a=effect_sizes,b=power_values))
}
## normal
power_normal=data.frame(x1= rep(seq(0,4, by = 0.1),5),y1=c(MW_right_tail(distn="normal",size1 =11,size2 = 12)$b,
                                                           MW_right_tail(distn="normal",size1 = 15,size2 = 17)$b,
                                                           MW_right_tail(distn="normal",size1 = 29,size2 = 42)$b,
                                                           MW_right_tail(distn="normal",size1 = 47,size2 = 55)$b,
                                                           MW_right_tail(distn="normal",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=power_normal,aes(x=x1,y=y1,col=size))+geom_line(lwd=0.9)+
  scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", N ( theta,1)," for ",theta >0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## cauchy
power_cauchy=data.frame(x1= rep(seq(0,4, by = 0.1),5),y1=c(MW_right_tail(distn="cauchy",size1 =11,size2 = 12)$b,
                                                           MW_right_tail(distn="cauchy",size1 = 15,size2 = 17)$b,
                                                           MW_right_tail(distn="cauchy",size1 = 29,size2 = 42)$b,
                                                           MW_right_tail(distn="cauchy",size1 = 47,size2 = 55)$b,
                                                           MW_right_tail(distn="cauchy",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=power_cauchy,aes(x=x1,y=y1,col=size))+geom_line(lwd=0.9)+
  scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", C ( theta,1)," for ",theta >0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## laplace
power_laplace=data.frame(x1= rep(seq(0,4, by = 0.1),5),y1=c(MW_right_tail(distn="laplace",size1 =11,size2 = 12)$b,
                                                           MW_right_tail(distn="laplace",size1 = 15,size2 = 17)$b,
                                                           MW_right_tail(distn="laplace",size1 = 29,size2 = 42)$b,
                                                           MW_right_tail(distn="laplace",size1 = 47,size2 = 55)$b,
                                                           MW_right_tail(distn="laplace",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=power_laplace,aes(x=x1,y=y1,col=size))+geom_line(lwd=0.9)+
  scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", Laplace ( theta,1)," for ",theta >0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## logistic
power_logistic=data.frame(x1= rep(seq(0,4, by = 0.1),5),y1=c(MW_right_tail(distn="logistic",size1 =11,size2 = 12)$b,
                                                            MW_right_tail(distn="logistic",size1 = 15,size2 = 17)$b,
                                                            MW_right_tail(distn="logistic",size1 = 29,size2 = 42)$b,
                                                            MW_right_tail(distn="logistic",size1 = 47,size2 = 55)$b,
                                                            MW_right_tail(distn="logistic",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=power_logistic,aes(x=x1,y=y1,col=size))+geom_line(lwd=0.9)+
  scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", Logistic ( theta,1)," for ",theta >0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
######################################################
# left tailed tests
MW_left_tail<-function(distn,size1,size2,alpha=0.025,iter=1000,rge=-4,...){
  effect_sizes <- seq(rge,0, by = 0.1);power_values=numeric(length(effect_sizes))
  for(i in 1:length(effect_sizes)){
    reject<-0
    for(j in 1:iter){
      sample1=draw_sample(dist=distn,par1=0,par2=1,size=size1);sample2=draw_sample(dist=distn,par1=(0+effect_sizes[i]),par2=1,size=size1)
      if(wilcox.test(sample1,sample2)$p.value <alpha) reject=reject+1
    }
    power_values[i]=reject/iter
  }
  return(list(a=effect_sizes,b=power_values))
}
## normal
pov=data.frame(x1= rep(seq(-4,0, by = 0.1),5),y1=c(MW_left_tail(distn="normal",size1 = 11,size2 = 12)$b,
                                                   MW_left_tail(distn="normal",size1 = 15,size2 = 17)$b,
                                                   MW_left_tail(distn="normal",size1 = 29,size2 = 42)$b,
                                                   MW_left_tail(distn="normal",size1 = 47,size2 = 55)$b,
                                                   MW_left_tail(distn="normal",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", N ( theta,1)," for ",theta <0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)

## cauchy
pov=data.frame(x1= rep(seq(-4,0, by = 0.1),5),y1=c(MW_left_tail(distn="cauchy",size1 = 11,size2 = 12)$b,
                                                   MW_left_tail(distn="cauchy",size1 = 15,size2 = 17)$b,
                                                   MW_left_tail(distn="cauchy",size1 = 29,size2 = 42)$b,
                                                   MW_left_tail(distn="cauchy",size1 = 47,size2 = 55)$b,
                                                   MW_left_tail(distn="cauchy",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", C ( theta,1)," for ",theta <0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)

## laplace
pov=data.frame(x1= rep(seq(-4,0, by = 0.1),5),y1=c(MW_left_tail(distn="laplace",size1 = 11,size2 = 12)$b,
                                                   MW_left_tail(distn="laplace",size1 = 15,size2 = 17)$b,
                                                   MW_left_tail(distn="laplace",size1 = 29,size2 = 42)$b,
                                                   MW_left_tail(distn="laplace",size1 = 47,size2 = 55)$b,
                                                   MW_left_tail(distn="laplace",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", Laplace ( theta,1)," for ",theta <0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## logistic
pov=data.frame(x1= rep(seq(-4,0, by = 0.1),5),y1=c(MW_left_tail(distn="logistic",size1 = 11,size2 = 12)$b,
                                                   MW_left_tail(distn="logistic",size1 = 15,size2 = 17)$b,
                                                   MW_left_tail(distn="logistic",size1 = 29,size2 = 42)$b,
                                                   MW_left_tail(distn="logistic",size1 = 47,size2 = 55)$b,
                                                   MW_left_tail(distn="logistic",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(0,4, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", logistic ( theta,1)," for ",theta <0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## Both sided
MW_both_sided<-function(distn,size1,size2,alpha=0.025,iter=1000,rge1=-3,rge2=3,...){
  effect_sizes <- seq(rge1,rge2, by = 0.1);power_values=numeric(length(effect_sizes))
  for(i in 1:length(effect_sizes)){
    reject<-0
    for(j in 1:iter){
      sample1=draw_sample(dist=distn,par1=0,par2=1,size=size1);sample2=draw_sample(dist=distn,par1=(0+effect_sizes[i]),par2=1,size=size1)
      if(wilcox.test(sample1,sample2)$p.value <alpha) reject=reject+1
    }
    power_values[i]=reject/iter
  }
  return(list(a=effect_sizes,b=power_values))
}
## normal
pov=data.frame(x1= rep(seq(-3,3, by = 0.1),5),y1=c(MW_both_sided(distn="normal",size1 = 11,size2 = 12)$b,
                                                   MW_both_sided(distn="normal",size1 = 15,size2 = 17)$b,
                                                   MW_both_sided(distn="normal",size1 = 29,size2 = 42)$b,
                                                   MW_both_sided(distn="normal",size1 = 47,size2 = 55)$b,
                                                   MW_both_sided(distn="normal",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(-3,3, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", N ( theta,1)," for ",theta != 0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## cauchy
pov=data.frame(x1= rep(seq(-3,3, by = 0.1),5),y1=c(MW_both_sided(distn="cauchy",size1 = 11,size2 = 12)$b,
                                                   MW_both_sided(distn="cauchy",size1 = 15,size2 = 17)$b,
                                                   MW_both_sided(distn="cauchy",size1 = 29,size2 = 42)$b,
                                                   MW_both_sided(distn="cauchy",size1 = 47,size2 = 55)$b,
                                                   MW_both_sided(distn="cauchy",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(-3,3, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", C ( theta,1)," for ",theta != 0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## laplace
pov=data.frame(x1= rep(seq(-3,3, by = 0.1),5),y1=c(MW_both_sided(distn="laplace",size1 = 11,size2 = 12)$b,
                                                   MW_both_sided(distn="laplace",size1 = 15,size2 = 17)$b,
                                                   MW_both_sided(distn="laplace",size1 = 29,size2 = 42)$b,
                                                   MW_both_sided(distn="laplace",size1 = 47,size2 = 55)$b,
                                                   MW_both_sided(distn="laplace",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(-3,3, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", laplace ( theta,1)," for ",theta != 0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)
## logistic
pov=data.frame(x1= rep(seq(-3,3, by = 0.1),5),y1=c(MW_both_sided(distn="logistic",size1 = 11,size2 = 12)$b,
                                                   MW_both_sided(distn="logistic",size1 = 15,size2 = 17)$b,
                                                   MW_both_sided(distn="logistic",size1 = 29,size2 = 42)$b,
                                                   MW_both_sided(distn="logistic",size1 = 47,size2 = 55)$b,
                                                   MW_both_sided(distn="logistic",size1 = 67,size2 = 85)$b),size=rep(c("m=11, n=12","m=15, n=17","m=29, n=35","m=42, n=55", "m=67, n=85"),each=length(seq(-3,3, by = 0.1))))
ggplot(data=pov,aes(x=x1,y=y1,col=size))+geom_line(lwd=.9)+scale_color_manual(values = c("cyan","violet","green","blue", "red"))+
  labs(x=expression(Theta),y="Power")+ggtitle(expression(paste("Power Function for alternative distribution  ", logistic ( theta,1)," for ",theta != 0)))+
  theme_stata()+geom_hline(yintercept = 0.05,linetype="dashed",lwd=0.5)

####################################################
# unbiasedness
set.seed(12345)
# Generate from baseline distribution F
rbase <- function(n, dist) {
  switch(
    dist,
    normal  = rnorm(n, mean = 0, sd = 1),
    cauchy  = rcauchy(n, location = 0, scale = 1),
    poisson = rpois(n, lambda = 4),
    stop("Unknown distribution")
  )
}

# One simulation run: two-sided Mann-Whitney test
mw_reject <- function(m, n, theta, dist, alpha) {
  x <- rbase(m, dist)
  y <- rbase(n, dist) + theta
  
  pval <- suppressWarnings(
    wilcox.test(
      x, y,
      alternative = "two.sided",
      exact = FALSE,
      correct = FALSE
    )$p.value
  )
  
  as.integer(pval <= alpha)
}

# Estimated power at one theta
mw_power <- function(m, n, theta, dist, alpha, nsim = 3000) {
  mean(replicate(nsim, mw_reject(m, n, theta, dist, alpha)))
}

# Full power curve over a grid of theta values
power_curve <- function(m, n, dist, alpha, theta_grid, nsim = 3000) {
  pow <- sapply(theta_grid, function(th) mw_power(m, n, th, dist, alpha, nsim))
  
  data.frame(
    dist = dist,
    m = m,
    n = n,
    alpha = alpha,
    theta = theta_grid,
    power = pow
  )
}

# Settings
alpha_vals <- c(0.025, 0.05, 0.1)
sample_sizes <- list(c(12, 17), c(60, 35))
dists <- c("normal", "cauchy", "poisson")

# theta grid for power curves
theta_grid <- seq(-2.5, 2.5, by = 0.05)

# Increase nsim if you want smoother curves
nsim <- 1000

# Run all simulations
all_results <- do.call(
  rbind,
  lapply(dists, function(dist) {
    do.call(
      rbind,
      lapply(sample_sizes, function(sz) {
        m <- sz[1]
        n <- sz[2]
        do.call(
          rbind,
          lapply(alpha_vals, function(a) {
            cat("Running:", dist, "m =", m, "n =", n, "alpha =", a, "\n")
            power_curve(m, n, dist, a, theta_grid, nsim = nsim)
          })
        )
      })
    )
  })
)
# Draw power curves
ggplot(all_results, aes(x = theta, y = power, color = alpha)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_grid(dist ~ paste0("(m,n)=(", m, ",", n, ")")) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  labs(
    x = expression(theta),
    y = "Power",
    color = expression(alpha),
    title = "Power curves of the two-sided Mann-Whitney test"
  ) +
  theme_bw(base_size = 12)
######################################################################

# Power Comparison: Z-test vs Mann-Whitney
####################################################
compare_power_t_vs_mannwhitney <- function(m, n, theta_vals, alpha = 0.025, iter = 1000) {
  power_t <- numeric(length(theta_vals))
  power_mw <- numeric(length(theta_vals))
  
  for (t in seq_along(theta_vals)) {
    theta <- theta_vals[t]
    reject_t <- 0
    reject_mw <- 0
    
    for (i in 1:iter) {
      # Generate samples: N(theta, 1) vs N(0, 1)
      sample1 <- rnorm(m, mean = theta, sd = 1)
      sample2 <- rnorm(n, mean = 0, sd = 1)
      
      # ---------- Parametric T-test ----------
      # Perform a one-sided t-test (H1: mu1 > mu2)
      t_res <- t.test(sample1, sample2, var.equal = TRUE, alternative = "greater")
      p_value_t <- t_res$p.value
      if (p_value_t < alpha) {
        reject_t <- reject_t + 1
      }
      
      # ---------- Mann-Whitney U-test (custom) ----------
      mw_stat <- m_whitney_ts(sample1, sample2)
      expected_mw <- m * n / 2
      sd_mw <- sqrt(m * n * (m + n + 1) / 12)
      z_mw <- (mw_stat - expected_mw) / sd_mw
      p_value_mw <- 1 - pnorm(z_mw)  # One-sided test
      
      if (p_value_mw < alpha) {
        reject_mw <- reject_mw + 1
      }
    }
    
    # Calculate empirical power for both tests
    power_t[t] <- reject_t / iter
    power_mw[t] <- reject_mw / iter
    
    cat("Completed for theta =", theta, "\n")
  }
  
  return(data.frame(theta = theta_vals, T_test_power = power_t, Mann_Whitney_power = power_mw))
}
##############
compare_power_z_vs_mannwhitney_less <- function(m, n, theta_vals, alpha = 0.025, iter = 1000) {
  
  power_z  <- numeric(length(theta_vals))
  power_mw <- numeric(length(theta_vals))
  
  for (t in seq_along(theta_vals)) {
    theta <- theta_vals[t]
    
    reject_z  <- 0
    reject_mw <- 0
    
    for (i in 1:iter) {
      
      # Generate samples: N(theta, 1) vs N(0, 1)
      sample1 <- rnorm(m, mean = theta, sd = 1)
      sample2 <- rnorm(n, mean = 0, sd = 1)
      
      # -----------------------------------
      # Z-test (known variance = 1)
      # H1: mu1 < mu2  ⇔ theta < 0
      # -----------------------------------
      xbar <- mean(sample1)
      ybar <- mean(sample2)
      
      z_stat <- (xbar - ybar) / sqrt(1/m + 1/n)
      
      # Left-tailed test
      p_value_z <- pnorm(z_stat)
      
      if (p_value_z < alpha) {
        reject_z <- reject_z + 1
      }
      
      # -----------------------------------
      # Mann–Whitney Test
      # -----------------------------------
      mw_stat <- m_whitney_ts(sample1, sample2)
      
      mu_mw <- m * n / 2
      sd_mw <- sqrt(m * n * (m + n + 1) / 12)
      
      z_mw <- (mw_stat - mu_mw) / sd_mw
      
      # LEFT-tailed now (important change)
      p_value_mw <- pnorm(z_mw)
      
      if (p_value_mw < alpha) {
        reject_mw <- reject_mw + 1
      }
    }
    
    power_z[t]  <- reject_z / iter
    power_mw[t] <- reject_mw / iter
    
    cat("Completed for theta =", theta, "\n")
  }
  
  return(data.frame(
    theta = theta_vals,
    Z_test_power = power_z,
    Mann_Whitney_power = power_mw
  ))
}
######################
compare_power_z_vs_mannwhitney_twosided <- function(m, n, theta_vals, alpha = 0.05, iter = 1000) {
  
  power_z  <- numeric(length(theta_vals))
  power_mw <- numeric(length(theta_vals))
  
  for (t in seq_along(theta_vals)) {
    theta <- theta_vals[t]
    
    reject_z  <- 0
    reject_mw <- 0
    
    for (i in 1:iter) {
      
      # Generate samples
      sample1 <- rnorm(m, mean = theta, sd = 1)
      sample2 <- rnorm(n, mean = 0, sd = 1)
      
      # -----------------------------------
      # Z-test (two-sided)
      # -----------------------------------
      xbar <- mean(sample1)
      ybar <- mean(sample2)
      
      z_stat <- (xbar - ybar) / sqrt(1/m + 1/n)
      
      # Two-sided p-value
      p_value_z <- 2 * (1 - pnorm(abs(z_stat)))
      
      if (p_value_z < alpha) {
        reject_z <- reject_z + 1
      }
      
      # -----------------------------------
      # Mann–Whitney test (two-sided)
      # -----------------------------------
      mw_stat <- m_whitney_ts(sample1, sample2)
      
      mu_mw <- m * n / 2
      sd_mw <- sqrt(m * n * (m + n + 1) / 12)
      
      z_mw <- (mw_stat - mu_mw) / sd_mw
      
      # Two-sided p-value
      p_value_mw <- 2 * (1 - pnorm(abs(z_mw)))
      
      if (p_value_mw < alpha) {
        reject_mw <- reject_mw + 1
      }
    }
    
    power_z[t]  <- reject_z / iter
    power_mw[t] <- reject_mw / iter
    
    cat("Completed for theta =", theta, "\n")
  }
  
  return(data.frame(
    theta = theta_vals,
    Z_test_power = power_z,
    Mann_Whitney_power = power_mw
  ))
}
# Plotting - 1
theta_vals <- seq(0, 3, by = 0.1)  # Range of effect sizes
m <- 11   # Sample size for group 1
n <- 19   # Sample size for group 2
iter <- 1000  # Number of iterations per theta value

# Run the power comparison simulation
power_results <- compare_power_t_vs_mannwhitney(m, n, theta_vals, alpha = 0.025, iter = iter)

# Prepare the data for plotting
power_long <- melt(power_results, id.vars = "theta", variable.name = "Test", value.name = "Power")

# Plot the power curves
ggplot(power_long, aes(x = theta, y = Power, color = Test, group = Test)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = bquote("Power Comparison for N(" * theta * ",1) vs N(0,1)"),
    subtitle = "Sample sizes: m = 11, n = 19; Alpha = 0.025" ,
    x = expression(theta),
    y = "Empirical Power"
  ) +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray")
# Plotting - 2
theta_vals <- seq(-3, 0, by = 0.1)  # Range of effect sizes
m <- 11   # Sample size for group 1
n <- 19   # Sample size for group 2
iter <- 1000  # Number of iterations per theta value

# Run the power comparison simulation
power_results <- compare_power_z_vs_mannwhitney_less(m, n, theta_vals, alpha = 0.025, iter = iter)
power_long <- melt(power_results, id.vars = "theta", variable.name = "Test", value.name = "Power")

# Plot the power curves
ggplot(power_long, aes(x = theta, y = Power, color = Test, group = Test)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = bquote("Power Comparison for N(" * theta * ",1) vs N(0,1)"),
    subtitle = "Sample sizes: m = 11, n = 19; Alpha = 0.025",
    x = expression(theta),
    y = "Empirical Power"
  ) +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray")
# Plotting - 3
theta_vals <- seq(-3, 3, by = 0.1)  # Range of effect sizes
m <- 11   # Sample size for group 1
n <- 19   # Sample size for group 2
iter <- 1000  # Number of iterations per theta value

# Run the power comparison simulation
power_results <- compare_power_z_vs_mannwhitney_twosided(m, n, theta_vals, alpha = 0.05, iter = iter)
power_long <- melt(power_results, id.vars = "theta", variable.name = "Test", value.name = "Power")
# Plot the power curves
ggplot(power_long, aes(x = theta, y = Power, color = Test, group = Test)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = bquote("Power Comparison for N(" * theta * ",1) vs N(0,1)"),
    subtitle = "Sample sizes: m = 11, n = 19; Alpha = 0.05",
    x = expression(theta),
    y = "Empirical Power"
  ) +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray")
##################################
# normal sample but unequal variance
theta_vals <- seq(-3, 3, by = 0.1)  # Range of effect sizes
m <- 30   # Sample size for group 1
n <- 30   # Sample size for group 2
iter <- 1000  # Number of iterations per theta value

compare_power_z_vs_mannwhitney_unequal_normal <- function(m, n, theta_vals, alpha = 0.05, iter = 1000) {
  
  power_z  <- numeric(length(theta_vals))
  power_mw <- numeric(length(theta_vals))
  
  for (t in seq_along(theta_vals)) {
    theta <- theta_vals[t]
    
    p_welch <- numeric(iter)
    p_mw <- numeric(iter)
    
    for (j in 1:iter) {
      # Generate samples
      x <- rnorm(m, mean = 0, sd = 1)
      y <- rnorm(n, mean = theta, sd = 9)
      
      # Perform Welch's t-test (Default in R is var.equal = FALSE)
      p_welch[j] <- t.test(x, y, var.equal = FALSE)$p.value
      
      # Perform Mann-Whitney U test
      p_mw[j] <- wilcox.test(x, y)$p.value
    }
    
    power_z[t]  <- mean(p_welch < alpha)
    power_mw[t] <- mean(p_mw < alpha)
    
    cat("Completed for theta =", theta, "\n")
  }
  
  return(data.frame(
    theta = theta_vals,
    t_test_power = power_z,
    Mann_Whitney_power = power_mw
  ))
}

# Run the power comparison simulation
power_results <- compare_power_z_vs_mannwhitney_unequal_normal(m, n, theta_vals, alpha = 0.05, iter = iter)
power_long <- melt(power_results, id.vars = "theta", variable.name = "Test", value.name = "Power")
# Plot the power curves
ggplot(power_long, aes(x = theta, y = Power, color = Test, group = Test)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = bquote("Power Comparison between Welch T-test and Mann–Whitney Test"),
    subtitle = "Sample sizes: m = 30, n = 30; Alpha = 0.05",
    x = expression(theta),
    y = "Empirical Power"
  ) +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray")

##########################################
# General Simulation Function
sim_power <- function(dist_type, n1=5, n2=7, theta_vals=seq(0, 3, 0.1), n_sim=1000) {
  results <- data.frame()
  
  for (theta in theta_vals) {
    p_welch <- numeric(n_sim)
    p_mw <- numeric(n_sim)
    
    for (j in 1:n_sim) {
      if (dist_type == "cauchy") {
        x <- rcauchy(n1, location = 0)
        y <- rcauchy(n2, location = theta)
      } else if (dist_type == "exponential") {
        # Using a rate shift; scale = 1/rate
        x <- rexp(n1, rate = 1)
        y <- rexp(n2, rate = 1 / (1 + theta)) 
      } else if (dist_type == "laplace") {
        x <- rlaplace(n1, mu = 0, sigma = 1)
        y <- rlaplace(n2, mu = theta_vals, sigma = 1)
      }
      
      # Suppress warnings for ties in Mann-Whitney (common in Poisson/Discrete)
      p_welch[j] <- suppressWarnings(t.test(x, y)$p.value)
      p_mw[j] <- suppressWarnings(wilcox.test(x, y)$p.value)
    }
    
    results <- rbind(results, data.frame(
      theta = theta,
      Welch = mean(p_welch < 0.05),
      Mann_Whitney = mean(p_mw < 0.05),
      Distribution = dist_type
    ))
  }
  return(results)
}

# Run simulations
res_cauchy <- sim_power("cauchy")
res_expo   <- sim_power("exponential")
res_laplace   <- sim_power("laplace")

# Combine and Plot
all_results <- rbind(res_cauchy, res_expo, res_laplace)
all_long <- pivot_longer(all_results, cols = c("Welch", "Mann_Whitney"), 
                         names_to = "Test", values_to = "Power")

ggplot(all_long, aes(x = theta, y = Power, color = Test, linetype = Test)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~Distribution, scales = "free_x") +
  theme_minimal() +
  labs(title = "Power Comparison for Non-Normal Distributions",
       subtitle = "Comparing Welch t-test vs Mann-Whitney U",
       x = expression(Delta * " (Shift Parameter)"), y = "Power (Rejection Rate)")

###################################################################
# High-Dimensional 2-Sample Location Problem Simulation

# Simulation Setup
set.seed(42) # For reproducibility

p <- 1000          # Total number of features (Dimensions)
n1 <- 25           # Sample size of Group 1
n2 <- 25           # Sample size of Group 2
n_signal <- 5      # Number of features with a true difference (Sparse signal)
delta <- 2.5       # The magnitude of the location shift

# Generate heavy-tailed data (t-distribution with 2 degrees of freedom)
# Group 1: All 1000 features are noise
X <- matrix(rt(n1 * p, df = 2), nrow = n1, ncol = p)

# Group 2: First 5 features have the signal, remaining 995 are noise
Y <- matrix(rt(n2 * p, df = 2), nrow = n2, ncol = p)
Y[, 1:n_signal] <- Y[, 1:n_signal] + delta 

# Compute Statistics Feature-by-Feature
pvals_t <- numeric(p)
pvals_mw <- numeric(p)

# Loop through all 1000 dimensions
for(k in 1:p) {
  # Classical t-test
  pvals_t[k] <- t.test(X[, k], Y[, k])$p.value
  
  # Mann-Whitney U test (Wilcoxon Rank Sum)
  # exact=FALSE to use asymptotic approximation, typical for large p simulations
  pvals_mw[k] <- wilcox.test(X[, k], Y[, k], exact = FALSE)$p.value 
}

# 3. Prepare Data for Visualization
# Calculate Bonferroni threshold (-log10 scale)
alpha <- 0.05
bonferroni_limit <- -log10(alpha / p)

# Create a data frame for plotting
df_plot <- data.frame(
  Feature = rep(1:p, 2),
  NegLogP = c(-log10(pvals_t), -log10(pvals_mw)),
  Method = rep(c("Classical T-Test", "Mann-Whitney (Rank-Based)"), each = p),
  Type = rep(ifelse(1:p <= n_signal, "True Signal", "Noise"), 2)
)

# 4. Create the Plot
ggplot(df_plot, aes(x = Feature, y = NegLogP, color = Type)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = bonferroni_limit, linetype = "dashed", color = "red", size = 0.8) +
  facet_wrap(~ Method, ncol = 1) +
  scale_color_manual(values = c("Noise" = "gray70", "True Signal" = "#E69F00")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "High-Dimensional Simulation: T-Test vs. Mann-Whitney",
    subtitle = "Comparing performance under heavy-tailed (df=2) noise",
    x = "Feature Index (Dimensions 1 to 1000)",
    y = expression(-log[10](p-value)),
    color = "Feature Category"
  ) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
