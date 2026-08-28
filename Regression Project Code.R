# installing required packages

install.packages("leaps")
install.packages("olsrr")
install.packages("lmtest")
install.packages("readxl")
install.packages("dplyr")
install.packages("zoo")
install.packages("leaps")
install.packages("car")
install.packages("ggplot2")
install.packages("ggcorrplot")
install.packages("MASS")
install.packages("lattice")

# Loading packages
library(leaps)
library(olsrr)
library(lmtest)
library(readxl)
library(dplyr)
library(zoo)
library(leaps)
library(car)
library(ggplot2)
library(ggcorrplot)
library(MASS)
library(lattice)

#########################################################

# Loading the data
data <- read.table("chicago.txt")

# summary of data
summary(data)

########################################################

# Creating histogram for each attributes
numeric_data <- data[, sapply(data, is.numeric)]
 
# Plot histograms for each numeric variable
par(mfrow = c(2, 4))  # Arrange plots in 2x4 grid (adjust based on number of columns)
for(col in colnames(numeric_data)) {
  hist(numeric_data[[col]],
       main = paste(col),
       xlab = col,
       col = "orange",
       border = "black")
}

# Plot histogram for log(income)
par(mfrow = c(1,1))
hist(log(data$income),
     main = "log(income)",
     xlab = "log(income)",
     col = "orange",
     border = "black")

#######################################################

# Log transforming income
data$log_income <- log(data$income)
data <- subset(data, select = -income)

#######################################################

# Pair Plots
plot(data)

# Selected Pair Plots
par(mfrow = c(2, 4))

# 1. Fire vs Theft
plot(data$fire, data$theft, main = "Fire vs Theft", 
     xlab = "Fire", ylab = "Theft", pch = 16, col = "red")
abline(lm(theft ~ fire, data = data), col = "darkred", lwd = 2)

# 2. Race vs Involact
plot(data$race, data$involact, main = "Race vs Involact", 
     xlab = "Race", ylab = "Involact", pch = 16, col = "blue")
abline(lm(involact ~ race, data = data), col = "darkblue", lwd = 2)

# 3. Age vs Volact
plot(data$age, data$volact, main = "Age vs Volact", 
     xlab = "Age", ylab = "Volact", pch = 16, col = "green")
abline(lm(volact ~ age, data = data), col = "darkgreen", lwd = 2)

# 4. Volact vs Log-income
plot(data$volact, data$log_income, main = "Volact vs Log-income", 
     xlab = "Volact", ylab = "Log-income", pch = 16, col = "purple")
abline(lm(log_income ~ volact, data = data), col = "darkviolet", lwd = 2)

# 5. Log-income vs Race
plot(data$race, data$log_income, main = "Log-income vs Race", 
     xlab = "Race", ylab = "Log-income", pch = 16, col = "orange")
abline(lm(log_income ~ race, data = data), col = "darkorange", lwd = 2)

# 6. Log-income vs Fire
plot(data$fire, data$log_income, main = "Log-income vs Fire", 
     xlab = "Fire", ylab = "Log-income", pch = 16, col = "orange")
abline(lm(log_income ~ fire, data = data), col = "darkorange", lwd = 2)

# 7. Log-income vs Theft
plot(data$theft, data$log_income, main = "Log-income vs Theft", 
     xlab = "Theft", ylab = "Log-income", pch = 16, col = "orange")
abline(lm(log_income ~ theft, data = data), col = "darkorange", lwd = 2)

# 8. Log-income vs Involact
plot(data$involact, data$log_income, main = "Log-income vs Involact", 
     xlab = "Involact", ylab = "Log-income", pch = 16, col = "orange")
abline(lm(log_income ~ involact, data = data), col = "darkorange", lwd = 2)

#######################################################

# creating boxplot
par(mfrow=c(2,4))
for(i in 1:7)
  boxplot(data[,i],main=names(data)[i])

#######################################################

# Compute correlation matrix
par(mfrow = c(1,1))
cor_matrix <- cor(data)

# Round and view
round(cor_matrix, 3)

# Simple correlation plot
ggcorrplot(cor_matrix, lab = TRUE, title = "Correlation Matrix", colors = c( "red", "white","green"))

#####################################################

# Initial model building
# Model for volact
model_volact <- lm(volact ~ race + fire + theft + age + log_income + involact, data=data)
summary(model_volact)
anova(model_volact)
par(mfrow=c(2,2)); plot(model_volact)

# Model for involact
model_involact <- lm(involact ~ race + fire + theft + age + log_income + volact, data=data)
summary(model_involact)
anova(model_involact)
par(mfrow=c(2,2)); plot(model_involact)

########################################################

## Normality test: QQ plot for the residuals

# For model volact
par(mfrow = c(1,1))
residuals_volact <- model_volact$residuals
qqnorm(residuals_volact, main = "QQ Plot for volact residuals")
qqline(residuals_volact, col = "red")

# For model involact
residuals_involact <- model_involact$residuals
qqnorm(residuals_involact, main = "QQ Plot for involact residuals")
qqline(residuals_involact, col = "red")

## Normality test: Sapiro-Wilk Test

# For model volact
shapiro.test(residuals_volact)
# For model involact
shapiro.test(residuals_involact)

#########################################################

# Box-cox transformation for volact
# Find the minimum value of the response
min_volact <- min(data$volact)

# Shift if needed
if (min_volact <= 0) {
  shift <- abs(min_volact) + 1
  data$volact_pos <- data$volact + shift
  cat("Shifted response by:", shift, "\n")
} else {
  data$volact_pos <- data$volact
}

# Refit model using the shifted response
model_volact_pos <- lm(volact_pos ~ race + fire + theft + age + involact + log_income, data = data)

# Perform Box–Cox transformation
bc <- boxcox(model_volact_pos, lambda = seq(-2, 2, 0.1), title = "For volact")
title(main = "Box-Cox Transformation Plot for volact",
      xlab = expression(lambda))

# Find best lambda
lambda_opt <- bc$x[which.max(bc$y)]
lambda_opt

## transformation
if (abs(lambda_opt) < 1e-4) {
  data$volact_bc <- log(data$volact_pos)
} else {
  data$volact_bc <- (data$volact_pos^lambda_opt - 1) / lambda_opt
}

# Fit model again
model_volact_bc <- lm(volact_bc ~ race + fire + theft + age + involact + log_income, data = data)
summary(model_volact_bc)

# Box-cox transformation for involact
# Find the minimum value of the response
min_involact <- min(data$involact)

# Shift if needed
if (min_involact <= 0) {
  shift <- abs(min_involact) + 1
  data$involact_pos <- data$involact + shift
  cat("Shifted response by:", shift, "\n")
} else {
  data$involact_pos <- data$involact
}

# Refit model using the shifted response
model_involact_pos <- lm(involact_pos ~ race + fire + theft + age + volact + log_income, data = data)

# Perform Box–Cox transformation 
bc <- boxcox(model_involact_pos, lambda = seq(-2, 2, 0.1))
title(main = "Box-Cox Transformation Plot for involact",
      xlab = expression(lambda))

# Find best lambda
lambda_opt <- bc$x[which.max(bc$y)]
lambda_opt

# transformation
if (abs(lambda_opt) < 1e-4) {
  data$involact_bc <- log(data$involact_pos)
} else {
  data$involact_bc <- (data$involact_pos^lambda_opt - 1) / lambda_opt
}

# Fit model again
model_involact_bc <- lm(involact_bc ~ race + fire + theft + age + volact + log_income, data = data)
summary(model_involact_bc)

########################################################

# QQ Plot and Shapiro Wilk test after 
# the Box-Cox transformation on both models
# current models are model_volact_bc and model_involact_bc
# Do the testing and Plotting using the previous codes

########################################################

# Testing for Serial correlation

# Durbin-Watson test

# Serial correlation for the residuals of volact
dwtest(model_volact_bc) # Durbin-Watson test (AR1)

# Serial correlation for the residuals of involact
dwtest(model_involact_bc) # Durbin-Watson test (AR1)

# Breusch-Godfrey test for up to lag 2
bgtest(model_volact_bc, order = 2) # model for volact
bgtest(model_involact_bc, order = 2) # model for involact

# Breusch-Godfrey test for up to lag 3
bgtest(model_volact_bc, order = 3) # model for volact
bgtest(model_involact_bc, order = 3) # model for involact

# Breusch-Godfrey test for up to lag 4
bgtest(model_volact_bc, order = 4) # model for volact
bgtest(model_involact_bc, order = 4) # model for involact

# Residual autocorrelation plot

# Autocorrelation function plot (ACF)
acf(residuals(model_volact_bc), main = "Autocorrelation of Residuals for model_volact_bc") # model for volact
acf(residuals(model_involact_bc), main = "Autocorrelation of Residuals for model_involact_bc") # model for involact

# Partial Autocorrelation function plot (PACF)
pacf(residuals(model_volact_bc), main = "Partial Autocorrelation of Residuals for model_volact_bc") # PACF for volact model
pacf(residuals(model_involact_bc), main = "Partial Autocorrelation of Residuals model_involact_bc") # PACF for involact model

#########################################################

# Multicollinearity check
# 1. VIF
vif(model_volact_bc)
vif(model_involact_bc)

# 2. Barplots for VIF per attribute
barplot(vif(model_volact_bc),ylim=c(0,6), col = "skyblue", main = "model_volact_bc")
abline(h=5)

barplot(vif(model_involact_bc),ylim=c(0,6), col = "skyblue", main = "model_involact_bc")
abline(h=5)

# 3. Condition number
kappa(model_volact_bc)
kappa(model_involact_bc)

#######################################################

# Heteroscedasticity check

# Bruesch Pagan Test
bptest(model_volact_bc)
bptest(model_involact_bc)

# Non-constant Variance Score Test 
ncvTest(model_volact_bc)
ncvTest(model_involact_bc)

######################################################

# Generalized Least Square Estimation
# Function for one GLS iteration
gls_iter <- function(model, data, i) {
  sq_resi <- residuals(model)^2
  fitted_vals <- fitted(model)
  ord <- order(fitted_vals)
  
  y_ma <- rollmean(sq_resi[ord], k = 3, fill = NA)
  weights <- ifelse(is.na(y_ma), mean(y_ma, na.rm = TRUE), y_ma)
  idx <- paste("Iteration", as.character(i))
  plot(fitted_vals[ord], sq_resi[ord],
       main = idx,
       xlab = "Fitted Values", ylab = "Squared Residuals",
       col = "blue")
  lines(fitted_vals[ord], y_ma, col = "red", lwd = 2)
  
  lm(volact_bc ~ race + fire + theft + age + involact + log_income,
     data = data, weights = weights)
}

# --- Iterations ---
par(mfrow = c(3,3))
model_1 <- model_volact_bc
model_2 <- gls_iter(model_1, data, 1)
model_3 <- gls_iter(model_2, data, 2)
model_4 <- gls_iter(model_3, data, 3)
model_5 <- gls_iter(model_4, data, 4)
model_6 <- gls_iter(model_5, data, 5)
model_7 <- gls_iter(model_6, data, 6)
model_8 <- gls_iter(model_7, data, 7)
model_9 <- gls_iter(model_8, data, 8)
model_10 <- gls_iter(model_9, data, 9)

# --- Final Weights ---
sq_resi_10 <- residuals(model_10)^2
fitted_vals_10 <- fitted(model_10)
ord10<- order(fitted_vals_10)
y_ma_10 <- rollmean(sq_resi_10[ord10], k = 10, fill = NA)
final_weights <- ifelse(is.na(y_ma_10), mean(y_ma_10, na.rm = TRUE), y_ma_10)
# dataframe with the weights
weights_df <- data.frame(
  zipcodes = names(final_weights),
  weights = as.numeric(final_weights)
)
print(weights_df)

# fitting WLSE for volact
model_volact_wls <- lm(volact_bc ~ race + fire + theft + age + involact + log_income,
                       data = data,
                       weights = final_weights)

summary(model_volact_wls)

# running again BP test
bptest(model_volact_wls)

#########################################################

par(mfrow = c(1,1))

# Partial residual plot
crPlots(model_volact_wls,main='Partial Residual Plot for Volact')
crPlots(model_involact_bc,main='Partial Residual Plot for Involact')

# Added variable plot
avPlots(model_volact_wls, main = "Added Variable Plot for Volact")
avPlots(model_involact_bc, main = "Added Variable Plot for Involact")

########################################################

# Outlier detection
par(mfrow = c(1,1))
# hat matrix diagonal - high leverage points for model_volact_wls

# Compute average leverage
n <- nrow(data)
p <- length(coef(model_volact_wls))
avg_h_volact <- p / n
plot(hatvalues(model_volact_wls),
     ylab="Hat Diagonals",
     xlab="Zip.code",
     main="Comparing hat matrix diagonals for model volact_wls",
     pch=19,col="blue"
)
text(hatvalues(model_volact_wls), pos = 4, cex = 0.8, col = "black")
# Add average leverage line
abline(h = avg_h_volact, col = "red", lwd = 2, lty = 2)
abline(h = 2 * avg_h_volact, col = "orange", lwd = 2, lty = 3)
abline(h = 3 * avg_h_volact, col = "darkred", lwd = 2, lty = 3)

# hat matrix diagonal - high leverage points for model_involact_bc

# Compute average leverage
n <- nrow(data)
p <- length(coef(model_involact_bc))
avg_h_involact <- p / n
plot(hatvalues(model_involact_bc),
     ylab="Hat Diagonals",
     xlab="Zip.code",
     main="Comparing hat matrix diagonals for model involact_bc",
     pch=19,col="blue"
)
text(hatvalues(model_involact_bc), pos = 4, cex = 0.8, col = "black")
# Add average leverage line
abline(h = avg_h_involact, col = "red", lwd = 2, lty = 2)
abline(h = 2 * avg_h_involact, col = "orange", lwd = 2, lty = 3)
abline(h = 3 * avg_h_involact, col = "darkred", lwd = 2, lty = 3)

#######################################################

# DFFITS

plot_dffits_with_labels <- function(model, data, model_name, label.col = "zip") {
  # Compute DFFITS
  dffit_vals <- dffits(model)
  
  # Indices corresponding to dffit_vals
  idx_all <- seq_along(dffit_vals)
  
  # Construct labels aligned with dffit_vals
  if (!is.null(rownames(data)) && length(rownames(data)) >= length(dffit_vals)) {
    labels_all <- rownames(data)[idx_all]
  } else if (label.col %in% names(data) && length(data[[label.col]]) >= length(dffit_vals)) {
    labels_all <- as.character(data[[label.col]][idx_all])
  }
  
  # Get p and n consistent with dffit vector
  p <- length(coef(model))
  n <- length(dffit_vals)
  threshold <- 2 * sqrt(p / n)
  
  # Plot
  plot(
    dffit_vals,
    main = paste("DFFITS Plot -", model_name),
    ylab = "DFFITS values",
    xlab = "Index",
    pch = 19,
    col = "blue"
  )
  abline(h = c(-threshold, threshold), col = "red", lty = 2)
  
  # Find indices exceeding threshold
  above_threshold <- which(abs(dffit_vals) > threshold)
  
  # Add labels only if there are any
  if (length(above_threshold) > 0) {
    text(
      x = above_threshold,
      y = dffit_vals[above_threshold],
      labels = labels_all[above_threshold],
      pos = 2,
      cex = 0.8,
      col = "black"
    )
  }
  
  # Return data frame (empty if none)
  return(data.frame(ZIP = labels_all[above_threshold],
                    DFFITS = as.numeric(dffit_vals[above_threshold]),
                    row.names = NULL))
}

influential_volact <- plot_dffits_with_labels(model_volact_wls, data, "Volact Model", label.col = "zip")
influential_involact <- plot_dffits_with_labels(model_involact_bc, data, "Involact Model", label.col = "zip")


# Inspect results
influential_volact
influential_involact

##########################################################

# Externally studentized residuals

# Compute externally studentized residuals
studentized_residuals_volact <- rstudent(model_volact_wls)
studentized_residuals_involact <- rstudent(model_involact_bc)

# Set threshold for outliers
threshold <- 2

# Plot for volact model
plot(
  studentized_residuals_volact,
  main = "Externally Studentized Residuals (volact_wls model)",
  ylab = "Studentized Residuals",
  xlab = "ZIP Code",
  pch = 19, col = "steelblue"
)

# Add reference lines
abline(h = c(-threshold, threshold), col = "red", lty = 2)

# Identify influential points
outliers_volact <- which(abs(studentized_residuals_volact) > threshold)

# Label the outliers
text(
  outliers_volact,
  studentized_residuals_volact[outliers_volact],
  labels = rownames(data)[outliers_volact],
  pos = 3, col = "black", cex = 0.8
)


# Plot for involact model
plot(
  studentized_residuals_involact,
  main = "Externally Studentized Residuals (involact_bc model)",
  ylab = "Studentized Residuals",
  xlab = "ZIP Code",
  pch = 19, col = "darkgreen"
)

# Add reference lines
abline(h = c(-threshold, threshold), col = "red", lty = 2)

# Identify influential points
outliers_involact <- which(abs(studentized_residuals_involact) > threshold)

# Label the outliers
text(
  outliers_involact,
  studentized_residuals_involact[outliers_involact],
  labels = rownames(data)[outliers_involact],
  pos = 3, col = "black", cex = 0.8
)

######################################################

# Cook's distance
# Compute Cook's distance
cooks_d_volact <- cooks.distance(model_volact_wls)
cooks_d_involact <- cooks.distance(model_involact_bc)

# plot volact
outliers_volact <- which(cooks_d_volact > 4 / length(cooks_d_volact))
ol_cook_v <- outliers_volact
plot(cooks_d_volact,main = "Cook's Distance for model volact_wls", ylab = "Cook's Distance",ylim=c(0,3),
     col="blue",)
abline(h = 4 / length(cooks_d_volact), col = "red",lty = 2)
text(outliers_volact, cooks_d_volact[outliers_volact],
     labels = outliers_volact, pos = 4, cex = 0.7, col = "darkred")
legend("topright",
       legend = c("Threshold = 4/n"),
       col = "red", lwd = 2, lty = 2, bty = "n")

# plot involact
outliers_involact <- which(cooks_d_involact > 4 / length(cooks_d_involact))
ol_cook_inv <- outliers_involact
plot(cooks_d_involact,main = "Cook's Distance for model involact_bc", ylab = "Cook's Distance",ylim=c(0,3),
     col="blue",)
abline(h = 4 / length(cooks_d_involact), col = "red",lty = 2)
text(outliers_involact, cooks_d_involact[outliers_involact],
     labels = outliers_involact, pos = 4, cex = 0.7, col = "darkred")
legend("topright",
       legend = c("Threshold = 4/n"),
       col = "red", lwd = 2, lty = 2, bty = "n")

######################################################

# Influential plot for model_volact_wls
influencePlot(model_volact_wls, main = "Influence Plot for model volact_wls",
              sub = "Circle size ∝ Cook’s Distance")

# Influential plot for model_involact_bc
influencePlot(model_involact_bc, main = "Influence Plot for model involact_bc",
              sub = "Circle size ∝ Cook’s Distance")

#######################################################

# removing the outliers: Zip 60610, 60607, 60611
newdata <- data[!(rownames(data) %in%
                    c("60610", "60607", "60611")), ]

# Run all the previous codes for 
# Box-Cox transformations
# on newdata and get the new models
# check model summary and plots
# check normality
# check heteroscedasticity
# check serial correlation
# check multicollinearity

########################################################

# Model Selection

# APR based on Mallow's Cp

# For model_volact_bc without outliers on newdata

# Run all possible subset regressions (APR, best model per subset size)

fit <- regsubsets(volact_bc ~ age + log_income + fire + theft + race + involact,
                  data = newdata, nbest = 1)

# Get summary output
fit_summary <- summary(fit)

# Plot Cp values

plot(fit_summary$cp,
     xlab = "Number of Predictors",
     ylab = "Mallows' Cp",
     main = "Mallows' Cp Plot",
     pch = 19, col = "blue")
abline(a = 0, b = 1, col = "red", lty = 2)

# Highlight minimum Cp

min_cp <- which.min(fit_summary$cp)
points(min_cp, fit_summary$cp[min_cp], col = "red", pch = 19, cex = 1.2)
text(min_cp, fit_summary$cp[min_cp],
     labels = paste("Min Cp =", round(fit_summary$cp[min_cp], 2)),
     pos = 3, col = "red")

#Best Model for Each Subset Size (1 to 6) 

cat("\nBest Model for Each Number of Predictors (p = 1 to 6):\n")

num_pred <- apply(fit_summary$which[, -1, drop = FALSE], 1, sum)
pred_names <- colnames(fit_summary$which)[-1]

for (p in 1:6) {
  idx <- which(num_pred == p)[1]
  if (!is.na(idx)) {
    sel_vars <- pred_names[fit_summary$which[idx, -1, drop = FALSE]]
    sel_vars <- sel_vars[sel_vars != "(Intercept)"]  # remove intercept if any
    cat("\nModel with", p, "predictor(s):\n")
    cat("Predictors:", paste(sel_vars, collapse = ", "), "\n")
    cat("Cp =", round(fit_summary$cp[idx], 3), "\n")
  }
}

# Highlight model with minimum Cp

min_cp <- which.min(fit_summary$cp)
points(min_cp, fit_summary$cp[min_cp], col = "red", pch = 19, cex = 1.2)
text(min_cp, fit_summary$cp[min_cp],
     labels = paste("Min Cp =", round(fit_summary$cp[min_cp], 2)),
     pos = 3, col = "red")

# Show which predictors are in the best model
fit_summary$which[min_cp, ]

# For model_involact_bc without outliers on newdata

# Run all possible subset regressions (APR, best model per subset size)

fit <- regsubsets(involact_bc ~ age + log_income + fire + theft + race + volact,
                  data = newdata, nbest = 1)

# Get summary output
fit_summary <- summary(fit)

# Plot Cp values

plot(fit_summary$cp,
     xlab = "Number of Predictors",
     ylab = "Mallows' Cp",
     main = "Mallows' Cp Plot",
     pch = 19, col = "blue")
abline(a = 0, b = 1, col = "red", lty = 2)

# Highlight minimum Cp

min_cp <- which.min(fit_summary$cp)
points(min_cp, fit_summary$cp[min_cp], col = "red", pch = 19, cex = 1.2)
text(min_cp, fit_summary$cp[min_cp],
     labels = paste("Min Cp =", round(fit_summary$cp[min_cp], 2)),
     pos = 3, col = "red")

#Best Model for Each Subset Size (1 to 6) 

cat("\nBest Model for Each Number of Predictors (p = 1 to 6):\n")

num_pred <- apply(fit_summary$which[, -1, drop = FALSE], 1, sum)
pred_names <- colnames(fit_summary$which)[-1]

for (p in 1:6) {
  idx <- which(num_pred == p)[1]
  if (!is.na(idx)) {
    sel_vars <- pred_names[fit_summary$which[idx, -1, drop = FALSE]]
    sel_vars <- sel_vars[sel_vars != "(Intercept)"]  # remove intercept if any
    cat("\nModel with", p, "predictor(s):\n")
    cat("Predictors:", paste(sel_vars, collapse = ", "), "\n")
    cat("Cp =", round(fit_summary$cp[idx], 3), "\n")
  }
}

# Highlight model with minimum Cp

min_cp <- which.min(fit_summary$cp)
points(min_cp, fit_summary$cp[min_cp], col = "red", pch = 19, cex = 1.2)
text(min_cp, fit_summary$cp[min_cp],
     labels = paste("Min Cp =", round(fit_summary$cp[min_cp], 2)),
     pos = 3, col = "red")

# Show which predictors are in the best model
fit_summary$which[min_cp, ]

#######################################################

# APR based on AIC/BIC

# For model volact

# Run all possible subset regressions (APR, best model per subset size)
fit <- regsubsets(volact_bc ~ age + log_income + fire + theft + race + involact,
                  data = newdata, nbest = 1)

fit_summary <- summary(fit)

# Prepare to compute AIC and BIC for each model

num_models <- nrow(fit_summary$which)
num_pred <- apply(fit_summary$which[, -1, drop = FALSE], 1, sum)
pred_names <- colnames(fit_summary$which)[-1]

# Create vectors to store AIC and BIC
aic_vals <- numeric(num_models)
bic_vals <- numeric(num_models)

# Compute AIC and BIC for each model
for (i in 1:num_models) {
  sel_vars <- pred_names[fit_summary$which[i, -1, drop = FALSE]]
  sel_vars <- sel_vars[sel_vars != "(Intercept)"]
  
  # Fit the model
  formula_str <- paste("volact_bc ~", paste(sel_vars, collapse = " + "))
  model <- lm(as.formula(formula_str), data = newdata)
  
  # Store AIC and BIC
  aic_vals[i] <- AIC(model)
  bic_vals[i] <- BIC(model)
}

# Plot AIC values
plot(aic_vals,
     xlab = "Model Index",
     ylab = "AIC",
     main = "AIC for All Possible Subset Models",
     pch = 19, col = "darkgreen")
min_aic <- which.min(aic_vals)
points(min_aic, aic_vals[min_aic], col = "red", pch = 19, cex = 1.2)
text(min_aic, aic_vals[min_aic],
     labels = paste("Min AIC =", round(aic_vals[min_aic], 2)),
     pos = 3, col = "red")

# Plot BIC values
plot(bic_vals,
     xlab = "Model Index",
     ylab = "BIC",
     main = "BIC for All Possible Subset Models",
     pch = 19, col = "purple")
min_bic <- which.min(bic_vals)
points(min_bic, bic_vals[min_bic], col = "red", pch = 19, cex = 1.2)
text(min_bic, bic_vals[min_bic],
     labels = paste("Min BIC =", round(bic_vals[min_bic], 2)),
     pos = 3, col = "red")

# Display best models based on AIC and BIC
cat("\nBest Model (Minimum AIC):\n")
best_aic_vars <- pred_names[fit_summary$which[min_aic, -1, drop = FALSE]]
best_aic_vars <- best_aic_vars[best_aic_vars != "(Intercept)"]
cat("Predictors:", paste(best_aic_vars, collapse = ", "), "\n")
cat("AIC =", round(aic_vals[min_aic], 3), "\n\n")

cat("Best Model (Minimum BIC):\n")
best_bic_vars <- pred_names[fit_summary$which[min_bic, -1, drop = FALSE]]
best_bic_vars <- best_bic_vars[best_bic_vars != "(Intercept)"]
cat("Predictors:", paste(best_bic_vars, collapse = ", "), "\n")
cat("BIC =", round(bic_vals[min_bic], 3), "\n\n")

# logical inclusion vector for best models
cat("\nPredictor inclusion (AIC-min model):\n")
print(fit_summary$which[min_aic, ])

cat("\nPredictor inclusion (BIC-min model):\n")
print(fit_summary$which[min_bic, ])

# for model involact

# Run all possible subset regressions (APR, best model per subset size)
fit <- regsubsets(involact_bc ~ age + log_income + fire + theft + race + volact,
                  data = newdata, nbest = 1)

fit_summary <- summary(fit)

# Prepare to compute AIC and BIC for each model

num_models <- nrow(fit_summary$which)
num_pred <- apply(fit_summary$which[, -1, drop = FALSE], 1, sum)
pred_names <- colnames(fit_summary$which)[-1]

# Create vectors to store AIC and BIC
aic_vals <- numeric(num_models)
bic_vals <- numeric(num_models)

# Compute AIC and BIC for each model
for (i in 1:num_models) {
  sel_vars <- pred_names[fit_summary$which[i, -1, drop = FALSE]]
  sel_vars <- sel_vars[sel_vars != "(Intercept)"]
  
  # Fit the model
  formula_str <- paste("involact_bc ~", paste(sel_vars, collapse = " + "))
  model <- lm(as.formula(formula_str), data = newdata)
  
  # Store AIC and BIC
  aic_vals[i] <- AIC(model)
  bic_vals[i] <- BIC(model)
}

# Plot AIC values
plot(aic_vals,
     xlab = "Model Index",
     ylab = "AIC",
     main = "AIC for All Possible Subset Models",
     pch = 19, col = "darkgreen")
min_aic <- which.min(aic_vals)
points(min_aic, aic_vals[min_aic], col = "red", pch = 19, cex = 1.2)
text(min_aic, aic_vals[min_aic],
     labels = paste("Min AIC =", round(aic_vals[min_aic], 2)),
     pos = 3, col = "red")

# Plot BIC values
plot(bic_vals,
     xlab = "Model Index",
     ylab = "BIC",
     main = "BIC for All Possible Subset Models",
     pch = 19, col = "purple")
min_bic <- which.min(bic_vals)
points(min_bic, bic_vals[min_bic], col = "red", pch = 19, cex = 1.2)
text(min_bic, bic_vals[min_bic],
     labels = paste("Min BIC =", round(bic_vals[min_bic], 2)),
     pos = 3, col = "red")

# Display best models based on AIC and BIC
cat("\nBest Model (Minimum AIC):\n")
best_aic_vars <- pred_names[fit_summary$which[min_aic, -1, drop = FALSE]]
best_aic_vars <- best_aic_vars[best_aic_vars != "(Intercept)"]
cat("Predictors:", paste(best_aic_vars, collapse = ", "), "\n")
cat("AIC =", round(aic_vals[min_aic], 3), "\n\n")

cat("Best Model (Minimum BIC):\n")
best_bic_vars <- pred_names[fit_summary$which[min_bic, -1, drop = FALSE]]
best_bic_vars <- best_bic_vars[best_bic_vars != "(Intercept)"]
cat("Predictors:", paste(best_bic_vars, collapse = ", "), "\n")
cat("BIC =", round(bic_vals[min_bic], 3), "\n\n")

# logical inclusion vector for best models
cat("\nPredictor inclusion (AIC-min model):\n")
print(fit_summary$which[min_aic, ])

cat("\nPredictor inclusion (BIC-min model):\n")
print(fit_summary$which[min_bic, ])

#######################################################

# Stepwise Regression

# For model volact

stepwise_model <- stepAIC(model_volact_bc, direction = "both")

# For model involact
stepwise_model <- stepAIC(model_involact_bc, direction = "both")

########################################################

# Priciple Component Analysis

# select only predictors
predictors <- newdata[, c(
  "race", "fire", "theft",
  "age", "involact", "log_income"
)]
pca_model <- prcomp(predictors, scale. = TRUE)
pca_model$rotation
pca_scores <- pca_model$x

plot(pca_model, type = "l", main = "Scree Plot")

pca_df <- as.data.frame(pca_scores)
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "blue") +
  ggtitle("PCA - First Two Components")

newdata_pca <- data.frame(
  volact_bc = newdata$volact_bc,
  involact_bc = newdata$involact_bc, pca_scores[, 1:3]
)

model_volact_pca <- lm(volact_bc ~ PC1 + PC2 + PC3, data = newdata_pca)
model_involact_pca <- lm(involact_bc ~ PC1 + PC2 + PC3, data = newdata_pca)

# Check then normality, heteroscedasticity, serial correlation, multicollinearity
# using the previous code
# on model_volact_pca and model_involact_pca
#######################################################

# Fitted values vs Actual data plot

# lambdas
lambda_volact <- 0.22222
lambda_involact <- -0.42424

# inverse Box–Cox
inv_boxcox <- function(z, lambda)
  if (lambda == 0) exp(z) else (lambda * z + 1)^(1 / lambda)

# get PCs for original data using saved pca_model
pca_vars <- rownames(pca_model$rotation)
pc_data <- as.data.frame(predict(pca_model, newdata = data[, pca_vars]))
pred_df <- pc_data[, c("PC1", "PC2", "PC3")]

# predict on Box–Cox scale
pred_volact_bc <- predict(model_volact_pca, newdata = pred_df)
pred_involact_bc <- predict(model_involact_pca, newdata = pred_df)

# back-transform
pred_volact <- inv_boxcox(pred_volact_bc, lambda_volact)
pred_involact <- (inv_boxcox(pred_involact_bc, lambda_involact) - 1)

# plot for volact
index <- seq_len(nrow(data))
plot(index, data$volact, type="l", col="black", main="Volact: Actual vs Predicted",
     xlab="Index", ylab="Volact")
lines(index, pred_volact, col="red", lwd=2)
legend("topleft", legend=c("Actual","Predicted"), col=c("black","red"), lty=1, bty="n")

# plot for involact
plot(index, data$involact, type="l", col="black", main="Involact: Actual vs Predicted",
     xlab="Index", ylab="Involact")
lines(index, pred_involact, col="red", lwd=2)
legend("topleft", legend=c("Actual","Predicted"), col=c("black","red"), lty=1, bty="n")

