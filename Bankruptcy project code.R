# installing required packages
install.packages("readr")
install.packages("GGally")
install.packages("GPArotation")
install.packages("ICSNP")
install.packages("nortest")
install.packages("candisc")
install.packages("bestNormalize")
install.packages("e1071")
install.packages("corrplot")
install.packages("caret")
install.packages("pROC")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("psych")
install.packages("biotools")
install.packages("MVN")
install.packages("car")
install.packages("MASS")
install.packages("dplyr")
install.packages("Hotelling")

# reqiured packages
library(readr)
library(GGally)
library(dplyr)
library(MASS)
library(car)
library(MVN)
library(biotools)
library(psych)
library(ggplot2)
library(reshape2)
library(pROC)
library(caret)
library(corrplot)
library(e1071)
library(bestNormalize)
library(candisc)
library(nortest)
library(ICSNP)
library(GPArotation)
library(Hotelling)

# data processing
data <- read.table("bankruptcy.txt")
colnames(data) <- c("x1", "x2", "x3", "x4", "y")
head(data)
print(sum(is.na(data))) # no missing values
data$y <- as.factor(data$y)
df_vars <- data[, 1:4]
# -------------------------------------------------------------------
# 1. EXPLORATORY DATA ANALYSIS (EDA)
# -------------------------------------------------------------------
# pairplot
par(mfrow = c(2, 2))
stripchart(x1 ~ y, data = data, method = "jitter", pch = 19,
           main = "x1 vs y", xlab = "y", ylab = "x1", col = "blue")

stripchart(x2 ~ y, data = data, method = "jitter", pch = 19,
           main = "x2 vs y", xlab = "y", ylab = "x2", col = "blue")

stripchart(x3 ~ y, data = data, method = "jitter", pch = 19,
           main = "x3 vs y", xlab = "y", ylab = "x3", col = "blue")

stripchart(x4 ~ y, data = data, method = "jitter", pch = 19,
           main = "x4 vs y", xlab = "y", ylab = "x4", col = "blue")
ggpairs(data, columns = 1:4, mapping = aes(color = y, alpha = 0.5))
pairs(data[, 1:4], col = data$y, main = "Scatterplot Matrix by Category")
# boxplot
for(i in 1:4) {
  boxplot(df_vars[,i] ~ data$y, main = colnames(df_vars)[i], 
          xlab = "Group", ylab = "", col = c("lightblue", "lightgreen"))
}
# histograms groupwise
par(mfrow = c(2,2))
for(i in 1:4){
  
  # Group 0
  x0 <- data[data$y == 0, i]
  
  hist(x0,
       probability = TRUE,
       main = paste("X", i, "- Group 0"),
       xlab = paste("X", i),
       col = "lightgreen",
       border = "black")
  
  lines(density(x0), col = "blue", lwd = 2.5)
  
  
  # Group 1
  x1 <- data[data$y == 1, i]
  
  hist(x1,
       probability = TRUE,
       main = paste("X", i, "- Group 1"),
       xlab = paste("X", i),
       col = "salmon",
       border = "black")
  
  lines(density(x1), col = "blue", lwd = 2.5)
}
# histograms for each covariate
for(i in 1:4) {
  hist(df_vars[,i], main = paste("Histogram of", colnames(df_vars)[i]), 
       xlab = colnames(df_vars)[i], col = "steelblue", breaks = 10)
}
# Correlation heatmap
cor_matrix <- cor(df_vars)
corrplot(cor_matrix, method = "color", addCoef.col = "black", 
         tl.col = "black", tl.srt = 45, main = "Correlation Heatmap")
# corr heatmap with the group variable
data$y <- as.numeric(data$y)
cor_mat <- cor(data[,1:5])
print(cor_mat)
# Convert matrix to long format
cor_df <- melt(cor_mat)
# Plot heatmap
ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2(low = "cyan", high = "red", mid = "yellow",
                       midpoint = 0, limit = c(-1,1)) +
  theme_minimal() +
  labs(title = "Correlation Heatmap",
       x = "", y = "")
# -------------------------------------------------------------------
# 2. CHECKING SKEWNESS
# -------------------------------------------------------------------
# skewness -> -1 → highly left-skewed, 1 → highly skewed, ~0 → symmetric
apply(df_vars, 2, skewness)
# -------------------------------------------------------------------
# 3. BOX'S M TEST (Homogeneity of Covariance Matrices)
# -------------------------------------------------------------------
# Box's M test for -> H0(null hypothesis) : Homoscedasticity
box_m_result <- boxM(df_vars, data$y)
print(box_m_result) # test result : not homoscedastic
plot(box_m_result) # plot the box M test result
# -------------------------------------------------------------------
# 4. NORMALITY TESTS
# -------------------------------------------------------------------
# checking normality
#QQ plot
par(mfrow=c(1,1))
for(i in 1:4){
  qqnorm(df_vars[,i], main = paste("QQ Plot X", i))
  qqline(df_vars[,i], col = "red", lwd = 2)
}
# Statistical tests for each variable
for(i in 1:4){
  var_name <- colnames(df_vars)[i]
  data_col <- df_vars[,i]
  
  cat(paste("\nVariable:", var_name, "\n"))
  
  # Shapiro-Wilk
  print(shapiro.test(data_col))
  
  # Anderson-Darling
  print(ad.test(data_col))
  
  # Kolmogorov-Smirnov (Comparing against a normal distribution with the sample mean/sd)
  print(ks.test(data_col, "pnorm", mean(data_col), sd(data_col)))
}
# Multivariate normality tests
mardia_res  <- MVN::mvn(data = df_vars, mvn_test = "mardia")
hz_res      <- MVN::mvn(data = df_vars, mvn_test = "hz")
royston_res <- MVN::mvn(data = df_vars, mvn_test = "royston")
# Printing the results of the Multivariate normality tests
mardia_res$multivariate_normality
hz_res$multivariate_normality
royston_res$multivariate_normality

# -------------------------------------------------------------------
# 5. VARIANCE INFLATION FACTOR (VIF)
# -------------------------------------------------------------------
# We fit a logistic regression model to check VIF among predictors
vif_model <- glm(y ~ x1 + x2 + x3 + x4, data = data, family = binomial)
print(vif(vif_model))
# -------------------------------------------------------------------
# 8. GROUP DIFFERENCES (MANOVA & Hotelling's T^2)
# -------------------------------------------------------------------
# MANOVA model
manova_model <- manova(cbind(x1, x2, x3, x4) ~ y, data = data)
summary(manova_model, test = "Wilks")
summary(manova_model, test = "Pillai")
summary(manova_model, test = "Hotelling-Lawley")
summary(manova_model, test = "Roy")
# group data are very much different from each other
# visualize MANOVA results
plot(candisc(manova_model))

# Hotelling T^2 test
# Split into the two groups
group0 <- data[data$y == 0, c("x1", "x2", "x3", "x4")]
group1 <- data[data$y == 1, c("x1", "x2", "x3", "x4")]
res <- hotelling.test(group0, group1)
print(res)

# ---------------------------------------------------
# Plotting the models using x1, x3 keeping x2, x4 fixed at their mean
# ---------------------------------------------------
#  Fit models on full data
lda_fit <- lda(y ~ x1 + x2 + x3 + x4, data = data)
qda_fit <- qda(y ~ x1 + x2 + x3 + x4, data = data)
log_fit <- glm(y ~ x1 + x2 + x3 + x4, data = data, family = binomial)
## -----------------------------
## 3. Create a grid for plotting
##    We plot over x1 and x3,
##    holding x2 and x4 fixed
## -----------------------------
x1_seq <- seq(min(data$x1) - 0.1, max(data$x1) + 0.1, length.out = 250)
x3_seq <- seq(min(data$x3) - 0.1, max(data$x3) + 0.1, length.out = 250)

grid <- expand.grid(
  x1 = x1_seq,
  x3 = x3_seq
)

grid$x2 <- mean(data$x2)
grid$x4 <- mean(data$x4)

## -----------------------------
## 4. Predict on grid
## -----------------------------
# LDA
lda_class <- predict(lda_fit, newdata = grid)$class
lda_prob1 <- predict(lda_fit, newdata = grid)$posterior[, "1"]

# QDA
qda_class <- predict(qda_fit, newdata = grid)$class
qda_prob1 <- predict(qda_fit, newdata = grid)$posterior[, "1"]

# Logistic regression
log_prob1 <- predict(log_fit, newdata = grid, type = "response")
log_class <- ifelse(log_prob1 > 0.5, "1", "0")

## Convert to matrices for contour/image plotting
zmat_lda_class <- matrix(as.numeric(lda_class), nrow = length(x1_seq), ncol = length(x3_seq))
zmat_lda_prob  <- matrix(lda_prob1, nrow = length(x1_seq), ncol = length(x3_seq))

zmat_qda_class <- matrix(as.numeric(qda_class), nrow = length(x1_seq), ncol = length(x3_seq))
zmat_qda_prob  <- matrix(qda_prob1, nrow = length(x1_seq), ncol = length(x3_seq))

zmat_log_class <- matrix(as.numeric(factor(log_class, levels = c("0","1"))), nrow = length(x1_seq), ncol = length(x3_seq))
zmat_log_prob  <- matrix(log_prob1, nrow = length(x1_seq), ncol = length(x3_seq))

## -----------------------------
## 5. Plot function
## -----------------------------
plot_region <- function(zclass, zprob, title_text, boundary_col = "black") {
  image(
    x1_seq, x3_seq, zclass,
    col = c("#f4a3a3", "#9ecae1"),
    xlab = "x1", ylab = "x3",
    main = title_text
  )
  
  points(
    data$x1, data$x3,
    pch = 19,
    col = ifelse(data$y == "0", "white", "blue4")
  )
  
  # Decision boundary
  contour(
    x1_seq, x3_seq, zprob,
    levels = 0.5,
    add = TRUE,
    drawlabels = FALSE,
    lwd = 2,
    col = boundary_col
  )
  
  legend(
    "topright",
    legend = c("Group 0", "Group 1", "Boundary"),
    col = c("white", "blue4", boundary_col),
    pch = c(19, 19, NA),
    lty = c(NA, NA, 1),
    lwd = c(NA, NA, 2),
    bty = "n"
  )
}

## -----------------------------
## 6. Plot all three models
## -----------------------------
par(mfrow = c(1, 3))

plot_region(
  zclass = zmat_lda_class,
  zprob  = zmat_lda_prob,
  title_text = "LDA"
)

plot_region(
  zclass = zmat_qda_class,
  zprob  = zmat_qda_prob,
  title_text = "QDA"
)

plot_region(
  zclass = zmat_log_class,
  zprob  = zmat_log_prob,
  title_text = "Logistic Regression"
)

# -------------------------------------------------------------------
# 9. CLASSIFICATION & CROSS VALIDATION
# -------------------------------------------------------------------
set.seed(123)
# Define cross-validation method
train_control <- trainControl(method = "cv", number = 5)
# CV for logistic model
model_logit_cv <- train(
  y ~ x1 + x2 + x3 + x4,
  data = data,
  method = "glm",
  family = "binomial",
  trControl = train_control
)
model_logit_cv
model_logit_cv$results$Accuracy

# LDA method
model_lda_cv <- train(
  y ~ x1 + x2 + x3 + x4,
  data = data,
  method = "lda",
  trControl = train_control
)
model_lda_cv
model_lda_cv$results$Accuracy

# QDA method
model_qda_cv <- train(
  y ~ x1 + x2 + x3 + x4,
  data = data,
  method = "qda",
  trControl = train_control
)
model_qda_cv
model_qda_cv$results$Accuracy

# comparing between the models
results <- data.frame(
  Model = c("Logistic", "LDA", "QDA"),
  Accuracy = c(
    model_logit_cv$results$Accuracy,
    model_lda_cv$results$Accuracy,
    model_qda_cv$results$Accuracy
  )
)

print(results)

# Compare Models
results <- resamples(list(Logistic = model_logit_cv, LDA = model_lda_cv, QDA = model_qda_cv))
cat("\nModel Comparison (Accuracy & Kappa):\n")
summary(results)

# Plot model comparisons
bwplot(results)

# Leave 10 out method for 2000 iterations
# 1. Setup Parameters
n_iterations <- 2000  # Number of random subsets to test
n_leave_out <- 10     # Leave 10 out for the test set
models <- c("Logistic", "LDA", "QDA")
metrics <- c("Accuracy", "Precision", "Recall", "F1")
# Store results
results <- data.frame(
  Iteration = integer(),
  Model = character(),
  Accuracy = numeric(),
  Precision = numeric(),
  Recall = numeric(),
  F1 = numeric(),
  stringsAsFactors = FALSE
)
# Function to compute metrics
calc_metrics <- function(actual, pred, positive_class = "1") {
  actual <- as.character(actual)
  pred <- as.character(pred)
  
  TP <- sum(pred == positive_class & actual == positive_class)
  TN <- sum(pred != positive_class & actual != positive_class)
  FP <- sum(pred == positive_class & actual != positive_class)
  FN <- sum(pred != positive_class & actual == positive_class)
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  recall <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
               NA,
               2 * precision * recall / (precision + recall))
  
  c(Accuracy = accuracy, Precision = precision, Recall = recall, F1 = f1)
}
# Vectors to store the accuracy of each iteration
#acc_log <- numeric(n_iterations)
#acc_lda <- numeric(n_iterations)
#acc_qda <- numeric(n_iterations)

# Set seed for reproducibility
set.seed(123)

# 2. The Cross-Validation Loop
for (i in 1:n_iterations) {
  
  # Randomly sample 10 row indices to be the test set
  test_indices <- sample(1:nrow(data), n_leave_out)
  
  # Split the data
  train_data <- data[-test_indices, ]
  test_data <- data[test_indices, ]
  
  # -----------------------------------------
  # Model 1: Logistic Regression
  # -----------------------------------------
  # suppressWarnings is used here because small sample sizes (n=36) 
  # can sometimes cause "perfect separation" warnings in logistic regression
  log_model <- suppressWarnings(glm(y ~ ., data = train_data, family = "binomial"))
  log_prob <- predict(log_model, test_data, type = "response")
  
  # Convert probabilities to class predictions (0 or 1)
  log_pred <- ifelse(log_prob >= 0.5, 1, 0)
  m <- calc_metrics(test_data$y, log_pred)
  results <- rbind(results, data.frame(
    Iteration = i, Model = "Logistic",
    Accuracy = m["Accuracy"], Precision = m["Precision"],
    Recall = m["Recall"], F1 = m["F1"]
  ))
  #acc_log[i] <- mean(log_pred == test_data$y)
  
  # -----------------------------------------
  # Model 2: Linear Discriminant Analysis (LDA)
  # -----------------------------------------
  p1 <- mean(train_data$y == "0")
  p2 <- mean(train_data$y == "1")
  lda_model <- lda(y ~ ., data = train_data, prior = c(p1, p2))
  lda_pred <- predict(lda_model, test_data)$class
  m <- calc_metrics(test_data$y, lda_pred)
  results <- rbind(results, data.frame(
    Iteration = i, Model = "LDA",
    Accuracy = m["Accuracy"], Precision = m["Precision"],
    Recall = m["Recall"], F1 = m["F1"]
  ))
  #acc_lda[i] <- mean(lda_pred == test_data$y)
  
  # -----------------------------------------
  # Model 3: Quadratic Discriminant Analysis (QDA)
  # -----------------------------------------
  qda_model <- qda(y ~ ., data = train_data)
  qda_pred <- predict(qda_model, test_data)$class
  m <- calc_metrics(test_data$y, qda_pred)
  results <- rbind(results, data.frame(
    Iteration = i, Model = "QDA",
    Accuracy = m["Accuracy"], Precision = m["Precision"],
    Recall = m["Recall"], F1 = m["F1"]
  ))
  #acc_qda[i] <- mean(qda_pred == test_data$y)
}
# Average metrics by model
avg_results <- aggregate(cbind(Accuracy, Precision, Recall, F1) ~ Model,
                         data = results,
                         FUN = function(x) mean(x, na.rm = TRUE))

cat("\nAverage metrics over 2000 iterations:\n")
print(avg_results)
# Plot the values
par(mfrow = c(2, 2))

boxplot(Accuracy ~ Model, data = results,
        main = "Accuracy", col = "lightblue", ylim = c(0, 1))

boxplot(Precision ~ Model, data = results,
        main = "Precision", col = "lightgreen", ylim = c(0, 1))

boxplot(Recall ~ Model, data = results,
        main = "Recall", col = "lightpink", ylim = c(0, 1))

boxplot(F1 ~ Model, data = results,
        main = "F1 Score", col = "yellow", ylim = c(0, 1))

# ---------------------------------------------------------
# Transformation
# ---------------------------------------------------------
min_val <- min(data$x1, na.rm = TRUE)
data$x1 <- if(min_val <= 0) data$x1 + abs(min_val) else data$x1
min_val <- min(data$x2, na.rm = TRUE)
x2_shifted <- if(min_val <= 0) data$x2 + abs(min_val) + 1 else data$x2
# log transformations
x2_log <- log(x2_shifted)
x3_log <- log(data$x3)
# sqrt transformations
x2_sqrt <- sqrt(x2_shifted)
x3_sqrt <- sqrt(data$x3)
# boxcox transformations
# For X2
bc_mod <- powerTransform(x2_shifted ~ 1, family="bcPower")
x2_bc <- bcPower(x2_shifted, lambda = bc_mod$lambda)
lambda_opt <- bc_mod$lambda
# For X3
bc_mod <- powerTransform(data$x3 ~ 1, family="bcPower")
x3_bc <- bcPower(data$x3, lambda = bc_mod$lambda)
lambda_opt <- bc_mod$lambda

# Yeo-Johnson Transformation
# X2
yj_x2 <- yeojohnson(data$x2)
x2_yj <- yj_x2$x.t
lambda_x2 <- yj_x2$lambda
# X3
yj_x3 <- yeojohnson(data$x3)
x3_yj <- yj_x3$x.t
lambda_x3 <- yj_x3$lambda
df_transformed <- data.frame(
  Log_X2 = x2_log, 
  Log_X3 = x3_log,
  Sqrt_X2 = x2_sqrt, 
  Sqrt_X3 = x3_sqrt,
  BoxCox_X2 = x2_bc, 
  BoxCox_X3 = x3_bc,
  YJ_X2 = x2_yj, 
  YJ_X3 = x3_yj
)
# Loop through the data frame by column index
for (i in 1:ncol(df_transformed)) {
  
  # Extract the current column data
  current_data <- df_transformed[, i]
  
  # Extract the column name for the title
  current_name <- colnames(df_transformed)[i]
  
  # Generate the QQ plot
  qqnorm(current_data, 
         main = paste("QQ Plot:", current_name), 
         pch = 19, 
         col = "darkblue")
  
  # Add the normal reference line
  qqline(current_data, col = "red", lwd = 2)
}
# Statistical tests for each variable
for(i in 1:ncol(df_transformed)) {
  var_name <- colnames(df_transformed)[i]
  data_col <- df_transformed[,i]
  
  cat(paste("\nVariable:", var_name, "\n"))
  
  # Shapiro-Wilk
  print(shapiro.test(data_col))
  
  # Anderson-Darling
  print(ad.test(data_col))
  
  # Kolmogorov-Smirnov (Comparing against a normal distribution with the sample mean/sd)
  print(ks.test(data_col, "pnorm", mean(data_col), sd(data_col)))
}
# Multivariate normality check
data$x2 <- x2_bc
data$x3 <- x3_bc
df_vars <- data[, 1:4]
# -------------------------------------------------------------------
# 6. CHI-SQUARE QQ PLOT (Mahalanobis Distance for Outliers)
# -------------------------------------------------------------------
# Chi square QQ plot
md <- mahalanobis(df_vars, center = colMeans(df_vars), cov = cov(df_vars))
p <- ncol(df_vars)
x_theory <- qchisq(ppoints(nrow(df_vars)), df = p)
y_sample <- sort(md)
ord <- order(md)
qqplot(
  x_theory,
  y_sample,
  main = expression("Chi-square Q-Q Plot"),
  xlab = expression("Theoretical " * chi^2 * " Quantiles"),
  ylab = "Ordered Mahalanobis Distances",
  pch = 19,
  col = "darkblue"
)
abline(0, 1, col = "green", lwd = 2)
cutoff <- qchisq(0.975, df = p) # alpha = 0.025
abline(h = cutoff, col = "red", lty = 2)
text(x_theory, y_sample, labels = ord, pos = 4, cex = 0.7)
outliers <- which(md > cutoff)
outliers
data[outliers, ]
# read the data again
# deleting the outlier
data <- data[-46, ]
# Do the boxcox on x2, x3
# Do the multivariate normality test
mardia_res  <- MVN::mvn(data = df_vars, mvn_test = "mardia")
hz_res      <- MVN::mvn(data = df_vars, mvn_test = "hz")
royston_res <- MVN::mvn(data = df_vars, mvn_test = "royston")
# Printing the results of the Multivariate normality tests
mardia_res$multivariate_normality
hz_res$multivariate_normality
royston_res$multivariate_normality


# 1. Perform PCA on the 4 continuous variables
# center = TRUE and scale. = TRUE standardizes the variables (mean = 0, variance = 1)
pca_result <- prcomp(df_vars, center = TRUE, scale. = TRUE)
print(summary(pca_result))
print(pca_result$rotation)
# 2. Extract Variance for the Scree Plot
# Calculate the variance explained by each principal component
pca_var <- pca_result$sdev^2

# Calculate the proportion of variance explained
prop_var <- pca_var / sum(pca_var)
cumulative_var <- cumsum(prop_var)

# 3. Create the Scree Plot
# Plot the individual proportion of variance explained
plot(1:4, prop_var, type = "b", 
     xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained",
     main = "PCA Scree Plot",
     ylim = c(0, 1),
     pch = 19, col = "steelblue", lwd = 2, xaxt = "n")

# Add the X-axis labels properly (PC1, PC2, etc.)
axis(1, at = 1:4, labels = paste0("PC", 1:4))

# Add a line showing the cumulative variance explained
lines(1:4, cumulative_var, type = "b", pch = 17, col = "darkorange", lwd = 2)

# Add text labels to the cumulative line to show exact percentages
text(1:4, cumulative_var, labels = paste0(round(cumulative_var * 100, 1), "%"), 
     pos = 3, col = "darkorange", cex = 0.8)

# Add a legend
legend("right", legend = c("Individual Variance", "Cumulative Variance"), 
       col = c("steelblue", "darkorange"), pch = c(19, 17), lwd = 2, bty = "n")

# Normality check for 2 PCs
pca_result <- prcomp(df_vars, center = TRUE, scale. = TRUE)
df_vars <- as.data.frame(predict(pca_result, newdata = df_vars)[, 1:2])
colnames(df_vars) <- c("PC1", "PC2")
df_vars$y <- data$y
df0 <- df_vars[df_vars$y == 0, 1:2]
df1 <- df_vars[df_vars$y == 1, 1:2]
# univariate and multivariate normality are satisfied
#------------------------------------
# Factoring
#------------------------------------
# --------------------------------------------------
# 2. KMO test for sampling adequacy
# --------------------------------------------------
kmo_res <- KMO(cor(df_vars))
print(kmo_res)

# Interpretation:
# KMO >= 0.90 : marvelous
# KMO >= 0.80 : meritorious
# KMO >= 0.70 : middling
# KMO >= 0.60 : mediocre
# KMO >= 0.50 : miserable
# KMO <  0.50 : unacceptable

if (kmo_res$MSA >= 0.5) {
  cat("\nConclusion: The data is adequate for factor analysis.\n")
} else {
  cat("\nConclusion: The data is NOT adequate for factor analysis.\n")
}
# --------------------------------------------------
# 3. Factor analysis with m = 1
# --------------------------------------------------
cat("\n============================\n")
cat("Factor Analysis: m = 1\n")
cat("============================\n")

fa_1_none <- factanal(df_vars, factors = 1, rotation = "none", scores = "regression")
print(fa_1_none, digits = 3, cutoff = 0)

# Note:
# For 1 factor, rotation is not meaningful, so only "none" is used.

# --------------------------------------------------
# 4. Factor analysis with m = 2
# --------------------------------------------------
cat("\n============================\n")
cat("Factor Analysis: m = 2, Varimax\n")
cat("============================\n")
fa_2_varimax <- factanal(df_vars, factors = 2, rotation = "varimax", scores = "regression")
print(fa_2_varimax, digits = 3, cutoff = 0)
kmo_res$MSA


# -----------------------------
# Naive Bayes classification
# -----------------------------
nb_model <- naiveBayes(y ~ ., data = data)

# Class predictions
nb_pred <- predict(nb_model, newdata = data)

# Posterior probabilities
nb_prob <- predict(nb_model, newdata = data, type = "raw")

# Confusion matrix
table(True = data$y, Predicted = nb_pred)

# Accuracy
mean(nb_pred == data$y)

# See probabilities
nb_pred

# -----------------------------
# 2000 iterations
# -----------------------------
set.seed(123)
n_iter <- 2000
test_size <- 10

all_metrics <- matrix(NA, nrow = n_iter, ncol = 4)
colnames(all_metrics) <- c("Accuracy", "Precision", "Recall", "F1")

for (i in 1:n_iter) {
  test_idx <- sample(1:nrow(data), test_size)
  train_data <- data[-test_idx, ]
  test_data  <- data[test_idx, ]
  
  nb_model <- naiveBayes(y ~ ., data = train_data)
  pred <- predict(nb_model, newdata = test_data)
  
  all_metrics[i, ] <- calc_metrics(test_data$y, pred, positive_class = "1")
}

colMeans(all_metrics, na.rm = TRUE)

# --------------------------
# Clustering
# --------------------------
x <- data[, c("x1", "x2", "x3", "x4")]
## Usually better to scale before K-means
x_scaled <- scale(x)
## PCA only for plotting
pca <- prcomp(x_scaled)
plot_df <- as.data.frame(pca$x[, 1:2])
names(plot_df) <- c("PC1", "PC2")
# -----------------------------
# Initial plot of observations
# -----------------------------
plot(plot_df$PC1, plot_df$PC2,
     pch = 19, col = "black",
     xlab = "PC1", ylab = "PC2",
     main = "Original observations")

text(plot_df$PC1, plot_df$PC2, labels = 1:nrow(plot_df), pos = 3, cex = 0.7)
# -----------------------------
# Functions for manual K-means
# -----------------------------
plot_clusters <- function(plot_df, clusters = NULL, centers = NULL, main_txt = "") {
  if (is.null(clusters)) {
    cols <- "black"
  } else {
    cols <- c("tomato", "steelblue")[clusters]
  }
  
  plot(plot_df$PC1, plot_df$PC2,
       pch = 19, col = cols,
       xlab = "PC1", ylab = "PC2",
       main = main_txt)
  
  text(plot_df$PC1, plot_df$PC2, labels = 1:nrow(plot_df), pos = 3, cex = 0.7)
  
  if (!is.null(centers)) {
    centers_pc <- predict(pca, newdata = centers)[, 1:2]
    points(centers_pc[,1], centers_pc[,2],
           pch = 8, cex = 2, lwd = 2, col = c("darkred", "darkblue"))
  }
}
assign_clusters <- function(x, centers) {
  x <- as.matrix(x)
  centers <- as.matrix(centers)
  d <- sapply(1:nrow(centers), function(k) {
    rowSums((x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2)
  })
  max.col(-d)
}
update_centers <- function(x, clusters, k) {
  x <- as.matrix(x)
  centers <- matrix(NA, nrow = k, ncol = ncol(x))
  for (j in 1:k) {
    centers[j, ] <- colMeans(x[clusters == j, , drop = FALSE])
  }
  colnames(centers) <- colnames(x)
  as.data.frame(centers)
}
# -----------------------------
# Manual K-means step by step
# -----------------------------
set.seed(123)
k <- 2
# choose 2 random initial centers from the observations
init_idx <- sample(1:nrow(x_scaled), k)
centers <- as.data.frame(x_scaled[init_idx, ])
plot_clusters(plot_df, centers = centers,
              main_txt = "Step 0: Initial centers")
max_iter <- 20
old_clusters <- NULL
for (iter in 1:max_iter) {
  
  ## Step A: assign each point to nearest center
  clusters <- assign_clusters(x_scaled, centers)
  
  plot_clusters(plot_df, clusters = clusters, centers = centers,
                main_txt = paste("Iteration", iter, "- After assignment"))
  
  ## Step B: recompute centers
  new_centers <- update_centers(x_scaled, clusters, k)
  
  plot_clusters(plot_df, clusters = clusters, centers = new_centers,
                main_txt = paste("Iteration", iter, "- After center update"))
  
  ## Stop if cluster assignments do not change
  if (!is.null(old_clusters) && all(clusters == old_clusters)) {
    cat("K-means converged at iteration", iter, "\n")
    break
  }
  
  centers <- new_centers
  old_clusters <- clusters
}

## final cluster labels
final_clusters <- clusters

# -----------------------------
# Compare clusters with true class y
# -----------------------------
# contingency table
tab <- table(Cluster = final_clusters, TrueClass = data$y)
print(tab)
predicted_class <- final_clusters %% 2
correct <- predicted_class == data$y
# plot
plot(plot_df$PC1, plot_df$PC2,
     col = ifelse(correct, "darkgreen", "red"),
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "Clustering vs true class")
text(plot_df$PC1, plot_df$PC2, labels = 1:nrow(plot_df), pos = 3, cex = 0.7)

legend("topright",
       legend = c("Correct", "Incorrect"),
       col = c("darkgreen", "red"),
       pch = 19)

x <- data[, c("x1", "x2", "x3", "x4")]
x_scaled <- scale(x)

# class-based initial centers
center0 <- colMeans(x_scaled[data$y == 0, ])
center1 <- colMeans(x_scaled[data$y == 1, ])
init_centers <- rbind(center0, center1)
km <- kmeans(x_scaled, centers = init_centers)
table(Cluster = km$cluster, TrueClass = data$y) # better
predicted_class <- km$cluster - 1
correct <- predicted_class == data$y
# plot
plot(plot_df$PC1, plot_df$PC2,
     col = ifelse(correct, "darkgreen", "red"),
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "Clustering vs true class")
text(plot_df$PC1, plot_df$PC2, labels = 1:nrow(plot_df), pos = 3, cex = 0.7)

legend("topright",
       legend = c("Correct", "Incorrect"),
       col = c("darkgreen", "red"),
       pch = 19)

# -------------------------------
# LDA + KMeans
# -------------------------------
lda_fit <- lda(y ~ x1 + x2 + x3 + x4, data = data)
lda_scores <- predict(lda_fit)$x
lda_df <- data.frame(LD1 = lda_scores[,1], y = data$y)
# Plot observations on LDA axis
plot(lda_df$LD1, rep(0, nrow(lda_df)),
     col = c("tomato", "steelblue")[lda_df$y],
     pch = 19,
     xlab = "LD1",
     ylab = "",
     main = "Data projected onto LDA direction",
     yaxt = "n")
text(lda_df$LD1, rep(0, nrow(lda_df)),
     labels = 1:nrow(lda_df), pos = 3, cex = 0.7)
legend("topright", legend = c("Class 0", "Class 1"),
       col = c("tomato", "steelblue"), pch = 19)
set.seed(123)
km_lda <- kmeans(lda_df["LD1"], centers = 2, nstart = 25)
lda_df$cluster <- km_lda$cluster
## Plot clusters on LDA axis
plot(lda_df$LD1, rep(0, nrow(lda_df)),
     col = c("darkgreen", "purple")[lda_df$cluster],
     pch = 19,
     xlab = "LD1",
     ylab = "",
     main = "K-means clustering on LDA scores",
     yaxt = "n")
text(lda_df$LD1, rep(0, nrow(lda_df)),
     labels = 1:nrow(lda_df), pos = 3, cex = 0.7)
tab <- table(Cluster = lda_df$cluster, TrueClass = lda_df$y)
print(tab)
predicted_class <- km_lda$cluster - 1
correct <- predicted_class == data$y
# plot
plot(plot_df$PC1, plot_df$PC2,
     col = ifelse(correct, "darkgreen", "red"),
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "Clustering vs true class")
text(plot_df$PC1, plot_df$PC2, labels = 1:nrow(plot_df), pos = 3, cex = 0.7)

legend("topright",
       legend = c("Correct", "Incorrect"),
       col = c("darkgreen", "red"),
       pch = 19)

# -------------------------
# Checking Accuracy
# -------------------------
# build confusion matrix
# -----------------------------------
get_conf_mat <- function(true_y, pred_y, y_levels) {
  table(
    Predicted = factor(pred_y, levels = y_levels),
    TrueClass = factor(true_y, levels = y_levels)
  )
}
lda_kmeans_result <- function(train_data, test_data, y_levels) {
  lda_fit <- lda(y ~ x1 + x2 + x3 + x4, data = train_data)
  
  train_ld <- predict(lda_fit, newdata = train_data)$x[, 1]
  test_ld  <- predict(lda_fit, newdata = test_data)$x[, 1]
  
  train_ld <- matrix(train_ld, ncol = 1)
  test_ld  <- matrix(test_ld,  ncol = 1)
  
  km <- kmeans(train_ld, centers = 2, nstart = 25)
  train_cluster <- km$cluster
  
  cluster_to_class <- rep(NA, 2)
  for (j in 1:2) {
    y_in_cluster <- train_data$y[train_cluster == j]
    cluster_to_class[j] <- names(which.max(table(y_in_cluster)))
  }
  
  centers <- km$centers[, 1]
  test_cluster <- sapply(test_ld[, 1], function(z) {
    which.min((z - centers)^2)
  })
  
  pred <- factor(cluster_to_class[test_cluster], levels = y_levels)
  conf_mat <- get_conf_mat(test_data$y, pred, y_levels)
  
  list(
    accuracy = mean(pred == test_data$y),
    conf_mat = conf_mat
  )
}
kmeans_result <- function(train_data, test_data, y_levels) {
  x_train <- train_data[, c("x1", "x2", "x3", "x4")]
  x_test  <- test_data[, c("x1", "x2", "x3", "x4")]
  
  train_center <- colMeans(x_train)
  train_scale  <- apply(x_train, 2, sd)
  
  x_train_scaled <- scale(x_train, center = train_center, scale = train_scale)
  x_test_scaled  <- scale(x_test,  center = train_center, scale = train_scale)
  
  km <- kmeans(x_train_scaled, centers = 2, nstart = 25)
  train_cluster <- km$cluster
  
  cluster_to_class <- rep(NA, 2)
  for (j in 1:2) {
    y_in_cluster <- train_data$y[train_cluster == j]
    cluster_to_class[j] <- names(which.max(table(y_in_cluster)))
  }
  
  centers <- km$centers
  test_cluster <- apply(x_test_scaled, 1, function(z) {
    dists <- apply(centers, 1, function(center) sum((z - center)^2))
    which.min(dists)
  })
  
  pred <- factor(cluster_to_class[test_cluster], levels = y_levels)
  conf_mat <- get_conf_mat(test_data$y, pred, y_levels)
  
  list(
    accuracy = mean(pred == test_data$y),
    conf_mat = conf_mat
  )
}
kmeans_groupmean_result <- function(train_data, test_data, y_levels) {
  x_train <- train_data[, c("x1", "x2", "x3", "x4")]
  x_test  <- test_data[, c("x1", "x2", "x3", "x4")]
  
  train_center <- colMeans(x_train)
  train_scale  <- apply(x_train, 2, sd)
  
  x_train_scaled <- scale(x_train, center = train_center, scale = train_scale)
  x_test_scaled  <- scale(x_test,  center = train_center, scale = train_scale)
  
  center0 <- colMeans(x_train_scaled[train_data$y == y_levels[1], , drop = FALSE])
  center1 <- colMeans(x_train_scaled[train_data$y == y_levels[2], , drop = FALSE])
  init_centers <- rbind(center0, center1)
  
  km <- kmeans(x_train_scaled, centers = init_centers)
  train_cluster <- km$cluster
  
  cluster_to_class <- rep(NA, 2)
  for (j in 1:2) {
    y_in_cluster <- train_data$y[train_cluster == j]
    cluster_to_class[j] <- names(which.max(table(y_in_cluster)))
  }
  
  centers <- km$centers
  test_cluster <- apply(x_test_scaled, 1, function(z) {
    dists <- apply(centers, 1, function(center) sum((z - center)^2))
    which.min(dists)
  })
  
  pred <- factor(cluster_to_class[test_cluster], levels = y_levels)
  conf_mat <- get_conf_mat(test_data$y, pred, y_levels)
  
  list(
    accuracy = mean(pred == test_data$y),
    conf_mat = conf_mat
  )
}
set.seed(123)
n_iter <- 2000
y_levels <- levels(data$y)

acc_lda <- numeric(n_iter)
acc_kmeans <- numeric(n_iter)
acc_kmeans_groupmean <- numeric(n_iter)

sum_conf_lda <- matrix(0, nrow = 2, ncol = 2,
                       dimnames = list(Predicted = y_levels, TrueClass = y_levels))
sum_conf_kmeans <- matrix(0, nrow = 2, ncol = 2,
                          dimnames = list(Predicted = y_levels, TrueClass = y_levels))
sum_conf_kmeans_groupmean <- matrix(0, nrow = 2, ncol = 2,
                                    dimnames = list(Predicted = y_levels, TrueClass = y_levels))

for (i in 1:n_iter) {
  test_idx <- sample(seq_len(nrow(data)), size = 10)
  train_data <- data[-test_idx, ]
  test_data  <- data[test_idx, ]
  
  res_lda <- lda_kmeans_result(train_data, test_data, y_levels)
  res_km  <- kmeans_result(train_data, test_data, y_levels)
  res_kmg <- kmeans_groupmean_result(train_data, test_data, y_levels)
  
  acc_lda[i] <- res_lda$accuracy
  acc_kmeans[i] <- res_km$accuracy
  acc_kmeans_groupmean[i] <- res_kmg$accuracy
  
  sum_conf_lda <- sum_conf_lda + res_lda$conf_mat
  sum_conf_kmeans <- sum_conf_kmeans + res_km$conf_mat
  sum_conf_kmeans_groupmean <- sum_conf_kmeans_groupmean + res_kmg$conf_mat
}

mean_conf_lda <- sum_conf_lda / n_iter
mean_conf_kmeans <- sum_conf_kmeans / n_iter
mean_conf_kmeans_groupmean <- sum_conf_kmeans_groupmean / n_iter
cat("LDA + KMeans\n")
cat("Mean accuracy =", mean(acc_lda), "\n")
cat("SD of accuracy =", sd(acc_lda), "\n")
print(round(mean_conf_lda, 3))

cat("\nPlain KMeans\n")
cat("Mean accuracy =", mean(acc_kmeans), "\n")
cat("SD of accuracy =", sd(acc_kmeans), "\n")
print(round(mean_conf_kmeans, 3))

cat("\nKMeans with group-mean initial centers\n")
cat("Mean accuracy =", mean(acc_kmeans_groupmean), "\n")
cat("SD of accuracy =", sd(acc_kmeans_groupmean), "\n")
print(round(mean_conf_kmeans_groupmean, 3))