#-------------------------------------------------------------------------------
# Toy example: Identical marginal distributions and correlations,
# but different distributions
# 
# last change: 2026-05-08
#-------------------------------------------------------------------------------

# Random seed for reproducibility ----------------------------------------------
set.seed(2026-04-02)

# Packages ---------------------------------------------------------------------
library(randomForest)
library(tidyverse)
library(here)
library(slicer)

# Custom functions -------------------------------------------------------------
# Generating empirical distributions
generate_data <- function(n = 1000, a = 2, angle = 0, Sigma = diag(1, 2)) {
  mus <- cbind(c(a, 0), c(0, a), c(-a, 0), c(0, -a))
  R <- rbind(
    c(cos(angle), -sin(angle)),
    c(sin(angle), cos(angle))
  )
  cluster <- sample(1:4, size = n, replace = TRUE)
  d <- matrix(0, nrow = n, ncol = 2)
  for (i in 1:4) {
    d[which(cluster == i), ] <- MASS::mvrnorm(
      sum(cluster == i), mu = mus[, i], Sigma = Sigma
    )
  }
  d %*% t(R)
}

# Currently, the outcome is just the identity on the angle. Change this here.
my_f <- function(x) {
  x
}

# Simulation constants ---------------------------------------------------------
sample_size <- 200
a <- 5
N_train <- 50
N_test  <- 50
training_idx <- 1:N_train
test_idx <- (N_train + 1):(N_train + N_test)
runs <- 20
L <- 50 # 50 thetas for sliced Wasserstein
n_dim <- 2
shutdown <- FALSE # shutdown computer 20' after running script?

# Scenario 1: Generate angles from Unif([0, pi/4]) -----------------------------
angle_range <- c(0, pi/4)

all_data_list <- vector("list", runs)
for (run in seq_len(runs)) {
  train_list <- vector("list", N_train)
  test_list <- vector("list", N_test)
  angles_train <- runif(N_train, angle_range[1], angle_range[2])
  angles_test <- runif(N_test, angle_range[1], angle_range[2])
  y_train <- my_f(angles_train)
  y_test <- my_f(angles_test)
  
  for (n in seq_len(N_train)) {
    train_list[[n]] <- generate_data(sample_size, a = a, 
      angle = angles_train[[n]])
  }
  for (n in seq_len(N_test)) {
    test_list[[n]] <- generate_data(sample_size, a = a, 
      angle = angles_test[[n]])
  }
  
  current_data <- list(
    y_train = y_train,
    y_test = y_test,
    train_list = train_list,
    test_list = test_list
  )
  
  all_data_list[[run]] <- current_data
}

# Analyse each run -------------------------------------------------------------
results <- data.frame(
  run = numeric(),
  method = character(),
  true = numeric(),
  prediction = numeric(),
  stringsAsFactors = FALSE
)

first_axis_models <- vector("list", runs)
second_axis_models <- vector("list", runs)
sliced_models <- vector("list", runs)
combined_models <- vector("list", runs)

for (run in seq_len(runs)) {
  message(paste0("Run ", run, " of ", runs, "."))
  current_data_list <- c(all_data_list[[run]]$train_list, 
                         all_data_list[[run]]$test_list)
  y_train <- all_data_list[[run]]$y_train
  y_test <- all_data_list[[run]]$y_test
  
  # Compute moment-based representation ----------------------------------------
  M <- matrix(0, nrow = N_train + N_test, ncol = 5)
  for (n in 1:(N_train + N_test)) {
    d <- current_data_list[[n]]
    M[n, 1:2] <- colMeans(d) |> as.vector()
    M[n, 3:4] <- apply(d, 2, var) |> as.vector()
    M[n, 5] <- cor(d)[1, 2]
  }
  
  M_train <- M[1:N_train, ]
  M_test <- M[(N_train + 1):(N_train + N_test), ]
  
  # Use mean of training set as prediction -------------------------------------
  my_results <- data.frame(
    run = run,
    method = "training set mean",
    true = y_test,
    prediction = mean(y_train),
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Vector-based model: random forests -----------------------------------------
  # Tune mtry
  best_mse <- Inf
  baseline.rf <- NULL
  for (m in seq_len(5)) {
    mod.rf <- randomForest(M_train, y_train, mtry = 2)
    if (mean(mod.rf$mse) < best_mse) {
      best_mse <- mean(mod.rf$mse)
      baseline.rf <- mod.rf
    }
  }
  predictions.rf <- predict(baseline.rf, M_test)
  my_results <- data.frame(
    run = run,
    method = "random forest",
    true = y_test,
    prediction = predictions.rf,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Compute Wasserstein distances along cardinal axes --------------------------
  cardinal_distances <- compute_all_distances(
    current_data_list, diag(1, n_dim), verbose = FALSE)
  
  # Wasserstein-based method: first axis ---------------------------------------
  first_axis <- fit_gpr(cardinal_distances[[1]], training_idx, test_idx, 
    y_train, y_test, verbose = FALSE)
  first_axis_models[[run]] <- first_axis
  my_results <- data.frame(
    run = run,
    method = "Wasserstein: first axis",
    true = y_test,
    prediction = first_axis$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Wasserstein-based method: first axis ---------------------------------------
  second_axis <- fit_gpr(cardinal_distances[[2]], training_idx, test_idx, 
    y_train, y_test, verbose = FALSE)
  second_axis_models[[run]] <- second_axis
  my_results <- data.frame(
    run = run,
    method = "Wasserstein: second axis",
    true = y_test,
    prediction = second_axis$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Sliced Wasserstein distances------------------------------------------------
  thetas <- generate_directions(L = L, d = n_dim)
  data_list <- c(train_list, test_list)
  distances <- compute_all_distances(current_data_list, thetas, 
    verbose = TRUE, keep_projections = FALSE)
  sliced_ws <- fit_gpr(distances, training_idx, test_idx, y_train, y_test)
  sliced_models[[run]] <- sliced_ws
  my_results <- data.frame(
    run = run,
    method = "sliced Wasserstein",
    true = y_test,
    prediction = sliced_ws$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Kernel learnt from cardinal axes. ------------------------------------------
  cardinal_axes <- fit_gpr_multiple(cardinal_distances, training_idx, test_idx, 
    y_train, y_test)
  combined_models[[run]] <- cardinal_axes
  my_results <- data.frame(
    run = run,
    method = "combined kernel (cardinal axes)",
    true = y_test,
    prediction = cardinal_axes$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
}

models45 <- list(
  first_axis_models,
  second_axis_models,
  sliced_models,
  combined_models
)
results45 <- results
saveRDS(models45, here("results", "toy_example", "models45.Rda"))
saveRDS(results45, here("results", "toy_example", "results45.Rda"))

# Scenario 2: Generate angles from Unif([0, pi/2]) -----------------------------
angle_range <- c(0, pi/2)

all_data_list <- vector("list", runs)
for (run in seq_len(runs)) {
  train_list <- vector("list", N_train)
  test_list <- vector("list", N_test)
  angles_train <- runif(N_train, angle_range[1], angle_range[2])
  angles_test <- runif(N_test, angle_range[1], angle_range[2])
  y_train <- my_f(angles_train)
  y_test <- my_f(angles_test)
  
  for (n in seq_len(N_train)) {
    train_list[[n]] <- generate_data(sample_size, a = a, 
      angle = angles_train[[n]])
  }
  for (n in seq_len(N_test)) {
    test_list[[n]] <- generate_data(sample_size, a = a,
      angle = angles_test[[n]])
  }
  
  current_data <- list(
    y_train = y_train,
    y_test = y_test,
    train_list = train_list,
    test_list = test_list
  )
  
  all_data_list[[run]] <- current_data
}

# Analyse each run -------------------------------------------------------------
results <- data.frame(
  run = numeric(),
  method = character(),
  true = numeric(),
  prediction = numeric(),
  stringsAsFactors = FALSE
)

first_axis_models <- vector("list", runs)
second_axis_models <- vector("list", runs)
sliced_models <- vector("list", runs)
combined_models <- vector("list", runs)

for (run in seq_len(runs)) {
  message(paste0("Run ", run, " of ", runs, "."))
  current_data_list <- c(all_data_list[[run]]$train_list, 
                         all_data_list[[run]]$test_list)
  y_train <- all_data_list[[run]]$y_train
  y_test <- all_data_list[[run]]$y_test
  
  # Compute moment-based representation ----------------------------------------
  M <- matrix(0, nrow = N_train + N_test, ncol = 5)
  for (n in 1:(N_train + N_test)) {
    d <- current_data_list[[n]]
    M[n, 1:2] <- colMeans(d) |> as.vector()
    M[n, 3:4] <- apply(d, 2, var) |> as.vector()
    M[n, 5] <- cor(d)[1, 2]
  }
  
  M_train <- M[1:N_train, ]
  M_test <- M[(N_train + 1):(N_train + N_test), ]
  
  # Use mean of training set as prediction -------------------------------------
  my_results <- data.frame(
    run = run,
    method = "training set mean",
    true = y_test,
    prediction = mean(y_train),
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Vector-based model: random forests -----------------------------------------
  # Tune mtry
  best_mse <- Inf
  baseline.rf <- NULL
  for (m in seq_len(5)) {
    mod.rf <- randomForest(M_train, y_train, mtry = 2)
    if (mean(mod.rf$mse) < best_mse) {
      best_mse <- mean(mod.rf$mse)
      baseline.rf <- mod.rf
    }
  }
  predictions.rf <- predict(baseline.rf, M_test)
  my_results <- data.frame(
    run = run,
    method = "random forest",
    true = y_test,
    prediction = predictions.rf,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Compute Wasserstein distances along cardinal axes --------------------------
  cardinal_distances <- compute_all_distances(
    current_data_list, diag(1, 2), verbose = FALSE)
  
  # Wasserstein-based method: first axis ---------------------------------------
  first_axis <- fit_gpr(cardinal_distances[[1]], training_idx, test_idx, 
    y_train, y_test, verbose = FALSE)
  first_axis_models[[run]] <- first_axis
  my_results <- data.frame(
    run = run,
    method = "Wasserstein: first axis",
    true = y_test,
    prediction = first_axis$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Wasserstein-based method: first axis ---------------------------------------
  second_axis <- fit_gpr(cardinal_distances[[2]], training_idx, test_idx, 
    y_train, y_test, verbose = FALSE)
  second_axis_models[[run]] <- second_axis
  my_results <- data.frame(
    run = run,
    method = "Wasserstein: second axis",
    true = y_test,
    prediction = second_axis$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Sliced Wasserstein distances -----------------------------------------------
  thetas <- generate_directions(L = L, d = n_dim)
  data_list <- c(train_list, test_list)
  distances <- compute_all_distances(current_data_list, thetas, 
    verbose = TRUE, keep_projections = FALSE)
  sliced_ws <- fit_gpr(distances, training_idx, test_idx, y_train, y_test)
  sliced_models[[run]] <- sliced_ws
  my_results <- data.frame(
    run = run,
    method = "sliced Wasserstein",
    true = y_test,
    prediction = sliced_ws$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
  
  # Kernel learnt from cardinal axes. ------------------------------------------
  cardinal_axes <- fit_gpr_multiple(cardinal_distances, training_idx, 
    test_idx, y_train, y_test)
  combined_models[[run]] <- cardinal_axes
  my_results <- data.frame(
    run = run,
    method = "combined kernel (cardinal axes)",
    true = y_test,
    prediction = cardinal_axes$test_predictions,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, my_results) 
}

models90 <- list(
  first_axis_models,
  second_axis_models,
  sliced_models,
  combined_models
)
results90 <- results
saveRDS(models90, here("results", "toy_example", "models90.Rda"))
saveRDS(results90, here("results", "toy_example", "results90.Rda"))

# Shut down computer after 20' (after running overnight) -----------------------
if (shutdown) system("shutdown /s /t 1200 /f")

# Software versions ------------------------------------------------------------
devtools::session_info("attached")
# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.5.0 (2025-04-11 ucrt)
# os       Windows 11 x64 (build 26200)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  English_United Kingdom.utf8
# ctype    English_United Kingdom.utf8
# tz       Europe/Zurich
# date     2026-05-11
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# (...)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────
# package      * version    date (UTC) lib source
# dplyr        * 1.2.1      2026-04-03 [1] CRAN (R 4.5.3)
# forcats      * 1.0.1      2025-09-25 [1] CRAN (R 4.5.3)
# ggplot2      * 4.0.3      2026-04-22 [1] CRAN (R 4.5.3)
# here         * 1.0.2      2025-09-15 [1] CRAN (R 4.5.3)
# lubridate    * 1.9.5      2026-02-04 [1] CRAN (R 4.5.3)
# purrr        * 1.0.4      2025-02-05 [1] CRAN (R 4.5.0)
# randomForest * 4.7-1.2    2024-09-22 [1] CRAN (R 4.5.2)
# readr        * 2.2.0      2026-02-19 [1] CRAN (R 4.5.3)
# slicer       * 0.0.0.9000 2026-05-11 [1] Github (janhove/slicer@83f3f04)
# stringr      * 1.6.0      2025-11-04 [1] CRAN (R 4.5.3)
# tibble       * 3.2.1      2023-03-20 [1] CRAN (R 4.5.0)
# tidyr        * 1.3.2      2025-12-19 [1] CRAN (R 4.5.3)
# tidyverse    * 2.0.0      2023-02-22 [1] CRAN (R 4.5.3)
# (...)
# END --------------------------------------------------------------------------
################################################################################