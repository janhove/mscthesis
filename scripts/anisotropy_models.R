#-------------------------------------------------------------------------------
# Anisotropy: GP models and figures
# 
# last change: 2026-05-08
#-------------------------------------------------------------------------------

# Random seed for reproducibility ----------------------------------------------
set.seed(2026-05-08)

# Packages ---------------------------------------------------------------------
library(here)
library(slicer)

# Custum functions -------------------------------------------------------------
rmse <- function(x, y) {
  sqrt(mean((x - y)^2))
}

# Outcome function
g <- function(x) ifelse(x == 0, 0, (cos(2*pi*x) - cos(pi * x))/(pi * x))

# Simulation settings ----------------------------------------------------------
# See anisotropy_distances.R
R <- 10   # number of simulation runs
N1 <- 50  # training sample size
N2 <- 50  # test sample size
n <- 200  # size of each empirical distribution
L <- 50   # number of projection directions

alphas <- 10^seq(-6, 6)

fixed_nugget <- 1e-6
fixed_ls <- 1
fixed_s2 <- 1

# Load data --------------------------------------------------------------------
mu_list <- readRDS(here("results", "anisotropy", "theoretical_means.Rda"))
data_list <- readRDS(here("results", "anisotropy", "data.Rda"))
theta_list <- readRDS(here("results", "anisotropy", "thetas.Rda"))
training_idx <- seq_len(N1)
test_idx <- N1 + seq_len(N2)

################################################################################
# Scenario 1: All signal contained in first dimension                          #
################################################################################

f <- function(x, y) 1.5*g(x)
outcome_list <- vector("list", R)
for (r in seq_len(R)) {
  my_mus <- mu_list[[r]]
  my_outcomes <- f(my_mus[, 1], my_mus[, 2])
  y_train <- my_outcomes[training_idx]
  y_test <- my_outcomes[test_idx]
  outcome_list[[r]] <- list(y_train = y_train, y_test = y_test)
}

# Analyse data -----------------------------------------------------------------
# (1) No tuning of hyperparameters.
# (2) With tuning of hyperparameters.
# (a) Raw transformation matrix.
# (b) Normalised transformation matrix.
# For each combination of (1/2) and (a/b), 
# transform the projection direction using the transformation matrix, 
# compute distances, run through both scenarios.
# When tuning hyperparameters, store them.

# Baseline for each run --------------------------------------------------------
baseline_rmse <- vector("numeric", length = R)
for (r in 1:R) {
  current_data <- data_list[[r]]
  current_outcomes <- outcome_list[[r]]
  baseline_rmse[r] <- rmse(current_outcomes$y_test, mean(current_outcomes$y_train))
}
my_results <- data.frame(
  scenario = 1,
  run = 1:R,
  method = "training set mean",
  normalised = NA,
  tuned = NA,
  alpha = NA,
  rmse = baseline_rmse,
  length_scale = NA,
  variance = NA,
  nugget = NA
)

# Raw transformation matrix, no tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_raw_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    K <- rbf(distances, length_scale = fixed_ls, variance = fixed_s2)
    Kxx <- K[training_idx, training_idx] 
    Kxstar <- K[test_idx, training_idx]
    predictions <- gpr_predict(Kxx, Kxstar, y_train, lambda2 = fixed_nugget, centre = TRUE)
    rmses[i] <- rmse(y_test, predictions)
    i <- i + 1
  }
}

results_raw_no_tuning <- expand.grid(
  scenario = 1,
  method = "unnormalised, untuned",
  normalised = FALSE,
  tuned = FALSE,
  length_scale = fixed_ls,
  variance = fixed_s2,
  nugget = fixed_nugget,
  alpha = alphas,
  run = 1:R
)
results_raw_no_tuning$rmse <- rmses
my_results <- my_results |> 
  bind_rows(results_raw_no_tuning)

# Raw transformation matrix, tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
est_ls <- est_s2 <- est_nugget <- rmses
i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_raw_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    my_model <- fit_gpr(distances, training_idx, test_idx, y_train, y_test)
    rmses[i] <- my_model$RMSE
    est_ls[i] <- my_model$length_scale
    est_s2[i] <- my_model$variance
    est_nugget[i] <- my_model$lambda2
    i <- i + 1
  }
}

results_raw_with_tuning <- expand.grid(
  scenario = 1,
  method = "unnormalised, tuned",
  normalised = FALSE,
  tuned = TRUE,
  alpha = alphas,
  run = 1:R
)
results_raw_with_tuning$rmse <- rmses
results_raw_with_tuning$length_scale <- est_ls
results_raw_with_tuning$nugget <- est_nugget
results_raw_with_tuning$variance <- est_s2
my_results <- my_results |> 
  bind_rows(results_raw_with_tuning)

# Normalised transformation matrix, no tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_normalised_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    K <- rbf(distances, length_scale = fixed_ls, variance = fixed_s2)
    Kxx <- K[training_idx, training_idx] 
    Kxstar <- K[test_idx, training_idx]
    predictions <- gpr_predict(Kxx, Kxstar, y_train, lambda2 = fixed_nugget, centre = TRUE)
    rmses[i] <- rmse(y_test, predictions)
    i <- i + 1
  }
}
results_normalised_no_tuning <- expand.grid(
  scenario = 1,
  method = "normalised, untuned",
  normalised = TRUE,
  tuned = FALSE,
  length_scale = fixed_ls,
  variance = fixed_s2,
  nugget = fixed_nugget,
  alpha = alphas,
  run = 1:R
)
results_normalised_no_tuning$rmse <- rmses
my_results <- my_results |> 
  bind_rows(results_normalised_no_tuning)


# Normalised transformation matrix, tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
est_ls <- est_s2 <- est_nugget <- rmses

i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_normalised_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    my_model <- fit_gpr(distances, training_idx, test_idx, y_train, y_test)
    rmses[i] <- my_model$RMSE
    est_ls[i] <- my_model$length_scale
    est_s2[i] <- my_model$variance
    est_nugget[i] <- my_model$lambda2
    i <- i + 1
  }
}

results_normalised_with_tuning <- expand.grid(
  scenario = 1,
  method = "normalised, tuned",
  normalised = TRUE,
  tuned = TRUE,
  alpha = alphas,
  run = 1:R
)
results_normalised_with_tuning$rmse <- rmses
results_normalised_with_tuning$length_scale <- est_ls
results_normalised_with_tuning$nugget <- est_nugget
results_normalised_with_tuning$variance <- est_s2
my_results <- my_results |> 
  bind_rows(results_normalised_with_tuning)

# Using cardinal distances -----------------------------------------------------
rmses <- vector("numeric", length = R)

i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  current_path <- paste0("cardinal_distances_run_", r, ".Rda")
  distances <- readRDS(here("results", "anisotropy", current_path))
  my_model <- fit_gpr_multiple(distances, training_idx, test_idx, y_train, y_test)
  rmses[i] <- my_model$RMSE
  i <- i + 1
}
multiple_1_rmse <- rmses

################################################################################
# Scenario 2: Global signal contained in first dimension, local in second      #
################################################################################

f <- function(x, y) g(x) + 0.5*g(y)
outcome_list <- vector("list", R)
for (r in seq_len(R)) {
  my_mus <- mu_list[[r]]
  my_outcomes <- f(my_mus[, 1], my_mus[, 2])
  y_train <- my_outcomes[training_idx]
  y_test <- my_outcomes[test_idx]
  outcome_list[[r]] <- list(y_train = y_train, y_test = y_test)
}

# Analyse data -----------------------------------------------------------------
# (1) No tuning of hyperparameters.
# (2) With tuning of hyperparameters.
# (a) Raw transformation matrix.
# (b) Normalised transformation matrix.
# For each combination of (1/2) and (a/b), 
# transform the projection direction using the transformation matrix, 
# compute distances, run through both scenarios.
# When tuning hyperparameters, store them.

# Baseline for each run --------------------------------------------------------
baseline_rmse <- vector("numeric", length = R)
for (r in 1:R) {
  current_data <- data_list[[r]]
  current_outcomes <- outcome_list[[r]]
  baseline_rmse[r] <- rmse(current_outcomes$y_test, mean(current_outcomes$y_train))
}
my_results <- my_results |> 
  bind_rows(data.frame(
    scenario = 2,
    run = 1:R,
    method = "training set mean",
    normalised = NA,
    tuned = NA,
    alpha = NA,
    rmse = baseline_rmse,
    length_scale = NA,
    variance = NA,
    nugget = NA
  ))

# Raw transformation matrix, no tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_raw_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    K <- rbf(distances, length_scale = fixed_ls, variance = fixed_s2)
    Kxx <- K[training_idx, training_idx] 
    Kxstar <- K[test_idx, training_idx]
    predictions <- gpr_predict(Kxx, Kxstar, y_train, lambda2 = fixed_nugget, centre = TRUE)
    rmses[i] <- rmse(y_test, predictions)
    i <- i + 1
  }
}

results_raw_no_tuning <- expand.grid(
  scenario = 2,
  method = "unnormalised, untuned",
  normalised = FALSE,
  tuned = FALSE,
  length_scale = fixed_ls,
  variance = fixed_s2,
  nugget = fixed_nugget,
  alpha = alphas,
  run = 1:R
)
results_raw_no_tuning$rmse <- rmses
my_results <- my_results |> 
  bind_rows(results_raw_no_tuning)

# Raw transformation matrix, tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
est_ls <- est_s2 <- est_nugget <- rmses
i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_raw_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    my_model <- fit_gpr(distances, training_idx, test_idx, y_train, y_test)
    rmses[i] <- my_model$RMSE
    est_ls[i] <- my_model$length_scale
    est_s2[i] <- my_model$variance
    est_nugget[i] <- my_model$lambda2
    i <- i + 1
  }
}

results_raw_with_tuning <- expand.grid(
  scenario = 2,
  method = "unnormalised, tuned",
  normalised = FALSE,
  tuned = TRUE,
  alpha = alphas,
  run = 1:R
)
results_raw_with_tuning$rmse <- rmses
results_raw_with_tuning$length_scale <- est_ls
results_raw_with_tuning$nugget <- est_nugget
results_raw_with_tuning$variance <- est_s2
my_results <- my_results |> 
  bind_rows(results_raw_with_tuning)

# Normalised transformation matrix, no tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_normalised_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    K <- rbf(distances, length_scale = fixed_ls, variance = fixed_s2)
    Kxx <- K[training_idx, training_idx] 
    Kxstar <- K[test_idx, training_idx]
    predictions <- gpr_predict(Kxx, Kxstar, y_train, lambda2 = fixed_nugget, centre = TRUE)
    rmses[i] <- rmse(y_test, predictions)
    i <- i + 1
  }
}
results_normalised_no_tuning <- expand.grid(
  scenario = 2,
  method = "normalised, untuned",
  normalised = TRUE,
  tuned = FALSE,
  length_scale = fixed_ls,
  variance = fixed_s2,
  nugget = fixed_nugget,
  alpha = alphas,
  run = 1:R
)
results_normalised_no_tuning$rmse <- rmses
my_results <- my_results |> 
  bind_rows(results_normalised_no_tuning)


# Normalised transformation matrix, tuning -----------------------------------------
rmses <- vector("numeric", length = R * length(alphas))
est_ls <- est_s2 <- est_nugget <- rmses

i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  for (alpha in alphas) {
    message(paste0("Current alpha: ", alpha, "."))
    current_path <- paste0("distances_normalised_alpha_", alpha, "_run_", r, ".Rda")
    distances <- readRDS(here("results", "anisotropy", current_path))
    my_model <- fit_gpr(distances, training_idx, test_idx, y_train, y_test)
    rmses[i] <- my_model$RMSE
    est_ls[i] <- my_model$length_scale
    est_s2[i] <- my_model$variance
    est_nugget[i] <- my_model$lambda2
    i <- i + 1
  }
}

results_normalised_with_tuning <- expand.grid(
  scenario = 2,
  method = "normalised, tuned",
  normalised = TRUE,
  tuned = TRUE,
  alpha = alphas,
  run = 1:R
)
results_normalised_with_tuning$rmse <- rmses
results_normalised_with_tuning$length_scale <- est_ls
results_normalised_with_tuning$nugget <- est_nugget
results_normalised_with_tuning$variance <- est_s2
my_results <- my_results |> 
  bind_rows(results_normalised_with_tuning)

# Using cardinal distances -----------------------------------------------------
rmses <- vector("numeric", length = R)

i <- 1
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  y_train <- outcome_list[[r]]$y_train
  y_test <- outcome_list[[r]]$y_test
  current_path <- paste0("cardinal_distances_run_", r, ".Rda")
  distances <- readRDS(here("results", "anisotropy", current_path))
  my_model <- fit_gpr_multiple(distances, training_idx, test_idx, y_train, y_test)
  predictions <- gpr_predict(Kxx, Kxstar, y_train, lambda2 = fixed_nugget, centre = TRUE)
  rmses[i] <- my_model$RMSE
  i <- i + 1
}
multiple_2_rmse <- rmses

################################################################################
# Graphs and summaries                                                         #
################################################################################

my_results$tuning <- ifelse(my_results$tuned, "tuned hyperparameters", "fixed hyperparameters")
my_results$trafo <- ifelse(my_results$normalised, 
                           "normalised transformation matrix",
                           "unnormalised transformation matrix") |> 
  factor() |> 
  relevel("unnormalised transformation matrix")

baseline_1 <- my_results |> 
  filter(scenario == 1) |> 
  filter(method == "training set mean") |> 
  select(rmse)
baseline_2 <- my_results |> 
  filter(scenario == 2) |> 
  filter(method == "training set mean") |> 
  select(rmse)

multiple_1 <- mean(multiple_1_rmse)
multiple_2 <- mean(multiple_2_rmse)

my_results |> 
  filter(scenario == 1) |> 
  filter(method != "training set mean") |> 
  group_by(method, alpha, trafo, tuning) |> 
  summarise(mean_rmse = mean(rmse), 
            sd_rmse = sd(rmse)) |> 
  ggplot(aes(x = alpha, y = mean_rmse,
             ymin = mean_rmse - sd_rmse,
             ymax = mean_rmse + sd_rmse)) +
  scale_x_log10(
    breaks = (10^seq(-3, 3))^2,
    labels = scales::label_math(10^.x, format = log10)
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  geom_line() +
  geom_linerange() +
  geom_hline(yintercept = baseline_1$rmse |> mean(), linetype = "dotted") +
  geom_hline(yintercept = multiple_1, linetype = "dashed") +
  xlab(expression(alpha)) +
  ylab("mean RMSE ± SD") +
  facet_grid(cols = vars(tuning),
             rows = vars(trafo)) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position = "bottom") +
  labs(
    title = "Predictive accuracy in Scenario 1",
    subtitle = "all signal contained in first dimension"
  )
ggsave(here("figures", "results_anisotropy_1.pdf"), width = 1.3*6, height = 1.3*4.5)

my_results |> 
  filter(scenario == 2) |> 
  filter(method != "training set mean") |> 
  group_by(method, alpha, trafo, tuning) |> 
  summarise(mean_rmse = mean(rmse), 
            sd_rmse = sd(rmse)) |> 
  ggplot(aes(x = alpha, y = mean_rmse,
             ymin = mean_rmse - sd_rmse,
             ymax = mean_rmse + sd_rmse)) +
  scale_x_log10(
    breaks = (10^seq(-3, 3))^2,
    labels = scales::label_math(10^.x, format = log10)
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  geom_line() +
  geom_linerange() +
  geom_hline(yintercept = baseline_2$rmse |> mean(), linetype = "dotted") +
  geom_hline(yintercept = multiple_2, linetype = "dashed") +
  xlab(expression(alpha)) +
  ylab("mean RMSE ± SD") +
  facet_grid(cols = vars(tuning),
             rows = vars(trafo)) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position = "bottom") +
  labs(
    title = "Predictive accuracy in Scenario 2",
    subtitle = "first dimension dominant"
  )
ggsave(here("figures", "results_anisotropy_2.pdf"), width = 1.3*6, height = 1.3*4.5)


my_results |> 
  group_by(scenario, method, alpha, normalised, tuned) |> 
  summarise(mean_rmse = mean(rmse),
            sd_rmse = sd(rmse))
