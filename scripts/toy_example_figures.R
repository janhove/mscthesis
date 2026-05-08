#-------------------------------------------------------------------------------
# Toy example: Identical marginal distributions and correlations,
# but different distributions: Figures
# 
# last change: 2026-05-08
#-------------------------------------------------------------------------------

# Packages and custom functions ------------------------------------------------
library(tidyverse)
library(here)

# Load data --------------------------------------------------------------------
models90 <- readRDS(here("results", "toy_example", "models90.Rda"))
results90 <- readRDS(here("results", "toy_example", "results90.Rda"))
models45 <- readRDS(here("results", "toy_example", "models45.Rda"))
results45 <- readRDS(here("results", "toy_example", "results45.Rda"))

# Summarise findings in terms of RMSE ------------------------------------------
results45 |>
  group_by(run, method) |>
  summarise(RMSE = sqrt(mean((prediction - true)^2)),
            .groups = "drop") |>
  ggplot(
    aes(y = method, x = RMSE)
  ) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 1, position = position_jitter(width = 0, height = 0.1)) +
  scale_y_discrete(
    limits = c("training set mean", "random forest", "Wasserstein: first axis",
               "Wasserstein: second axis", "combined kernel (cardinal axes)",
               "sliced Wasserstein") |> rev(),
    labels = c("Baseline: Training set mean", "Baseline: Random forest",
               "Wasserstein kernel: first dimension", "Wasserstein kernel: second dimension",
               "Sum of Wasserstein kernels", "Sliced Wasserstein kernel") |> rev()
  ) + 
  ylab(NULL) +
  labs(title = "Predictive accuracy in Scenario 1",
       subtitle = expression(paste(omega, " from Uniform(", 0, ",", pi/4, ")"))) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))

ggsave(here("figures", "toy_example45_rmse.pdf"), width = 6, height = 3.8)

results45 |>
  group_by(run, method) |>
  summarise(RMSE = sqrt(mean((prediction - true)^2)),
            .groups = "drop") |>
  group_by(method) |>
  summarise(mean(RMSE), sd(RMSE))

results90 |>
  group_by(run, method) |>
  summarise(RMSE = sqrt(mean((prediction - true)^2)),
            .groups = "drop") |>
  ggplot(
    aes(y = method, x = RMSE)
  ) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 1, position = position_jitter(width = 0, height = 0.1)) +
  scale_y_discrete(
    limits = c("training set mean", "random forest", "Wasserstein: first axis",
               "Wasserstein: second axis", "combined kernel (cardinal axes)",
               "sliced Wasserstein") |> rev(),
    labels = c("Baseline: Training set mean", "Baseline: Random forest",
               "Wasserstein kernel: first dimension", "Wasserstein kernel: second dimension",
               "Sum of Wasserstein kernels", "Sliced Wasserstein kernel") |> rev()
  ) + 
  ylab(NULL) +
  labs(title = "Predictive accuracy in Scenario 2",
       subtitle = expression(paste(omega, " from Uniform(", 0, ",", pi/2, ")"))) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))

ggsave(here("figures", "toy_example90_rmse.pdf"), width = 6, height = 3.8)

results90 |>
  group_by(run, method) |>
  summarise(RMSE = sqrt(mean((prediction - true)^2)),
            .groups = "drop") |>
  group_by(method) |>
  summarise(mean(RMSE), sd(RMSE))

# Summarise findings in terms of R^2 -------------------------------------------
results45 <- results45 |> mutate(scenario = "Scenario 1")
results90 <- results90 |> mutate(scenario = "Scenario 2")
results <- results45 |> bind_rows(results90)
r2_results <- results |> 
  group_by(scenario, method, run) |> 
  summarise(
    mse = mean((prediction - true)^2),
    .groups = "drop"
  ) |> 
  pivot_wider(names_from = method, values_from = mse) |> 
  mutate(
    R2_wasserstein_first = 1 - `Wasserstein: first axis`/`training set mean`,
    R2_wasserstein_second = 1 - `Wasserstein: second axis`/`training set mean`,
    R2_random_forest = 1 - `random forest`/`training set mean`,
    R2_combined_kernel = 1 - `combined kernel (cardinal axes)`/`training set mean`,
    R2_sliced_wasserstein = 1 - `sliced Wasserstein`/`training set mean`
  ) |> 
  pivot_longer(
    cols = starts_with("R2_"),
    values_to = "R2",
    names_to = "method",
    names_prefix = "R2_"
  ) |> 
  select(scenario, run, method, R2)

r2_results |>
  ggplot(
    aes(y = method, x = R2)
  ) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 1, position = position_jitter(width = 0, height = 0.1)) +
  scale_y_discrete(
    limits = c("random_forest", "wasserstein_first",
               "wasserstein_second", "combined_kernel",
               "sliced_wasserstein") |> rev(),
    labels = c("Baseline: Random forest",
               "Wasserstein kernel: first dimension", "Wasserstein kernel: second dimension",
               "Sum of Wasserstein kernels", "Sliced Wasserstein kernel") |> rev()
  ) +
  facet_wrap(vars(scenario), ncol = 1) +
  xlab(expression(R^2)) +
  ylab(NULL) +
  labs(title = expression(paste(R^2, " in both scenarios"))) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))

ggsave(here("figures", "toy_exampleR2.pdf"), width = 6, height = 5)

r2_results |>
  group_by(method, scenario) |>
  summarise(mean(R2), sd(R2))

# END --------------------------------------------------------------------------
################################################################################