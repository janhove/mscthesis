#-------------------------------------------------------------------------------
# Anisotropy: Distance computations
# 
# last change: 2026-05-11
#-------------------------------------------------------------------------------

# Random seed for reproducibility ----------------------------------------------
set.seed(2026-04-17)

# Packages ---------------------------------------------------------------------
library(here)
library(slicer)

# Simulation settings ----------------------------------------------------------
shutdown <- FALSE # shutdown computer 20' after running script
R <- 10           # number of simulation runs
N1 <- 50          # training sample size
N2 <- 50          # test sample size
n <- 200          # size of each empirical distribution
L <- 50           # number of projection directions for sliced WS
mu_range <- pi
alphas <- 10^seq(-6, 6)

# Generate data ----------------------------------------------------------------
mu_list <- vector("list", R)
for (r in 1:R) {
  mu_list[[r]] <- cbind(runif(N1 + N2, -mu_range, mu_range), 
                        runif(N1 + N2, -mu_range, mu_range))
}
saveRDS(mu_list, here("results", "anisotropy", "theoretical_means.Rda"))

data_list <- vector("list", R)
for (r in 1:R) {
  distribution_list <- vector("list", N1 + N2)
  outcomes <- vector("numeric", N1 + N2)
  for (i in 1:(N1 + N2)) {
    current_mu <- mu_list[[r]][i, ] |> as.vector()
    distribution_list[[i]] <- MASS::mvrnorm(n = n, mu = current_mu, Sigma = diag(1, 2))
  }
  data_list[[r]] <- list(
    train_list = distribution_list[1:N1], 
    test_list = distribution_list[(N1+1):(N1+N2)]
  )
}

saveRDS(data_list, here("results", "anisotropy", "data.Rda"))

# Generate projection directions -----------------------------------------------
theta_list <- vector("list", R)
for (r in 1:R) {
  thetas <- generate_directions(L, 2)
  theta_list[[r]] <- thetas
}

saveRDS(theta_list, here("results", "anisotropy", "thetas.Rda"))

# Compute distances ------------------------------------------------------------
# (a) Raw transformation matrix.
# (b) Normalised transformation matrix.
# (c) Distances along cardinal directions.

# Distance computation, raw transformation matrix ------------------------------
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  my_distributions <- c(data_list[[r]]$train_list, data_list[[r]]$test_list)
  my_thetas <- theta_list[[r]]
  for (alpha in alphas) {
    message(paste("Current alpha:", alpha, "..."))
    current_path <- paste0("distances_raw_alpha_", alpha, "_run_", r, ".Rda")
    A <- diag(c(alpha, 1))
    distances <- compute_all_distances(my_distributions, my_thetas, A,
      verbose = FALSE, keep_projections = FALSE)
    saveRDS(distances, here("results", "anisotropy", current_path))
  }
}

# Distance computation, normalised transformation matrix -----------------------
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  my_distributions <- c(data_list[[r]]$train_list, data_list[[r]]$test_list)
  my_thetas <- theta_list[[r]]
  for (alpha in alphas) {
    message(paste("Current alpha:", alpha, "..."))
    current_path <- paste0("distances_normalised_alpha_", alpha, "_run_", r, ".Rda")
    A <- diag(c(alpha, 1)) / sqrt(alpha)
    distances <- compute_all_distances(my_distributions, my_thetas, A,
      verbose = FALSE, keep_projections = FALSE)
    saveRDS(distances, here("results", "anisotropy", current_path))
  }
}

# Distances along cardinal directions ------------------------------------------
for (r in 1:R) {
  message(paste0("Run ", r, " of ", R, "."))
  my_distributions <- c(data_list[[r]]$train_list, data_list[[r]]$test_list)
  cardinal_distances <- compute_all_distances(my_distributions, 
    diag(1, 2), verbose = TRUE)
  current_path <- paste0("cardinal_distances_run_", r, ".Rda")
  saveRDS(cardinal_distances, here("results", "anisotropy", current_path))
}

# Shut down 20' after finishing the script -------------------------------------
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