#-------------------------------------------------------------------------------
# Analysis of brain data
# 
# last change: 2026-05-11
#-------------------------------------------------------------------------------

# Random seed ------------------------------------------------------------------
set.seed(2026-05-11)

# Packages and custom functions ------------------------------------------------
library(tidyverse)
library(purrr)
library(readxl)
library(here)
library(slicer)
# You'll also need the whitening package.

# Scaling / whitening functions ------------------------------------------------
# Note: centring doesn't matter for either Wasserstein or sliced Wasserstein
get_trafo <- function(distributions, type) {
  sigmas <- map(distributions, cov)
  joint_sigma <- Reduce(`+`, sigmas) / length(sigmas)
  
  if (type == "scale") return(diag(sqrt(1/diag(joint_sigma))))
  whitening::whiteningMatrix(joint_sigma, method = type)
}

# Magical constants ------------------------------------------------------------
L <- 200  # number of projection directions for sliced Wasserstein
d <- 4    # inputs are in 4 dimensions
thetas <- generate_directions(L, d)

# Data -------------------------------------------------------------------------
files <- list.files(here("data", "eckermann"), pattern = "*.xlsx")
ids <- str_split_i(files, "_", -2)
classes <- str_split_i(files, "_", -3)
tibbles <- map(here("data", "eckermann", files), 
               ~ read_xlsx(.x) |> 
                 mutate(id = str_split_i(.x, "_", -2)) |> 
                 mutate(class = str_split_i(.x, "_", -3))) |> 
  list_rbind()

# Eckermann et al. applied a couple of filtering and transformation steps to
# their data; see their Analyse_features_sep_PB_public.R. Below I carry out
# these same steps, referencing the corresponding lines in Eckermann et al.’s R
# script.

# Filtering
tibbles <- tibbles |> 
  filter(cellSphericity <= 1) |>                                      # l. 36-37
  filter(NNdist_1_um >= 0.6 * 8.8) |>                                 # l. 38-39
  filter(NNdist_1_um <= 8.8 + 3*2.3212) |>                            # l. 40-41
  filter(cellVolume_um3 >= 35)                                        # l. 42-43

# Compute relative heterogeneity in electron density                  # l. 45-46
tibbles <- tibbles |> 
  group_by(id) |> 
  mutate(edensityVar_rel = edensityVar_nm3 / mean(edensityMean_nm3)^2) |> 
  ungroup()

# Compute number of neighbouring objects within 13.5 µm.
# The original code is pretty strange (l. 51-79), but I'm pretty confident
# that it boils down to this:
tibbles$NumNeighs_withinR14um <- (tibbles[paste0("NNdist_", 1:20, "_um")] <= 13.5) |> 
  rowSums(na.rm = TRUE)

tibbles <- tibbles |> 
  select(id, class, 
         compactness = edensityMedian_nm3,
         heterogeneity = edensityVar_rel,
         size = cellVolume_um3,
         form = cellSphericity,
         packing = NumNeighs_withinR14um)

log_trafo <- TRUE
tibbles <- tibbles |>
  filter(compactness > 0)
if (log_trafo) {
  tibbles <- tibbles |> 
    mutate(size = log10(size))
}

saveRDS(tibbles, here("results", "brains", "tibbles.Rda"))

tibbles |> 
  group_by(id) |> 
  summarise(n = n()) |> 
  arrange(desc(n))

# id <- "04"
# n <- tibbles |> filter(id == id) |> nrow()
# tibbles |> 
#   filter(id == id) |> 
#   select(packing, compactness, heterogeneity, size, form) |> 
#   scatterplot_matrix(
#     labels = c("Neighbours", "Compactness", "Norm. variance", "Volume", "Sphericity"),
#     top = panel.density,
#     main = paste0("Distribution for brain ", id, ", (n = ", n, ")")
#   )

# Goal -------------------------------------------------------------------------
# We'll use the first four dimensions (compactness, heterogeneity, size, form)
# to predict a property of the fourth dimension.
get_matrix <- function(.x, data, cols = NULL) {
  d <- data |> 
    filter(id == .x) |> 
    select(-id, -class) |> 
    as.matrix()
  if (is.null(cols)) return(d)
  d[, cols]
}
input_distributions <- map(ids, get_matrix, data = tibbles, cols = 1:4)
output_distributions <- map(ids, get_matrix, data = tibbles, cols = 5)

# Distance computations without transformations --------------------------------
# These distances would be the same even in LOOCV.

# Sliced Wasserstein
sw_dist <- compute_all_distances(input_distributions, thetas, verbose = TRUE,
                                 keep_projections = FALSE)
saveRDS(sw_dist, here("results", "brains", "no_trafo_sliced.Rda"))

# Marginal Wasserstein
marginal_dist <- compute_all_distances(input_distributions, diag(1, d), verbose = TRUE,
                                       keep_projections = TRUE)
saveRDS(marginal_dist, here("results", "brains", "no_trafo_margins.Rda"))

# Distance computations after scaling ------------------------------------------
# The scaling will vary depending on the input distribution that's been left out.
# For the scaling, we're going to use a pooled estimate.
loocv_sw_dists_scaled <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "scale")
  loocv_sw_dists_scaled[[i]] <- compute_all_distances(input_distributions, 
                                                      thetas, A = W, verbose = TRUE, keep_projections = FALSE)
  cat("\n")
}
saveRDS(loocv_sw_dists_scaled, here("results", "brains", "scaled_sliced.Rda"))

# The marginal distances will only be affected by scaling up to a factor,
# so no need to recompute them. But to make sure...
loocv_marginal_dists_scaled <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "scale")
  loocv_marginal_dists_scaled[[i]] <- compute_all_distances(input_distributions, 
                                                            diag(1, d), A = W, verbose = TRUE, keep_projections = TRUE)
  cat("\n")
}
saveRDS(loocv_marginal_dists_scaled, here("results", "brains", "scaled_margins.Rda"))

# Distance computations after whitening ----------------------------------------
# The scaling will vary depending on the input distribution that's been left out.
# For the whitening, we're going to use a pooled estimate.
# The precise whitening procedure doesn't matter.
loocv_sw_dists_whitened <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "ZCA")
  loocv_sw_dists_whitened[[i]] <- compute_all_distances(input_distributions, 
                                                        thetas, A = W, verbose = TRUE, keep_projections = FALSE)
  cat("\n")
}
saveRDS(loocv_sw_dists_whitened, here("results", "brains", "whitened_sliced.Rda"))

# This time, we do need to recompute the marginal distances.
# The output will be a list of length 20, each element of which is a list of 4 matrices.
loocv_marginal_dists_zca <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "ZCA")
  loocv_marginal_dists_zca[[i]] <- compute_all_distances(input_distributions, 
                                                         diag(1, d), A = W, verbose = TRUE, keep_projections = TRUE)
  cat("\n")
}
saveRDS(loocv_marginal_dists_zca, here("results", "brains", "zca_margins.Rda"))

loocv_marginal_dists_zca_cor <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "ZCA-cor")
  loocv_marginal_dists_zca_cor[[i]] <- compute_all_distances(input_distributions, 
                                                             diag(1, d), A = W, verbose = TRUE, keep_projections = TRUE)
  cat("\n")
}
saveRDS(loocv_marginal_dists_zca_cor, here("results", "brains", "zca_cor_margins.Rda"))

loocv_marginal_dists_pca <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "PCA")
  loocv_marginal_dists_pca[[i]] <- compute_all_distances(input_distributions, 
                                                         diag(1, d), A = W, verbose = TRUE, keep_projections = TRUE)
  cat("\n")
}
saveRDS(loocv_marginal_dists_pca, here("results", "brains", "pca_margins.Rda"))

loocv_marginal_dists_pca_cor <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "PCA-cor")
  loocv_marginal_dists_pca_cor[[i]] <- compute_all_distances(input_distributions, 
                                                             diag(1, d), A = W, verbose = TRUE, keep_projections = TRUE)
  cat("\n")
}
saveRDS(loocv_marginal_dists_pca_cor, here("results", "brains", "pca_cor_margins.Rda"))

loocv_marginal_dists_cholesky <- vector("list", length = length(input_distributions))
for (i in seq_along(input_distributions)) {
  message(paste0("Fold ", i, " of ", length(input_distributions), "..."))
  W <- get_trafo(input_distributions[-i], type = "Cholesky")
  loocv_marginal_dists_cholesky[[i]] <- compute_all_distances(input_distributions, 
                                                              diag(1, d), A = W, verbose = TRUE, keep_projections = TRUE)
  cat("\n")
}
saveRDS(loocv_marginal_dists_cholesky, here("results", "brains", "cholesky_margins.Rda"))

# Note: The *marginal* Wasserstein distances needn't be the same regardless of whitening procedure!

# Predicting the 10% trimmed mean ----------------------------------------------
# We'll predict the 10% trimmed mean in LOOCV using the different distances:
# - raw sliced Wasserstein
# - raw marginal Wasserstein
# - scaled sliced Wasserstein
# - whitened sliced Wasserstein
# - whitened marginal Wasserstein (different versions)
sw_dist <- here("results", "brains", "no_trafo_sliced.Rda") |> readRDS()
marginal_dist <- here("results", "brains", "no_trafo_margins.Rda") |> readRDS()
loocv_sw_dists_scaled <- here("results", "brains", "scaled_sliced.Rda") |> readRDS()
loocv_sw_dists_whitened <- here("results", "brains", "whitened_sliced.Rda") |> readRDS()
loocv_marginal_dists_zca <- here("results", "brains", "zca_margins.Rda") |> readRDS()
loocv_marginal_dists_zca_cor <- here("results", "brains", "zca_cor_margins.Rda") |> readRDS()
loocv_marginal_dists_pca <- here("results", "brains", "pca_margins.Rda") |> readRDS()
loocv_marginal_dists_pca_cor <- here("results", "brains", "pca_cor_margins.Rda") |> readRDS()
loocv_marginal_dists_cholesky <- here("results", "brains", "cholesky_margins.Rda") |> readRDS()
loocv_marginal_dists_scaled <- here("results", "brains", "scaled_margins.Rda") |> readRDS()

outcomes <- sapply(output_distributions, mean, trim = 0.1)
predictions <- matrix(NA, nrow = 20, ncol = 11)
pb <- txtProgressBar(0, 20, style = 3)
for (i in seq_along(input_distributions)) {
  training_idx <- seq_along(input_distributions)[-i]
  outcome_train <- outcomes[training_idx]
  outcome_test <- outcomes[i]
  
  predictions[i, 1] <- fit_gpr(sw_dist, training_idx, i, 
                                 outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 2] <- fit_gpr_multiple(marginal_dist, training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 3] <- fit_gpr(loocv_sw_dists_scaled[[i]], training_idx, i,
                                 outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 4] <- fit_gpr(loocv_sw_dists_whitened[[i]], training_idx, i,
                                 outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 5] <- fit_gpr_multiple(loocv_marginal_dists_zca[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 6] <- fit_gpr_multiple(loocv_marginal_dists_zca_cor[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 7] <- fit_gpr_multiple(loocv_marginal_dists_pca[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 8] <- fit_gpr_multiple(loocv_marginal_dists_pca_cor[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 9] <- fit_gpr_multiple(loocv_marginal_dists_cholesky[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 10] <- mean(outcome_train)
  predictions[i, 11] <- fit_gpr_multiple(loocv_marginal_dists_scaled[[i]], training_idx, i,
                                           outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  setTxtProgressBar(pb, i)
}
cat("\n")

df_predictions <- cbind(outcomes, predictions) |> 
  as.data.frame()
colnames(df_predictions) <- c("outcome", "sliced Wasserstein: raw", "marginal Wasserstein: raw",
                              "sliced Wasserstein: scaled", "sliced Wasserstein: whitened",
                              "marginal Wasserstein: ZCA-whitened", "marginal Wasserstein: ZCA-cor-whitened",
                              "marginal Wasserstein: PCA-whitened", "marginal Wasserstein: PCA-cor-whitened",
                              "marginal Wasserstein: Cholesky-whitened", "baseline (training mean imputation)",
                              "marginal Wasserstein: scaled")
df_predictions <- df_predictions |> 
  pivot_longer(cols = -outcome,
               names_to = "Method", values_to = "Prediction")
rmses <- df_predictions |> 
  group_by(Method) |> 
  summarise(RMSE = sqrt(mean((outcome - Prediction)^2)) |> round(3))
rmses
df_predictions <- df_predictions |> 
  filter(Method != "marginal Wasserstein: scaled") |> 
  filter(Method != "baseline (training mean imputation)") |> 
  left_join(rmses) |> 
  mutate(Method = paste0(Method, "\n(RMSE: ", RMSE, ")"))
df_predictions |> 
  ggplot(aes(x = outcome, y = Prediction)) +
  geom_point(shape = 1) + 
  facet_wrap(facets = vars(reorder(Method, -RMSE))) +
  xlab("10% trimmed mean of number of neighbours") +
  ylab("out-of-fold prediction") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))
ggsave(here("figures", "brain_trimmed_mean.pdf"), width = 1.3*6.3, height = 1.3*4.7)

entropy <- function(x) {
  counts <- table(x)
  p <- counts / sum(counts)
  -sum(p * log(p))
}
outcomes <- sapply(output_distributions, entropy)
predictions <- matrix(NA, nrow = 20, ncol = 11)
pb <- txtProgressBar(0, 20, style = 3)
for (i in seq_along(input_distributions)) {
  training_idx <- seq_along(input_distributions)[-i]
  outcome_train <- outcomes[training_idx]
  outcome_test <- outcomes[i]
  
  predictions[i, 1] <- fit_gpr(sw_dist, training_idx, i, 
                                 outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 2] <- fit_gpr_multiple(marginal_dist, training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 3] <- fit_gpr(loocv_sw_dists_scaled[[i]], training_idx, i,
                                 outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 4] <- fit_gpr(loocv_sw_dists_whitened[[i]], training_idx, i,
                                 outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 5] <- fit_gpr_multiple(loocv_marginal_dists_zca[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 6] <- fit_gpr_multiple(loocv_marginal_dists_zca_cor[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 7] <- fit_gpr_multiple(loocv_marginal_dists_pca[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 8] <- fit_gpr_multiple(loocv_marginal_dists_pca_cor[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 9] <- fit_gpr_multiple(loocv_marginal_dists_cholesky[[i]], training_idx, i,
                                          outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  predictions[i, 10] <- mean(outcome_train)
  predictions[i, 11] <- fit_gpr_multiple(loocv_marginal_dists_scaled[[i]], training_idx, i,
                                           outcome_train, outcome_test, verbose = FALSE, runs = 50)$test_predictions
  setTxtProgressBar(pb, i)
}
cat("\n")

df_predictions <- cbind(outcomes, predictions) |> 
  as.data.frame()
colnames(df_predictions) <- c("outcome", "sliced Wasserstein: raw", "marginal Wasserstein: raw",
                              "sliced Wasserstein: scaled", "sliced Wasserstein: whitened",
                              "marginal Wasserstein: ZCA-whitened", "marginal Wasserstein: ZCA-cor-whitened",
                              "marginal Wasserstein: PCA-whitened", "marginal Wasserstein: PCA-cor-whitened",
                              "marginal Wasserstein: Cholesky-whitened", "baseline (training mean imputation)",
                              "marginal Wasserstein: scaled")
df_predictions <- df_predictions |> 
  pivot_longer(cols = -outcome,
               names_to = "Method", values_to = "Prediction")
rmses <- df_predictions |> 
  group_by(Method) |> 
  summarise(RMSE = sqrt(mean((outcome - Prediction)^2)) |> round(3))
rmses
df_predictions <- df_predictions |> 
  filter(Method != "marginal Wasserstein: scaled") |> 
  filter(Method != "baseline (training mean imputation)") |> 
  left_join(rmses) |> 
  mutate(Method = paste0(Method, "\n(RMSE: ", RMSE, ")"))
df_predictions |> 
  ggplot(aes(x = outcome, y = Prediction)) +
  geom_point(shape = 1) + 
  facet_wrap(facets = vars(reorder(Method, -RMSE))) +
  xlab("entropy of number of neighbours") +
  ylab("out-of-fold prediction") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))
ggsave(here("figures", "brain_entropy.pdf"), width = 1.3*6.3, height = 1.3*4.7)
