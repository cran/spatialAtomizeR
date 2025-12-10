## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = FALSE  # Set to TRUE if you want examples to run during build
)

## ----install, eval = FALSE----------------------------------------------------
# # Install from GitHub
# devtools::install_github("bellayqian/spatialAtomizeR")

## ----libraries----------------------------------------------------------------
# library(spatialAtomizeR)
# library(nimble)  # Required for ABRM models

## ----simulate-----------------------------------------------------------------
# sim_data <- simulate_misaligned_data(
#   seed = 42,
#   res1 = c(5, 5),    # Y grid resolution (coarser)
#   res2 = c(10, 10),  # X grid resolution (finer)
# 
#   # Distribution specifications
#   dist_covariates_x = c('normal', 'poisson', 'binomial'),
#   dist_covariates_y = c('normal', 'poisson', 'binomial'),
#   dist_y = 'poisson',
# 
#   # Intercepts (critical parameters)
#   x_intercepts = c(4, -1, -1),
#   y_intercepts = c(4, -1, -1),
#   beta0_y = -1,
# 
#   # Spatial correlation
#   x_correlation = 0.5,
#   y_correlation = 0.5,
# 
#   # True effect sizes
#   beta_x = c(-0.03, 0.1, -0.2),
#   beta_y = c(0.03, -0.1, 0.2)
# )

## ----model--------------------------------------------------------------------
# model_code <- get_abrm_model()

## ----run_abrm-----------------------------------------------------------------
# results <- run_abrm(
#   sim_data = sim_data,
#   model_code = model_code,
# 
#   # Map distribution indices to positions
#   norm_idx_x = 1,   # 'normal' is 1st in dist_covariates_x
#   pois_idx_x = 2,   # 'poisson' is 2nd
#   binom_idx_x = 3,  # 'binomial' is 3rd
#   norm_idx_y = 1,
#   pois_idx_y = 2,
#   binom_idx_y = 3,
# 
#   # Outcome distribution: 1=normal, 2=poisson, 3=binomial
#   dist_y = 2,
# 
#   # MCMC parameters
#   niter = 50000,
#   nburnin = 30000,
#   nchains = 2
# )

## ----results------------------------------------------------------------------
# # Parameter estimates
# print(results$parameter_estimates)
# 
# # Convergence diagnostics
# print(results$mcmc_results$convergence)

## ----advanced-----------------------------------------------------------------
# # Prepare spatial structure
# bookkeeping <- prepare_spatial_bookkeeping(sim_data)
# 
# # Create adjacency matrices
# adjacency <- prepare_adjacency_matrices(
#   bookkeeping$gridy_yorder,
#   bookkeeping$gridx_xorder
# )
# 
# # Prepare NIMBLE inputs
# nimble_inputs <- prepare_nimble_inputs(
#   bookkeeping, adjacency, sim_data,
#   norm_idx_x = 1, pois_idx_x = 2, binom_idx_x = 3,
#   norm_idx_y = 1, pois_idx_y = 2, binom_idx_y = 3,
#   dist_y = 2
# )
# 
# # Run NIMBLE MCMC
# sim_metadata <- list(
#   sim_number = 1,
#   x_correlation = 0.5,
#   y_correlation = 0.5
# )
# 
# mcmc_results <- run_nimble_model(
#   constants = nimble_inputs$constants,
#   data = nimble_inputs$data,
#   inits = nimble_inputs$inits,
#   sim_metadata = sim_metadata,
#   model_code = model_code,
#   niter = 50000,
#   nburnin = 30000,
#   nchains = 2,
#   thin = 10,
#   save_plots = TRUE
# )

## ----sensitivity--------------------------------------------------------------
# # Define base parameters
# base_params <- list(
#   res1 = c(5, 5),
#   res2 = c(10, 10),
#   dist_covariates_x = c('normal','poisson','binomial'),
#   dist_covariates_y = c('normal','poisson','binomial'),
#   dist_y = 'poisson',
#   x_intercepts = c(4, -1, -1),
#   y_intercepts = c(4, -1, -1),
#   beta0_y = -1,
#   beta_x = c(-0.03, 0.1, -0.2),
#   beta_y = c(0.03, -0.1, 0.2)
# )
# 
# # Run sensitivity analysis
# sensitivity_results <- run_sensitivity_analysis(
#   correlation_grid = c(0.2, 0.6),
#   n_sims_per_setting = 3,
#   base_params = base_params,
#   model_code = model_code,
#   base_seed = 123
# )
# 
# # View summary
# print(sensitivity_results$summary_by_correlation)

## ----citation, eval = FALSE---------------------------------------------------
# citation("spatialAtomizeR")

