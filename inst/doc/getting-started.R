## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = FALSE  # Set to TRUE to run code during build (requires external data for Example 2)
)

## ----install, eval = FALSE----------------------------------------------------
# # Install from CRAN
# install.packages("spatialAtomizeR")
# 
# # Or install development version from GitHub
# devtools::install_github("bellayqian/spatialAtomizeR")

## ----libraries, message = FALSE-----------------------------------------------
# library(spatialAtomizeR)
# library(nimble)  # Required for ABRM models

## ----simulate-----------------------------------------------------------------
# sim_data <- simulate_misaligned_data(
#   seed = 42,
# 
#   # Distribution specifications for covariates
#   dist_covariates_x = c('normal', 'poisson', 'binomial'),
#   dist_covariates_y = c('normal', 'poisson', 'binomial'),
#   dist_y = 'poisson',  # Outcome distribution
# 
#   # Intercepts for data generation (REQUIRED)
#   x_intercepts = c(4, -1, -1),      # One per X covariate
#   y_intercepts = c(4, -1, -1),      # One per Y covariate
#   beta0_y = -1,                     # Outcome model intercept
# 
#   # Spatial correlation parameters
#   x_correlation = 0.5,  # Correlation between X covariates
#   y_correlation = 0.5,  # Correlation between Y covariates
# 
#   # True effect sizes for outcome model
#   beta_x = c(-0.03, 0.1, -0.2),    # Effects of X covariates
#   beta_y = c(0.03, -0.1, 0.2)      # Effects of Y covariates
# )

## ----examine_data-------------------------------------------------------------
# # Check the class
# class(sim_data)
# 
# # Print method shows basic information
# print(sim_data)
# 
# # Summary method shows more details
# summary(sim_data)

## ----model--------------------------------------------------------------------
# model_code <- get_abrm_model()

## ----run_abrm-----------------------------------------------------------------
# abrm_results <- run_abrm(
#   gridx = sim_data$gridx,
#   gridy = sim_data$gridy,
#   atoms = sim_data$atoms,
#   model_code = model_code,
#   true_params = sim_data$true_params,  # Optional: for validation
# 
#   # Map distribution indices to positions in dist_covariates_x/y
#   norm_idx_x = 1,   # 'normal' is 1st in dist_covariates_x
#   pois_idx_x = 2,   # 'poisson' is 2nd
#   binom_idx_x = 3,  # 'binomial' is 3rd
#   norm_idx_y = 1,   # Same for Y covariates
#   pois_idx_y = 2,
#   binom_idx_y = 3,
# 
#   # Outcome distribution: 1=normal, 2=poisson, 3=binomial
#   dist_y = 2,
# 
#   # MCMC parameters
#   niter = 50000,    # Total iterations per chain
#   nburnin = 30000,  # Burn-in iterations
#   nchains = 2       # Number of chains
# )

## ----examine_results----------------------------------------------------------
# # Check the class
# class(abrm_results)
# 
# # Print method shows summary statistics
# print(abrm_results)
# 
# # Summary method shows detailed parameter estimates
# summary(abrm_results)
# 
# # Plot method shows MCMC diagnostics (if available)
# plot(abrm_results)

## ----access_results-----------------------------------------------------------
# # Parameter estimates table
# abrm_results$parameter_estimates
# 
# # MCMC convergence diagnostics
# abrm_results$mcmc_results$convergence
# 
# # Full MCMC samples
# # abrm_results$mcmc_results$samples

## ----load_spatial_packages, message = FALSE-----------------------------------
# library(tigris)    # For US Census shapefiles
# library(sf)        # Spatial features
# library(sp)        # Spatial objects
# library(spdep)     # Spatial dependence
# library(raster)    # Spatial data manipulation
# library(dplyr)     # Data manipulation
# library(ggplot2)   # Plotting

## ----load_utah_data-----------------------------------------------------------
# set.seed(500)
# 
# # Load Utah counties (Y-grid)
# cnty <- counties(state = 'UT')
# gridy <- cnty %>%
#   rename(ID = GEOID) %>%
#   mutate(ID = as.numeric(ID))  # ID must be numeric
# 
# # Load Utah school districts (X-grid)
# scd <- school_districts(state = 'UT')
# gridx <- scd %>%
#   rename(ID = GEOID) %>%
#   mutate(ID = as.numeric(ID))

## ----plot_misalignment, fig.width = 7, fig.height = 5-------------------------
# ggplot() +
#   geom_sf(data = gridx, fill = NA, color = "orange", linewidth = 1.2) +
#   geom_sf(data = gridy, fill = NA, color = "blue",
#           linetype = 'dashed', linewidth = 0.5) +
#   labs(title = "Spatial Misalignment in Utah",
#        subtitle = "Orange: School Districts | Blue: Counties") +
#   theme_void()

## ----create_atoms-------------------------------------------------------------
# # Intersect the two grids to create atoms
# atoms <- raster::intersect(as(gridy, 'Spatial'), as(gridx, 'Spatial'))
# atoms <- sf::st_as_sf(atoms)
# 
# # Rename ID variables to match expected names
# atoms <- atoms %>%
#   rename(ID_y = ID_1,    # Y-grid (county) ID
#          ID_x = ID_2)     # X-grid (school district) ID

## ----load_population----------------------------------------------------------
# # NOTE: You must download LandScan data separately
# # Available at: https://landscan.ornl.gov/
# # This example assumes the file is in your working directory
# 
# pop_raster <- raster("landscan-global-2022.tif")
# 
# # Extract population for each atom
# pop_atoms <- raster::extract(
#   pop_raster,
#   st_transform(atoms, crs(pop_raster)),
#   fun = sum,
#   na.rm = TRUE
# )
# 
# atoms$population <- pop_atoms

## ----generate_synthetic_data--------------------------------------------------
# # Create atom-level spatial adjacency matrix
# W_a <- nb2mat(
#   poly2nb(as(atoms, "Spatial"), queen = TRUE),
#   style = "B",
#   zero.policy = TRUE
# )
# 
# # Generate spatial random effects
# atoms <- atoms %>%
#   mutate(
#     re_x = gen_correlated_spat(W = W_a, n_vars = 1, correlation = 1, rho = 0.8),
#     re_z = gen_correlated_spat(W = W_a, n_vars = 1, correlation = 1, rho = 0.8),
#     re_y = gen_correlated_spat(W = W_a, n_vars = 1, correlation = 1, rho = 0.8)
#   )
# 
# # Generate atom-level covariates and outcome
# atoms <- atoms %>%
#   mutate(
#     # X-grid covariate (Poisson)
#     covariate_x_1 = rpois(
#       n = nrow(atoms),
#       lambda = population * exp(-1 + re_x)
#     ),
#     # Y-grid covariate (Normal)
#     covariate_y_1 = rnorm(
#       n = nrow(atoms),
#       mean = population * (3 + re_z),
#       sd = 1 * population
#     )
#   ) %>%
#   mutate(
#     # Outcome (Poisson)
#     y = rpois(
#       n = nrow(atoms),
#       lambda = population * exp(
#         -5 +
#         1 * (covariate_x_1 / population) -
#         0.3 * (covariate_y_1 / population) +
#         re_y
#       )
#     )
#   )

## ----aggregate_data-----------------------------------------------------------
# # Aggregate X covariates to X-grid
# gridx[["covariate_x_1"]] <- sapply(gridx$ID, function(j) {
#   atom_indices <- which(atoms$ID_x == j)
#   sum(atoms[["covariate_x_1"]][atom_indices])
# })
# 
# # Aggregate Y covariates to Y-grid
# gridy[["covariate_y_1"]] <- sapply(gridy$ID, function(j) {
#   atom_indices <- which(atoms$ID_y == j)
#   sum(atoms[["covariate_y_1"]][atom_indices])
# })
# 
# # Aggregate outcome to Y-grid
# gridy$y <- sapply(gridy$ID, function(j) {
#   atom_indices <- which(atoms$ID_y == j)
#   sum(atoms$y[atom_indices])
# })

## ----run_abrm_utah------------------------------------------------------------
# model_code <- get_abrm_model()
# 
# mcmc_results <- run_abrm(
#   gridx = gridx,
#   gridy = gridy,
#   atoms = atoms,
#   model_code = model_code,
# 
#   # Specify covariate distributions
#   pois_idx_x = 1,  # X covariate is Poisson
#   norm_idx_y = 1,  # Y covariate is Normal
#   dist_y = 2,      # Outcome is Poisson
# 
#   # MCMC settings (increase for real analysis)
#   niter = 300000,
#   nburnin = 200000,
#   nchains = 2
# )
# 
# # View results
# summary(mcmc_results)

## ----comparison, eval = FALSE-------------------------------------------------
# comparison <- run_both_methods(
#   sim_data = sim_data,
#   sim_metadata = list(sim_number = 1, x_correlation = 0.5, y_correlation = 0.5),
#   model_code = model_code,
#   nimble_params = list(niter = 1000, nburnin = 500, thin = 2, nchains = 2),
#   output_dir = tempdir(),
#   norm_idx_x = c(1),
#   pois_idx_x = c(2),
#   binom_idx_x = c(3),
#   norm_idx_y = c(1),
#   pois_idx_y = c(2),
#   binom_idx_y = c(3),
#   dist_y = 2,
#   outcome_type = 'poisson'
# )
# 
# # View comparison results
# print(comparison)
# summary(comparison)

## ----citation, eval = FALSE---------------------------------------------------
# citation("spatialAtomizeR")

