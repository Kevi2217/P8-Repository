########################## NON-PARALLEL APPLICATION ############################
set.seed(1)

final_results <- list()
# DGP1 SVAR estimation
for(n in 1:length(Ns)){
  # Put these function here since they only change per N
  VAR_coefs <- make_VAR_coefs(Ns[n])
  irf_1to1 <- irf_from_VAR(VAR_coefs, hmax = hmax)
  
  for(t in 1:length(Ts)){
    # Define case name
    case_name <- paste0("N", Ns[n], "_T", Ts[t])
    
    # Lists for storing results in replication loop
    svar_point_estimates_list <- vector("list", M)
    svar_lower_ci_list <- vector("list", M)
    svar_upper_ci_list <- vector("list", M)
    
    for(rep in 1:M){
      # --- Simulate Data ---
      sim_output <- tryCatch({
        simulate_VAR_DGP(
          VAR_coefficients = VAR_coefs,
          N = Ns[n],
          T_sample = Ts[t],
          Sigma_u = NULL,
          intercept = NULL,
          burnin = 100
        )
      }, error = function(e) {
        warning(paste("Replication", rep, "DGP generation failed:", e$message))
        return(NULL) # Return NULL if DGP fails
      })
      
      # --- Estimate SVAR IRF ---
      svar_analysis_result <- tryCatch({
        estimate_svar_irf(sim_data = sim_output$X, # Use the simulated data matrix
                          p = 4,          # Assuming lags = 4
                          hmax = hmax,
                          impulse_var = 1, 
                          response_var = 1,
                          ci_level = 0.95,
                          n_boot = 200      # Used for bootstrapping
        )
      }, error = function(e) {
        warning(paste("Replication", rep, "SVAR estimation failed:", e$message))
        return(NULL) # Return NULL if estimation fails
      })
      
      
      # --- Store Results ---
      if (!is.null(svar_analysis_result)) {
        svar_point_estimates_list[[rep]] <- svar_analysis_result$point_estimate
        svar_lower_ci_list[[rep]]      <- svar_analysis_result$ci_lower
        svar_upper_ci_list[[rep]]      <- svar_analysis_result$ci_upper
      } else {
        # Store NAs if estimation failed
        svar_point_estimates_list[[rep]] <- rep(NA_real_, hmax + 1)
        svar_lower_ci_list[[rep]]      <- rep(NA_real_, hmax + 1)
        svar_upper_ci_list[[rep]]      <- rep(NA_real_, hmax + 1)
      }
      
    } # --- End Replication Loop ---
    
    # --- Combine Results into Matrices ---
    # cbind combines the list elements (vectors) into columns of a matrix
    # Resulting matrices will have (hmax+1) rows and M columns
    svar_points_matrix <- do.call(cbind, svar_point_estimates_list)
    svar_lowers_matrix <- do.call(cbind, svar_lower_ci_list)
    svar_uppers_matrix <- do.call(cbind, svar_upper_ci_list)
    
    # --- AGGREGATE RESULTS ---
    # List to store aggregated results for this N, T case
    svar_results_agg <- list()
    
    
    # Calculating estimate/upper/lower
    svar_results_agg$estimate <- rowMeans(svar_points_matrix, na.rm = TRUE)
    svar_results_agg$upper <- rowMeans(svar_uppers_matrix, na.rm = TRUE)
    svar_results_agg$lower <- rowMeans(svar_lowers_matrix, na.rm = TRUE)
    
    # Calculate Bias: Mean estimate across replications minus true value
    # Assumes svar_points_matrix is (hmax+1) x M and irf_1to1 is length (hmax+1)
    svar_results_agg$bias <- rowMeans(svar_points_matrix, na.rm = TRUE) - irf_1to1
    
    # Calculate RMSE: Root Mean Squared Error
    svar_errors_sq <- (sweep(svar_points_matrix, 1, irf_1to1, "-"))^2
    svar_results_agg$rmse <- sqrt(rowMeans(svar_errors_sq, na.rm = TRUE))
    
    # Calculate Coverage Rate: Proportion of replications where CI contains the true value
    # Assumes svar_lowers_matrix and svar_uppers_matrix are (hmax+1) x M
    svar_results_agg$coverage <- rowMeans((irf_1to1 >= svar_lowers_matrix) & (irf_1to1 <= svar_uppers_matrix), na.rm = TRUE)
    
    # Calculate Average Width: Mean width of confidence intervals
    svar_results_agg$avg_width <- rowMeans(svar_uppers_matrix - svar_lowers_matrix, na.rm = TRUE)
    
    # --- Store Aggregated Results ---
    # Store the results for this case (N, T) in the main results list
    final_results[[case_name]] <- list(SVAR = svar_results_agg)
    
    print(paste("Finished aggregation for case:", case_name))
    
  } # --- End T Loop ---
}   # --- End N Loop ---