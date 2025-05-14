############################# PARALLEL SETUP ###################################
# --- Define Function for a Single Replication ---
# This function will be executed in parallel for each replication index (1 to M)
run_application_replications <- function(rep_index) {
  # Note: rep_index is passed by parLapply but not explicitly used here
  # More advanced seeding could use rep_index
  
  # --- Simulate Data ---
  sim_output <- tryCatch({
    simulate_VAR_DGP(
      VAR_coefficients = VAR_coeffs,
      N = current_N,
      T_sample = current_T,
      Sigma_u = NULL, # Or pass specific Sigma_u if needed
      intercept = NULL,
      burnin = 100
    )
  }, error = function(e) {
    # Return NULL or an indicator of failure if DGP fails
    # warning(paste("Rep", rep_index, "DGP generation failed:", e$message)) # Optional warning
    return(NULL)
  })
  
  # If DGP failed, return NULL for this replication
  if (is.null(sim_output)) {
    return(NULL)
  }
  
  # --- Estimate SVAR IRF ---
  svar_analysis_result <- tryCatch({
    estimate_svar_irf(sim_data = sim_output$X,
                      p = p_lags,
                      hmax = hmax,
                      impulse_var = 1,
                      response_var = 1,
                      ci_level = 0.95,
                      n_boot = 200 # Bootstrap runs within each replication
    )
  }, error = function(e) {
    # Return NULL or an indicator of failure if estimation fails
    # warning(paste("Rep", rep_index, "SVAR estimation failed:", e$message)) # Optional warning
    return(NULL)
  })
  
  # Return the result (list with point_estimate, ci_lower, ci_upper or NULL)
  return(svar_analysis_result)
}
############################# PARALLEL SETUP ###################################

for(n in 1:length(Ns)){
  current_Ts_for_N <- Ts[[n]]
  current_N <- Ns[n]
  # --- Generate Coefficients and True IRF for this N ---
  # These are calculated once per N, outside the parallel part
  VAR_coeffs <- make_VAR_coefs2(current_N)
  true_irf_11 <- irf_from_VAR(VAR_coeffs, hmax = hmax)
  p_lags <- dim(VAR_coeffs)[3] # Determine lag order p
  
  for(t in 1:length(current_Ts_for_N)){
    current_T <- current_Ts_for_N[t]
    case_name <- paste0("N", current_N, "_T", current_T)
    print(paste("Time: ", Sys.time(), "Running case:", case_name, "with", M, "replications on", threads, "threads..."))
    
    # --- Setup Parallel Backend ---
    cl <- makeCluster(threads)
    
    # --- Export Necessary Objects and Load Packages on Workers ---
    # Export variables needed within the parallel task
    clusterExport(cl, varlist = c("simulate_VAR_DGP", "estimate_svar_irf",
                                  "VAR_coeffs", "current_N", "current_T",
                                  "p_lags", "hmax"),
                  envir = environment())
    # Load the 'vars' package on each worker node
    clusterEvalQ(cl, library(vars))
    
    # --- RUNNING REPLICATIONS in Parallel ---
    # parLapply runs the function 'run_application_replications' for each number from 1 to M
    # distributing the tasks across the cluster 'cl'
    parallel_results_list <- parLapply(cl, 1:M, run_application_replications)
    
    # --- Stop Parallel Backend ---
    stopCluster(cl)
    print(paste("  Parallel execution finished for:", case_name))
    
    # --- Process Results from Parallel Execution ---
    # Initialize lists to store extracted results
    svar_point_estimates_list <- vector("list", M)
    svar_lower_ci_list      <- vector("list", M)
    svar_upper_ci_list      <- vector("list", M)
    
    # Extract results, handling NULLs from failed replications
    for (rep in 1:M) {
      result <- parallel_results_list[[rep]]
      if (!is.null(result) && is.list(result) &&
          !is.null(result$point_estimate) && !is.null(result$ci_lower) && !is.null(result$ci_upper)) {
        # Check if lengths match hmax+1, otherwise fill with NA
        expected_len <- hmax + 1
        if(length(result$point_estimate) == expected_len &&
           length(result$ci_lower) == expected_len &&
           length(result$ci_upper) == expected_len) {
          
          svar_point_estimates_list[[rep]] <- result$point_estimate
          svar_lower_ci_list[[rep]]      <- result$ci_lower
          svar_upper_ci_list[[rep]]      <- result$ci_upper
        } else {
          warning(paste("Case", case_name, "Rep", rep, "result has unexpected length. Storing NAs."))
          svar_point_estimates_list[[rep]] <- rep(NA_real_, expected_len)
          svar_lower_ci_list[[rep]]      <- rep(NA_real_, expected_len)
          svar_upper_ci_list[[rep]]      <- rep(NA_real_, expected_len)
        }
      } else {
        # Store NAs if replication failed (result is NULL or malformed)
        svar_point_estimates_list[[rep]] <- rep(NA_real_, hmax + 1)
        svar_lower_ci_list[[rep]]      <- rep(NA_real_, hmax + 1)
        svar_upper_ci_list[[rep]]      <- rep(NA_real_, hmax + 1)
      }
    } # --- End rep loop ---
    
    # --- Combine Results into Matrices ---
    svar_points_matrix <- do.call(cbind, svar_point_estimates_list)
    svar_lowers_matrix <- do.call(cbind, svar_lower_ci_list)
    svar_uppers_matrix <- do.call(cbind, svar_upper_ci_list)
    
    # --- AGGREGATE RESULTS ---
    svar_results_agg <- list()
    
    
    # Calculating estimate/upper/lower
    svar_results_agg$estimate <- rowMeans(svar_points_matrix, na.rm = TRUE)
    svar_results_agg$upper <- rowMeans(svar_uppers_matrix, na.rm = TRUE)
    svar_results_agg$lower <- rowMeans(svar_lowers_matrix, na.rm = TRUE)
    
    # Calculate Bias: Mean estimate across replications minus true value
    # Assumes svar_points_matrix is (hmax+1) x M and irf_1to1 is length (hmax+1)
    svar_results_agg$bias <- rowMeans(svar_points_matrix, na.rm = TRUE) - true_irf_11
    
    # Calculate RMSE: Root Mean Squared Error
    svar_errors_sq <- (sweep(svar_points_matrix, 1, true_irf_11, "-"))^2
    svar_results_agg$rmse <- sqrt(rowMeans(svar_errors_sq, na.rm = TRUE))
    
    # Calculate Coverage Rate: Proportion of replications where CI contains the true value
    # Assumes svar_lowers_matrix and svar_uppers_matrix are (hmax+1) x M
    svar_results_agg$coverage <- rowMeans((true_irf_11 >= svar_lowers_matrix) & (true_irf_11 <= svar_uppers_matrix), na.rm = TRUE)
    
    # Calculate Average Width: Mean width of confidence intervals
    svar_results_agg$avg_width <- rowMeans(svar_uppers_matrix - svar_lowers_matrix, na.rm = TRUE)
    
    # --- Store Aggregated Results ---
    # Store the results for this case (N, T) in the main results list
    final_results_switching[[case_name]] <- list(SVAR = svar_results_agg)
    
  } # --- End T Loop ---
}   # --- End N Loop ---