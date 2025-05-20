# Makes DGP1 (VAR) coefficients
make_VAR_coefs<-function(N,rhos=c(0.2,0.15,0.1,0.05)){
  VAR_lags=length(rhos)
  VAR_coefficients<-array(dim=c(N,N,VAR_lags))
  for(l in 1:VAR_lags){
    for(i in 1:N){
      for(j in 1:N){
        if(abs(i-j)>=N/2){
          VAR_coefficients[i,j,l]=0
        }else{
          VAR_coefficients[i,j,l]=rhos[l]^(abs(i-j)+1)
        }
      }
    }
  }
  return(VAR_coefficients)
}

# Makes DGP2 (VAR) coefficients
make_VAR_coefs2<-function(N,rhos=c(0.2,0.15,0.1,0.05)){
  VAR_lags=length(rhos)
  VAR_coefficients<-array(dim=c(N,N,VAR_lags))
  for(l in 1:VAR_lags){
    for(i in 1:N){
      for(j in 1:N){
        if(abs(i-j)>=N/2){
          VAR_coefficients[i,j,l]=0
        }else{
          VAR_coefficients[i,j,l]=rhos[l]^(abs(i-j)+1)
        }
      }
    }
  }
  VAR_coefficients[,,2]= -VAR_coefficients[,,2]
  VAR_coefficients[,,4]= -VAR_coefficients[,,4]
  return(VAR_coefficients)
}

# Recusively obatins VMA coefficients from DGP1,2 coefficients
# Meaning the true IRF
irf_from_VAR<-function(VAR_coefficients, hmax=10){
  # Inverting the VAR polynomial to the VMA ---------------------------------
  # Find math at slide 19 of https://kevinkotze.github.io/ts-7-slide/#19
  #B will be the VMA coefficients, A will be the VAR coefficients, with 0
  #matrices for lags past 4
  #the indexing of B will start at 0
  N<-dim(VAR_coefficients)[1]
  VAR_lags<-dim(VAR_coefficients)[3]
  B<-array(0,dim=c(N,N,hmax+1))
  B[,,1]<-diag(N)
  
  A<-array(0,dim=c(N,N,hmax))
  A[,,1:VAR_lags]<-VAR_coefficients
  
  for(i in 2:(hmax+1)){
    sum=matrix(0,N,N)
    for(j in 1:(i-1)){
      sum<-sum+B[,,i-j]%*%A[,,j]
    }
    B[,,i]<-sum
  }
  irf_1to1<-B[1,1,]
  return(irf_1to1)
}


# Simulation function 1
simulate_VAR_DGP <- function(
    VAR_coefficients, # Changed input: now takes coefficients directly
    N,
    T_sample,
    Sigma_u = NULL,
    intercept = NULL,
    burnin = 100
) {
  
  # --- Input Validation ---
  if (!is.array(VAR_coefficients) || length(dim(VAR_coefficients)) != 3) {
    stop("VAR_coefficients must be a 3D array (N x N x p).")
  }
  if (dim(VAR_coefficients)[1] != N || dim(VAR_coefficients)[2] != N) {
    stop("First two dimensions of VAR_coefficients must be N x N.")
  }
  p <- dim(VAR_coefficients)[3] # Number of lags determined from input
  
  # --- 2. Define Intercept ---
  if (is.null(intercept)) {
    intercept <- matrix(0.0, nrow = N, ncol = 1) # Default zero intercept
  } else {
    if (length(intercept) != N || !is.numeric(intercept)) {
      stop("Intercept must be a numeric vector of length N.")
    }
    intercept <- matrix(intercept, nrow = N, ncol = 1)
  }
  
  # --- 3. Define Error Covariance Matrix Sigma_u ---
  if (is.null(Sigma_u)) {
    Sigma_u_internal <- diag(N) # Default: u_t ~ N(0, I)
  } else {
    if (!is.matrix(Sigma_u) || !all(dim(Sigma_u) == c(N, N))) {
      stop("Sigma_u must be an N x N matrix.")
    }
    # Check if Sigma_u is positive definite (required for mvrnorm)
    eigenvalues <- eigen(Sigma_u, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigenvalues <= 1e-8)) { # Use a small tolerance
      warning("Sigma_u might not be positive definite.")
    }
    Sigma_u_internal <- Sigma_u
  }
  
  # --- 4. Initialize Simulation ---
  # Total length needed = p initial values + burnin + T_sample observations
  T_total <- T_sample + burnin + p # Use the burnin argument here
  # Matrix to store simulated data: N rows (variables) x T_total columns (time)
  z_sim <- matrix(0.0, nrow = N, ncol = T_total)
  # Matrix to store errors: N rows x T_total columns
  u_sim <- matrix(0.0, nrow = N, ncol = T_total) # Initialize error storage
  
  # --- 5. Run Simulation Loop ---
  # Generate all errors at once for efficiency
  # Note: all_errors will be T_total x N
  all_errors <- MASS::mvrnorm(n = T_total, mu = rep(0, N),
                              Sigma = Sigma_u_internal)
  
  # Loop starts after the initial p zero columns, runs for burnin + T_sample steps
  for (t in (p + 1):T_total) {
    # Calculate deterministic part: c + sum(A_k * z_{t-k})
    z_deterministic <- intercept
    for (l in 1:p) {
      # Access columns t-l which contain previous values
      z_deterministic <- z_deterministic + VAR_coefficients[,, l] %*% z_sim[, t - l, drop = FALSE]
    }
    
    # Get the pre-generated error u_t for this time step
    # Note: all_errors is T_total x N, so we take row t and transpose it to N x 1
    u_t <- matrix(all_errors[t, ], ncol = 1)
    u_sim[, t] <- u_t # Store the error
    
    # Calculate z_t and store it in the current column t
    z_sim[, t] <- z_deterministic + u_t
  }
  
  # --- 6. Select Final Sample (Discard Burn-in) ---
  # Get the T_sample observations generated "after" the initial p zeros and burnin period
  # These are the last T_sample columns
  start_col <- burnin + p + 1 # First column of the final sample
  end_col <- T_total         # Last column of the final sample
  z_intermediate <- z_sim[, start_col:end_col, drop = FALSE] # Still N x T_sample here
  u_intermediate <- u_sim[, start_col:end_col, drop = FALSE] # Get corresponding errors
  
  # --- 7. Transpose and Return Result ---
  # Transpose the results to be T_sample x N
  X_final <- t(z_intermediate)
  u_final <- t(u_intermediate)
  
  # Returns a list containing the T_sample x N data matrix and error matrix
  return(list(X = X_final, u = u_final))
}

# Estimate SVAR and Calculate Structural IRF with Bootstrap CIs
# Fits a VAR, identifies using Cholesky, calculates the specified structural IRF,
# and computes bootstrap confidence intervals.
estimate_svar_irf <- function(sim_data,
                              p = 4,
                              hmax = 10,
                              impulse_var = 1,
                              response_var = 1,
                              ci_level = 0.95,
                              n_boot = 200
) {
  
  # --- 1. Prepare Data ---
  if (nrow(sim_data) < ncol(sim_data)) {
    data_for_var <- t(sim_data)
  } else {
    data_for_var <- sim_data # Assume it's already T rows x N cols
  }
  # Ensure variable names are generic if none exist, helps irf function access
  if (is.null(colnames(data_for_var))) {
    colnames(data_for_var) <- paste0("V", 1:ncol(data_for_var))
    # Adjust impulse/response if they were integers to match new names
    if(is.numeric(impulse_var)) impulse_var <- paste0("V", impulse_var)
    if(is.numeric(response_var)) response_var <- paste0("V", response_var)
  }
  
  # --- 2. Estimate Reduced-Form VAR ---
  var_model <- tryCatch({
    vars::VAR(y = data_for_var,
              p = p,
              type = "none") # Assuming no intercept/trend in DGP based on image
  }, error = function(e) {
    warning(paste("VAR estimation failed:", e$message))
    return(NULL)
  })
  
  # If VAR estimation failed, return NULL
  if (is.null(var_model)) {
    return(NULL)
  }
  
  # --- 3. Calculate Structural IRF and Bootstrap CIs ---
  irf_results <- tryCatch({
    vars::irf(var_model,
              impulse = impulse_var,   # impulse variable
              response = response_var, # response variable
              n.ahead = hmax,          # Horizon
              ortho = TRUE,            # Get orthogonalized (structural) shocks via Cholesky
              ci = ci_level,           # CI
              boot = TRUE,
              runs = n_boot,           # Bootstrap runs
              seed = NULL)             # RNG/seed is specified before function runs
  }, error = function(e) {
    warning(paste("IRF calculation failed:", e$message))
    return(NULL)
  })
  
  # If IRF calculation failed, return NULL
  if (is.null(irf_results)) {
    # Clean up VAR model object if needed (optional)
    # rm(var_model)
    return(NULL)
  }
  
  # --- 4. Extract Results ---
  # The irf function returns a list containing irf, lower, upper for the specified variables
  # The results include horizon 0
  point_estimate_vec <- irf_results$irf[[1]][, 1]
  ci_lower_vec       <- irf_results$Lower[[1]][, 1]
  ci_upper_vec       <- irf_results$Upper[[1]][, 1]
  
  # --- 5. Return ---
  return(list(
    point_estimate = point_estimate_vec,
    ci_lower = ci_lower_vec,
    ci_upper = ci_upper_vec
  ))
}



analyze_svar_failure_boundary <- function(Ts, Ns, p_fixed = 4, rhos = c(0.2, 0.15, 0.1, 0.05), hmax = 10, estimate_function = estimate_svar_irf) {
  # Ensure Ns and Ts are sorted
  if (is.unsorted(Ns)) {
    warning("Input vector 'Ns' was not sorted. Sorting ascendingly.")
    Ns <- sort(Ns)
  }
  if (is.unsorted(Ts)) {
    warning("Input vector 'Ts' was not sorted. Sorting ascendingly.")
    Ts <- sort(Ts)
  }
  max_n_tested <- max(Ns) # Store the highest N value
  n_count <- length(Ns)   # Number of N values to check
  
  results_list <- list()
  result_counter <- 1
  start_n_index <- 1 # Start checking from the first N for the first T
  
  # --- Outer Loop for T ---
  for (t_val in Ts) {
    print(paste("Running T =", t_val, "- Starting check from N =", Ns[start_n_index]))
    failure_occurred_for_T <- FALSE # Flag for this T
    
    # --- Add Success Records for Ns skipped by optimization ---
    if (start_n_index > 1) {
      print(paste("  Recording assumed success for N =", paste(Ns[1:(start_n_index - 1)], collapse=", ")))
      for (n_success_idx in 1:(start_n_index - 1)) {
        results_list[[result_counter]] <- data.frame(N = Ns[n_success_idx], T = t_val, Status = 1)
        result_counter <- result_counter + 1
      }
    }
    
    # --- Inner Loop for N (starting from optimized index) ---
    # Only loop if start_n_index is within bounds
    if (start_n_index <= n_count) {
      for (n_idx in start_n_index:n_count) {
        n_val <- Ns[n_idx]
        case_name <- paste0("N", n_val, "_T", t_val)
        print(paste("  Checking N =", n_val))
        
        # --- Generate Coefficients for this N ---
        VAR_coeffs <- tryCatch({
          make_VAR_coefs(N = n_val, rhos = rhos)
        }, error = function(e) {
          warning(paste("Case", case_name, "Failed to generate VAR coefficients:", e$message))
          return(NULL)
        })
        if(is.null(VAR_coeffs)) {
          print(paste("    -> Failure (Coefficient Generation Error)"))
          results_list[[result_counter]] <- data.frame(N = n_val, T = t_val, Status = 0)
          result_counter <- result_counter + 1
          failure_occurred_for_T <- TRUE
          start_n_index <- n_idx # Next T starts checking from this failed N index
          break # Stop checking higher N for this T
        }
        # Check if p_fixed is valid
        current_p <- p_fixed
        if (p_fixed > dim(VAR_coeffs)[3]) {
          warning(paste("Case", case_name, "Fixed p =", p_fixed, "exceeds VAR order from rhos. Adjusting p."))
          current_p <- dim(VAR_coeffs)[3]
          print(paste("    -> Adjusted p to", current_p))
        }
        
        # --- Simulate Data ---
        sim_output <- tryCatch({
          simulate_VAR_DGP(
            VAR_coefficients = VAR_coeffs,
            N = n_val,
            T_sample = t_val,
            ... # Pass additional arguments like Sigma_u, intercept, burnin
          )
        }, error = function(e) {
          warning(paste("Case", case_name, "DGP generation failed:", e$message))
          return(NULL) # Treat DGP failure as estimation failure
        })
        
        if (is.null(sim_output)) {
          print(paste("    -> Failure (DGP Error)"))
          results_list[[result_counter]] <- data.frame(N = n_val, T = t_val, Status = 0)
          result_counter <- result_counter + 1
          failure_occurred_for_T <- TRUE
          start_n_index <- n_idx # Next T starts checking from this failed N index
          break # Stop checking higher N for this T
        }
        
        # --- Estimate SVAR IRF ---
        svar_result <- tryCatch({
          estimate_function( # Use the specified estimation function
            sim_data = sim_output$X,
            p = current_p, # Use potentially adjusted p
            hmax = hmax,
            impulse_var = 1,
            response_var = 1,
            n_boot = 1 # Only need 1 bootstrap run to check estimation success
          )
        }, error = function(e) {
          # Catch errors within the estimation function itself
          warning(paste("Case", case_name, "Estimation error:", e$message))
          return(NULL)
        })
        
        # --- Check Success and Store/Break ---
        if (is.null(svar_result)) {
          # First failure for this T (at this N index or higher)
          print(paste("    -> Failure (Estimation Error / Cholesky Failed)"))
          results_list[[result_counter]] <- data.frame(N = n_val, T = t_val, Status = 0)
          result_counter <- result_counter + 1
          failure_occurred_for_T <- TRUE
          start_n_index <- n_idx # Next T starts checking from this failed N index
          break # Stop checking higher N for this T
        } else {
          # Success for this N, T
          print(paste("    -> Success"))
          results_list[[result_counter]] <- data.frame(N = n_val, T = t_val, Status = 1)
          result_counter <- result_counter + 1
          # Update start_n_index ONLY if this was the last N and it succeeded
          if (n_idx == n_count) {
            start_n_index <- n_idx # If last N succeeded, next T starts checking last N again
          }
          # Continue to the next N for this T
        }
      } # --- End N Loop ---
    } else {
      # This happens if start_n_index was already beyond the last N index
      print(paste("  Skipping N checks as previous T failed at or before N =", Ns[start_n_index-1]))
      # We still need to record the failure boundary for this T at the same N as the previous T
      results_list[[result_counter]] <- data.frame(N = Ns[start_n_index-1], T = t_val, Status = 0)
      result_counter <- result_counter + 1
      failure_occurred_for_T <- TRUE # Mark as failed for consistency, although loop did not run
    }
    
    
    # --- Handle case where the inner N loop finished without failure ---
    # This implies all N from start_n_index to max_n_tested succeeded
    # The start_n_index for the next T should remain at max_n_tested index
    if (!failure_occurred_for_T) {
      print(paste("  All checked N values succeeded for T =", t_val))
      # No failure occurred, so the next T should start checking from the last N index again
      start_n_index <- n_count
    }
    
  } # --- End T Loop ---
  
  # Combine results into a single data frame
  if (length(results_list) == 0) {
    warning("No results were generated.")
    return(data.frame(N = integer(), T = integer(), Status = integer()))
  }
  results_df <- dplyr::bind_rows(results_list)
  
  # Ensure unique entries per N, T combination, keeping the minimum status
  # Handles the boundary point (Status=0) if it occurred
  results_df <- results_df %>%
    group_by(N, T) %>%
    summarise(Status = min(Status), .groups = 'drop')
  
  return(results_df)
}


plot_svar_failure_boundary <- function(boundary_data) {
  # --- Input Checks ---
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'dplyr' are required. Please install them.", call. = FALSE)
  }
  if (!all(c("N", "T", "Status") %in% names(boundary_data))) {
    stop("Input data frame must contain columns 'N', 'T', and 'Status'.", call. = FALSE)
  }
  if (nrow(boundary_data) == 0) {
    warning("Input data frame is empty. Returning an empty plot.")
    return(ggplot() + labs(title="SVAR Failure Boundary (No Data)") + theme_minimal())
  }
  
  # Separate successful and first failure points
  success_points <- boundary_data %>% dplyr::filter(Status == 1)
  failure_points <- boundary_data %>%
    dplyr::filter(Status == 0) %>%
    dplyr::arrange(T) # Arrange by T for line plotting
  
  # --- Create Plot ---
  p <- ggplot() +
    # Plot successful points
    geom_point(data = success_points, aes(x = T, y = N),
               color = "darkgreen", size = 2, shape = 1) + # e.g., green circles
    
    # Plot the first failure points
    geom_point(data = failure_points, aes(x = T, y = N),
               color = "red", size = 3, shape = 4) + # e.g., red 'x'
    
    # Draw line connecting the first failure points
    geom_line(data = failure_points, aes(x = T, y = N),
              color = "red", linetype = "dashed", size = 1) +
    
    # Use all tested T and N values for axis breaks
    scale_x_continuous(breaks = unique(boundary_data$T)) +
    scale_y_continuous(breaks = unique(boundary_data$N)) +
    
    labs(
      title = "SVAR Estimation Failure Boundary",
      subtitle = "Red 'x' marks the first N where estimation failed for a given T.",
      x = "Sample Size (T)",
      y = "Number of Variables (N)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Add annotation if there are T values with no failures
  max_n_tested <- max(boundary_data$N)
  successful_t <- boundary_data %>%
    group_by(T) %>%
    filter(all(Status == 1)) %>%
    distinct(T)
  
  if(nrow(successful_t) > 0) {
    p <- p + annotate("text", x = successful_t$T, y = max_n_tested + (max_n_tested * 0.05), # Position above max N
                      label = "✔", color = "darkgreen", size = 5, vjust = 0)
  }
  
  return(p)
}


sims_to_coverage<-function(path, M, Ns, Ts, PIs, hmax, partial, switching) {
  coverage<-array(0,dim = c(length(Ns), length(Ts), length(PIs), hmax+1), dimnames = list(N=Ns, T=Ts, PI=PIs, horizon=0:hmax))
  if(partial==TRUE){
    prefix<-"partial"
  }else{
    prefix<-"regular"
  }
  
  if(switching==TRUE){
    suffix<-"_switching_signs"
  }else{
    suffix<-NULL
  }
  
  for(n in 1:length(Ns)){
    current_Ts_for_N <- Ts[[n]]
    for(t in 1:length(current_Ts_for_N)){
      for(pi in 1:length(PIs)){
        N=Ns[n]
        T_=current_Ts_for_N[t]
        PIconstant=PIs[pi]
        if(switching==TRUE){
          VAR_coef<-make_VAR_coefs2(N)
        }else{
          VAR_coef<-make_VAR_coefs(N)
        }
        true_irf_1to1<-irf_from_VAR(VAR_coef)
        sim<-readRDS(file=paste0(path,"/", prefix,"_N",N,"_T",T_,"_PI",PIconstant,suffix,".RDS"))
        for(h in 1:(hmax+1)){
          for(i in 1:M){
            if(sim$intervals[h,1,i]<=true_irf_1to1[h] && true_irf_1to1[h]<=sim$intervals[h,3,i]){
              coverage[n,t,pi,h]<- coverage[n,t,pi,h]+1/M
            }
          }
        }
      }
    }
  }
  return(coverage)
}

sims_to_width<-function(path, M, Ns, Ts, PIs, hmax, partial, switching){
  width<-array(0,dim = c(length(Ns), length(Ts), length(PIs), hmax+1), dimnames = list(N=Ns, T=Ts, PI=PIs, horizon=0:hmax))
  if(partial==TRUE){
    prefix<-"partial"
  }else{
    prefix<-"regular"
  }
  
  if(switching==TRUE){
    suffix<-"_switching_signs"
  }else{
    suffix<-NULL
  }
  
  for(n in 1:length(Ns)){
    current_Ts_for_N <- Ts[[n]]
    for(t in 1:length(current_Ts_for_N)){
      for(pi in 1:length(PIs)){
        N=Ns[n]
        T_=current_Ts_for_N[t]
        PIconstant=PIs[pi]
        sim<-readRDS(file=paste0(path,"/", prefix,"_N",N,"_T",T_,"_PI",PIconstant,suffix,".RDS"))
        for(h in 1:(hmax+1)){
          for(i in 1:M){
            width[n,t,pi,h]<- width[n,t,pi,h]+(sim$intervals[h,3,i]-sim$intervals[h,1,i])/M
          }
        }
      }
    }
  }
  return(width)
}


sims_to_bias <- function(path, M, Ns, Ts, PIs, hmax, partial, switching) {
  # --- Initialize Output Array ---
  bias <- array(0,
                dim = c(length(Ns), length(Ts), length(PIs), hmax + 1),
                dimnames = list(N = Ns, T = Ts, PI = PIs, horizon = 0:hmax))
  
  # --- Determine File Name Prefix/Suffix ---
  if (partial == TRUE) {
    prefix <- "partial"
    } else {
    prefix <- "regular"
  }
  
  if (switching == TRUE) {
    suffix <- "_switching_signs"
    } else {
    suffix <- NULL
  }
  
  # --- Loop Through Simulation Parameters ---
  for(n in 1:length(Ns)){
    current_Ts_for_N <- Ts[[n]]
    for(t in 1:length(current_Ts_for_N)){
      for(pi in 1:length(PIs)){
        N=Ns[n]
        T_=current_Ts_for_N[t]
        PIconstant <- PIs[pi]
        
        if(switching==TRUE){
          VAR_coef<-make_VAR_coefs2(N)
        }else{
          VAR_coef<-make_VAR_coefs(N)
        }
        true_irf_1to1<-irf_from_VAR(VAR_coef)
        
        # Construct file path
        file_path <- paste0(path, "/", prefix, "_N", N, "_T", T_, "_PI", PIconstant, suffix, ".RDS")
        
        # Read simulation results file
        sim <- tryCatch({
          readRDS(file = file_path)
        }, error = function(e) {
          warning(paste("Could not read file:", file_path, "-", e$message))
          return(NULL) # Return NULL if file reading fails
        })
        
        # Skip this iteration if file reading failed
        if (is.null(sim) || !("intervals" %in% names(sim)) || dim(sim$intervals)[3] != M) {
          warning(paste("Skipping combination N=", N, "T=", T_, "PI=", PIconstant, "due to missing/invalid file or incorrect M dimension."))
          # Fill bias with NA for this combination
          bias[n, t, pi, ] <- NA_real_
          next
        }
        
        # Loop through horizons
        for (h in 0:hmax) {
          horizon_index <- h + 1
          
          # Calculate sum of point estimates across replications
          # Assumes point estimate is the 2nd element in the interval dimension
          sum_of_estimates <- sum(sim$intervals[horizon_index, 2, 1:M], na.rm = TRUE) # Sum estimates for horizon h across M reps
          
          # Calculate average estimate
          num_valid_reps <- sum(!is.na(sim$intervals[horizon_index, 2, 1:M]))
          if (num_valid_reps > 0) {
            average_estimate <- sum_of_estimates / num_valid_reps
            # Calculate and store average bias
            bias[n, t, pi, horizon_index] <- average_estimate - true_irf_1to1[horizon_index]
          } else {
            # If no valid estimates for this horizon, bias is NA
            bias[n, t, pi, horizon_index] <- NA_real_
            warning(paste("No valid estimates found for N=", N, "T=", T_, "PI=", PIconstant, "Horizon=", h))
          }
        } # End horizon loop
      } # End PI loop
    } # End T loop
  } # End N loop
  
  return(bias)
}

# simulations to RMSE
sims_to_rmse <- function(path, M, Ns, Ts, PIs, hmax, partial, switching) {
  # --- Initialize Output Array ---
  rmse <- array(0,
                dim = c(length(Ns), length(Ts), length(PIs), hmax + 1),
                dimnames = list(N = Ns, T = Ts, PI = PIs, horizon = 0:hmax))
  
  # --- Determine File Name Prefix/Suffix ---
  if (partial == TRUE) {
    prefix <- "partial"
  } else {
    prefix <- "regular"
  }
  
  if (switching == TRUE) {
    suffix <- "_switching_signs"
  } else {
    suffix <- NULL # Or "" if you prefer empty string
  }
  
  # --- Loop Through Simulation Parameters ---
  for(n in 1:length(Ns)){
    current_Ts_for_N <- Ts[[n]]
    for(t in 1:length(current_Ts_for_N)){
      for(pi in 1:length(PIs)){
        N=Ns[n]
        T_=current_Ts_for_N[t]
        PIconstant <- PIs[pi]
        
        if(switching==TRUE){
          VAR_coef<-make_VAR_coefs2(N)
        }else{
          VAR_coef<-make_VAR_coefs(N)
        }
        true_irf_1to1<-irf_from_VAR(VAR_coef)
        
        # Constructing file path
        file_path <- paste0(path, "/", prefix, "_N", N, "_T", T_, "_PI", PIconstant, suffix, ".RDS")
        
        # Read simulation results file
        sim <- tryCatch({
          readRDS(file = file_path)
        }, error = function(e) {
          warning(paste("Could not read file:", file_path, "-", e$message))
          return(NULL) # Return NULL if file reading fails
        })
        
        # Skip this iteration if file reading failed
        if (is.null(sim) || !("intervals" %in% names(sim)) || dim(sim$intervals)[3] != M) {
          warning(paste("Skipping combination N=", N, "T=", T_, "PI=", PIconstant, "due to missing/invalid file or incorrect M dimension."))
          # Fill bias with NA for this combination
          bias[n, t, pi, ] <- NA_real_
          next
        }
        
        # Loop through horizons (0 to hmax)
        for (h in 0:hmax) {
          horizon_index <- h + 1 # R index is 1-based
          
          # Calculate sum of point estimates across replications
          rmse[n, t, pi, horizon_index] <- sqrt(mean((sim$intervals[horizon_index, 2, 1:M] - true_irf_1to1[horizon_index])^2, na.rm = TRUE)) # Sum estimates for horizon h across M reps
        } # End horizon loop
      } # End PI loop
    } # End T loop
  } # End N loop
  
  return(rmse)
}