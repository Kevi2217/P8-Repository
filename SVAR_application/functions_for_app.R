estimate_fixed_post_lasso_svar_structured <- function(structured_data, p_lags, 
                                                target_CPI_name = "CPIAUCSL", 
                                                target_IP_name = "INDPRO",
                                                n_ahead_irf = 48, 
                                                n_bootstrap_irf = 199, 
                                                ci_level = 0.95,
                                                num_cores = NULL,
                                                n_bootstrap_paths_to_return = 0,
                                                fixed_nfolds_cv = 10) {
  
  # --- 0. Data Preparation ---
  slow_vars_df <- as.data.frame(structured_data$slow_data)
  fast_vars_df <- as.data.frame(structured_data$fast_data)
  
  if (!is.data.frame(structured_data$FFR) && !is.matrix(structured_data$FFR)) {
    ffr_df <- data.frame(FFR = as.vector(structured_data$FFR))
  } else if (ncol(as.data.frame(structured_data$FFR)) > 1) { 
    stop("FFR data (structured_data$FFR) should have only one column.")
  } else {
    ffr_df <- as.data.frame(structured_data$FFR)
    if(is.null(colnames(ffr_df)) || colnames(ffr_df)[1] == "" || length(colnames(ffr_df))==0) {
      colnames(ffr_df) <- "FFR" 
    }
  }
  ffr_col_name_actual <- colnames(ffr_df)[1]
  
  Y_data_ordered <- cbind(slow_vars_df, ffr_df, fast_vars_df)
  Y_data_ordered <- as.matrix(Y_data_ordered) 
  
  K_vars_total <- ncol(Y_data_ordered)
  var_names <- colnames(Y_data_ordered)
  
  if(is.null(var_names) || any(var_names == "") || length(var_names) != K_vars_total) {
    num_slow <- if(!is.null(slow_vars_df)) ncol(slow_vars_df) else 0
    num_fast <- if(!is.null(fast_vars_df)) ncol(fast_vars_df) else 0
    default_slow_names <- if(num_slow > 0) paste0("Slow", 1:num_slow) else character(0)
    default_ffr_name <- ffr_col_name_actual 
    default_fast_names <- if(num_fast > 0) paste0("Fast", 1:num_fast) else character(0)
    var_names <- c(default_slow_names, default_ffr_name, default_fast_names)
    if(length(var_names) != K_vars_total) { var_names <- paste0("V", 1:K_vars_total) }
    colnames(Y_data_ordered) <- var_names
  }
  
  cat("Total variables in ordered system:", K_vars_total, "\n")
  cat("Observations provided:", nrow(Y_data_ordered), "\n")
  cat("Variable order for VAR and Cholesky:", paste(var_names, collapse=", "), "\n")
  
  ffr_shock_idx <- which(var_names == ffr_col_name_actual)
  if(length(ffr_shock_idx) == 0) stop(paste("FFR column '", ffr_col_name_actual, "' not found.", sep=""))
  
  target_response_vars_map <- list(FFR = ffr_col_name_actual, CPI = target_CPI_name, IP = target_IP_name)
  target_response_indices <- sapply(target_response_vars_map, function(name) {
    idx <- which(var_names == name)
    if(length(idx) == 0) {
      warning(paste("Target response '", name, "' not found. IRF will be NA.", sep=""))
      return(NA_integer_)
    }
    return(idx[1]) 
  })
  names(target_response_indices) <- names(target_response_vars_map)
  
  # --- 1. Lagged Regressors for VAR ---
  create_VAR_lags_glmnet <- function(data, p) { 
    K <- ncol(data)
    data_colnames <- colnames(data) 
    if(is.null(data_colnames) || length(data_colnames) != K) data_colnames <- paste0("V", 1:K)
    T_orig <- nrow(data)
    if (T_orig <= p) {
      stop(paste0("Obs (", T_orig, ") must be > lags (", p, ")."))
    }
    X_lagged <- matrix(NA, nrow = T_orig - p, ncol = K * p)
    Y_target <- data[(p + 1):T_orig, , drop = FALSE]
    col_names_X <- character(K * p) 
    for (lag_val in 1:p) {
      lag_data_block <- data[(p + 1 - lag_val):(T_orig - lag_val), , drop = FALSE]
      start_col_idx <- (lag_val - 1) * K + 1; end_col_idx <- lag_val * K
      X_lagged[, start_col_idx:end_col_idx] <- lag_data_block
      col_names_X[start_col_idx:end_col_idx] <- paste0(data_colnames, ".L", lag_val)
    }
    colnames(X_lagged) <- col_names_X
    return(list(Y_target = Y_target, X_lagged = X_lagged))
  }
  
  prepared_data <- create_VAR_lags_glmnet(Y_data_ordered, p_lags)
  Y_target_matrix <- prepared_data$Y_target
  X_predictors_matrix <- prepared_data$X_lagged
  T_eff <- nrow(Y_target_matrix) 
  
  current_nfolds <- ifelse(T_eff >= fixed_nfolds_cv, fixed_nfolds_cv, max(3, T_eff %/% 3))
  if (T_eff < 3) current_nfolds <- T_eff 
  
  # --- 2. Equation-by-Equation Post-LASSO OLS Estimation (Original Data) ---
  num_predictors <- ncol(X_predictors_matrix)
  B_hat_post_lasso_orig <- matrix(0, nrow = K_vars_total, ncol = num_predictors) 
  intercepts_post_lasso_orig <- numeric(K_vars_total)
  selected_vars_indices_orig_list <- list() # To store selected var indices for each equation
  
  cat("\nStarting Post-LASSO OLS estimation for each VAR equation (Original Data)...\n")
  for (k in 1:K_vars_total) {
    y_k <- Y_target_matrix[, k]
    cv_fit <- cv.glmnet(X_predictors_matrix, y_k, alpha = 1, nfolds = current_nfolds, standardize = TRUE)
    optimal_lambda <- cv_fit$lambda.min 
    lasso_model <- glmnet(X_predictors_matrix, y_k, alpha = 1, lambda = optimal_lambda, standardize = TRUE)
    lasso_coeffs_k <- coef(lasso_model)
    selected_vars_indices_orig_list[[k]] <- which(lasso_coeffs_k[-1, 1] != 0) 
    
    if (length(selected_vars_indices_orig_list[[k]]) > 0) {
      X_selected_k <- X_predictors_matrix[, selected_vars_indices_orig_list[[k]], drop = FALSE]
      if (qr(X_selected_k)$rank < ncol(X_selected_k)) {
        warning(paste("Rank deficiency (Original Data), eq:", var_names[k], ". Using LASSO coeffs."))
        intercepts_post_lasso_orig[k] <- lasso_coeffs_k[1,1]
        B_hat_post_lasso_orig[k,] <- as.vector(lasso_coeffs_k[-1,1]) # Fallback to LASSO
      } else {
        tryCatch({
          ols_fit_k <- lm(y_k ~ X_selected_k)
          intercepts_post_lasso_orig[k] <- coef(ols_fit_k)[1]
          B_hat_post_lasso_orig[k, selected_vars_indices_orig_list[[k]]] <- coef(ols_fit_k)[-1]
        }, error = function(e) {
          warning(paste("OLS failed (Original Data), eq:", var_names[k], ". Using LASSO. Error:", e$message))
          intercepts_post_lasso_orig[k] <- lasso_coeffs_k[1,1]
          B_hat_post_lasso_orig[k,] <- as.vector(lasso_coeffs_k[-1,1]) 
        })
      }
    } else {
      intercepts_post_lasso_orig[k] <- mean(y_k) 
      cat("  No variables selected by LASSO for equation:", var_names[k], "(Original Data).\n")
    }
  }
  names(selected_vars_indices_orig_list) <- var_names
  rownames(B_hat_post_lasso_orig) <- var_names; colnames(B_hat_post_lasso_orig) <- colnames(X_predictors_matrix)
  
  estimated_A_coeffs_post_lasso_orig <- list()
  for (i in 1:p_lags) {
    A_i_matrix <- B_hat_post_lasso_orig[, ((i - 1) * K_vars_total + 1):(i * K_vars_total)]
    colnames(A_i_matrix) <- var_names; rownames(A_i_matrix) <- var_names 
    estimated_A_coeffs_post_lasso_orig[[i]] <- A_i_matrix
  }
  names(estimated_A_coeffs_post_lasso_orig) <- paste0("A", 1:p_lags)
  
  # --- 3. Calculate Reduced-Form Residuals (from Original Data Post-LASSO OLS) ---
  U_hat_post_lasso_orig <- matrix(NA, nrow = T_eff, ncol = K_vars_total); colnames(U_hat_post_lasso_orig) <- var_names
  for (t_loop in 1:T_eff) { 
    y_lagged_vec_t <- X_predictors_matrix[t_loop, ] 
    fitted_val_t <- intercepts_post_lasso_orig + B_hat_post_lasso_orig %*% y_lagged_vec_t
    U_hat_post_lasso_orig[t_loop, ] <- Y_target_matrix[t_loop, ] - fitted_val_t
  }
  
  # --- 4. Estimate Residual Covariance Matrix ---
  Sigma_U_hat_post_lasso_orig <- (t(U_hat_post_lasso_orig) %*% U_hat_post_lasso_orig) / (T_eff - 1) 
  
  # --- 5. Structural Identification ---
  eigen_Sigma_orig <- eigen(Sigma_U_hat_post_lasso_orig, symmetric = TRUE, only.values = TRUE)$values
  if(any(eigen_Sigma_orig <= 1e-8)){ warning("Sigma_U (Original) not positive definite.") }
  B0_hat_orig <- tryCatch({ t(chol(Sigma_U_hat_post_lasso_orig)) }, error = function(e) { stop(paste("Cholesky failed (Original):", e$message)) })
  colnames(B0_hat_orig) <- var_names; rownames(B0_hat_orig) <- var_names 
  cat("\nStructural identification via Cholesky complete (Original Data). FFR ('", ffr_col_name_actual, "') is var #", ffr_shock_idx, "\n")
  
  # --- 6. Compute Structural IRFs (Point Estimates), Normalize, and Cumulate ---
  compute_VMA_coeffs_from_VAR <- function(A_coeffs_list_func, p_lags_var, K_vars_var, n_ahead_vma) { 
    Psi_matrices <- list(); Psi_matrices[[1]] <- diag(K_vars_var) 
    A_full <- list()
    for(i_vma in 1:p_lags_var) { 
      if(i_vma <= length(A_coeffs_list_func)) { A_full[[i_vma]] <- A_coeffs_list_func[[i_vma]] } 
      else { A_full[[i_vma]] <- matrix(0, nrow=K_vars_var, ncol=K_vars_var) }
    }
    for (h_vma in 1:n_ahead_vma) { 
      Psi_h <- matrix(0, nrow = K_vars_var, ncol = K_vars_var)
      for (j_vma in 1:min(h_vma, p_lags_var)) { Psi_h <- Psi_h + (Psi_matrices[[h_vma - j_vma + 1]] %*% A_full[[j_vma]]) }
      Psi_matrices[[h_vma + 1]] <- Psi_h
    }
    return(Psi_matrices)
  }
  
  Psi_matrices_post_lasso_orig <- compute_VMA_coeffs_from_VAR(estimated_A_coeffs_post_lasso_orig, p_lags, K_vars_total, n_ahead_irf)
  structural_IRFs_raw_orig <- array(NA, dim = c(K_vars_total, K_vars_total, n_ahead_irf + 1))
  for (h_loop in 0:n_ahead_irf) { structural_IRFs_raw_orig[,,h_loop+1] <- Psi_matrices_post_lasso_orig[[h_loop+1]] %*% B0_hat_orig }
  
  ffr_on_ffr_h0_orig <- structural_IRFs_raw_orig[ffr_shock_idx, ffr_shock_idx, 1] 
  if (abs(ffr_on_ffr_h0_orig) < 1e-9) {
    warning("Initial FFR response (Original) is near zero. Normalization might be unstable.")
    ffr_on_ffr_h0_scale_factor_orig <- 1 
  } else { ffr_on_ffr_h0_scale_factor_orig <- ffr_on_ffr_h0_orig }
  structural_IRFs_normalized_ffr_shock_orig <- structural_IRFs_raw_orig[, ffr_shock_idx, ] / ffr_on_ffr_h0_scale_factor_orig
  
  point_estimates_targeted_cumulative <- list(); point_estimates_targeted_noncumulative <- list()
  for(resp_target_name in names(target_response_vars_map)){
    actual_resp_idx <- target_response_indices[resp_target_name]
    if(!is.na(actual_resp_idx)){
      irf_vec_normalized <- structural_IRFs_normalized_ffr_shock_orig[actual_resp_idx, ]
      point_estimates_targeted_noncumulative[[resp_target_name]] <- irf_vec_normalized 
      if(resp_target_name == "CPI" || resp_target_name == "IP"){
        point_estimates_targeted_cumulative[[resp_target_name]] <- cumsum(irf_vec_normalized) 
      } else { point_estimates_targeted_cumulative[[resp_target_name]] <- irf_vec_normalized }
    } else {
      point_estimates_targeted_noncumulative[[resp_target_name]] <- rep(NA, n_ahead_irf + 1)
      point_estimates_targeted_cumulative[[resp_target_name]] <- rep(NA, n_ahead_irf + 1)
    }
  }
  cat("Point estimates of IRFs computed (Original Data).\n")
  
  # --- 7. Bootstrap Confidence Intervals for IRFs (Fixed Model Bootstrap) ---
  if (is.null(num_cores)) { num_cores <- detectCores() - 1; if (num_cores < 1) num_cores <- 1 }
  cat("\nSetting up parallel backend with", num_cores, "cores for FIXED MODEL bootstrap...\n")
  cl <- makeCluster(num_cores); registerDoParallel(cl)
  cat("Starting parallel bootstrap for IRF CIs (", n_bootstrap_irf, "replications)...\n")
  
  vars_to_export <- c("U_hat_post_lasso_orig", "T_eff", "p_lags", "K_vars_total", "var_names", 
                      "Y_data_ordered", "intercepts_post_lasso_orig", "B_hat_post_lasso_orig", 
                      "selected_vars_indices_orig_list", # Export indices of originally selected vars
                      "num_predictors", "ffr_shock_idx", "n_ahead_irf",
                      "target_response_indices", "create_VAR_lags_glmnet", "compute_VMA_coeffs_from_VAR")
  
  bootstrap_results_list <- foreach(
    b_loop = 1:n_bootstrap_irf, .packages = c("glmnet"), .export = vars_to_export,
    .combine = 'list', .multicombine = TRUE, .errorhandling = 'pass' 
  ) %dopar% {
    U_hat_centered_worker <- scale(U_hat_post_lasso_orig, center = TRUE, scale = FALSE) 
    boot_indices_worker <- sample(1:T_eff, T_eff, replace = TRUE)
    U_boot_worker <- U_hat_centered_worker[boot_indices_worker, , drop = FALSE]
    Y_data_boot_worker <- matrix(NA, nrow = T_eff + p_lags, ncol = K_vars_total)
    colnames(Y_data_boot_worker) <- var_names 
    Y_data_boot_worker[1:p_lags, ] <- Y_data_ordered[1:p_lags, ] 
    
    for (t_boot_gen in (p_lags + 1):(T_eff + p_lags)) { 
      x_pred_t_boot_worker <- numeric(K_vars_total * p_lags)
      for(lag_b in 1:p_lags){ x_pred_t_boot_worker[((lag_b-1)*K_vars_total + 1) : (lag_b*K_vars_total)] <- Y_data_boot_worker[t_boot_gen-lag_b, ] }
      Y_fitted_t_boot_worker <- intercepts_post_lasso_orig + B_hat_post_lasso_orig %*% x_pred_t_boot_worker 
      Y_data_boot_worker[t_boot_gen, ] <- Y_fitted_t_boot_worker + U_boot_worker[t_boot_gen - p_lags, ]
    }
    Y_data_boot_for_var_worker <- Y_data_boot_worker
    
    prepared_data_boot_worker <- create_VAR_lags_glmnet(Y_data_boot_for_var_worker, p_lags)
    Y_target_matrix_boot_worker <- prepared_data_boot_worker$Y_target
    X_predictors_matrix_boot_worker <- prepared_data_boot_worker$X_lagged 
    
    # Re-estimate OLS VAR on bootstrap sample USING FIXED SET OF PREDICTORS
    B_hat_ols_boot_worker <- matrix(0, nrow = K_vars_total, ncol = num_predictors)
    intercepts_ols_boot_worker <- numeric(K_vars_total)
    
    for (k_boot in 1:K_vars_total) {
      y_k_boot_worker <- Y_target_matrix_boot_worker[, k_boot]
      # Use the selected variables from the ORIGINAL LASSO fit for this equation
      selected_vars_for_this_eq_boot <- selected_vars_indices_orig_list[[k_boot]] 
      
      if (length(selected_vars_for_this_eq_boot) > 0) {
        X_selected_k_boot <- X_predictors_matrix_boot_worker[, selected_vars_for_this_eq_boot, drop = FALSE]
        if (qr(X_selected_k_boot)$rank < ncol(X_selected_k_boot)) { # Check rank deficiency
   
          intercepts_ols_boot_worker[k_boot] <- NA 
          next 
        }
        tryCatch({
          ols_fit_k_boot <- lm(y_k_boot_worker ~ X_selected_k_boot)
          intercepts_ols_boot_worker[k_boot] <- coef(ols_fit_k_boot)[1]
          B_hat_ols_boot_worker[k_boot, selected_vars_for_this_eq_boot] <- coef(ols_fit_k_boot)[-1]
        }, error = function(e_boot) { 
          intercepts_ols_boot_worker[k_boot] <- NA # Mark as failed
        })
      } else { # No variables were selected in the original fit for this equation
        intercepts_ols_boot_worker[k_boot] <- mean(y_k_boot_worker)
      }
    }
    if(any(is.na(intercepts_ols_boot_worker))) { return(list(error="OLS_boot_failed_for_some_eq"))} # Skip if any OLS failed
    
    estimated_A_coeffs_ols_boot_worker <- list()
    for (i_boot in 1:p_lags) { estimated_A_coeffs_ols_boot_worker[[i_boot]] <- B_hat_ols_boot_worker[, ((i_boot - 1) * K_vars_total + 1):(i_boot * K_vars_total)] }
    
    U_hat_ols_boot_worker <- matrix(NA, nrow = T_eff, ncol = K_vars_total)
    for (t_boot_res in 1:T_eff) { 
      y_lagged_vec_t_boot_worker <- X_predictors_matrix_boot_worker[t_boot_res, ] 
      fitted_val_t_boot <- intercepts_ols_boot_worker + B_hat_ols_boot_worker %*% y_lagged_vec_t_boot_worker
      U_hat_ols_boot_worker[t_boot_res, ] <- Y_target_matrix_boot_worker[t_boot_res, ] - fitted_val_t_boot
    }
    Sigma_U_hat_ols_boot_worker <- (t(U_hat_ols_boot_worker) %*% U_hat_ols_boot_worker) / (T_eff - 1)
    
    eigen_Sigma_boot_worker <- eigen(Sigma_U_hat_ols_boot_worker, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigen_Sigma_boot_worker <= 1e-8)){ return(list(error = "Sigma_U_boot_OLS not PD")) } 
    B0_hat_boot_worker <- tryCatch(t(chol(Sigma_U_hat_ols_boot_worker)), error = function(e) NULL)
    if (is.null(B0_hat_boot_worker)) { return(list(error = "Cholesky_boot_OLS failed")) } 
    
    Psi_matrices_ols_boot_worker <- compute_VMA_coeffs_from_VAR(estimated_A_coeffs_ols_boot_worker, p_lags, K_vars_total, n_ahead_irf)
    current_boot_irfs_raw_worker <- array(NA, dim = c(K_vars_total, K_vars_total, n_ahead_irf + 1))
    for (h_boot_irf in 0:n_ahead_irf) { current_boot_irfs_raw_worker[,,h_boot_irf+1] <- Psi_matrices_ols_boot_worker[[h_boot_irf+1]] %*% B0_hat_boot_worker }
    
    ffr_on_ffr_h0_boot_scale_worker <- current_boot_irfs_raw_worker[ffr_shock_idx, ffr_shock_idx, 1]
    if (abs(ffr_on_ffr_h0_boot_scale_worker) < 1e-9) ffr_on_ffr_h0_boot_scale_worker <- 1 
    current_boot_irfs_normalized_ffr_shock_worker <- current_boot_irfs_raw_worker[, ffr_shock_idx, ] / ffr_on_ffr_h0_boot_scale_worker
    
    list_for_this_bootstrap_rep <- list()
    for(resp_target_idx_boot in 1:length(target_response_indices)){ 
      actual_resp_idx_boot <- target_response_indices[resp_target_idx_boot] 
      resp_target_name_boot <- names(target_response_indices)[resp_target_idx_boot] 
      if(!is.na(actual_resp_idx_boot)){
        irf_vec_boot_normalized_worker <- current_boot_irfs_normalized_ffr_shock_worker[actual_resp_idx_boot, ]
        list_for_this_bootstrap_rep[[paste0(resp_target_name_boot, "_noncum")]] <- irf_vec_boot_normalized_worker
        if(resp_target_name_boot == "CPI" || resp_target_name_boot == "IP"){
          list_for_this_bootstrap_rep[[paste0(resp_target_name_boot, "_cum")]] <- cumsum(irf_vec_boot_normalized_worker) 
        } else { list_for_this_bootstrap_rep[[paste0(resp_target_name_boot, "_cum")]] <- irf_vec_boot_normalized_worker }
      }
    }
    return(list_for_this_bootstrap_rep) 
  } 
  
  stopCluster(cl)
  cat("\nParallel bootstrap finished.\n")
  
  bootstrap_IRFs_targeted_cumulative_storage <- array(NA, dim = c(length(target_response_indices), n_ahead_irf + 1, n_bootstrap_irf))
  bootstrap_IRFs_targeted_noncumulative_storage <- array(NA, dim = c(length(target_response_indices), n_ahead_irf + 1, n_bootstrap_irf))
  dimnames(bootstrap_IRFs_targeted_cumulative_storage) <- list(Response = names(target_response_vars_map), Horizon = 0:n_ahead_irf, Replication = 1:n_bootstrap_irf)
  dimnames(bootstrap_IRFs_targeted_noncumulative_storage) <- list(Response = names(target_response_vars_map), Horizon = 0:n_ahead_irf, Replication = 1:n_bootstrap_irf)
  
  valid_reps_count <- 0
  for (b_idx in 1:length(bootstrap_results_list)) {
    rep_result <- bootstrap_results_list[[b_idx]]
    if (!is.null(rep_result) && is.null(rep_result$error) && length(rep_result) > 0) { 
      valid_reps_count <- valid_reps_count + 1
      for (resp_target_idx_proc in 1:length(target_response_indices)) {
        resp_name_proc <- names(target_response_indices)[resp_target_idx_proc]
        noncum_name <- paste0(resp_name_proc, "_noncum")
        cum_name <- paste0(resp_name_proc, "_cum")
        if (!is.null(rep_result[[noncum_name]])) {
          bootstrap_IRFs_targeted_noncumulative_storage[resp_target_idx_proc, , valid_reps_count] <- rep_result[[noncum_name]]
        }
        if (!is.null(rep_result[[cum_name]])) {
          bootstrap_IRFs_targeted_cumulative_storage[resp_target_idx_proc, , valid_reps_count] <- rep_result[[cum_name]]
        }
      }
    }
  }
  
  if (valid_reps_count < n_bootstrap_irf) {
    warning(paste(n_bootstrap_irf - valid_reps_count, "bootstrap replication(s) failed or were empty and discarded."))
    if (valid_reps_count == 0) stop("All bootstrap replications failed. Cannot compute CIs.")
    bootstrap_IRFs_targeted_cumulative_storage <- bootstrap_IRFs_targeted_cumulative_storage[,,1:valid_reps_count, drop=FALSE]
    bootstrap_IRFs_targeted_noncumulative_storage <- bootstrap_IRFs_targeted_noncumulative_storage[,,1:valid_reps_count, drop=FALSE]
  }
  
  results_list <- list()
  if (n_bootstrap_paths_to_return > 0 && valid_reps_count > 0) {
    results_list$bootstrap_paths_cumulative <- list()
    results_list$bootstrap_paths_noncumulative <- list()
    num_paths_to_sample <- min(n_bootstrap_paths_to_return, valid_reps_count)
    sampled_boot_indices <- sample(1:valid_reps_count, num_paths_to_sample)
    for(resp_target_name_path in names(target_response_vars_map)){
      storage_idx_path <- which(names(target_response_vars_map) == resp_target_name_path)
      if(!is.na(target_response_indices[resp_target_name_path])) { 
        results_list$bootstrap_paths_cumulative[[resp_target_name_path]] <- 
          bootstrap_IRFs_targeted_cumulative_storage[storage_idx_path, , sampled_boot_indices, drop=FALSE]
        results_list$bootstrap_paths_noncumulative[[resp_target_name_path]] <- 
          bootstrap_IRFs_targeted_noncumulative_storage[storage_idx_path, , sampled_boot_indices, drop=FALSE]
      }
    }
  }
  
  lower_quantile <- (1 - ci_level) / 2; upper_quantile <- 1 - lower_quantile
  for(resp_target_name_out in names(target_response_vars_map)){ 
    actual_resp_idx_out <- target_response_indices[resp_target_name_out] 
    if(!is.na(actual_resp_idx_out)){
      point_irf_vec_cum <- point_estimates_targeted_cumulative[[resp_target_name_out]]
      lower_ci_vec_cum <- numeric(n_ahead_irf + 1); upper_ci_vec_cum <- numeric(n_ahead_irf + 1)
      point_irf_vec_noncum <- point_estimates_targeted_noncumulative[[resp_target_name_out]]
      lower_ci_vec_noncum <- numeric(n_ahead_irf + 1); upper_ci_vec_noncum <- numeric(n_ahead_irf + 1)
      for (h_val in 1:(n_ahead_irf + 1)) {
        storage_idx_for_resp_out <- which(names(target_response_vars_map) == resp_target_name_out) 
        valid_boot_irfs_cum <- bootstrap_IRFs_targeted_cumulative_storage[storage_idx_for_resp_out, h_val, ] 
        valid_boot_irfs_cum <- valid_boot_irfs_cum[!is.na(valid_boot_irfs_cum)] 
        if(length(valid_boot_irfs_cum) > max(1, valid_reps_count * 0.05)){ # Adjusted threshold for quantile
          lower_ci_vec_cum[h_val] <- quantile(valid_boot_irfs_cum, lower_quantile, type = 4, na.rm = TRUE)
          upper_ci_vec_cum[h_val] <- quantile(valid_boot_irfs_cum, upper_quantile, type = 4, na.rm = TRUE)
        } else { lower_ci_vec_cum[h_val] <- NA; upper_ci_vec_cum[h_val] <- NA }
        valid_boot_irfs_noncum <- bootstrap_IRFs_targeted_noncumulative_storage[storage_idx_for_resp_out, h_val, ] 
        valid_boot_irfs_noncum <- valid_boot_irfs_noncum[!is.na(valid_boot_irfs_noncum)]
        if(length(valid_boot_irfs_noncum) > max(1, valid_reps_count * 0.05)){ # Adjusted threshold
          lower_ci_vec_noncum[h_val] <- quantile(valid_boot_irfs_noncum, lower_quantile, type = 4, na.rm = TRUE)
          upper_ci_vec_noncum[h_val] <- quantile(valid_boot_irfs_noncum, upper_quantile, type = 4, na.rm = TRUE)
        } else { lower_ci_vec_noncum[h_val] <- NA; upper_ci_vec_noncum[h_val] <- NA }
      }
      if (resp_target_name_out == "FFR") {
        results_list[[paste0("IRF_", resp_target_name_out)]] <- data.frame(
          Horizon = 0:n_ahead_irf, PointEstimate = point_irf_vec_noncum,
          LowerCI = lower_ci_vec_noncum, UpperCI = upper_ci_vec_noncum)
      } else {
        results_list[[paste0("IRF_", resp_target_name_out, "_cumulative")]] <- data.frame(
          Horizon = 0:n_ahead_irf, PointEstimate = point_irf_vec_cum,
          LowerCI = lower_ci_vec_cum, UpperCI = upper_ci_vec_cum)
        results_list[[paste0("IRF_", resp_target_name_out, "_noncumulative")]] <- data.frame(
          Horizon = 0:n_ahead_irf, PointEstimate = point_irf_vec_noncum,
          LowerCI = lower_ci_vec_noncum, UpperCI = upper_ci_vec_noncum)}} 
    else {
      results_list[[paste0("IRF_", resp_target_name_out, "_cumulative")]] <- paste("Var '", target_response_vars_map[[resp_target_name_out]], "' not found.", sep="")
      results_list[[paste0("IRF_", resp_target_name_out, "_noncumulative")]] <- paste("Var '", target_response_vars_map[[resp_target_name_out]], "' not found.", sep="")
    }
  }
  
  cat("\nTargeted, (Pre-)Transformed, Normalized, Cumulative and Non-Cumulative Bootstrap CIs computed.\n")
  return(results_list)
}