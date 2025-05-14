library(glmnet)
library(vars)

estimate_lasso_svar_structured <- function(structured_data, p_lags, n_ahead_irf = 48, n_bootstrap_irf = 199, ci_level = 0.95) {
  
  # --- 0. Data Preparation: Combine and Order Data ---
  # Combine data in the specified "Slow -> FFR -> Fast" order
  Y_data_ordered <- as.matrix(cbind(structured_data$slow_data, 
                          structured_data$FFR, 
                          structured_data$fast_data))
  
  K_vars_total <- ncol(Y_data_ordered)
  var_names <- colnames(Y_data_ordered)
  if(is.null(var_names)) var_names <- paste0("V", 1:K_vars_total) # Default names if none
  
  # cat("Total variables in ordered system:", K_vars_total, "\n")
  # cat("Variable order for VAR and Cholesky:", var_names, "\n")
  
  # --- 1. Create Lagged Regressors for VAR ---
  # This helper function creates the Y_target (LHS) and X_lagged (RHS) matrices for VAR estimation
  create_VAR_lags_glmnet <- function(data, p) {
    K <- ncol(data)
    T_orig <- nrow(data)
    
    if (T_orig <= p) {
      stop("Number of observations must be greater than the number of lags.")
    }
    
    # Matrix to hold lagged predictors (RHS)
    # Each block of K columns corresponds to a lag
    X_lagged <- matrix(NA, nrow = T_orig - p, ncol = K * p)
    
    # Target variables (LHS, unlagged, starting from p+1)
    Y_target <- data[(p + 1):T_orig, , drop = FALSE]
    
    col_names_X <- c()
    for (lag_val in 1:p) {
      lag_data_block <- data[(p + 1 - lag_val):(T_orig - lag_val), , drop = FALSE]
      X_lagged[, ((lag_val - 1) * K + 1):(lag_val * K)] <- lag_data_block
      col_names_X <- c(col_names_X, paste0(colnames(data), ".L", lag_val))
    }
    colnames(X_lagged) <- col_names_X
    
    return(list(Y_target = Y_target, X_lagged = X_lagged))
  }
  
  prepared_data <- create_VAR_lags_glmnet(Y_data_ordered, p_lags)
  Y_target_matrix <- prepared_data$Y_target
  X_predictors_matrix <- prepared_data$X_lagged
  
  T_eff <- nrow(Y_target_matrix) # Effective number of observations
  
  # --- 2. Equation-by-Equation LASSO Estimation (Phase 1 of LASSO-SVAR) ---
  num_predictors <- ncol(X_predictors_matrix)
  
  # B_hat_lasso stores the (K x Kp) matrix of slope coefficients
  B_hat_lasso <- matrix(0, nrow = K_vars_total, ncol = num_predictors) 
  intercepts_lasso <- numeric(K_vars_total)
  
  cat("\nStarting LASSO estimation for each VAR equation...\n")
  for (k in 1:K_vars_total) {
    y_k <- Y_target_matrix[, k]
    
    # Find optimal lambda using cross-validation (nfolds can be adjusted)
    # For time series, specialized CV like rolling window might be better, but 10-fold is common.
    cv_fit <- cv.glmnet(X_predictors_matrix, y_k, alpha = 1, nfolds = min(10, T_eff %/% 3), standardize = TRUE)
    optimal_lambda <- cv_fit$lambda.min # Or cv_fit$lambda.1se for a more parsimonious model
    
    # Fit final LASSO model with optimal lambda
    lasso_model <- glmnet(X_predictors_matrix, y_k, alpha = 1, lambda = optimal_lambda, standardize = TRUE)
    
    coeffs_k <- coef(lasso_model) # This is a sparse matrix object
    intercepts_lasso[k] <- coeffs_k[1, 1]
    B_hat_lasso[k, ] <- as.vector(coeffs_k[-1, 1]) 
    
    cat("  Estimated LASSO for equation:", var_names[k], "(lambda_min:", round(optimal_lambda,5) ,")\n")
  }
  rownames(B_hat_lasso) <- var_names
  colnames(B_hat_lasso) <- colnames(X_predictors_matrix)
  
  # Reshape B_hat_lasso into list of A1, A2, ..., Ap matrices (K x K)
  estimated_A_coeffs_lasso <- list()
  for (i in 1:p_lags) {
    A_i_matrix <- B_hat_lasso[, ((i - 1) * K_vars_total + 1):(i * K_vars_total)]
    colnames(A_i_matrix) <- var_names # Columns are the variables whose lags are taken
    rownames(A_i_matrix) <- var_names # Rows are the dependent variables in VAR equations
    estimated_A_coeffs_lasso[[i]] <- A_i_matrix
  }
  names(estimated_A_coeffs_lasso) <- paste0("A", 1:p_lags)
  
  # --- 3. Calculate Reduced-Form Residuals (Phase 2 of LASSO-SVAR) ---
  U_hat_lasso <- matrix(NA, nrow = T_eff, ncol = K_vars_total)
  colnames(U_hat_lasso) <- var_names
  
  Y_fitted_lasso <- matrix(NA, nrow = T_eff, ncol = K_vars_total)
  for (t in 1:T_eff) {
    # y_lagged_vec_t is the t-th row of X_predictors_matrix, representing [Y_{t-1}', Y_{t-2}', ..., Y_{t-p}']'
    y_lagged_vec_t <- X_predictors_matrix[t, ] 
    Y_fitted_lasso[t, ] <- intercepts_lasso + B_hat_lasso %*% y_lagged_vec_t
  }
  U_hat_lasso <- Y_target_matrix - Y_fitted_lasso
  
  # --- 4. Estimate Residual Covariance Matrix (Sigma_U) ---
  # Degrees of freedom adjustment is complex for LASSO.
  # A common simplification is T_eff or T_eff - 1, or sum of non-zero coeffs.
  # Using T_eff - 1 for simplicity here.
  Sigma_U_hat_lasso <- (t(U_hat_lasso) %*% U_hat_lasso) / (T_eff - 1) 
  
  # --- 5. Structural Identification using Cholesky Decomposition ---
  # The ordering "Slow -> FFR -> Fast" is already embedded in Y_data_ordered,
  # and thus in U_hat_lasso and Sigma_U_hat_lasso.
  # B0_hat is the contemporaneous impact matrix (U_t = B0_hat * epsilon_t)
  # For Cholesky, if U = B0 * eps, and Sigma_U = B0 * B0', then B0 = t(chol(Sigma_U))
  # Or B0 = chol(Sigma_U) depending on how you define B0 and structural shocks.
  # Standardly, for impulse = var_j, shock_j is the j-th column of B0.
  # Let's assume U_t = P * eps_t where P is lower triangular from Cholesky(Sigma_U)
  # Then P = t(chol(Sigma_U_hat_lasso)) if Sigma_U = P P'
  
  # Check for positive definiteness before Cholesky
  eigen_Sigma <- eigen(Sigma_U_hat_lasso, symmetric = TRUE, only.values = TRUE)$values
  if(any(eigen_Sigma <= 1e-8)){ # Check for very small or negative eigenvalues
    warning("Estimated Sigma_U is not positive definite or near singular. Cholesky may fail or be unreliable. Consider adding a small diagonal component (ridge) or re-evaluating LASSO.")
    # Example of adding a small ridge: Sigma_U_hat_lasso <- Sigma_U_hat_lasso + diag(1e-7, K_vars_total)
  }
  
  B0_hat <- tryCatch({
    t(chol(Sigma_U_hat_lasso)) # B0 such that Sigma_U = B0 * B0'
  }, error = function(e) {
    stop(paste("Cholesky decomposition failed. Sigma_U may not be positive definite. Error:", e$message))
  })
  colnames(B0_hat) <- var_names # Shocks are named after the variable they primarily affect in Cholesky
  rownames(B0_hat) <- var_names # Response of variables
  
  cat("\nStructural identification via Cholesky complete.\n")
  cat("FFR is variable number:", which(var_names == colnames(ffr_df)[1]), "in the ordering.\n")
  
  # --- 6. Compute Structural Impulse Responses (Phase 3 of LASSO-SVAR) ---
  # Helper function to get VMA (Psi) matrices from VAR (A) matrices
  compute_VMA_coeffs_from_VAR <- function(A_coeffs_list, p_lags_var, K_vars_var, n_ahead_vma) {
    Psi_matrices <- list()
    Psi_matrices[[1]] <- diag(K_vars_var) # Psi_0 = I
    
    # Need full A matrices including those that might be zero if p_lags_var > length(A_coeffs_list)
    # This means A_coeffs_list should contain p_lags_var elements, some can be zero matrices
    A_full <- list()
    for(i in 1:p_lags_var) {
      if(i <= length(A_coeffs_list)) {
        A_full[[i]] <- A_coeffs_list[[i]]
      } else {
        A_full[[i]] <- matrix(0, nrow=K_vars_var, ncol=K_vars_var)
      }
    }
    
    for (h in 1:n_ahead_vma) { # Psi_1 to Psi_n_ahead
      Psi_h <- matrix(0, nrow = K_vars_var, ncol = K_vars_var)
      for (j in 1:min(h, p_lags_var)) {
        Psi_h <- Psi_h + (Psi_matrices[[h - j + 1]] %*% A_full[[j]])
      }
      Psi_matrices[[h + 1]] <- Psi_h
    }
    return(Psi_matrices)
  }
  
  Psi_matrices_lasso <- compute_VMA_coeffs_from_VAR(estimated_A_coeffs_lasso, p_lags, K_vars_total, n_ahead_irf)
  
  # structural_IRFs_point: Array [response_var, shock_var, horizon]
  structural_IRFs_point <- array(NA, dim = c(K_vars_total, K_vars_total, n_ahead_irf + 1),
                                 dimnames = list(Response = var_names, Shock = var_names, Horizon = 0:n_ahead_irf))
  
  for (h in 0:n_ahead_irf) {
    structural_IRFs_point[,,h+1] <- Psi_matrices_lasso[[h+1]] %*% B0_hat
  }
  
  cat("Point estimates of structural IRFs computed.\n")
  
  # --- 7. Bootstrap Confidence Intervals for IRFs ---
  # This is a complex part. The bootstrap should ideally account for LASSO selection uncertainty.
  # A common (though potentially imperfect) approach is residual bootstrap, re-running LASSO.
  # For simplicity, a basic residual bootstrap structure is shown.
  
  # Store IRFs from each bootstrap replication: [response_var, shock_var, horizon, boot_rep]
  bootstrap_IRFs_storage <- array(NA, dim = c(K_vars_total, K_vars_total, n_ahead_irf + 1, n_bootstrap_irf))
  
  cat("\nStarting bootstrap for IRF confidence intervals (", n_bootstrap_irf, "replications)...\n")
  pb <- utils::txtProgressBar(min = 0, max = n_bootstrap_irf, style = 3)
  
  for (b in 1:n_bootstrap_irf) {
    # 1. Resample residuals (simple version: sample rows with replacement)
    # Ensure residuals are centered if not already
    U_hat_centered <- scale(U_hat_lasso, center = TRUE, scale = FALSE)
    boot_indices <- sample(1:T_eff, T_eff, replace = TRUE)
    U_boot <- U_hat_centered[boot_indices, , drop = FALSE]
    
    # 2. Generate bootstrap Y_data_boot using estimated LASSO parameters and U_boot
    Y_data_boot <- matrix(NA, nrow = T_eff + p_lags, ncol = K_vars_total)
    Y_data_boot[1:p_lags, ] <- Y_data_ordered[1:p_lags, ] # Initial conditions
    
    current_lags_matrix_boot <- X_predictors_matrix[1, , drop=FALSE] # Placeholder, will be updated
    
    for (t in (p_lags + 1):(T_eff + p_lags)) {
      # Prepare lagged Y for prediction
      # This needs to construct the X_predictors_matrix row for Y_data_boot at time t-1, ..., t-p
      # Simplified: assume we can reconstruct the predictor vector
      # This is a tricky part for dynamic simulation in bootstrap
      
      # For dynamic simulation, we need to build the X vector from Y_data_boot itself
      x_pred_t_boot <- numeric(K_vars_total * p_lags)
      for(lag_b in 1:p_lags){
        x_pred_t_boot[((lag_b-1)*K_vars_total + 1) : (lag_b*K_vars_total)] <- Y_data_boot[t-lag_b, ]
      }
      
      Y_fitted_t_boot <- intercepts_lasso + B_hat_lasso %*% x_pred_t_boot
      Y_data_boot[t, ] <- Y_fitted_t_boot + U_boot[t - p_lags, ]
    }
    Y_data_boot_for_var <- Y_data_boot # Full series for re-estimation
    
    # 3. Re-estimate LASSO-VAR on Y_data_boot_for_var
    # This is the most computationally intensive part and crucial for selection uncertainty
    prepared_data_boot <- create_VAR_lags_glmnet(Y_data_boot_for_var, p_lags)
    Y_target_matrix_boot <- prepared_data_boot$Y_target
    X_predictors_matrix_boot <- prepared_data_X_lagged # Using original for simplicity or re-calculate
    
    B_hat_lasso_boot <- matrix(0, nrow = K_vars_total, ncol = num_predictors)
    intercepts_lasso_boot <- numeric(K_vars_total)
    
    for (k_boot in 1:K_vars_total) {
      y_k_boot <- Y_target_matrix_boot[, k_boot]
      # For speed, could use the original optimal_lambda or re-estimate (better)
      cv_fit_boot <- cv.glmnet(X_predictors_matrix_boot, y_k_boot, alpha = 1, nfolds = min(5, T_eff %/% 5), standardize = TRUE) # Fewer folds for speed
      optimal_lambda_boot <- cv_fit_boot$lambda.min
      lasso_model_boot <- glmnet(X_predictors_matrix_boot, y_k_boot, alpha = 1, lambda = optimal_lambda_boot, standardize = TRUE)
      coeffs_k_boot <- coef(lasso_model_boot)
      intercepts_lasso_boot[k_boot] <- coeffs_k_boot[1, 1]
      B_hat_lasso_boot[k_boot, ] <- as.vector(coeffs_k_boot[-1, 1])
    }
    
    estimated_A_coeffs_lasso_boot <- list()
    for (i_boot in 1:p_lags) {
      estimated_A_coeffs_lasso_boot[[i_boot]] <- B_hat_lasso_boot[, ((i_boot - 1) * K_vars_total + 1):(i_boot * K_vars_total)]
    }
    
    # 4. Re-calculate residuals and Sigma_U_boot
    U_hat_lasso_boot <- matrix(NA, nrow = T_eff, ncol = K_vars_total)
    Y_fitted_lasso_boot <- matrix(NA, nrow = T_eff, ncol = K_vars_total)
    for (t_boot in 1:T_eff) {
      y_lagged_vec_t_boot <- X_predictors_matrix_boot[t_boot, ] 
      Y_fitted_lasso_boot[t_boot, ] <- intercepts_lasso_boot + B_hat_lasso_boot %*% y_lagged_vec_t_boot
    }
    U_hat_lasso_boot <- Y_target_matrix_boot - Y_fitted_lasso_boot
    Sigma_U_hat_lasso_boot <- (t(U_hat_lasso_boot) %*% U_hat_lasso_boot) / (T_eff - 1)
    
    # 5. Re-do Cholesky and compute IRFs for this bootstrap sample
    eigen_Sigma_boot <- eigen(Sigma_U_hat_lasso_boot, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigen_Sigma_boot <= 1e-8)){
      # Skip this bootstrap rep or add ridge
      # For simplicity, skipping if problematic
      utils::setTxtProgressBar(pb, b)
      next 
    }
    B0_hat_boot <- tryCatch(t(chol(Sigma_U_hat_lasso_boot)), error = function(e) NULL)
    if (is.null(B0_hat_boot)) { # Cholesky failed
      utils::setTxtProgressBar(pb, b)
      next
    }
    
    Psi_matrices_lasso_boot <- compute_VMA_coeffs_from_VAR(estimated_A_coeffs_lasso_boot, p_lags, K_vars_total, n_ahead_irf)
    for (h_boot in 0:n_ahead_irf) {
      bootstrap_IRFs_storage[,,h_boot+1,b] <- Psi_matrices_lasso_boot[[h_boot+1]] %*% B0_hat_boot
    }
    utils::setTxtProgressBar(pb, b)
  } # End bootstrap loop
  close(pb)
  
  # Calculate confidence intervals
  lower_quantile <- (1 - ci_level) / 2
  upper_quantile <- 1 - lower_quantile
  
  IRF_lower_ci <- array(NA, dim = dim(structural_IRFs_point))
  IRF_upper_ci <- array(NA, dim = dim(structural_IRFs_point))
  dimnames(IRF_lower_ci) <- dimnames(IRF_upper_ci) <- dimnames(structural_IRFs_point)
  
  for (resp_v in 1:K_vars_total) {
    for (shock_v in 1:K_vars_total) {
      for (h_val in 1:(n_ahead_irf + 1)) {
        valid_boot_irfs <- bootstrap_IRFs_storage[resp_v, shock_v, h_val, !is.na(bootstrap_IRFs_storage[resp_v, shock_v, h_val, ])]
        if(length(valid_boot_irfs) > n_bootstrap_irf * 0.5){ # Require at least 50% valid bootstrap reps
          IRF_lower_ci[resp_v, shock_v, h_val] <- quantile(valid_boot_irfs, lower_quantile, type = 4, na.rm = TRUE)
          IRF_upper_ci[resp_v, shock_v, h_val] <- quantile(valid_boot_irfs, upper_quantile, type = 4, na.rm = TRUE)
        } else {
          IRF_lower_ci[resp_v, shock_v, h_val] <- NA
          IRF_upper_ci[resp_v, shock_v, h_val] <- NA
        }
      }
    }
  }
  cat("\nBootstrap CIs computed.\n")
  
  # --- Return Results ---
  return(list(
    ordered_variable_names = var_names,
    ffr_column_index = which(var_names == colnames(ffr_df)[1]),
    estimated_A_coeffs_lasso = estimated_A_coeffs_lasso,
    intercepts_lasso = intercepts_lasso,
    B_hat_lasso_Kp_form = B_hat_lasso, # K x Kp matrix
    Sigma_U_hat_lasso = Sigma_U_hat_lasso,
    B0_hat_cholesky = B0_hat,
    structural_IRFs_point = structural_IRFs_point,
    structural_IRFs_lower_ci = IRF_lower_ci,
    structural_IRFs_upper_ci = IRF_upper_ci,
    X_predictors_matrix = X_predictors_matrix, # For reference
    Y_target_matrix = Y_target_matrix # For reference
  ))
}