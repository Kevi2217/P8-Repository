library(HDLPrepro) 
library(desla)
library(MASS)
library(vars)
library(ggplot2)
library(reshape2)
library(dplyr)
library(vars)
library(foreach)
library(doParallel)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(parallel)

# Sources functions used in analysis
source("/Users/frederikdyrup/Desktop/8.semester/projekt_kode/egen_kode/DGP_bias_analysis/functions_for3_1.R")

# For measuring how long script took
# start_time <- Sys.time()

################################ settings ######################################
seed_value <- 1
Ns <- c(10,20,40) # numbers of variables
Ts <- Ts_list_example <- list(# sample size
    c(100, 200, 500),   # T values for N=10
    c(200, 300, 500),   # T values for N=20
    c(300, 400, 500))   # T values for N=40
PIs <- c(0.8) # values of the plug-in constant used for the data-dependent selection of the lasso penalization. Generally, higher value gives stronger penalization. For details, see Algorithm 1 in the supplementary appendix C.5 of https://doi.org/10.1016/j.jeconom.2022.08.008. 
run_only_PI0.8 <- TRUE #  When true, the code only runs and saves the simulation for the plug-in value of 0.8, which is what we present in the paper. To run the simulation for all plug-in values in the PIs vector, set this to true. Note that to actually plot the results for other plug-in values, changes need to be made to the "plotting_simulation3_1.R" script. Running the full simulation may cause errors on non-Windows operating systems. 
M <- 100 # number of replications in the simulation
hmax <- 10 #  maximum horizon - the impulse response function is evaluated from horizon 0 to hmax
VAR_lags <- 4 # number of VAR lags in the DGP
LP_lags <- 4 # number of lags included in the local projection estimation equations
selection <- 4 # type of selection of the penalization parameter. 4 denotes the data-dependent selection which uses the plug-in constant above. 1-3 denote AIC, BIC and EBIC respectively
progress_bar <- TRUE # setting this to TRUE will show a progress bar each time simulate_LP is ran, so showing progress over simulation replications
OLS <- FALSE # setting whether OLS should be used in estimating the local projections
threads <- parallel::detectCores()-2 # the number of cores used in parallel computation 
alphas <- 0.05 # desired level of the test; equivalently, (1-alpha)% confidence intervals 
z_quantiles <- qnorm(alphas/2.0,0.0,1.0,FALSE,FALSE); # quantiles of the Normal distribution associated with alpha
chi2_quantiles <- qchisq(alphas,1,FALSE,FALSE); # quantiles of the Chi-squared distribution associated with alpha
################################################################################
############################ Initial analysis ##################################
# set.seed(seed_value)
# succes_data <- analyze_svar_failure_boundary(Ns = seq(10, 40, 1),
#                                              Ts = seq(50, 200, 5))
# 
# plot_svar_failure_boundary(succes_data)
################################################################################
################################# SVAR IRF #####################################
final_results <- list()
set.seed(seed_value)
# Parallel Computing SVAR IRF (DGP1)
source("/Users/frederikdyrup/Desktop/8.semester/projekt_kode/egen_kode/DGP_bias_analysis/parallel_3_1.R")
saveRDS(final_results, file="final_results.RDS")

final_results_switching <- list()
set.seed(seed_value)
# Parallel Computing SVAR IRF (DGP2)
source("/Users/frederikdyrup/Desktop/8.semester/projekt_kode/egen_kode/DGP_bias_analysis/parallel_3_1_switching.R")
saveRDS(final_results, file="final_results.RDS")
  
  for (case_name in names(final_results)) {
    final_results[[case_name]]$SVAR$estimate[1]   <- 1.0
    final_results[[case_name]]$SVAR$lower[1]      <- 1.0
    final_results[[case_name]]$SVAR$upper[1]      <- 1.0
    final_results[[case_name]]$SVAR$coverage[1]   <- 1.0
    final_results[[case_name]]$SVAR$rmse[1]       <- 0
    final_results[[case_name]]$SVAR$avg_width[1]  <- 0
    final_results[[case_name]]$SVAR$bias[1]       <- 0
    # switching
    final_results_switching[[case_name]]$SVAR$estimate[1]   <- 1.0
    final_results_switching[[case_name]]$SVAR$lower[1]      <- 1.0
    final_results_switching[[case_name]]$SVAR$upper[1]      <- 1.0
    final_results_switching[[case_name]]$SVAR$coverage[1]   <- 1.0
    final_results_switching[[case_name]]$SVAR$rmse[1]       <- 0
    final_results_switching[[case_name]]$SVAR$avg_width[1]  <- 0
    final_results_switching[[case_name]]$SVAR$bias[1]       <- 0
  }
################################################################################
############################## Partial LP IRF ##################################
# DGP1 Proposed desla: partial DL ---------------------------------------------------
init_partial=TRUE # should the parameters of interest remain unpenalized in the first step?

set.seed(seed_value)

for(n in 1:length(Ns)){
  current_Ts_for_N <- Ts[[n]]
  for(t in 1:length(current_Ts_for_N)){
    for(pi in 1:length(PIs)){
      N=Ns[n]
      T_=current_Ts_for_N[t]
      PIconstant=PIs[pi]
      VAR_coefficients<-make_VAR_coefs(N)
      irf_1to1<-irf_from_VAR(VAR_coefficients)
      Sigma_epsilon<-diag(N)
      seeds_gen<-sample(1e8, M)
      seeds_DL<-array(data=sample(1e8, 2*(hmax+1)*M), dim=c(2, hmax+1, M))
      if(PIconstant==0.8 || run_only_PI0.8==FALSE){
sim<-simulate_LP(M, T_, LP_lags, hmax, VAR_coefficients, Sigma_epsilon, irf_1to1,
                 init_partial, z_quantiles, chi2_quantiles, selection, PIconstant,
                 progress_bar, OLS, threads, seeds_gen, seeds_DL)
        saveRDS(sim, file=paste0("partial_N",N,"_T",T_,"_PI",PIconstant,".RDS"))
      }
    }
  }
}
print("Partial DL done")

# DGP2 Proposed desla: partial DL switching ----------------------------------------
init_partial=TRUE # should the parameters of interest remain unpenalized in the first step?

set.seed(seed_value)
################################################################################
for(n in 1:length(Ns)){
  current_Ts_for_N <- Ts[[n]]
  for(t in 1:length(current_Ts_for_N)){
    for(pi in 1:length(PIs)){
      N=Ns[n]
      T_=current_Ts_for_N[t]
      PIconstant=PIs[pi]
      VAR_coefficients<-make_VAR_coefs2(N) # uses VAR coefficients with switched signs
      irf_1to1<-irf_from_VAR(VAR_coefficients)
      Sigma_epsilon<-diag(N)
      seeds_gen<-sample(1e8, M)
      seeds_DL<-array(data=sample(1e8, 2*(hmax+1)*M), dim=c(2, hmax+1, M))
      if(PIconstant==0.8 || run_only_PI0.8==FALSE){
        sim<-simulate_LP(M, T_, LP_lags, hmax, VAR_coefficients, Sigma_epsilon, irf_1to1,
                         init_partial, z_quantiles, chi2_quantiles, selection, PIconstant,
                         progress_bar, OLS, threads, seeds_gen, seeds_DL)
        saveRDS(sim, file=paste0("partial_N",N,"_T",T_,"_PI",PIconstant,"_switching_signs.RDS"))
      }
    }
  }
}
print("Partial DL Switching done")
##############################################################################
#################################  Plotting  #################################

# Loading and processing simulations saved above
# coverage
svar_coverage<-sapply(final_results, function(item) item$SVAR$coverage)
coverage_partial_lp<-sims_to_coverage(path=getwd(), M, Ns,Ts,PIs,hmax,partial=TRUE, switching=FALSE)
svar_coverage_switching<-sapply(final_results_switching, function(item) item$SVAR$coverage)
coverage_partial_lp_switching<-sims_to_coverage(path=getwd(), M, Ns,Ts,PIs,hmax,partial=TRUE, switching=TRUE)

# width
svar_width<-sapply(final_results, function(item) item$SVAR$avg_width)
width_partial_lp<-sims_to_width(path=getwd(), M, Ns, Ts, PIs, hmax,partial=TRUE, switching=FALSE)
svar_width_switching<-sapply(final_results_switching, function(item) item$SVAR$avg_width)
width_partial_lp_switching<-sims_to_width(path=getwd(), M, Ns, Ts, PIs, hmax, partial=TRUE, switching=TRUE)

# bias
svar_bias<-sapply(final_results, function(item) item$SVAR$bias)
bias_partial_lp<-sims_to_bias(path=getwd(), M, Ns, Ts, PIs, hmax,partial=TRUE, switching=FALSE)
svar_bias_switching<-sapply(final_results_switching, function(item) item$SVAR$bias)
bias_partial_lp_switching<-sims_to_bias(path=getwd(), M, Ns, Ts, PIs, hmax, partial=TRUE, switching=TRUE)
# RMSE
svar_rmse<-sapply(final_results, function(item) item$SVAR$rmse)
rmse_partial_lp<-sims_to_rmse(path=getwd(), M, Ns, Ts, PIs, hmax,partial=TRUE, switching=FALSE)
svar_rmse_switching<-sapply(final_results_switching, function(item) item$SVAR$rmse)
rmse_partial_lp_switching<-sims_to_rmse(path=getwd(), M, Ns, Ts, PIs, hmax, partial=TRUE, switching=TRUE)


####################### Custom Legend used in all plots ########################
custom_legend <- get_legend(ggplot(data.frame(
  category = factor(c("SVAR", "LP"), # Descriptive names for your legend
                    levels = c("SVAR", "LP")),
  # Dummy x and y values, these won't actually be plotted in the final figure
  x_val = c(1, 1),
  y_val = c(1, 1)
), aes(x = x_val, y = y_val, color = category)) +
  geom_line(linewidth = 1) + # The geom_line here ensures the legend key is a line
  scale_color_manual(
    name = NULL, # Set to NULL if you don't want a title for the legend, or provide one e.g., "Method Type"
    values = c("SVAR" = "red", "LP" = "blue"),
    labels = c("SVAR", "LP") # These are the labels that will appear
  ) +
  guides(color = guide_legend(
    direction = "horizontal",  # Arrange legend items side-by-side
  )) +
  theme_void() + # Minimal theme, we only want the legend
  theme(
    legend.position = "bottom", # This is the position within the dummy plot, ggarrange controls final
    legend.direction = "horizontal", # Reinforce horizontal layout for the whole legend box
    legend.spacing.x = unit(0.5, "cm") 
  ))
################################################################################
# Plot Coverages ----------------------------------------------------------
# code for Figure S.1 in the supplementary appendix and Figure 1
pi <- 1 # pi=1 plots the simulation results for plug-in constant = 0.8, which is what we include in the paper. If more values were simulated, this needs to be manually adjusted
count <- 0
P <- list()
for(n in 1:length(Ns)){
  for(t in 1:length(Ts)){
    count<-count+1
    df<-data.frame(horizon=0:hmax,
                   svar_not_switching=svar_coverage[,count], svar_switching=svar_coverage_switching[,count],
                   not_switching_partial=coverage_partial_lp[n,t,pi,], switching_partial=coverage_partial_lp_switching[n,t,pi,])
    df_melt<-melt(df, id="horizon")
    P[[count]]<-ggplot(data=df,
                       aes(x=horizon, y=value, colour=variable)) +
      geom_hline(yintercept=0.95, color="black", linewidth=0.5)+
      geom_line(aes(x=horizon,y=svar_not_switching), colour="red", linewidth=0.75)+ #regular DL
      geom_line(aes(x=horizon,y=svar_switching), colour="red", linewidth=0.75, linetype="twodash")+ #regular DL switching
      geom_line(aes(x=horizon,y=not_switching_partial), colour="blue", linewidth=0.75)+ #partial DL
      geom_line(aes(x=horizon,y=switching_partial), colour="blue",linewidth=0.75,  linetype="twodash")+ #partial DL switching
      ylim(c(0,1.0001))+
      ylab(NULL)+
      xlab("Horizon")+
      scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      ggtitle(paste0("P=",Ns[n], ", T=",Ts[[n]][[t]]))+theme_bw()
  }
}
# This ensures that ggarrange doesn't try to use their original legends.
P <- lapply(P, function(plot_object) {
  plot_object + theme(legend.position = "none")
})
# Plot figure
annotate_figure(
  ggarrange(P[[1]],P[[2]],P[[3]],P[[4]],P[[5]],P[[6]],P[[7]],P[[8]],P[[9]],
            ncol=3, nrow=3,
            legend="bottom",common.legend = TRUE,
            legend.grob = custom_legend
  ),
  top = text_grob("Coverage",
                  color = "black",
                  face = "bold",
                  size = 12)
)
# Save Figure
ggsave(filename="coverage.png",device="png",width=18, height = 18, units="cm",dpi=1000)

# Plot Interval Widths ----------------------------------------------------
# code for Figure S.2 in the supplementary appendix
pi <- 1 # pi=1 plots the simulation results for plug-in constant = 0.8, which is what we include in the paper. If more values were simulated, this needs to be manually adjusted
count <- 0
P <- list()
for(n in 1:length(Ns)){
  for(t in 1:length(Ts)){
    count<-count+1
    df<-data.frame(horizon=0:hmax,
                   svar_not_switching=svar_width[, count], svar_switching=svar_width_switching[, count],
                   not_switching_partial=width_partial_lp[n,t,pi,], switching_partial=width_partial_lp_switching[n,t,pi,])
    df_melt<-melt(df, id="horizon")
    P[[count]]<-ggplot(data=df,
                       aes(x=horizon, y=value, colour=variable)) +
      geom_line(aes(x=horizon,y=svar_not_switching), colour="red", linewidth=0.75)+ #regular DL
      geom_line(aes(x=horizon,y=svar_switching), colour="red", linewidth=0.75,  linetype="twodash")+ #regular DL switching
      geom_line(aes(x=horizon,y=not_switching_partial), colour="blue", linewidth=0.75)+ #partial DL
      geom_line(aes(x=horizon,y=switching_partial), colour="blue", linewidth=0.75, linetype="twodash")+ #partial DL switching
      scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      ggtitle(paste0("P=",Ns[n], ", T=",Ts[[n]][[t]]))+
      ylab(NULL)+
      xlab("Horizon")+theme_bw()
  }
}
# This ensures that ggarrange doesn't try to use their original legends.
P <- lapply(P, function(plot_object) {
  plot_object + theme(legend.position = "none")
})
# Plot figure
annotate_figure(
  ggarrange(P[[1]],P[[2]],P[[3]],P[[4]],P[[5]],P[[6]],P[[7]],P[[8]],P[[9]],
            ncol=3, nrow=3,
            legend="bottom",common.legend = TRUE,
            legend.grob = custom_legend
  ),
  top = text_grob("Width",
                  color = "black",
                  face = "bold",
                  size = 12)
)
# Save Figure
ggsave(filename="width.png",device="png",width=18, height = 18, units="cm",dpi=1000)



# Plot Bias -----------------------------------------------------
pi <- 1 # pi=1 plots the simulation results for plug-in constant = 0.8, which is what we include in the paper. If more values were simulated, this needs to be manually adjusted
count <- 0
P <- list()
for(n in 1:length(Ns)){
  for(t in 1:length(Ts)){
    count<-count+1
    df<-data.frame(horizon=0:hmax,
                   svar_not_switching=svar_bias[, count], svar_switching=svar_bias_switching[, count],
                   not_switching_partial=bias_partial_lp[n,t,pi,], switching_partial=bias_partial_lp_switching[n,t,pi,])
    
    df_melt<-melt(df, id="horizon")
    P[[count]]<-ggplot(data=df,
                       aes(x=horizon, y=value, colour=variable)) +
      geom_line(aes(x=horizon,y=svar_not_switching), colour="red", linewidth=0.75)+ #regular DL
      geom_line(aes(x=horizon,y=svar_switching), colour="red", linewidth=0.75,  linetype="twodash")+ #regular DL switching
      geom_line(aes(x=horizon,y=not_switching_partial), colour="blue", linewidth=0.75)+ #partial DL
      geom_line(aes(x=horizon,y=switching_partial), colour="blue", linewidth=0.75, linetype="twodash")+ #partial DL switching
      scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      ggtitle(paste0("P=",Ns[n], ", T=",Ts[[n]][[t]]))+
      ylab(NULL)+
      xlab("Horizon")+
      theme_bw()
  }
}
# This ensures that ggarrange doesn't try to use their original legends.
P <- lapply(P, function(plot_object) {
  plot_object + theme(legend.position = "none")
})
# Plot figure
annotate_figure(
  ggarrange(P[[1]],P[[2]],P[[3]],P[[4]],P[[5]],P[[6]],P[[7]],P[[8]],P[[9]],
            ncol=3, nrow=3,
            legend="bottom",common.legend = TRUE,
            legend.grob = custom_legend
            ),
  top = text_grob("Bias",
                  color = "black",
                  face = "bold",
                  size = 12)
)
# Save Figure
ggsave(filename="bias.png",device="png",width=18, height = 18, units="cm",dpi=1000)


# Plot RMSE -----------------------------------------------------
# code for Figure S.2 in the supplementary appendix
pi <- 1 # pi=1 plots the simulation results for plug-in constant = 0.8, which is what we include in the paper. If more values were simulated, this needs to be manually adjusted
count <- 0
P <- list()
for(n in 1:length(Ns)){
  for(t in 1:length(Ts)){
    count<-count+1
    df<-data.frame(horizon=0:hmax,
                   svar_not_switching=svar_rmse[, count], svar_switching=svar_rmse_switching[, count],
                   not_switching_partial=rmse_partial_lp[n,t,pi,], switching_partial=rmse_partial_lp_switching[n,t,pi,])
    df_melt<-melt(df, id="horizon")
    P[[count]]<-ggplot(data=df,
                       aes(x=horizon, y=value, colour=variable)) +
      geom_line(aes(x=horizon,y=svar_not_switching), colour="red", linewidth=0.75)+ #regular DL
      geom_line(aes(x=horizon,y=svar_switching), colour="red", linewidth=0.75,  linetype="twodash")+ #regular DL switching
      geom_line(aes(x=horizon,y=not_switching_partial), colour="blue", linewidth=0.75)+ #partial DL
      geom_line(aes(x=horizon,y=switching_partial), colour="blue", linewidth=0.75, linetype="twodash")+ #partial DL switching
      scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      ggtitle(paste0("P=",Ns[n], ", T=",Ts[[n]][[t]]))+
      ylab(NULL)+
      xlab("Horizon")+theme_bw()
  }
}
# This ensures that ggarrange doesn't try to use their original legends.
P <- lapply(P, function(plot_object) {
  plot_object + theme(legend.position = "none")
})
# Plot figure
annotate_figure(
  ggarrange(P[[1]],P[[2]],P[[3]],P[[4]],P[[5]],P[[6]],P[[7]],P[[8]],P[[9]],
            ncol=3, nrow=3,
            legend="bottom",common.legend = TRUE,
            legend.grob = custom_legend
  ),
  top = text_grob("RMSE",
                  color = "black",
                  face = "bold",
                  size = 12)
)
# Save Figure
ggsave(filename="rmse.png",device="png",width=18, height = 18, units="cm",dpi=1000)

