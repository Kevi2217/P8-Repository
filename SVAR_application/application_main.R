library(glmnet)
library(vars)
library(HDLPrepro)
library(xtable)
library(ggplot2)
library(ggpubr)
library(foreach)
library(doParallel)
library(desla)
library(glmnet)
library(hdi)
library(MASS)


# Set working directory to current folder

# "Processing_FREDMD.R" is from Adamek's own repository
source(paste0(getwd(),"/processing_FREDMD.R"))
source(paste0(getwd(),"/functions_for_app.R"))

# Load CDbmedium if not running "processing_FREDMD.R"
# CDbmedium <- readRDS(file = "C:/Users/dyrup/Desktop/fred_applicaton/CDbmedium.RData")
################################ settings ######################################
hmax<-48 # maximum horizon
lags<-13 # number of lags
PI_constant<-0.4 # this is the plug-in constant used for the data-dependent selection of the lasso penalization.
threads<-parallel::detectCores()-2 # number of cores used in parallel computation 
seed_value <- 1 # seed for the random number generation for reproducibility.
################################################################################

set.seed(seed_value)
######################### Fixed boostrap post LASSO SVAR #######################
results_list <- estimate_fixed_post_lasso_svar_structured(
  structured_data = CDbmedium,
  p_lags = lags,
  target_CPI_name = "CPIAUCSL", 
  target_IP_name = "INDPRO",
  n_ahead_irf = hmax, 
  n_bootstrap_irf = 199,
  ci_level = 0.95)

# Save or load
saveRDS(results_list, file = paste0("results_list.RData"))
# load(paste0(getwd(),"/SVAR_results.RData")

################################################################################
#################################### HDLP ######################################
other_slow<-as.matrix(CDbmedium$slow_data[,-c(which(colnames(CDbmedium$slow_data)=="INDPRO"),
                                              which(colnames(CDbmedium$slow_data)=="CPIAUCSL"))])
fast<-as.matrix(CDbmedium$fast_data)
IP<-as.matrix(CDbmedium$data_all[,"INDPRO"])
CPI<-as.matrix(CDbmedium$data_all[,"CPIAUCSL"])
FFR<-as.matrix(CDbmedium$FFR)

# Set working directory to current folder
# Load or run HDLP() functions
HDLP_FFR_irf <- readRDS(paste0(getwd(),"/HDLP_FFR_irf.RData"))
HDLP_IP_irf <- readRDS(paste0(getwd(),"/HDLP_IP_irf.RData"))
HDLP_CPI_irf <- readRDS(paste0(getwd(),"/HDLP_CPI_irf.RData"))

# HDLP_FFR_irf<-HDLP(r=cbind(other_slow,IP,CPI), x=FFR, y=FFR, q=fast,
#                    y_predetermined = F, cumulate_y = F, hmax=hmax, lags=lags,
#                    threads = threads, PI_constant = PI_constant)
# 
# 
# HDLP_IP_irf<-HDLP(r=cbind(other_slow,CPI), x=FFR, y=IP, q=fast,
#                   y_predetermined = T, cumulate_y = T, hmax=hmax, lags=lags,
#                   threads=threads, PI_constant = PI_constant)
# 
# 
# HDLP_CPI_irf<-HDLP(r=cbind(other_slow,IP), x=FFR, y=CPI, q=fast,
#                    y_predetermined = T, cumulate_y = T, hmax=hmax, lags=lags,
#                    threads=threads, PI_constant = PI_constant)
################################################################################
################################## PLOTTING ####################################
p1<-ggplot(data=results_list$IRF_FFR, mapping=aes(x=Horizon, y=PointEstimate))+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_line()+
  geom_ribbon(mapping=aes(ymin=LowerCI, ymax=UpperCI), alpha=0.5)+
  xlab("Horizon")+ylab("FFR")

p2<-ggplot(data=results_list$IRF_IP_cumulative, mapping=aes(x=Horizon, y=PointEstimate))+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_line()+
  geom_ribbon(mapping=aes(ymin=LowerCI, ymax=UpperCI), alpha=0.5)+
  xlab("Horizon")+ylab("IP")

p3<-ggplot(data=results_list$IRF_CPI_cumulative, mapping=aes(x=Horizon, y=PointEstimate))+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_line()+
  geom_ribbon(mapping=aes(ymin=LowerCI, ymax=UpperCI), alpha=0.5)+
  xlab("Horizon")+ylab("CPI")

o1<-HDLP_FFR_irf
p4<-ggplot(mapping=aes(x=0:(dim(o1$intervals)[1]-1)))+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_line(mapping=aes(y=o1$intervals[,2,]), color="blue")+
  geom_ribbon(mapping=aes(ymin=o1$intervals[,1,], ymax=o1$intervals[,3,]), alpha=0.5, fill="blue")+
  xlab("Horizon")+ylab("FFR")

o2<-HDLP_IP_irf
p5<-ggplot(mapping=aes(x=0:(dim(o2$intervals)[1]-1)))+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_line(mapping=aes(y=o2$intervals[,2,]), color="blue")+
  geom_ribbon(mapping=aes(ymin=o2$intervals[,1,], ymax=o2$intervals[,3,]), alpha=0.5, fill="blue")+
  xlab("Horizon")+ylab("IP")


o3<-HDLP_CPI_irf
p6<-ggplot(mapping=aes(x=0:(dim(o3$intervals)[1]-1)))+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_line(mapping=aes(y=o3$intervals[,2,]), color="blue")+
  geom_ribbon(mapping=aes(ymin=o3$intervals[,1,], ymax=o3$intervals[,3,]), alpha=0.5, fill="blue")+
  xlab("Horizon")+ylab("CPI")

#################### CREATING PLOTS HDLP AND SVAR ########################
titled_row1 <- annotate_figure(
  ggarrange(
    p4, p5, p6,
    nrow = 1,
    ncol = 3
  ),
  top = text_grob("HDLP", 
                  color = "black", 
                  face = "bold", 
                  size = 14) 
)

titled_row2 <- annotate_figure(
  ggarrange(
    p1, p2, p3,
    nrow = 1,
    ncol = 3
  ),
  top = text_grob("SVAR", 
                  color = "black", 
                  face = "bold", 
                  size = 14) 
)

final_figure <- annotate_figure(ggarrange(
  titled_row1,
  titled_row2,
  nrow = 2
))

# Save plot
# png("hdlp_svar_application.png", width = 1200, height = 750, units = "px", res = 100)
