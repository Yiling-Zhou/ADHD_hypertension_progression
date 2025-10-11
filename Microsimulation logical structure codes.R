#-------------------------------------------------------------------------------
# Multistate Microsimulation Model for ADHD and Hypertension
#
# PURPOSE: This script outlines the R code used to run a continuous-time, 
# semi-Markov microsimulation model using the latent failure time approach.
# It simulates patient trajectories across multiple cardiovascular-related states
# based on previously fitted Cox proportional hazards models and a Poisson model.
#
# NOTE TO REVIEWERS: the original data and fitted model objects are not included. 
# objects are defined to ensure the script's structure and logic can be verified. 
# The script is NOT intended to be run directly to reproduce results.
#-------------------------------------------------------------------------------

# Load R packages
# data.table is loaded last to avoid masking functions from other packages.
library(ggplot2)
library(survminer)
library(survival)
library(tibble)
library(data.table)

#-------------------------------------------------------------------------------
# 1. Data and Model Import/Fitting
#-------------------------------------------------------------------------------

# NOTE: Original lines to load proprietary data and fitted models are commented 
# out and replaced with placeholders (e.g., loading a small dummy data frame).

# knitr::opts_chunk$set(root.dir = "H:/Multistate_ADHDandHTN")
# # cox models at each transition
# ls_model <- list()
# for(i in 1:14){
#   ls_model[[i]] <- readRDS(paste0("processData/adhd/model_transition", i , ".rds"))
# }
# # poisson model
# poisson_model <- readRDS("processData/adhd/fit_poisson.rds")
# # longitudinal co-variate data
# dt_cov <- fread("processData/adhd/dt_covariate.csv", colClasses = list(character = "pin"))

# Tidy data
dt_cov <- dt_cov[, .(pin, baseline_date, follow_end, tstart, tstop, exposure, on_adhd_med, sex, age, baseline_yr, baseline_yr_cat, t2d, hyperlipid)]
dt_cov[, sex:=factor(sex, levels = c("male", "female"))]
# population with ADHD before or at baseline
adhd_pop <- dt_cov[tstart==0&exposure==1,]
dt_cov <- dt_cov[pin%in%adhd_pop$pin, ]
dt_cov <- dt_cov[order(pin, tstart),]
rm(adhd_pop)


#-------------------------------------------------------------------------------
# 2. Multistate Structure Setting (Critical for Review)
#-------------------------------------------------------------------------------

tmat <- matrix(NA, 7, 7)
dimnames(tmat) <- list(from = c("hypertension", "CKD", "stroke","AMI", "HHF","CRdeath", "noCRdeath"),
                       to = c("hypertension", "CKD", "stroke", "AMI", "HHF", "CRdeath", "noCRdeath"))
tmat[1,] <- c(NA,1,2,3,4,5,6)
tmat[2,] <- c(NA,NA,NA,NA,NA,7,8)
tmat[3,] <- c(NA,NA,NA,NA,NA,9,10)
tmat[4,] <- c(NA,NA,NA,NA,NA,11,12)
tmat[5,] <- c(NA,NA,NA,NA,NA,13,14)
tmat[6,] <- c(NA,NA,NA,NA,NA,NA,NA)
tmat[7,] <- c(NA,NA,NA,NA,NA,NA,NA)

print(tmat)

# state space
state_space <- as.data.frame(tmat) 
setDT(state_space)
state_space <- melt(state_space, variable.name = "state", value.name = "trans", na.rm=TRUE)
state_space[, state:=as.character(state)]
# time horizon: 10 years (in days)
admin_end_time <- 10*365


#-------------------------------------------------------------------------------
# 3. Baseline Hazard ~10 years to 10 years
#-------------------------------------------------------------------------------
# NOTE: The code to import ls_basehaz.rds and the check loop are needed before.
ls_basehaz 

#-------------------------------------------------------------------------------
# 4. Core Simulation Functions
#-------------------------------------------------------------------------------

my_survival_fun <- function(t, bh, lp) {
  Ht <- subset(bh, time<=t)[, "hazard"]*exp(lp) 
  Ht_t <- Ht[length(Ht)] 
  return(exp(-(Ht_t)))
}


compute_event_time_vector_fun <- function(bh, lp_vector){
  
  u_vector <- runif(length(lp_vector)) 
  
  compute_event_time <- function(lp, u){
    
    event_time_try <- try(uniroot(\(x) my_survival_fun(t=x, bh, lp) - u, 
                                  interval = c(1, admin_end_time)), silent = TRUE)
    
    if(class(event_time_try) == "try-error") {event_time <- admin_end_time + 1} else {event_time <- event_time_try$root}
  }
  
  
  event_time_vector <- sapply(seq_along(lp_vector), function(k){
    compute_event_time(lp_vector[k], u_vector[k])
  })
  
  return(event_time_vector)
}


#-------------------------------------------------------------------------------
# 5. Simulation Data Setup
#-------------------------------------------------------------------------------

# Create simulation datasets by fixing exposure (ADHD status/medication use)

### Simulation for adhd in ADHD population
dt_cov_adhd <- copy(dt_cov)
dt_cov_adhd <- dt_cov_adhd[, exposure:=1] # Set exposure for ADHD group
dt_cov_adhd_last <- dt_cov_adhd[order(pin, tstart),][,.SD[.N], by=pin][tstop<admin_end_time, ][, tstart:=tstop][,tstop:=admin_end_time]
dt_cov_adhd <- rbind(dt_cov_adhd,dt_cov_adhd_last)
dt_cov_adhd <- dt_cov_adhd[order(pin, tstart),]
rm(dt_cov_adhd_last)


### Simulation for non adhd in ADHD population
dt_cov_noadhd <- copy(dt_cov)
dt_cov_noadhd <- dt_cov_noadhd[, exposure:=0][, on_adhd_med:=0] # Set exposure to 0 and no medication for non-ADHD scenario
dt_cov_noadhd_last <- dt_cov_noadhd[order(pin, tstart),][,.SD[.N], by=pin][tstop<admin_end_time, ][, tstart:=tstop][,tstop:=admin_end_time]
dt_cov_noadhd <- rbind(dt_cov_noadhd,dt_cov_noadhd_last)
dt_cov_noadhd <- dt_cov_noadhd[order(pin, tstart),]
rm(dt_cov_noadhd_last)
rm(dt_cov) 


#-------------------------------------------------------------------------------
# 6. Parameter Selection Preparation
#-------------------------------------------------------------------------------

########### poisson model #####################
coef_poisson <- poisson_model$coefficients 
# covariance_poisson <- vcov(poisson_model) # PSA
formula_poisson <- poisson_model$formula

############ cox model #####################
ls_coef_cox <- list()
#ls_covariance_cox <- list()
ls_formula_cox <- list()
for(m in 1:14) {ls_coef_cox[[m]] <- ls_model[[m]]$coefficients} 
#for(m in 1:14) {ls_covariance_cox[[m]] <- ls_model[[m]]$var2} # PSA
for(m in 1:14) {ls_formula_cox[[m]] <- as.formula(paste0("~", as.character(ls_model[[m]]$formula)[3]))}


#-------------------------------------------------------------------------------
# 7. Simulation Loop (Point Estimate)
#-------------------------------------------------------------------------------

# This section is the core of the methodology and is preserved entirely.
# Note: 'coef_poisson' and 'ls_coef_cox' are used directly (point estimates).
# parallel processing
# comuptational burden

library(foreach)
library(doParallel)
# cores is set to a placeholder value for demonstration
registerDoParallel(cores = 4) 
#system time
start_time <- Sys.time()
p <- 0 # Fixed value for point estimate simulation
# Use point estimates directly
coef_poisson_S <- coef_poisson
ls_coef_cox_S <- ls_coef_cox

sim_output <- foreach(S = 1:500, 
                      .combine = "rbind",
                      .packages = c("data.table", "survival", "purrr"))%dopar%{
                        
                        #####################Part 1: ADHD#################                                          
                        #################################################                    
                        sim_table_0 <- data.table(pin = unique(dt_cov_adhd$pin), time = 0, state = "hypertension",  record = 1, sim = S, adhd = "ADHD", PSA=p)                     
                        ####################stage 1#######################
                        ####1. prepare covariate data first transition to six states
                        dt_cov_1 <- dt_cov_adhd[tstart == 0, ]
                        dt_cov_1[, hypertension_duration:= tstart/365]
                        ## update age, baseline_yr alongside time_k
                        dt_cov_1[, age:=age+as.integer(tstart/365)][, baseline_yr:=baseline_yr+as.integer(tstart/365)]
                      
                        #####2. update hypertension medication use by predicting rate using poisson model
                        dt_cov_1[, time:=1] # 1 day for Poisson model
                        dt_cov_1[, trans_state:="hypertension"] # strata for Poisson model
                        dt_cov_1[, trans_state:=factor(trans_state, levels=c("AMI", "CKD", "HHF", "hypertension", "stroke"))]
                        
                        ## design matrix
                        designMatrix <- model.matrix(formula_poisson, data = dt_cov_1, contrasts.arg = poisson_model$contrasts)
                        ## predict using poisson model
                        ln_y <- designMatrix %*% coef_poisson
                        
                        dt_cov_1[, hyper_med_rate:= exp(ln_y[,1])]
                        
                        dt_cov_1[, hyper_med_p:=1-exp(-hyper_med_rate*1)]
                        
                        dt_cov_1[, hypertension_med:= rbinom(1, 1, hyper_med_p), by = pin] ### first-order error
                        
                        
                        dt_cov_1[,time:=NULL]
                        dt_cov_1[,trans_state:=NULL]
                        dt_cov_1[,hyper_med_rate:=NULL]
                        dt_cov_1[,hyper_med_p:=NULL]
                        rm(designMatrix)
                        rm(ln_y)
                        dt_cov_1[, tstart:=NULL]
                        dt_cov_1[, tstop:=NULL]
                        
                        ######3. predict LP using cox model
                        for(j in 1:6){
                          bh <- ls_basehaz[[j]]
                          
                          ## prepare matrix
                          designMatrix <- model.matrix(ls_formula_cox[[j]], data = dt_cov_1, contrasts.arg = ls_model[[j]]$contrasts)
                          ## predict using cox model
                          lp_matrix <- (designMatrix[, -1] - rep(ls_model[[j]]$means, each = nrow(designMatrix))) %*% ls_coef_cox[[j]]
                          
                          dt_cov_1[, lp := lp_matrix[,1]]
                          
                          rm(designMatrix)
                          rm(lp_matrix)
                          
                          event_time_vector <- compute_event_time_vector_fun(bh = bh, lp_vector = dt_cov_1$lp)
                          name_event <- state_space[trans==j,]$state
                          dt_cov_1[, (name_event):= event_time_vector]
                          
                          dt_cov_1[, lp:=NULL]
                          print(j)
                          rm(bh)
                          rm(event_time_vector)
                          rm(name_event)
                        }
                        
                        rm(j)
                        #####4. determine the minimum event time and next state
                        columns <- dput(names(na.omit(tmat["hypertension", ]))) 
                        
                        dt_cov_1[, event_time:= do.call(pmin, .SD), .SDcols = columns]
                        ### sim table
                        dt_cov_1[, state:= columns[max.col(-.SD)], .SDcols = columns]
                        dt_cov_1[event_time==admin_end_time + 1, `:=` (state="survival", event_time= admin_end_time)]
                        
                        #####5. extract the result table for the first transition
                        sim_table_1 <- dt_cov_1[, .(pin, state, event_time)]
                        sim_table_1[, record:=2]
                        
                        setnames(sim_table_1, old ="event_time", new ="time")
                        rm(dt_cov_1)
                        
                        ####################stage 2#######################
                        dt_cov_2 <- sim_table_1[!state%in%c("survival", "CRdeath", "noCRdeath"),][,.(pin, state, time)]
                        dt_cov_2 <- dt_cov_2[time<admin_end_time,]
                        dt_cov_2 <- dt_cov_adhd[dt_cov_2, on=.(pin, tstop > time, tstart <=time)]
                        dt_cov_2 <- dt_cov_2[order(pin, tstart),][,.SD[1], by = pin]
                        dt_cov_2[, hypertension_duration:= tstart/365]
                        ## update age, baseline_yr alongside time_k
                        dt_cov_2[, age:=age+as.integer(tstart/365)][, baseline_yr:=baseline_yr+as.integer(tstart/365)]
                
                        #####2. update hypertension medication use by predicting rate using poisson model
                        dt_cov_2[, time:=1] # 1 day for Poisson model
                        dt_cov_2[, trans_state:=state] # strata for Poisson model
                        dt_cov_2[, trans_state:=factor(trans_state, levels=c("AMI", "CKD", "HHF", "hypertension", "stroke"))]
                        ###similar loop until the endpoint see the supplement figure, and finally calculate the continuous time from baseline.
                        
                        ##############################################################################    
                        ##################### Part 2: non-ADHD ########################################                                
                        ###############################################################################            
                        ###similar loop as Part 1 ADHD, until the endpoint see the supplement figure
                       
                        
                        output <- rbind(sim_table_adhd, sim_table_noadhd)
                        rm(sim_table_adhd)
                        rm(sim_table_noadhd)
                        
                        output
                      }

# save output
