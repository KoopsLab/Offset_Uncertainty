################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# Efficacy of a quantitative method to define biodiversity offset multipliers to 
# account for uncertainty in complex projects. 
#
# R code exploring a method to estimate uncertainty compensation ratios for
# complex offsetting projects
# Based off of Bradford (2017)
# The method is extended to allow for multiple impacts and offsets and allow
# for CR multiplied to be applied differently (in different amounts) across
# offsets
# Simulations explore:
# - the efficacy of the approach across level of uncertainty,
#   correlation, number of projects, and equivalency levels
# - the impact of offset weights and uncertainty have on CRs
# - the impact of correlation among impacts and offsets
# - what are expected CRs across levels of uncertainty for data-limited 
#   scenarios

#-------------------------------------------------------------------------------

# clear workspace
rm(list=ls())

# load libraries
library(data.table)
library(ggplot2)
library(parallel)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#-------------------------------------------------------------------------------
# Functions & ggplot theme
#-------------------------------------------------------------------------------

# ggplot theme
theme_me <- theme_bw() +
  theme(axis.title = element_text(size = 14, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 12, family = "sans", colour ="black"),
        axis.text.y = element_text(size = 12, family = "sans", hjust = 0.6,
                                   angle = 90, colour = "black"),
        legend.title = element_text(size = 14, family = "sans"),
        legend.text = element_text(size = 14, family = "sans"),
        strip.text = element_text(size = 14, family = "sans"))

#-------------------------------------------------------------------------------

# draws correlated random variables - N(0,1) - normal distn with mean 0 and sd 1
var.cor.draw <- function(n, cor.mat){
  lapply(1:n, function(i){
    rawvars <- MASS::mvrnorm(n = 1, mu = rep(0,nrow(cor.mat)), Sigma = cor.mat)
    pvars <- pnorm(rawvars, mean = 0, sd = 1)
  })
}

# Function to generate random distributions for impact(s) and offset(s)
dist.f <- function(
    mean.i = 1,          # Impact mean(s) - arithmetic 
    mu.i = NA,           # impact geometric mean
    CV.i = 0.5,          # Impact CV - uncertainty in mean
    mean.o = 1,          # Offset mean
    mu.o = NA,           # offset geometric mean
    CV.o = 0.5,          # Offset CV - uncertainty in mean
    cor = NA,            # correlation among impact and offset - matrix, vector of length 2, or NA 
                         # - NA means no correlation among impact or offset
                         # - 2 element vector gives correlation among ipmacts and offsets separately, impact and offset are uncorrelated
                         # - matrix is the correlation matrix dim(ni+no, ni+no)
    n = 10000            # number of draws
) {
  
  n.i <- length(mean.i) # number of impacts
  n.o <- length(mean.o) # number of offsets
  
  # Find moments of distributions - for log-normal dist. 
  sigma.i <- sqrt(log(CV.i^2+1))      # log-sd of distributions
  if(is.na(mu.i)){
    mu.i <- log(mean.i) - sigma.i^2/2 # log-mean of distributions
  }
  
  sigma.o <- sqrt(log(CV.o^2+1))      # log-sd of distributions
  if(is.na(mu.o)){
    mu.o <- log(mean.o) - sigma.o^2/2 # log-mean of distributions
  }
  
  # Correlation matrix - among impacts and offsets
  # if correlation provided by user 
  if(is.matrix(cor)){           
    cor.mat = cor
  # if a vector of length 2 provided - values represent correlation among
  # impacts and correlation among offsets; no correlation bewteen impacts and
  # offsets
  } else if(length(cor) == 2) { 
    cor.mat <- matrix(0, nrow = (n.i+n.o), ncol = (n.i+n.o))
    cor.mat[1:n.i,1:n.i] = cor[1]
    cor.mat[(n.i + 1):(n.i+n.o),(n.i + 1):(n.i+n.o)] = cor[2]
    diag(cor.mat) <- 1 
  # if NA provided - no correlation
  } else if(is.na(cor)) {      
    cor.mat = diag(n.i+n.o)
  } 
  
  # Correlated random deviates 
  norm.vars <- var.cor.draw(n, cor.mat) # draw random deviates with correlations
  norm.vars.i <- lapply(norm.vars, function(x) x[1:n.i])               # impacts  
  norm.vars.o <- lapply(norm.vars, function(x) x[(n.i + 1):(n.o+n.i)]) # offsets
  
  # Covert to prob. dist'n selected
  # Requires calc of dist'n params based on provided means and sds
  d.i <- do.call(rbind, lapply(norm.vars.i, function(x){ # impacts
    qlnorm(x, meanlog = mu.i,  sdlog = sigma.i)
  }))
  d.o <- do.call(rbind, lapply(norm.vars.o, function(x){ # offsets
    qlnorm(x, meanlog = mu.o,  sdlog = sigma.o)
  }))
  
  # output
  list(impacts = d.i, # matrix of impacts
       offsets = d.o) # matrix of offsets
  
}

# Function to simulated and calculated Compensation Ratios
CRu.f <- function(
    d.i,     # impact distributions
    d.o,     # offset distributions
    p = 0.8, # Risk tolerance threshold,
    CR_weight# offset weights - the ability to apply CR to an offset (0:1)
) {
  
  # sum offsets across 
  di <- apply(d.i, 1, sum)
  
  # sum offsets across replicates
  do <- apply(d.o, 1, sum)
  
  # CR distributions 
  M = di/do   # CR assuming applied across all offsets
  # CR dist for weighted offsets
  M.weight = 1 + (di - do) / apply(d.o, 1, function(x) sum(CR_weight * x)) 
 
  # Exceedancy thresholds
  ET <- round(quantile(di, p),2)
  
  # Compensation Ratios
  # Assuming pplied across all offsets   
  CR <- round(quantile(M, p) ,2)[[1]]
  
  # Weighted CR by available offsets
  CR.adjusted <- round(quantile(M.weight, p) , 2)[[1]] # values
  CR.adjusted <- (CR.adjusted - 1) * CR_weight + 1     # turn into vector
  
  # Output
  list("EQ" = p,                      # Equivalency threshold
       "ET" =  ET,                    # Exceedance threshold
       "CR - All" = CR,               # CR with equal weights among offsets
       "CR - Adjusted" = CR.adjusted  # CR with defined weights
  )

}

# Function to compare the efficacy of offset relative to impact over years
CR_comp_f <- function(
  d.i,
  d.o
) {
  
  # Sum impacts across replicates
  di <- apply(d.i, 1, sum)
  
  # sum offsets across replicates
  do <- apply(d.o, 1, sum)
  
  # Output
  list(
    Impact = sum(di),        # Total impact
    Offset = sum(do),        # Total offset
    ratio = sum(do)/sum(di), # ratio of offset to impact
    success = (do) >= (di)   # Boolean is offset > impact
  )
}

#-------------------------------------------------------------------------------
# TEST - display distributions 
#-------------------------------------------------------------------------------

# Parameter 
exData  <- list(
  mean.i = 6,             # Mean of impact - one impact
  CV.i = 0.5,             # uncertainty (CV) of impact
  mean.o = c(2, 2, 2),    # Mean of Offsets - 3 offset projects
  CV.o = c(0.4, 0.1, 0.8),# uncertainty (CV) of offset projects
  cor = NA,               # correlation among impact and offset - uncorrelated
  n = 10000,              # number of Monte Carlo draws
  p = 0.8,                # Equivalency threshold - prob NNL
  CR_weight = c(1, 1, 1)  # Offset CR weights
)

# Initial Impact/offset distributions
dist.data <- dist.f(
  mean.i = exData$mean.i, # mean impact
  CV.i = exData$CV.i,     # CV impact
  mean.o = exData$mean.o, # mean offset   
  CV.o = exData$CV.o,     # CV offset
  cor = exData$cor,       # correlation parameters     
  n = exData$n            # number of simluaitons
)

# calculate compensation ratios (offset multipliers) for uncertainty
CR = CRu.f(
  d.i = dist.data$impacts,     # impact dist.
  d.o = dist.data$offset,      # offset dist.
  p = exData$p,                # risk level
  CR_weight = exData$CR_weight # offset weights
)

# multiplier values
CR[["CR - All"]]      # ws = 1   
CR[["CR - Adjusted"]] # by project

# PLOT distributions
plot.data <- data.table(
  x = c(apply(dist.data$impact, 1, sum),
        apply(dist.data$offset, 1, sum),
        apply(dist.data$offset * CR[["CR - All"]], 1, sum)
        ),
  type = c(rep("impact", 10000),
           rep("offset - no CR", 10000),
           rep("offset - CR all", 10000)
           )
  )
plot.data[, mean(x), by = type]
plot.data[, sd(x), by = type]

ggplot(plot.data, aes()) +
  geom_density(aes(x, fill = type), alpha = 0.5) +
  theme_me

# clean up
rm(exData, dist.data, CR, plot.data)

#-------------------------------------------------------------------------------
# Simulation to compare the probability of NNL for difference levels of 
# uncertain with a multiplier is and is not applied
#-------------------------------------------------------------------------------

# Parameters 
n.i <- 1                # Number of impacts
n.o <- 1:3              # number of offsets
CV.i <- c(0.25, 0.5, 1) # uncertainty of impacts
CV.o <- c(0.25, 0.5, 1) # uncertainty of offsets
cor <- c(0, 0.5, 0.95)  # correlation among impacts and offsets

# parameter combos
params <- tidyr::crossing(n.i, n.o,
                          CV.i, CV.o,
                          cor)

# Run loops - in parallel
cl <- makeCluster(9)    # make clusters
clusterExport(cl, ls()) # export everything to clusters - careful

# run sim
CR_comp_sim <- parLapply(cl, 1:nrow(params), function(i){
  
  require(data.table)
  
  # Extract parametesr
  mean.i = rep(1/params$n.i[i], params$n.i[i]) # Impact mean(s) - arithmetic 
  CV.i = rep(params$CV.i[i], params$n.i[i])    # Impact CV - uncertainty in mean
  mean.o = rep(1/params$n.o[i], params$n.o[i]) # Offset mean
  CV.o = rep(params$CV.o[i], params$n.o[i])    # Offset CV - uncertainty in mean
  cor = simstudy::genCorMat(params$n.i[i] + params$n.o[i], # correlation matrix
                            rho = params$cor[i])                 
  CR_weight = rep(1, params$n.o[i])            # offset weights
  
  n = 100000        # number of draws
  # Generate initial impact and offset distributions
  dist.data <- dist.f(
    mean.i = mean.i, # impact mean
    CV.i =  CV.i,    # impact uncertainty (CV)
    mean.o = mean.o, # Offset mean
    CV.o = CV.o,     # Offset uncertainty (CV)
    cor = cor,       # correlation matrix
    n = n            # number of draws
  )
  
  # Prob of no net loss with no CR
  # no. of offset > impacts
  Psuccess_noCR <- sum(apply(dist.data$offsets,1,sum) > 
                         apply(dist.data$impacts,1,sum)) / n
  
  # Identify CR - for different risk tolerances
  ps <- c(0.8, 0.9, 0.95)
  CR_by_p <- lapply(ps, function(p){
    CR_sim <- CRu.f( 
      d.i = dist.data$impacts, # impacts distributions
      d.o = dist.data$offsets, # offset distributions
      p = p,                   # equivalency threshold,
      CR_weight = CR_weight    # multipliers weights
    ) 
    CR_sim
  })
  
  # Prob of no net loss with CR
  Psuccess_CR <- lapply(CR_by_p, function(cr){# loop over p values
    # output data table
    data.table(
      EQ = cr$EQ,         # equivalency threshold - 1- risk tolerance
      ET = cr$ET,         # exceedance threshold
      CR = cr$`CR - All`, # multiplier with ws = 1
      # prob of achieving NNL
      Psuccess = sum(apply(sweep(dist.data$offsets, 2, 
                                 cr[["CR - Adjusted"]], "*"),1,sum) > 
                     apply(dist.data$impacts,1,sum)) / n,
      # prob of > ET
      Peq = sum(apply(sweep(dist.data$offsets, 2, cr[["CR - Adjusted"]], "*"), 
                      1, sum) > 
                  cr[["ET"]]) / n
    )
  })
    
  # Output 
  list(params = params[i,],          # sim parameter
       Psuccess_noCR = Psuccess_noCR,# NNL with no multiplier
       Psuccess_CR = Psuccess_CR)    # NNL with multiiplier

})
stopCluster(cl) # close clusters

# Extract Data
CR_comp_data <- do.call(rbind, lapply(CR_comp_sim, function(x){
  
  data.table(cbind(x$params,
        p = round(1-c(0.8, 0.9, 0.95),2),
        "Psuccess_noCR" = x$Psuccess_noCR,
        "Psuccess_CR" = sapply(x$Psuccess_CR, function(y)y$Psuccess),
        "CR" = sapply(x$Psuccess_CR, function(y) unlist(y$CR))
        ))
}))
CR_comp_data$n.o <- paste("# of Offsets:", CR_comp_data$n.o)

# PLOT
##png("Fig1.png", width = 6.5, height = 3.5, unit = "in", res = 300)
ggplot() +
  geom_line(data = CR_comp_data[p == 0.20 & CV.i == 0.25,], 
            aes(CV.o, Psuccess_noCR, linetype = as.factor(cor))) +
  geom_line(data = CR_comp_data[CV.i == 0.25], 
            aes(CV.o, Psuccess_CR, 
                colour = as.factor(p), 
                linetype = as.factor(cor))) +
  
  facet_grid(.~n.o) +
  labs(x = "Offset Uncertainty (CV)", y = "Probability of NNL", 
       colour = "Risk Tolerance Threshold",
       linetype = "Correlation") +
  theme_me + theme(legend.position = 'bottom',
                   legend.box="vertical",legend.margin=margin())
#dev.off()

# Plot - CR with uncertainty
ggplot() +
  geom_line(data = CR_comp_data[CV.i == 0.25], 
            aes(CV.o, CR, 
                colour = as.factor(p), 
                linetype = as.factor(cor))) +
  
  facet_grid(.~n.o) +
  labs(x = "Offset Uncertainty (CV)", y = "Offset Multiplier", 
       colour = "Risk Tolerance Threshold",
       linetype = "Correlation") +
  theme_me + theme(legend.position = 'bottom',
                   legend.box="vertical",legend.margin=margin())

#-------------------------------------------------------------------------------
# Simulation to compare the probability of NNL when different weights are 
# applied to offsets with various uncertainties
#-------------------------------------------------------------------------------

# Parameters
n.i <- 1               # Number of impacts
n.o <- 2               # Number of offsets
CV.i <- 0.25           # uncertainty of impact
cor <- c(0, 0.5, 1)    # correlation between impacts and offsets
weight.n <- 1:3        # CR weight option 3 scenarios
CV.o <- 1:2            # Uncertainty of impact option - 2 scenarios
p = c(0.8, 0.90, 0.95) # risk tolerance 

# CR weight scenarios
CR_weight.list = list(
  c(1, 1), # equal CRs
  c(1, 0), # CR only applied to offset 1
  c(0, 1)) # CR only applies to offset 2

# Offset uncertainty scenarios
CV.o.list = list(
  c(0.5, 0.5),  # Equal uncertainty
  c(1.0, 0.1))  # Offset 1 in more uncertain than offset 2

# Param combos
params <- tidyr::crossing(n.i, n.o,
                          CV.i, CV.o,
                          cor,
                          p,
                          weight.n)

# Run simulation
cl <- makeCluster(9)
clusterExport(cl, ls())

CR_comp_sim <- parLapply(cl, 1:nrow(params), function(i){
  
  require(data.table)
  
  mean.i = rep(1/params$n.i[i], params$n.i[i])   # Impact mean(s) - arithmetic 
  CV.i = rep(params$CV.i[i], params$n.i[i])      # Impact CV 
  mean.o = rep(1/params$n.o[i], params$n.o[i])   # Offset mean
  CV.o = CV.o.list[[params$CV.o[i]]]             # Offset CV 
  cor = simstudy::genCorMat(params$n.i[i] + params$n.o[i], # correlation matrix
                            rho = params$cor[i]) # correlation matrix
  CR_weight = rep(1, params$n.o[i])              # offset weights
  n = 100000                                     # number of draws
  p = params$p[i]                                # equivalency threshold 
  years = 1
  CR_weight <- CR_weight.list[[params$weight.n[i]]]
  
  # Generate initial impact and offset distributions
  dist.data <- dist.f(
    mean.i = mean.i, # impact mean
    CV.i =  CV.i,    # impact CV - uncertainty in mean
    mean.o = mean.o, # Offset mean
    CV.o = CV.o,     # Offset CV - uncertainty in mean
    cor = cor,       # correlation matrix
    n = n            # number of draws
  )
  
  # Prob of no net loss with no CR
  Psuccess_noCR <- sum(apply(dist.data$offsets,1,sum) > 
                         apply(dist.data$impacts,1,sum)) / n
  
  # Calc CR
  CR_sim <- CRu.f( 
    d.i = dist.data$impacts, # impact distributions(s)
    d.o = dist.data$offsets, # offset distribution(s)
    p = p,                   # equivalency threshold,
    CR_weight = CR_weight    # weights of CR multiplier by offset
  ) 
  
  # Prob of no net loss with CR
  Psuccess_CR = sum(
    apply(sweep(dist.data$offsets, 2, CR_sim[["CR - Adjusted"]], `*`),1,sum) > 
                   apply(dist.data$impacts,1,sum)) / n
  
    # Output
  list(params = params[i,],          # sim paramerters
       CR_sim = CR_sim,              # compendation multipiers
       Psuccess_noCR = Psuccess_noCR,# prob NNL without CR
       Psuccess_CR = Psuccess_CR)    # prob NNL with CR

})
stopCluster(cl)

# Organize into dataframe
CR_comp_data <- do.call(rbind, lapply(CR_comp_sim, function(x){
  
  data.table(cbind(x$params,                              # sim parameters
                   "CR" = x$CR_sim[["CR - All"]],         # CR for both offsets (weights equal)
                   "CR1" = x$CR_sim[["CR - Adjusted"]][1],# CR for offset 1
                   "CR2" = x$CR_sim[["CR - Adjusted"]][2],# CR for offset 2
                   "Psuccess_noCR" = x$Psuccess_noCR,     # prob NNL w/o CR
                   "Psuccess_CR" = x$Psuccess_CR          # prob NNL w CR
  ))
}))

# Set offset CV to catergorical indicator
CR_comp_data$CV.o <- ifelse(CR_comp_data$CV.o == 1,
                            "Equal",  # equal uncertinaty
                            "Unequal")# unequal uncertinaty O1 > O2

# Set CR weigths to categorical indicator
CR_comp_data$weight.n <- ifelse(CR_comp_data$weight.n  == 1,
                            "Both",                           # Both get CR
                            ifelse(CR_comp_data$weight.n  == 2, 
                                   "Offset 1",            # O1 gets CR
                                   "Offset 2")            # O2 gets CR
)

CR_comp_data[,"CR_total" := CR1+CR2]

# convert EQ to risk tolerance (1-p)
CR_comp_data$p <- round(1-CR_comp_data$p, 2)

# PLOt
#png("Fig2.png", width = 6.5, height = 3.25, unit = "in", res = 300)
ggplot() +
  geom_point(data = CR_comp_data[cor == 0], 
             aes(x = as.factor(CV.o),
                 y = CR_total, 
                 #shape = as.factor(p),
                 colour = as.factor(p))
  ) +
  facet_grid(.~weight.n)+
  labs(x = "Uncertainty", 
       y = "Offset Multiplier", 
       colour = "Risk Tolerance") +
  theme_me + theme(legend.position = 'bottom',
                   legend.box="vertical",legend.margin=margin())
#dev.off()

#-------------------------------------------------------------------------------
# Simulation to investigate the impact of correlation on CRs
#-------------------------------------------------------------------------------

# Parameters 
n.i <- 3                # Number of impacts
n.o <- 3                # number of offsets
CV.i <- c(0.1, 0.5, 1)  # uncertainty of impacts
CV.o <- c(0.1, 0.5, 1)  # uncertainty of offsets
cor <- seq(0.0, 1, 0.1) # correlation among impacts and offsets

# parameter combos
params <- tidyr::crossing(n.i, n.o,
                          CV.i, CV.o,
                          cor)

# Run loops
cl <- makeCluster(15)
clusterExport(cl, ls())

CR_comp_sim <- parLapply(cl, 1:nrow(params), function(i){
  
  require(data.table)
  
  mean.i = rep(1/params$n.i[i], params$n.i[i])   # Impact mean(s) - arithmetic 
  CV.i = rep(params$CV.i[i], params$n.i[i])      # Impact CV 
  mean.o = rep(1/params$n.o[i], params$n.o[i])   # Offset mean
  CV.o = rep(params$CV.o[i], params$n.o[i])      # Offset CV 
  cor =  simstudy::genCorMat(length(mean.i)+length(mean.o), 
                             rho = params$cor[i])# correlation value
  p = 0.8                                        # equivalency threshold 
  CR_weight <- rep(1, params$n.i[i])             # offset weights
  n = 100000                                     # number of draws
  
  # Generate initial impact and offset distributions
  dist.data <- dist.f(
    mean.i = mean.i, # impact mean
    CV.i =  CV.i,    # Impact CV - uncertainty in mean
    mean.o = mean.o, # Offset mean
    CV.o = CV.o,     # Offset CV - uncertainty in mean
    cor = cor,       # correlation matrix
    n = n            # number of draws
  )

  # Prob of no net loss with no CR
  Psuccess_noCR <- sum(apply(dist.data$offsets,1,sum) > 
                         apply(dist.data$impacts,1,sum)) / n
  
  # Identify CR - for different risk tolerances
  ps <- c(0.8, 0.9, 0.95)
  CR_by_p <- lapply(ps, function(p){
    CR_sim <- CRu.f( 
      d.i = dist.data$impacts, # impacts distributions
      d.o = dist.data$offsets, # offset distributions
      p = p,                   # equivalency threshold,
      CR_weight = CR_weight    # multipliers weights
    ) 
    CR_sim
  })
  
  # Prob of no net loss with CR
  Psuccess_CR <- lapply(CR_by_p, function(cr){
    list(
      EQ = sapply(CR_by_p, function(x) x$EQ),             # Equivalency threshold
      ET = sapply(CR_by_p, function(x) x$ET),             # Exceedance threshold
      CR = lapply(CR_by_p, function(x) x$`CR - Adjusted`),# compensation ratio
      Psuccess = sum(                                     # Prob NNL
        apply(sweep(dist.data$offsets, 2, cr[["CR - Adjusted"]], "*"),1,sum) > 
          apply(dist.data$impacts,1,sum)) / n
    )
  })
  
  # Output
  list(params = params[i,],          # parameters
       CR_sim = CR_by_p,             # compendation ratio
       Psuccess_noCR = Psuccess_noCR,# Prob NNL w/o CR
       Psuccess_CR = Psuccess_CR)    # Prob NNL w CR
  
})
stopCluster(cl)

# Extract results into data frame
CR_comp_data <- do.call(rbind, lapply(CR_comp_sim, function(x){
  
  data.table(cbind(
    x$params,                                                    # sim params
    "EQ" = sapply(x$CR_sim, function(cr) cr$EQ),                 # p value
    "CR" = sapply(x$CR_sim, function(cr) cr$`CR - All`),         # multiplier
    "Psuccess_noCR" = x$Psuccess_noCR,                           # prob NNL w/o CR
    "Psuccess_CR" = sapply(x$Psuccess_CR, function(P) P$Psuccess)# prob NNL w CR
  ))
}))

# Convert EQ to risk tolerance and assign as factor with risk tolerance decreasing
CR_comp_data$EQ <- paste("r =", 1-CR_comp_data$EQ)
CR_comp_data$EQ <- factor(CR_comp_data$EQ, 
                          levels = c("r = 0.2", "r = 0.1", "r = 0.05"),
                          labels = c(expression(italic(r) == 0.2),
                                     expression(italic(r) == 0.1),
                                     expression(italic(r) == 0.05)))

CR_comp_data[CV.i == 0.1 & CV.o == 1.0 & (cor == 0 | cor == 1.0)]
# Plot
#png("Fig3.png", width = 6.5, height = 3.5, unit = "in", res = 300)
ggplot() +
  geom_line(data = CR_comp_data, 
            aes(cor, CR, 
                linetype = as.factor(CV.i), 
                colour = as.factor(CV.o))) +
  facet_grid(.~EQ, labeller = label_parsed) +
  labs(x = "Correlation", y = "Offset Multiplier", 
       linetype = "Impact Uncertainty (CV)",
       colour = "Offset Uncertainty (CV)") +
  theme_me + theme(legend.position = 'bottom',
                   legend.box="vertical",legend.margin=margin())
#dev.off()

#png("Fig3a.png", width = 6.5, height = 6.5, unit = "in", res = 300)
ggplot() +
  geom_line(data = CR_comp_data, 
            aes(cor, CR, 
                colour = (EQ))) +
  facet_grid(CV.i~CV.o, labeller = label_parsed) +
  labs(x = "Correlation", y = "Offset Multiplier", 
       colour = "Offset Uncertainty (CV)") +
  scale_colour_discrete(labels = function(x) parse(text = x)) +
  theme_me + theme(legend.position = 'bottom',
                   legend.box = "vertical", legend.margin = margin())
#dev.off()

#-------------------------------------------------------------------------------
# Simulation to calculate CR for different levels or impact/offset uncertainty
# assume equivalency, no correlation and all offset weights = 1. 
#-------------------------------------------------------------------------------

# Parameters
EQ <- c(0.8, 0.9, 0.95)  # risk tolerances 5 - 20%
CV.i <- seq(0.00,1,0.05) # impact uncertainty (CV) - 0:1
CV.o <- seq(0.00,1,0.05) # offset uncertainty (CV) - 0:1

# param combos
CV_mat <- tidyr::crossing(EQ, CV.i,CV.o)

# Run Simulations
no_cores <- 14               # number of cores
cl <- makeCluster(no_cores)  # create clusters
clusterExport(cl, ls()) # send data to clusters

CV_mat$CR <- parSapply(cl, 1:nrow(CV_mat), function(i) {
  
  require(data.table)
  
  mean.i = 1            # Impact mean(s) - arithmetic 
  CV.i = CV_mat$CV.i[i] # Impact CV - uncertainty in mean
  mean.o = 1            # Offset mean
  CV.o = CV_mat$CV.o[i] # Offset CV - uncertainty in mean
  cor =  NA             # correlation value
  p = CV_mat$EQ[i]      # equivalency threshold 
  CR_weight <- 1        # Multiplier application
  
  # Generate initial impact and offset distributions
  dist.data <- dist.f(
    mean.i = mean.i, # Impact mean
    CV.i =  CV.i,    # Impact uncertainty
    mean.o = mean.o, # Offset mean
    CV.o = CV.o,     # Offset CV - uncertainty in mean
    cor = cor,       # correlation value - none
    n = 100000       # number of draws
  )

  # Identify CR - for different risk tolerances
  CR_sim <- CRu.f( 
    d.i = dist.data$impacts, # impacts distributions
    d.o = dist.data$offsets, # offset distributions
    p = p,                   # equivalency threshold,
    CR_weight = CR_weight    # multipliers weights
  ) 
  CR_sim$`CR - All`

})
stopCluster(cl)

# round CR values to 0.5 increments 
CV_mat$CR.round <- round(CV_mat$CR/0.5)*0.5

# Convert EQ to risk tolerance and assign as factor with risk tolerance decreasing
CV_mat$EQ <- paste("r =", 1-CV_mat$EQ)
CV_mat$EQ <- factor(CV_mat$EQ, 
                    levels = c("r = 0.2", "r = 0.1", "r = 0.05"),
                    labels = c(expression(italic(r) == 0.2),
                               expression(italic(r) == 0.1),
                               expression(italic(r) == 0.05)))



# PLOT
#png("Fig4.png", width = 4, height = 8, unit = "in", res = 300)
ggplot(CV_mat) +
  geom_tile(aes(x = CV.i, y = CV.o, fill = as.factor(CR.round)))+
  scale_fill_discrete() + 
  labs(x = "Impact Uncertainty (CV)", y = "Offset Uncertainty (CV)", 
       fill = "m")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.2)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))+
  facet_wrap(~EQ, ncol = 1, labeller = label_parsed) +
  theme_me + theme(legend.title = element_text(face = "italic"))
#dev.off()

#-------------------------------------------------------------------------------
# EXAMPLE CALCULATION
#-------------------------------------------------------------------------------
# example calculation for manuscript
# 5 ha of riparian habitat impacted
# mean = 5, CV = 0.1
# Two offset projects - riparian restoration % off channel pool
# riparian - PERT uncertainty
# mode = 0.25, min = 0, max = 0.75, cost = 17$/m2
# pool - log-normal with CV
# mean = 1, CV = 0.5, cost = 85$/m2
# optimize cost of offsetting by adjusting wj
#-------------------------------------------------------------------------------

# ---- PARAMETERS ----

n.sim <- 100000 # number of simulations
p <- 0.9        # equivalency threshold - 1- risk tolerance

# Impact
mean.i <- 5 # mean impact 5 ha of
cv.i <- 0.1 # uncertainty of impact

# offset 1 - riparian restoration
mode.o1 <- 0.2      # most likely value
min.o1 <- 0         # minimum value
max.o1 <- 1         # maximum value
cost.o1 <- 17*10000 # cost per hectare (Theis et al. 2022)

# offset 2 - off-channel pool
mean.o2 <- 1        # mean erected value
cv.o2 <- 0.5        # uncertainty 
cost.o2 <- 85*10000 # cost per hectare (Theis et al. 2022)

# ---- UNCERTAINTY DISTRBUTIONS ----
# generate value distributions

# Impact
sigma.i <- sqrt(log(cv.i^2+1))    # convert CV to sd
mu.i <- log(mean.i) - sigma.i^2/2 # covert mean to log-mean
d.i <- rlnorm(n.sim, meanlog = mu.i,  sdlog = sigma.i) # generate random draws

# offset 1 
# - PERT error
d.o1 = mc2d::rpert(n.sim,         # number of random values
                  mode = mode.o1, # mode
                  min = min.o1,   # min
                  max = max.o1,   # max
                  shape = 4)      # shape parameter - 4 default

# offset 2
sigma.o2 <- sqrt(log(cv.o2^2+1))     # convert CV to sd
mu.o2 <- log(mean.o2) - sigma.o2^2/2 # covert mean to log-mean
d.o2 <- rlnorm(n.sim, meanlog = mu.o2,  sdlog = sigma.o2)# generate random draws

# ---- MULTIPLIERS & COSTS ----
# calculated offset multipliers under different offsets weights and calculate
# total cost of offset, then optimize to identify weighting with minim cost

## - Equal Weight - 
CR = CRu.f(
  d.i = matrix(d.i),       # impact dist
  d.o = cbind(d.o1, d.o2), # offset dist
  p = p,                   # risk threshold
  CR_weight = c(1,1)       # Weights - equal
)

# multiplier values
CR$`CR - Adjusted` 
# 6.89 each

# Costs
sum(CR$`CR - Adjusted` * c(cost.o1, cost.o2))
# $7.03M

## - Only offset 1 used -
CR = CRu.f(
  d.i = matrix(d.i), # impact dist
  d.o = cbind(d.o1), # offset dist
  p = p,             # risk threshold
  CR_weight = c(1)   # Weights - must be 1 with i offset project
)

# multiplier values
CR$`CR - Adjusted`
# 55.71

# Costs
sum(CR$`CR - Adjusted` * c(cost.o1))
# $9.47M

## - only offset 2 -
CR = CRu.f(
  d.i = matrix(d.i), # impact dist
  d.o = cbind(d.o2), # offset dist
  p = p,             # risk threshold
  CR_weight = c(1)   # Weights - must be 1 with i offset project
)

# multiplier values
CR$`CR - Adjusted`
# 10.31

# Costs
sum(CR$`CR - Adjusted` * c(cost.o2))
# $8.8M

## - Optimize for min. costs -
# optimization function 
# find the wj that minimized total cost 
# optimize over only one w value while keeping the other = 1 (one projects need
# w = 1 always)
optim.f <- function(w,         # par for optimization
                    offset = 1 # identify which offset 
                    ) {
  
  # set w vector depending on offset project
  if(offset == 1) {
    ws <- c(w, 1) # offset 1
  } else if (offset == 2) {
    ws <- c(1, w) # offset 2
  }
  
  # calculate multiplier
  CR = CRu.f(
    d.i = matrix(d.i),      # impact dist. 
    d.o = cbind(d.o1, d.o2),# offset dist.
    p = p,                  # risk level
    CR_weight = ws          # offset weights
  )

  # calculate costs
  cost <- sum(CR$`CR - Adjusted` * c(cost.o1, cost.o2))

  # return costs
  return(cost)
  
}

# Offset 1 weights optimized
min.01 <- optim(   # optimizer
  par = c(1),      # initial value
  fn = optim.f,    # function
  method = "Brent",# method
  lower = c(0),    # lower limit = 0
  upper = c(1),    # upper limit = 1
  offset = 1       # offset project 1 
)
# w1 = 0.98

# calc what the multiplier is for optimized weights 
CR = CRu.f(
  d.i = matrix(d.i),           # impact dist.
  d.o = cbind(d.o1, d.o2),     # offset dist
  p = p,                       # risk levevl
  CR_weight = c(min.01$par, 1) # offset weight - w1 optimized value 
)

# multiplier values
CR$`CR - Adjusted`
# 6.84 and 6.98

# costs
sum(CR$`CR - Adjusted` * c(cost.o1, cost.o2))
# $7.02M

# offset 2 weights optimized
min.02 <- optim(   # optimizer
  par = c(1),      # initial value
  fn = optim.f,    # function
  method = "Brent",# M3thod
  lower = c(0),    # lower limit = 0
  upper = c(1),    # upper limit = 1
  offset = 2       # offset project 2
)
# w2 = 0.27
CR = CRu.f(
  d.i = matrix(d.i),
  d.o = cbind(d.o1, d.o2),
  p = p,
  CR_weight = c(1, min.02$par) # how CR will be applied across offset types
)
# w2 = 0.28

# multiplier values
CR$`CR - Adjusted`
# 14.4, 4.86

# costs
sum(CR$`CR - Adjusted` * c(cost.o1, cost.o2))
# 6.58M

################################################################################