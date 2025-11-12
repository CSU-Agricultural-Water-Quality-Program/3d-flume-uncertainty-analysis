# bayesflume.r
# Bayesian Calibration of 3D-Printed Trapezoidal Flume
# AJ Brown
# Agricultural Water Quality Program
# Colorado State University

# Load necessary libraries
library(rethinking)
library(ggplot2)
library(dplyr)
library(tidyr)

# set working directory
setwd("C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/data")

# Load the flume calibration data
df <- read.csv("../data/flume_data.csv")

# create data list for stan model
dat <- list(
  N = nrow(df),
  Q_obs = df$flow_cfs,
  H_obs = df$depth_ft,
  F = as.integer(as.factor(df$flume_type)),
  J = length(unique(df$flume_type))
)

# Define the Bayesian model
flume_model <- ulam(
    alist(
        # model for D* observed flow
        Q_obs ~ dnorm(Q_true, sigma_true[F]),
        vector[J]:sigma_true ~ dexp(1),
        
        # model for D (true flow; unobserved)
        vector[N]:Q_true ~ dnorm(mu, sigma_obs[F]),
        mu <- a[F]*pow(H_true[i], b[F]),
        
        vector[J]:a ~ dlnorm(log(1.55), 0.5),
        vector[J]:b ~ dlnorm(log(2.58), 0.5),
        vector[J]:sigma_obs ~ dexp(1),
        
        # model for H with measurement error
        vector[N]:H_true ~ dnorm(H_obs, tau_H),
        tau_H ~ dexp(1)
    ),
    data = dat,
    chains = 4,
    cores = 4,
    iter = 2000
)

