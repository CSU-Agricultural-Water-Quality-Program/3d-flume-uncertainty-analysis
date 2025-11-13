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
df_raw <- read.csv("../data/flume_data.csv")

# Drop flume calibration rows (i.e., where flume_type == EQ)
df <- df_raw %>% filter(flume_type != "EQ")

# create data list for Stan model
dat <- list(
  N       = nrow(df),
  Q_obs   = as.numeric(df$flow_cfs),                 # observed discharge (cfs)
  H_obs   = as.numeric(df$depth_ft),                 # observed head (ft)
  F       = as.numeric(as.factor(df$flume_type)),    # flume index 1..J
  J       = length(unique(df$flume_type)),           # number of flume groups
  sigma_H = 0.005                                    # ft, from point gauge accuracy
)


# Define the Bayesian models

# Model 1.0 
# - basic model, not looking at flume types separately
# - True flow Q is latent variable
flume_model_1 <- ulam(
  alist(
    # Observation model for discharge
    Q_obs ~ dnorm(Q_true, sigma_true),
    
    # Process model for true discharge (power-law)
    vector[N]:Q_true ~ dnorm(a * H_obs^b, sigma_proc),
    
    # Priors
    a ~ dlnorm(log(1.55) , 0.5),
    b ~ dlnorm(log(2.58) , 0.5),
    sigma_true ~ dexp(1),
    sigma_proc ~ dexp(1)
  ),
  data = dat, chains = 4, cores = 4, iter = 2000, log_lik = TRUE
)

# Summarize the model
precis(flume_model_1, depth = 2)
plot(precis(flume_model_1, depth = 2))


# Model 2.0
# - hierarchical model, looking at flume types separately
# - True flow Q is latent variable
flume_model_2 <- ulam(
  alist(
    # Observation model for discharge (sigma varies by flume via F)
    Q_obs ~ dnorm( Q_true , sigma_true[F] ),
    
    # Process model for latent true discharge (sigma varies by flume)
    vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
    
    # Vectorized mean, per row using flume-indexed a and b
    mu <- a[F] * H_obs^b[F],
    
    # Per-flume priors
    vector[J]:a ~ dlnorm( log(1.55) , 0.5 ),
    vector[J]:b ~ dlnorm( log(2.58) , 0.5 ),
    
    # Per-flume scatter
    vector[J]:sigma_true ~ dexp(1),
    vector[J]:sigma_proc ~ dexp(1)
  ),
  data = dat,
  chains = 4, cores = 4, iter = 2000, log_lik = TRUE
)


# Summarize the model
precis(flume_model_2, depth = 2)
plot(precis(flume_model_2, depth = 2))


# Model 3.0
# - hierarchical model, looking at flume types separately
# - True flow Q is observed with error
# - includes measurement error explicitly

# assumes: dat$sigma_H exists (scalar), and H_obs > 0

flume_model_3 <- ulam(
  alist(
    # Latent true discharge with per-flume process scatter (vectorized)
    vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
    
    # Vectorized mean by flume
    mu <- a[F] * H_obs^b[F],
    
    # Observation model: heteroscedastic SD
    # global percent-of-reading + per-flume additive infiltration
    Q_obs ~ dnorm(
      Q_true,
      sqrt(
        square( cv_Q * Q_true ) +          # global % of reading
          square( sigma_Q_add[F] )           # per-flume infiltration
      )
    ),
    
    # Per-flume rating-curve parameters
    vector[J]:a ~ dlnorm( log(1.55) , 0.5 ),
    vector[J]:b ~ dlnorm( log(2.58) , 0.5 ),
    
    # Per-flume process scatter
    vector[J]:sigma_proc ~ dexp(1),
    
    # Error hyperparameters
    cv_Q ~ dlnorm( log(2e-4) , 0.5 ),            # single scalar
    vector[J]:sigma_Q_add ~ dlnorm( log(0.02) , 0.5 )  # per flume
  ),
  data = dat, chains = 4, cores = 4, iter = 4000
)

# model convergence check
traceplot(flume_model_3)
trankplot(flume_model_3)

# Summarize the model and save output
p3 <- precis(flume_model_3, depth = 2)
write.csv(p3, "../output/flume_model_3_precis.csv", row.names = TRUE)
p3

# plot model summary and save as jpg
jpeg("../output/flume_model_3_precis.jpg", width = 6.5, height = 6.5, units = "in", res = 600)
plot(p3)
dev.off()
plot(p3)


post3 <- extract.samples(flume_model_3)

# Graphing Outputs

# Stage discharge curves by flume type with uncertainty intervals
#' Plot Bayesian stage–discharge curves with posterior ribbons
#'
#' @param post  list returned by rethinking::extract.samples(ulam_fit)
#'              expects components: a[J], b[J], sigma_proc[J], cv_Q, sigma_Q_add[J]
#' @param dat   list you passed to ulam, needs: Q_obs[N], H_obs[N], F[N], J
#' @param prob  probability mass for ribbons (default 0.89)
#' @param ndraws number of posterior draws to use for ribbons (default 800)
#' @param show_manufacturer optional data.frame with columns flume (factor or int 1..J),
#'              a_m (manufacturer a), b_m (manufacturer b) to add dotted curves
#' @return      ggplot object
# --- Dependencies ---
# install.packages(c("ggplot2","dplyr","tidyr"))

plot_flume_posterior <- function(post, dat, prob = 0.95, ndraws = 2000,
                                 show_manufacturer = NULL,
                                 flume_names = c("1" = "3D Printed 1",
                                                 "2" = "3D Printed 2",
                                                 "3" = "Fiberglass")) {
  stopifnot(all(c("Q_obs","H_obs","F","J") %in% names(dat)))
  Q    <- as.numeric(dat$Q_obs)
  H    <- pmax(as.numeric(dat$H_obs), 1e-8)
  Fidx <- as.integer(dat$F)
  J    <- as.integer(dat$J)
  
  # use a single global grid so all lines span the full x-range
  Hgrid <- seq(min(H), max(H), length.out = 200)
  
  # choose posterior draws
  Sfull <- nrow(post$a)
  S     <- min(ndraws, Sfull)
  sel   <- if (S < Sfull) sample.int(Sfull, S) else seq_len(S)
  
  a_draw  <- post$a[sel, , drop=FALSE]
  b_draw  <- post$b[sel, , drop=FALSE]
  sp_draw <- post$sigma_proc[sel, , drop=FALSE]
  cv_draw <- as.numeric(post$cv_Q[sel, 1])
  sq_draw <- post$sigma_Q_add[sel, , drop=FALSE]
  
  flumes_present <- sort(unique(Fidx))
  
  # ribbons on common grid
  rib_list <- lapply(flumes_present, function(j){
    logH   <- log(Hgrid)
    log_mu <- outer(b_draw[, j], logH, `*`) +
      matrix(log(a_draw[, j]), nrow = S, ncol = length(Hgrid))
    mu_mat <- exp(log_mu)
    sd_mat <- sqrt( (cv_draw^2) * (mu_mat^2) + (sq_draw[, j]^2) + (sp_draw[, j]^2) )
    yrep   <- mu_mat + sd_mat * matrix(rnorm(length(mu_mat)), nrow = S)
    data.frame(
      flume = j,
      H     = Hgrid,
      Q_lo  = apply(yrep, 2, quantile, probs = (1 - prob)/2),
      Q_hi  = apply(yrep, 2, quantile, probs = 1 - (1 - prob)/2),
      Q_med = apply(yrep, 2, median)
    )
  })
  ribbon_df <- do.call(rbind, rib_list)
  
  # observations with median error model
  cv_med   <- median(post$cv_Q)
  sq_med_j <- apply(post$sigma_Q_add, 2, median)
  sp_med_j <- apply(post$sigma_proc,   2, median)
  obs_df <- data.frame(H = H, Q = Q, flume = Fidx)
  obs_df$se_Q <- sqrt( (cv_med * obs_df$Q)^2 +
                         (sq_med_j[obs_df$flume])^2 +
                         (sp_med_j[obs_df$flume])^2 )
  
  # manufacturer curve on global grid
  manuf_df <- NULL
  if (!is.null(show_manufacturer)) {
    if (!all(c("a_m","b_m") %in% names(show_manufacturer)))
      stop("show_manufacturer must have columns a_m and b_m (and optional flume).")
    if ("flume" %in% names(show_manufacturer)) {
      tmp <- subset(show_manufacturer, flume %in% flumes_present)
      mlist <- lapply(tmp$flume, function(j){
        a_m <- tmp$a_m[tmp$flume == j][1]
        b_m <- tmp$b_m[tmp$flume == j][1]
        data.frame(flume = j, H = Hgrid, Qm = a_m * Hgrid^b_m)
      })
      manuf_df <- do.call(rbind, mlist)
    } else {
      # one global manufacturer curve (no flume split)
      manuf_df <- data.frame(H = Hgrid,
                             Qm = show_manufacturer$a_m[1] * Hgrid^(show_manufacturer$b_m[1]))
    }
  }
  
  # match Python colors: 1=blue, 2=green, 3=orange (edit if your indexing differs)
  color_map <- c("1" = "#1f77b4",  # 3D Printed 1
                 "2" = "#ff7f0e",  # 3D Printed 2
                 "3" = "#2ca02c")  # Fiberglass
  # keep only present flumes
  color_map <- color_map[intersect(names(color_map), as.character(flumes_present))]
  
  # shapes by flume, stable across figures
  shape_map <- c("1" = 16, "2" = 18, "3" = 15)
  
  # plotting
  ggplot() +
    geom_ribbon(data = ribbon_df,
                aes(H, ymin = Q_lo, ymax = Q_hi, fill = factor(flume)),
                alpha = 0.25, color = NA) +
    geom_line(data = ribbon_df,
              aes(H, Q_med, color = factor(flume)), linewidth = 1.2) +
    geom_errorbar(data = obs_df,
                  aes(H, ymin = pmax(Q - se_Q, 0), ymax = Q + se_Q,
                      color = factor(flume)),
                  width = 0, alpha = 0.6) +
    geom_point(data = obs_df,
               aes(H, Q, color = factor(flume), shape = factor(flume)),
               size = 2.2, stroke = 0.6) +
    { if (!is.null(manuf_df) && "flume" %in% names(manuf_df))
      geom_line(data = manuf_df, aes(H, Qm, group = flume),
                linewidth = 1.0, linetype = 3, color = "black")
      else if (!is.null(manuf_df))
        geom_line(data = manuf_df, aes(H, Qm),
                  linewidth = 1.0, linetype = 3, color = "black") } +
    scale_color_manual(values = color_map, name = "Flume", labels = flume_names[names(color_map)]) +
    scale_fill_manual(values = scales::alpha(color_map, 0.25), guide = "none") +
    scale_shape_manual(values = shape_map[names(color_map)], name = "Flume", labels = flume_names[names(color_map)]) +
    labs(x = "Stage (ft)", y = "Flow rate (cfs)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())
}


# Assuming you already have:
#   post3 <- extract.samples(flume_model_3)
#   dat   <- list(N=nrow(df), Q_obs=df$flow_cfs, H_obs=df$depth_ft,
#                 F=as.integer(as.factor(df$flume_type)), J=length(unique(df$flume_type)))

# Optional manufacturer coefficients by flume index (edit as needed)
# Equation: Q = a_m * H^b_m
manuf <- data.frame(
  flume = c(1),
  a_m   = c(1.55),
  b_m   = c(2.58)
)

g <- plot_flume_posterior(post3, dat, prob = 0.95, ndraws = 4000,
                          show_manufacturer = manuf)
print(g)
# save the plot
ggsave("../output/flume_posterior_curves.jpg", g, width = 8, height = 5.5, dpi = 300)








