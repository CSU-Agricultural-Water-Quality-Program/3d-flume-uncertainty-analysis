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
#df_raw <- read.csv("../data/flume_data.csv")
df_raw <- read.csv("../data/sim_data.csv") # for testing with simulated data

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
# flume_model_1 <- ulam(
#     alist(
#       # Observation model for discharge
#       Q_obs ~ dnorm(Q_true, sigma_true),
#       
#       # Process model for true discharge (power-law)
#       vector[N]:Q_true ~ dnorm(a * H_obs^b, sigma_proc),
#       
#       # Priors
#       a ~ dlnorm(log(1.55) , 0.5),
#       b ~ dlnorm(log(2.58) , 0.5),
#       sigma_true ~ dexp(1),
#       sigma_proc ~ dexp(1)
#     ),
#     data = dat, chains = 4, cores = 4, iter = 2000, log_lik = TRUE
#   )

# Summarize the model
# precis(flume_model_1, depth = 2)
# plot(precis(flume_model_1, depth = 2))


# Model 2.0
# - hierarchical model, looking at flume types separately
# - True flow Q is latent variable
# flume_model_2 <- ulam(
#     alist(
#       # Observation model for discharge (sigma varies by flume via F)
#       Q_obs ~ dnorm( Q_true , sigma_true[F] ),
#       
#       # Process model for latent true discharge (sigma varies by flume)
#       vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
#       
#       # Vectorized mean, per row using flume-indexed a and b
#       mu <- a[F] * H_obs^b[F],
#       
#       # Per-flume priors
#       vector[J]:a ~ dlnorm( log(1.55) , 0.5 ),
#       vector[J]:b ~ dlnorm( log(2.58) , 0.5 ),
#       
#       # Per-flume scatter
#       vector[J]:sigma_true ~ dexp(1),
#       vector[J]:sigma_proc ~ dexp(1)
#     ),
#     data = dat,
#     chains = 4, cores = 4, iter = 2000, log_lik = TRUE
#   )


# Summarize the model
# precis(flume_model_2, depth = 2)
# plot(precis(flume_model_2, depth = 2))


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
traceplot_ulam(flume_model_3, ask=F, window=c(2100,4000))
trankplot(flume_model_3, ask=F, window=c(2100,4000))

# Summarize the model and save output
p3 <- precis(flume_model_3, depth = 2)
write.csv(p3, "../output/flume_model_3_precis.csv", row.names = TRUE)
p3

# plot model summary and save as jpg
jpeg("../output/flume_model_3_precis.jpg", width = 6.5, height = 6.5, units = "in", res = 600)
plot(p3)
dev.off()
par(mfrow=c(1,1))
plot(p3)


post3 <- extract.samples(flume_model_3)


# Model 4: same as Model 3, but reparameterize a = exp(alpha)
# so the mean is computed via log space:
#   mu = exp(alpha[F] + b[F] * log(H))
#
# Likelihood and uncertainty structure remain on original Q scale:
#   Q_true ~ Normal(mu, sigma_proc[F])
#   Q_obs  ~ Normal(Q_true, sqrt((cv_Q*Q_true)^2 + sigma_Q_add[F]^2))

# Model 4.0
# - identical likelihood / uncertainty structure to Model 3
# - reparameterizes a[F] as exp(alpha[F]) so the power law mean is computed in log space
# - goal: improve sampling geometry for a/b without changing your manuscript’s uncertainty story

flume_model_4 <- ulam(
  alist(
    # Latent true discharge with per-flume process scatter (vectorized)
    vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
    
    # Mean computed in log space, then exponentiated back to original scale
    # (this is mathematically identical to mu <- a[F] * H_obs^b[F] when a[F] = exp(alpha[F]))
    mu <- exp( alpha[F] + b[F] * log(H_obs + 1e-9) ),
    
    # Observation model: SAME as Model 3 (heteroscedastic SD)
    Q_obs ~ dnorm(
      Q_true,
      sqrt(
        square( cv_Q * Q_true ) +
          square( sigma_Q_add[F] )
      )
    ),
    
    # Priors: alpha is the log of a, consistent with your original a prior
    # If a ~ LogNormal(log(1.55), 0.5), then alpha = log(a) ~ Normal(log(1.55), 0.5)
    vector[J]:alpha ~ dnorm( log(1.55) , 0.5 ),
    
    # Same b prior as Model 3
    vector[J]:b ~ dlnorm( log(2.58) , 0.5 ),
    
    # Same per-flume process scatter
    vector[J]:sigma_proc ~ dexp(1),
    
    # Same error hyperparameters
    cv_Q ~ dlnorm( log(2e-4) , 0.5 ),
    vector[J]:sigma_Q_add ~ dlnorm( log(0.02) , 0.5 )
  ),
  data = dat,
  chains = 4, cores = 4, iter = 4000
  #log_lik = TRUE
)


# model convergence check
traceplot_ulam(flume_model_4, ask=F, window=c(2100,4000))
trankplot(flume_model_4, ask=F, window=c(2100,4000))

# Summarize the model and save output
p4 <- precis(flume_model_4, depth = 2)
write.csv(p4, "../output/flume_model_4_precis.csv", row.names = TRUE)
p4

# plot model summary and save as jpg
jpeg("../output/flume_model_4_precis.jpg", width = 6.5, height = 6.5, units = "in", res = 600)
plot(p4)
dev.off()
par(mfrow=c(1,1))
plot(p4)


######################### Graphing Outputs #####################################

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

plot_flume_posterior <- function(post, dat, prob = 0.95, ndraws = 4000,
                                 show_manufacturer = NULL,
                                 ribbon_kind = c("process", "parameter", "predictive"),
                                 show_qobs_errors = FALSE,
                                 flume_names = c("1" = "3D Printed 1",
                                                 "2" = "3D Printed 2",
                                                 "3" = "Fiberglass")) {
  
  ribbon_kind <- match.arg(ribbon_kind)
  stopifnot(all(c("Q_obs","H_obs","F","J") %in% names(dat)))
  
  Q    <- as.numeric(dat$Q_obs)
  H    <- pmax(as.numeric(dat$H_obs), 1e-8)
  Fidx <- as.integer(dat$F)
  J    <- as.integer(dat$J)
  
  # global grid (smooth coverage across x-range)
  Hgrid <- seq(min(H), max(H), length.out = 200)
  logH  <- log(Hgrid)
  
  # subset posterior draws
  Sfull <- nrow(post$a)
  S     <- min(ndraws, Sfull)
  sel   <- if (S < Sfull) sample.int(Sfull, S) else seq_len(S)
  
  a_draw  <- post$a[sel, , drop=FALSE]
  b_draw  <- post$b[sel, , drop=FALSE]
  sp_draw <- post$sigma_proc[sel, , drop=FALSE]
  cv_draw <- as.numeric(post$cv_Q[sel])
  sq_draw <- post$sigma_Q_add[sel, , drop=FALSE]
  
  flumes_present <- sort(unique(Fidx))
  
  # central interval endpoints and z for prob
  p_lo  <- (1 - prob)/2
  p_hi  <- 1 - p_lo
  zcrit <- qnorm((1 + prob)/2)
  
  # ---------- Smooth ribbons for Q_true ----------
  # No per-point random noise; everything is built from smooth functions
  rib_list <- lapply(flumes_present, function(j){
    # posterior of mu(H) on the grid
    log_mu <- outer(b_draw[, j], logH, `*`) +
      matrix(log(a_draw[, j]), nrow=S, ncol=length(Hgrid))
    mu_mat <- exp(log_mu)
    
    mu_med <- apply(mu_mat, 2, median)
    mu_lo  <- apply(mu_mat, 2, quantile, p_lo)
    mu_hi  <- apply(mu_mat, 2, quantile, p_hi)
    
    if (ribbon_kind == "parameter") {
      Q_lo <- mu_lo
      Q_hi <- mu_hi
    } else if (ribbon_kind == "process") {
      # add smooth process width using the median sigma_proc
      sp_med <- median(sp_draw[, j])
      Q_lo <- mu_med - zcrit * sp_med
      Q_hi <- mu_med + zcrit * sp_med
    } else { # predictive: include measurement + process on top of mu_med
      cv_med  <- median(cv_draw)
      sq_med  <- median(sq_draw[, j])
      sp_med  <- median(sp_draw[, j])
      sd_pred <- sqrt((cv_med * mu_med)^2 + sq_med^2 + sp_med^2)
      Q_lo <- mu_med - zcrit * sd_pred
      Q_hi <- mu_med + zcrit * sd_pred
    }
    
    data.frame(
      flume = j,
      H     = Hgrid,
      Q_lo  = Q_lo,
      Q_hi  = Q_hi,
      Q_med = mu_med
    )
  })
  ribbon_df <- do.call(rbind, rib_list)
  
  # ---------- Q_true intervals at observed H ----------
  Qtrue_draws <- post$Q_true[sel, , drop=FALSE]
  Qtrue_df <- data.frame(
    H    = H,
    Q    = apply(Qtrue_draws, 2, median),
    Q_lo = apply(Qtrue_draws, 2, quantile, p_lo),
    Q_hi = apply(Qtrue_draws, 2, quantile, p_hi),
    flume = Fidx,
    ptype = "Q_true"
  )
  
  # ---------- Q_obs points (optional error bars controlled by prob) ----------
  cv_med   <- median(post$cv_Q)
  sq_med_j <- apply(post$sigma_Q_add, 2, median)
  sp_med_j <- apply(post$sigma_proc,   2, median)
  se_obs   <- sqrt( (cv_med * Q)^2 + (sq_med_j[Fidx])^2 + (sp_med_j[Fidx])^2 )
  Qobs_df  <- data.frame(
    H = H, Q = Q,
    Q_lo = pmax(Q - zcrit*se_obs, 0),
    Q_hi = Q + zcrit*se_obs,
    flume = Fidx,
    ptype = "Q_obs"
  )
  
  # ---------- Manufacturer curve ----------
  manuf_df <- NULL
  if (!is.null(show_manufacturer)) {
    if ("flume" %in% names(show_manufacturer)) {
      tmp <- subset(show_manufacturer, flume %in% flumes_present)
      manuf_df <- do.call(rbind, lapply(tmp$flume, function(j){
        a_m <- tmp$a_m[tmp$flume == j][1]; b_m <- tmp$b_m[tmp$flume == j][1]
        data.frame(flume = j, H = Hgrid, Qm = a_m * Hgrid^b_m)
      }))
    } else {
      manuf_df <- data.frame(H = Hgrid,
                             Qm = show_manufacturer$a_m[1] * Hgrid^(show_manufacturer$b_m[1]))
    }
  }
  
  # ---------- aesthetics ----------
  color_map <- c("1"="#1f77b4","2"="#ff7f0e","3"="#2ca02c")
  color_map <- color_map[intersect(names(color_map), as.character(flumes_present))]
  
  # We’ll keep a separate “Points” legend using shape only, with custom keys
  # Q_true: filled with black outline; Q_obs: white fill with colored outline (legend shows generic)
  shape_map <- c("Q_true" = 21, "Q_obs" = 21)
  
  library(ggplot2)
  ggplot() +
    # ribbons and median curve
    geom_ribbon(data = ribbon_df,
                aes(H, ymin = Q_lo, ymax = Q_hi, fill = factor(flume)),
                alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line(data = ribbon_df,
              aes(H, Q_med, color = factor(flume)), linewidth = 1.2) +
    
    # Q_true: filled with flume color, black outline
    geom_linerange(data = Qtrue_df,
                   aes(H, ymin = Q_lo, ymax = Q_hi, color = factor(flume)),
                   linewidth = 1.6, alpha = 0.6, show.legend = FALSE) +
    geom_point(data = Qtrue_df,
               aes(H, Q, fill = factor(flume), shape = ptype),
               size = 3.6, stroke = 0.9, color = "black") +
    
    # Q_obs: white fill, flume-colored outline
    { if (show_qobs_errors)
      geom_linerange(data = Qobs_df,
                     aes(H, ymin = Q_lo, ymax = Q_hi, color = factor(flume)),
                     linewidth = 0.9, alpha = 0.45, show.legend = FALSE) } +
    geom_point(data = Qobs_df,
               aes(H, Q, color = factor(flume), shape = ptype),
               fill = "white", size = 2.0, stroke = 1.2) +
    
    # manufacturer equation
    { if (!is.null(manuf_df))
      geom_line(data = manuf_df,
                aes(H, Qm, linetype = "Manufacturer equation"),
                linewidth = 1.0, color = "black") } +
    
    # scales and guides
    scale_color_manual(values = color_map,
                       name = "Flume",
                       labels = flume_names[names(color_map)]) +
    scale_fill_manual(values = color_map, guide = "none") +
    scale_shape_manual(values = shape_map,
                       name = "Points",
                       breaks = c("Q_true","Q_obs"),
                       labels = c("Estimated True Flow","Measured Flow")) +
    scale_linetype_manual(values = c("Manufacturer equation" = "dotted"), name = "") +
    guides(
      color   = guide_legend(order = 1),
      shape   = guide_legend(
        order = 2,
        override.aes = list(
          # two keys, in order of breaks: Q_true then Q_obs
          shape = c(21, 21),
          size  = c(3.6, 3.0),
          stroke= c(0.9, 1.2),
          color = c("black", "black"),  # outline in legend
          fill  = c("grey60", "white")  # filled vs open for legend clarity
        )
      ),
      linetype= guide_legend(order = 3)
    ) +
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
  a_m   = c(1.55),
  b_m   = c(2.58)
)

# Generate the plot
g <- plot_flume_posterior(
  post  = post3,              # posterior samples from ulam
  dat   = dat,                # data list used for the model
  prob  = 0.95,               # credible interval probability
  ndraws = 4000,              # number of posterior draws used
  show_manufacturer = manuf,  # add manufacturer equation line
  ribbon_kind = "parameter",  # options: "parameter", "process", or "predictive"
  show_qobs_errors = FALSE,   # measurement error bars for Q_obs off (recommended)
  flume_names = c("1" = "3D Printed 1",
                  "2" = "3D Printed 2",
                  "3" = "Fiberglass")
)

# Print to screen
print(g)

# Save a high-resolution version
ggsave("../output/flume_posterior_curves_final.jpg", g,
       width = 8, height = 5.5, dpi = 600)


# Generative predictions output
# Generate posterior summaries of Q(H) on a grid for each flume
# post  : extract.samples(ulam_fit)
# depths: numeric vector of H (ft) at which to report results
# prob  : central interval mass (e.g., 0.95)
# kind  : "parameter" (default), "process", or "predictive"
# manuf : data.frame with either global a_m,b_m columns, or per-flume (flume,a_m,b_m)
# flume_names: named character vector mapping "1","2","3",... -> labels to use
generate_flume_table <- function(post, depths = seq(0.1, 0.8, by = 0.1),
                                 prob = 0.95,
                                 kind = c("parameter", "process", "predictive"),
                                 ndraws = 4000,
                                 manuf = NULL,
                                 flume_names = c("1"="3D Printed 1","2"="3D Printed 2","3"="Fiberglass"),
                                 digits = 3) {
  
  kind <- match.arg(kind)
  stopifnot(is.matrix(post$a), is.matrix(post$b))
  
  # subset draws
  Sfull <- nrow(post$a)
  S     <- min(ndraws, Sfull)
  sel   <- if (S < Sfull) sample.int(Sfull, S) else seq_len(S)
  
  a  <- post$a[sel, , drop=FALSE]              # S x J
  b  <- post$b[sel, , drop=FALSE]              # S x J
  sp <- post$sigma_proc[sel, , drop=FALSE]     # S x J
  cv <- as.numeric(post$cv_Q[sel])             # S
  sq <- post$sigma_Q_add[sel, , drop=FALSE]    # S x J
  
  J <- ncol(a)
  H <- depths
  G <- length(H)
  p_lo <- (1 - prob)/2
  p_hi <- 1 - p_lo
  
  # compute per-flume summaries
  out_list <- lapply(seq_len(J), function(j){
    # S x G mu draws
    mu_mat <- exp( outer(b[,j], log(H), `*`) + matrix(log(a[,j]), nrow=S, ncol=G) )
    
    if (kind == "parameter") {
      lo <- apply(mu_mat, 2, quantile, p_lo)
      md <- apply(mu_mat, 2, median)
      hi <- apply(mu_mat, 2, quantile, p_hi)
    } else if (kind == "process") {
      # add process noise with draw-specific sigma_proc (smooth summary via simulation)
      yrep <- mu_mat + matrix(rnorm(S*G), nrow=S, ncol=G) * matrix(sp[,j], nrow=S, ncol=G)
      lo <- apply(yrep, 2, quantile, p_lo)
      md <- apply(mu_mat, 2, median)  # center is mu
      hi <- apply(yrep, 2, quantile, p_hi)
    } else { # predictive: add process + measurement
      sd_pred <- sqrt( (cv^2) %o% (rep(1,G)) * (mu_mat^2) + (sq[,j]^2) %o% rep(1,G) + (sp[,j]^2) %o% rep(1,G) )
      yrep <- mu_mat + matrix(rnorm(S*G), nrow=S, ncol=G) * sd_pred
      lo <- apply(yrep, 2, quantile, p_lo)
      md <- apply(mu_mat, 2, median)
      hi <- apply(yrep, 2, quantile, p_hi)
    }
    
    data.frame(Depth_ft = H,
               flume = names(flume_names)[j],
               median = md, lo = lo, hi = hi)
  })
  
  df_long <- do.call(rbind, out_list)
  
  # wide table with one set of columns per flume
  # value columns as "median [lo, hi]" strings for quick reporting
  library(dplyr); library(tidyr)
  df_fmt <- df_long %>%
    mutate(label = paste0(sprintf(paste0("%.",digits,"f"), median),
                          " [", sprintf(paste0("%.",digits,"f"), lo),
                          ", ", sprintf(paste0("%.",digits,"f"), hi), "]"),
           Flume = flume_names[flume]) %>%
    select(Depth_ft, Flume, label) %>%
    pivot_wider(names_from = Flume, values_from = label) %>%
    arrange(Depth_ft)
  
  # manufacturer column
  if (!is.null(manuf)) {
    if (all(c("a_m","b_m") %in% names(manuf)) && !("flume" %in% names(manuf))) {
      df_fmt$`Manufacturer Eq.` <- sprintf(paste0("%.",digits,"f"),
                                           manuf$a_m[1] * H^(manuf$b_m[1]))
    } else if (all(c("flume","a_m","b_m") %in% names(manuf))) {
      # if provided per flume, add a column per flume with "Manu." suffix
      for (j in seq_len(J)) {
        if (as.character(j) %in% manuf$flume) {
          a_m <- manuf$a_m[manuf$flume == j][1]
          b_m <- manuf$b_m[manuf$flume == j][1]
          colname <- paste0(flume_names[as.character(j)], " (Manu.)")
          df_fmt[[colname]] <- sprintf(paste0("%.",digits,"f"), a_m * H^(b_m))
        }
      }
    }
  }
  
  # also return a numeric version with separate median/lo/hi columns per flume (for CSV/analysis)
  df_num <- df_long %>%
    mutate(Flume = flume_names[flume]) %>%
    select(Depth_ft, Flume, median, lo, hi) %>%
    pivot_wider(names_from = Flume,
                values_from = c(median, lo, hi),
                names_sep = "_") %>%
    arrange(Depth_ft)
  
  list(display = df_fmt, numeric = df_num)
}

# grid and manufacturer equation
depth_grid <- seq(0.1, 0.8, by = 0.1)
manuf <- data.frame(a_m = 1.55, b_m = 2.58)  # global manufacturer curve

res <- generate_flume_table(
  post = post3,
  depths = depth_grid,
  prob = 0.95,
  kind = "parameter",     # "parameter" for comparing flume calibration curves
  ndraws = 4000,
  manuf = manuf,
  flume_names = c("1"="3D Printed 1","2"="3D Printed 2","3"="Fiberglass"),
  digits = 3
)

# Quick view for reporting (one column per flume with "median [lo, hi]")
res$display

# Numeric table for export or further analysis
res$numeric

# Save both to CSVs if you like
# write.csv(res$display, "../output/flume_generated_display.csv", row.names = FALSE)
write.csv(res$numeric, "../output/flume_generated_numeric.csv", row.names = FALSE)







