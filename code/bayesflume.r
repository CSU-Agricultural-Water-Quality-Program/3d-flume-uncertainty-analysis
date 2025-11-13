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
  data = dat, chains = 4, cores = 4, iter = 2000, log_lik = TRUE
)

# Summarize the model
precis(flume_model_3, depth = 1)
plot(precis(flume_model_3, depth = 1))
post3 <- extract.samples(flume_model_3)



# works but diverges too much - need to non center it.
flume_model_3nc <- ulam(
  alist(
    # Per-flume rating-curve parameters (directly positive)
    vector[J]:a ~ dlnorm( log(1.55) , 0.3 ),
    vector[J]:b ~ dlnorm( log(2.58) , 0.3 ),
    
    # Per-flume process scatter
    vector[J]:sigma_proc ~ dexp(2),
    
    # Observation-error hyperparameters
    cv_Q ~ dlnorm( log(2e-4) , 0.3 ),                 # global % of reading (scalar)
    vector[J]:sigma_Q_add ~ dlnorm( log(0.02) , 0.3 ),# per-flume infiltration
    
    # Vectorized mean
    mu <- a[F] * H_obs^b[F],
    
    # Latent true discharge (vectorized)
    vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
    
    # Observation model (element-wise => builds log_lik)
    Q_obs ~ dnorm(
      Q_true,
      sqrt( (cv_Q * Q_true)^2 + (sigma_Q_add[F])^2 )
    )
  ),
  data = dat,
  chains = 4, cores = 4, iter = 2000, log_lik = TRUE,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

# Summarize the model
precis(flume_model_3nc, depth = 2)


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
plot_flume_posterior <- function(post, dat, prob = 0.89, ndraws = 800,
                                 show_manufacturer = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
  if (!requireNamespace("dplyr", quietly = TRUE))   stop("Install dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE))   stop("Install tidyr")
  if (!requireNamespace("purrr", quietly = TRUE))   stop("Install purrr")
  
  library(ggplot2); library(dplyr); library(tidyr); library(purrr)
  
  N <- length(dat$Q_obs); J <- dat$J; Fidx <- dat$F
  H <- as.numeric(dat$H_obs); Q <- as.numeric(dat$Q_obs)
  
  # Posterior draws: thin for speed if needed
  Sfull <- nrow(post$a)
  S <- min(ndraws, Sfull)
  sel <- if (S < Sfull) sample.int(Sfull, S) else seq_len(S)
  
  a_draw  <- post$a[sel, , drop = FALSE]                # S x J
  b_draw  <- post$b[sel, , drop = FALSE]                # S x J
  sp_draw <- post$sigma_proc[sel, , drop = FALSE]       # S x J
  cv_draw <- as.numeric(post$cv_Q[sel, 1])              # S
  sq_draw <- post$sigma_Q_add[sel, , drop = FALSE]      # S x J
  
  # Per-flume grids over observed head range
  grids <- lapply(1:J, function(j){
    Hj <- H[Fidx == j]
    hseq <- seq(max(1e-6, min(Hj)), max(Hj), length.out = 100)
    data.frame(flume = j, H = hseq)
  })
  grid <- dplyr::bind_rows(grids)
  
  # Build ribbons: for each flume j and grid H, compute predictive draws
  # Draw-level predictive SD uses: sqrt( (cv * Q_true)^2 + sigma_Q_add_j^2 + sigma_proc_j^2 )
  # We simulate one predictive per draw (good mixture approx for ribbons)
  ribbon_df <- map_dfr(1:J, function(j){
    Hgrid <- grid$H[grid$flume == j]
    # matrices S x G
    mu_mat <- outer(1:S, Hgrid, function(s, k) a_draw[s, j] * Hgrid[k]^b_draw[s, j])
    sd_mat <- outer(1:S, Hgrid, function(s, k) {
      sqrt( (cv_draw[s] * mu_mat[s, k])^2 + sq_draw[s, j]^2 + sp_draw[s, j]^2 )
    })
    # predictive draws: one per draw (vectorized)
    yrep <- mu_mat + sd_mat * matrix(rnorm(length(mu_mat)), nrow = S)
    tibble(
      flume = j,
      H = rep(Hgrid, each = 1),
      # summarize across draws
      Q_lo = apply(yrep, 2, quantile, probs = (1 - prob)/2),
      Q_hi = apply(yrep, 2, quantile, probs = 1 - (1 - prob)/2),
      Q_med = apply(yrep, 2, median),
      Q_mean = apply(yrep, 2, mean)
    )
  })
  
  # Observations with posterior-based error bars (median params)
  cv_med   <- median(post$cv_Q)
  sq_med_j <- apply(post$sigma_Q_add, 2, median)        # length J
  sp_med_j <- apply(post$sigma_proc,   2, median)
  obs_df <- tibble(H = H, Q = Q, flume = Fidx) %>%
    mutate(se_Q = sqrt( (cv_med * Q)^2 + sq_med_j[flume]^2 + sp_med_j[flume]^2 ))
  
  # Optional manufacturer curves
  manuf_df <- NULL
  if (!is.null(show_manufacturer)) {
    tmp <- show_manufacturer %>%
      mutate(flume = as.integer(flume)) %>%
      inner_join(grid %>% distinct(flume), by = "flume")
    manuf_df <- tmp %>%
      group_by(flume) %>%
      group_modify(~{
        hseq <- grid$H[grid$flume == .y$flume[1]]
        tibble(flume = .y$flume[1], H = hseq,
               Qm = .x$a_m[1] * hseq^.x$b_m[1])
      }) %>% ungroup()
  }
  
  # Aesthetics
  palette <- scales::hue_pal()(J)
  flabs <- factor(1:J, labels = paste0("Flume ", 1:J))
  
  p <- ggplot() +
    # ribbons first
    geom_ribbon(data = ribbon_df,
                aes(x = H, ymin = Q_lo, ymax = Q_hi, fill = factor(flume)),
                alpha = 0.25, color = NA) +
    # posterior medians
    geom_line(data = ribbon_df,
              aes(x = H, y = Q_med, color = factor(flume)),
              linewidth = 1.2) +
    # observations with error bars
    geom_errorbar(data = obs_df,
                  aes(x = H, ymin = pmax(Q - se_Q, 0), ymax = Q + se_Q,
                      color = factor(flume)),
                  width = 0, alpha = 0.6) +
    geom_point(data = obs_df,
               aes(x = H, y = Q, color = factor(flume), shape = factor(flume)),
               size = 2.2, stroke = 0.6) +
    # optional manufacturer dotted curve
    { if (!is.null(manuf_df))
      geom_line(data = manuf_df,
                aes(x = H, y = Qm, group = flume),
                linetype = 3, linewidth = 1.0, color = "black") } +
    scale_color_manual(values = palette, name = "Flume") +
    scale_fill_manual(values = palette, guide = "none") +
    scale_shape_discrete(name = "Flume") +
    labs(x = "Stage (ft)", y = "Flow rate (cfs)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())
  
  return(p)
}

# Assuming you already have:
#   post3 <- extract.samples(flume_model_3)
#   dat   <- list(N=nrow(df), Q_obs=df$flow_cfs, H_obs=df$depth_ft,
#                 F=as.integer(as.factor(df$flume_type)), J=length(unique(df$flume_type)))

# Optional manufacturer coefficients by flume index (edit as needed)
manuf <- data.frame(
  flume = c(1,2,3),
  a_m   = c(3.01, 2.05, 1.96),
  b_m   = c(3.76, 3.01, 2.91)
)

g <- plot_flume_posterior(post3, dat, prob = 0.95, ndraws = 800,
                          show_manufacturer = manuf)
print(g)
# ggsave("flume_posterior_curves.png", g, width = 8, height = 5.5, dpi = 300)







