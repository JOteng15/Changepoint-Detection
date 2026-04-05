#Declaring libraries and important data
library(tidyverse)
library(dplyr)
library(nbastatR)
library(changepoint)
library(ggplot2)
library(hoopR)
library(FOCuS)
library(tibble)
library(tictoc)
library(data.table)
source("load_data.R")


#--------------------Gaussian model--------------------------------------------------------------------------------------
detect_change_cusum <- function(y, sigma2, c) {
  n <- length(y)
  C_max <- 0
  t_hat <- 0
  S_n <- sum(y)
  S <- 0
  
  C_values <- numeric(n - 1)  # store all C_t^2 values
  
  for (t in 1:(n - 1)) {
    S <- S + y[t]
    
    # Means for two segments
    y_bar_1_t <- S / t
    y_bar_t1_n <- (S_n - S) / (n - t)
    
    # Compute C_t^2
    C_t2 <- (as.numeric(t) * (as.numeric(n) - as.numeric(t)) / as.numeric(n)) * (y_bar_1_t - y_bar_t1_n)^2
    
    # Save C_t^2 value
    C_values[t] <- C_t2
    
    # Update maximum
    if (C_t2 > C_max) {
      C_max <- C_t2
      t_hat <- t
    }
  }
  
  # Decision rule
  if (C_max / sigma2 > c) {
    result <- list(
      change_point = t_hat,
      C_max = C_max,
      C_values = C_values,
      message = "Changepoint detected",
      c=c
    )
  } else {
    result <- list(
      change_point = NULL,
      C_max = C_max,
      C_values = C_values,
      message = "No changepoint detected",
      c = c
    )
  }
  
  return(result)
}
#-----------------------------------Poisson Model----------------------------------------------------------------------------------------------------------
detect_changepoint_poisson <- function(y, c) {
  n <- length(y)
  C_max <- 0
  t_hat <- 0
  S_n <- sum(y)
  S <- 0
  C_values <- numeric(n-1)
  y_bar <- S_n / n
  
  stat_output <- NULL
  
  for (t in 1:(n - 1)) {
    S <- S + y[t]
    y_bar_1_t <- S / t
    y_bar_t1_n <- (S_n - S) / (n - t)
    
    # Avoid log(0)
    if (y_bar_1_t == 0) y_bar_1_t <- 1e-10
    if (y_bar_t1_n == 0) y_bar_t1_n <- 1e-10
    
    C_t <- 2 * (t * y_bar_1_t * log(y_bar_1_t / y_bar) +
                  (n - t) * y_bar_t1_n * log(y_bar_t1_n / y_bar))
    C_values[t] <- C_t
    if (C_t > C_max) {
      C_max <- C_t
      t_hat <- t
    }
  }
  
  if (C_max > c) {
    return(list(changepoint = t_hat,C_values = C_values, C_max = C_max, detected = TRUE))
  } else {
    return(list(changepoint = NULL, C_values = C_values, C_max = C_max, detected = FALSE))
  }
}
#--------------------------------------------change in variance--------------------------------------------------------------------------------------------------------------
detect_changepoint_variance <- function(y, c){
  n <- length(y)
  y_bar <- mean(y)
  s2_total <- mean((y - y_bar)^2)
  
  C_max <- 0
  t_hat <- 0
  C_values <- numeric(n-1)
  
  for (t in 2:(n - 2)) {  
    s2_1 <- mean((y[1:t] - mean(y[1:t]))^2)
    s2_2 <- mean((y[(t + 1):n] - mean(y[(t + 1):n]))^2)
    
    # Avoid log(0)
    if (s2_1 <= 0) s2_1 <- 1e-10
    if (s2_2 <= 0) s2_2 <- 1e-10
    
    C_t <- n * log(s2_total) - as.numeric(t) * log(s2_1) - (n - as.numeric(t)) * log(s2_2)
    C_values[t] <- C_t
    
    if (C_t > C_max) {
      C_max <- C_t
      t_hat <- t
    }
  }
  
  if (C_max > c) {
    return(list(changepoint = t_hat, C_values = C_values, Cmax = C_max, detected = TRUE))
  } else {
    return(list(changepoint = NULL, Cmax = C_max, C_values = C_values, detected = FALSE))
  }
}




graphplot <- function(data1,data){
  par(mfrow=c(2,1))
  
  plot.ts(data1,
          main = paste("Plus Minus Over Time"),
          ylab = "Plus Minus Score",
          xlab = "t",
          lwd = 2)
  
  text(350,5000,labels = paste("CP: ",as.character(data$change_point)))
  abline(v = data$change_point, col = "red", lty = 2)

  plot.ts(data$C_values,
          main = paste("CUSUM Statistics for","Plus Minus Score"),
          ylab = expression(C[t]^2),
          xlab = "t",
          lwd = 2)
  

  text(350,5000,labels = paste("CP: ",as.character(data$change_point)))
  abline(v = data$change_point, col = "red", lty = 2)
}
# thres <- tune_cusum_penalty()
# ans <- detect_change_cusum(cavs$plusminusTeam,1,)
# 
# graphplot(result)
#---------------------------------------------Mean-variance changepoint detection-------------------------------------------------------------------------------------------------------------------------
detect_change_mean_variance <- function(y, c) {
  n <- length(y)
  
  # Precompute cumulative sums and cumulative squared sums
  cs <- cumsum(y)
  cs2 <- cumsum(y^2)
  
  # Function to compute segment mean and variance MLEs
  segment_stats <- function(S, S2, len) {
    mean_hat <- S / len
    var_hat <- (S2 - len * mean_hat^2) / len   # MLE variance (biased)
    var_hat <- max(var_hat, 1e-12)             # avoid zero division
    return(list(mean = mean_hat, var = var_hat))
  }
  
  # Storage
  LLR_values <- numeric(n - 1)

  LLR_max <- -Inf
  t_hat <- NULL
  
  for (t in 1:(n - 1)) {
    # Segment 1
    S1  <- cs[t]
    S21 <- cs2[t]
    seg1 <- segment_stats(S1, S21, t)
    
    # Segment 2
    S2  <- cs[n] - cs[t]
    S22 <- cs2[n] - cs2[t]
    seg2 <- segment_stats(S2, S22, n - t)
    
    # Full segment
    whole <- segment_stats(cs[n], cs2[n], n)
    
    # Log-likelihoods (Gaussian)
    LL_full <- -n/2 * log(whole$var) - sum((y - whole$mean)^2) / (2 * whole$var)
    
    LL_split <-
      -t/2 * log(seg1$var) - sum((y[1:t] - seg1$mean)^2) / (2 * seg1$var) +
      -(n - t)/2 * log(seg2$var) - sum((y[(t+1):n] - seg2$mean)^2) / (2 * seg2$var)
    
    # Likelihood ratio statistic
    LLR <- 2 * (LL_split - LL_full)
    LLR_values[t] <- LLR
    
    if (LLR > LLR_max) {
      LLR_max <- LLR
      t_hat <- t
    }
  }
  
  # Decision rule: LLR > threshold 'c'
  if (LLR_max > c) {
    result <- list(
      change_point = t_hat,
      LLR = LLR_max,
      LLR_values = LLR_values,
      message = "Change point detected",
      c = c
    )
  } else {
    result <- list(
      change_point = NULL,
      LLR = LLR_max,
      LLR_values = LLR_values,
      message = "No change point detected",
      c = c
    )
  }
  
  return(result)
}



#------------------------------------------Monte Carlo simulation to tune an empirical penalty ----------------------------------------------------------------------------------------------------------
tune_cusum_penalty <- function(n,sigma2 = 1, alpha = 0.05, replicates = 1000) {
  ratios <- numeric(replicates)
  
  for (r in 1:replicates) {
    y <- rnorm(n, mean = 0,sigma2 )
    res <- detect_change_cusum(y,1, c = Inf)
    ratios[r] <- res$C_max/sigma2
  }
  
  c_empirical <- quantile(ratios, probs = 1 - alpha)
  return(as.numeric(c_empirical))
}  


#------------------------------------Monte Carlo simulation to tune an empirical penalty for Poisson-----------------------------------------------------------
tune_cusum_penalty_poisson <- function(n, lambda0 = 5, alpha = 0.05, replicates = 1000) {
  ratios <- numeric(replicates)
  
  for (r in 1:replicates) {

    y <- rpois(n, lambda = lambda0)
    res <- detect_changepoint_poisson(y, c = Inf )

    ratios[r] <- res$C_max

  }
  
  # Empirical quantile (penalty threshold)
  c_empirical <- quantile(ratios, probs = 1 - alpha)
  return(as.numeric(c_empirical))
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tune_threshold_mean_variance <- function(y, M = 500, alpha = 0.95) {
  
  n <- length(y)
  mu_hat <- mean(y)
  sigma_hat <- sd(y)
  
  max_llr <- numeric(M)
  
  for (i in 1:M) {
    # simulate null Gaussian data
    sim <- rnorm(n, mu_hat, sigma_hat)
    
    # compute LLR and store the max statistic
    result <- detect_change_mean_variance(sim, c = -Inf)
    max_llr[i] <- max(result$LLR_values)
  }
  
  threshold <- quantile(max_llr, alpha)
  
  return(list(
    threshold = threshold,
    simulated_maxima = max_llr
  ))
}

#-----------------------------------------------------Comparing alpha values--------------------------------------------------------------------------------------  
# --- Theoretical threshold function (Gumbel quantile, scaled) -----------------
# theoretical_threshold <- function(n, alpha, scale_factor = 1) {
#   a_n <- 1 / sqrt(2 * log(log(n)))
#   b_n <- 2 * log(log(n)) + 0.5 * log(log(log(n))) - 0.5 * log(pi)
#   # Gumbel quantile for upper tail 1 - alpha
#   q_alpha <- -log(-0.5 * log(1 - alpha))
#   # Scale to match empirical statistic
#   return((b_n + a_n * q_alpha)^2)
# }
# 
# n_values <- c(100, 500, 1000, 10000)
# alpha_values <- c(0.1, 0.05, 0.01)
# sigma2 <- 1
# replicates <- 1000
# 
# results <- expand.grid(n = n_values, alpha = alpha_values)
# results$empirical_c <- NA
# results$theoretical_c <- NA
# 
# for (i in 1:nrow(results)) {
# n <- results$n[i]
# alpha <- results$alpha[i]
# 
# cat("Running simulation for n =", n, "alpha =", alpha, "\n")
# 
# # Empirical threshold
# results$empirical_c[i] <- tune_cusum_penalty(n, sigma2, alpha, replicates)
# 
# # Theoretical threshold (choose one form, e.g. 2 * log(log(n)))
# results$theoretical_c[i] <- theoretical_threshold(n, alpha)
# }
# 
# # --- Plot comparison-----------------------------------------------------------------------------------------------------------------
# 
# 
# ggplot(results, aes(x = n)) +
# geom_line(aes(y = empirical_c, color = as.factor(alpha), linetype = "Empirical"), linewidth = 1) +
# geom_line(aes(y = theoretical_c, color = as.factor(alpha), linetype = "Theoretical"), linewidth = 1, alpha = 0.7) +
# scale_x_log10() +
# labs(
#   title = "CUSUM Thresholds: Monte Carlo vs Theoretical",
#   x = "Signal length (n)",
#   y = "CUSUM threshold (c)",
#   color = "Alpha",
#   linetype = "Type"
# ) +
# theme_minimal(base_size = 14)
# 
# #----------Using the monte carlo threshold and Simpson data-------------------------------------------------------------------------------#
# sigma2_hat <- var(y)
# c_empirical <- tune_poisson_cusum_penalty(
#   n = length(y),
#   sigma2 = sigma2_hat,
#   alpha = 0.05,       # desired false alarm rate
#   replicates = 1000   # number of Monte Carlo runs
# )
# 
# result <- detect_change_cusum(y, sigma2 = sigma2_hat, c = c_empirical)
# plot(result$C_values, type = "l", lwd = 2,
#      main = "CUSUM Statistic (Change-in-Mean Detection)",
#      ylab = expression(C[t]^2), xlab = "t")
# abline(h = sigma2_hat * c_empirical, col = "red", lty = 2)
# if (!is.null(result$change_point)) {
#   abline(v = result$change_point, col = "blue", lwd = 2)
#   legend("topright", legend = c("CUSUM", "Threshold", "Detected Change"),
#          col = c("black", "red", "blue"), lty = c(1,2,1), lwd = 2)
# }
# 
# #-----------------------------NBA -----------------------------------------------------------------------------------------------------------#
# cavsoreb <- cavs$orebTeam
# 
# sigma2_hat <- var(cavsoreb)
# c_empirical <- tune_cusum_penalty_poisson(
#   n = length(cavsoreb),
#   sigma2_hat,
#   alpha = 0.05,       # desired false alarm rate
#   replicates = 1000   # number of Monte Carlo runs
# )
# 
# c_empirical
# 
# 
# #-------------At most one change (changepoint)--------------------------------------------------------------------------------------
# par(mfrow=c(3,2))
# for (i in 37:45){
#   x <- cavs[,i]
#   sigma2_hat <- var(x)
#   c_empirical <- tune_cusum_penalty_poisson(
#     n = length(x),
#     sigma2_hat,
#     alpha = 0.05,       # desired false alarm rate
#     replicates = 1000   # number of Monte Carlo runs
#   )
#   a <- detect_changepoint_poisson(x,c_empirical)
#   graphplot(a)
#   abline(v = a$changepoint)
#   plot(x)
#   
# }
# par(mfrow=c(1,1))
# 
# 
# 
# set.seed(123)
# y <- rpois(100, lambda = 5)
# res <- detect_changepoint_poisson(y, c = Inf)
# str(res)


#-------------------------------COMPARING MODELS BETWEEN GAUSSIAN AND MEAN-VARIANCE MODEL----------------------------------
compare_models <- function(y, cusum_sigma2, cusum_c, mv_c) {
  
  # Run Gaussian CUSUM
  cusum_res <- detect_change_cusum(
    y = y,
    sigma2 = cusum_sigma2,
    c = cusum_c
  )
  
  # Run Mean–Variance LRT
  mv_res <- detect_change_mean_variance(
    y = y,
    c = mv_c
  )
  
  list(
    cusum = cusum_res,
    mean_variance = mv_res
  )
}
plot_changepoints <- function(y, results) {
  n <- length(y)
  
  cusum_cp <- results$cusum$change_point
  mv_cp <- results$mean_variance$change_point
  
  plot(y, type = "l", lwd = 2,
       main = "Time Series with CUSUM + Mean-Variance Changepoints",
       ylab = "Value", xlab = "Time")
  
  if (!is.null(cusum_cp)) {
    abline(v = cusum_cp, col = "blue", lwd = 2, lty = 2)
  }
  
  if (!is.null(mv_cp)) {
    abline(v = mv_cp, col = "red", lwd = 2, lty = 2)
  }
  
  legend("topright",
         legend = c("CUSUM", "Mean-Variance"),
         col = c("blue", "red"),
         lwd = 2, lty = 2)
}
plot_statistics <- function(results) {
  
  par(mfrow=c(2,1))
  
  # Plot CUSUM statistic
  plot(results$cusum$C_values, type="l", col="blue", lwd=2,
       main="CUSUM Detection Statistic (C_t^2)",
       ylab="C_t^2", xlab="t")
  
  abline(h = as.numeric(results$cusum$c), col="darkgrey", lwd=2, lty=3)
  
  
  # Plot LLR values
  plot(results$mean_variance$LLR_values, type="l", col="red", lwd=2,
       main="Mean-Variance LRT Statistic (LLR)",
       ylab="LLR", xlab="t")
  
  abline(h = as.numeric(results$mean_variance$c), col="darkgrey", lwd=2, lty=3)
  
  par(mfrow=c(1,1))
}

# # Example: your plus-minus data vector
# y <- cavs_data$plusminusTeam
# 
# # Tuned penalties
# cusum_c <- c_emp1     
# cusum_sigma2 <- var(y)
# mv_c <- c_emp2$threshold          
# 
# results <- compare_models(y, cusum_sigma2, cusum_c, mv_c)
# 
# plot_changepoints(y, results)
# plot_statistics(results)
# 
# detect_change_mean_variance(
#   y = y,
#   c = mv_c
# )

#-----------------------------CUSUM_ONLINE-----------------------------------------------------------------------------------------------------
cusum_online <- function(x, mu0=0, mu1=5, h=6) {
  S <- 0
  cp <- NA  # changepoint detection time
  
  for (t in 1:length(x)) {
    # Log-likelihood ratio CUSUM (Gaussian)
    k <- (mu0 + mu1) / 2      # reference value
    S <- max(0, S + (x[t] - k))
    
    if (S > h && is.na(cp)) {
      cp <- t  # record first time the threshold is crossed
    }
  }
  
  return(list(cp=cp))
}




plot1 <- function(data){
  res <- cusum_online(data)
  pen <- tune_online(data)
  d <- online_cumsum(data,mu_1 = mean(data),pen)
  
  par(mfrow = c(2,1))
  plot(data, type="l", col="black", lwd=2,
       main=paste("Time Series with Online CUSUM Detection"),
       xlab="Game", ylab=var)
  
  # Add changepoint line if detected
  if (!is.na(res$cp)) {
    abline(v = res$cp, col="red", lwd=3)
    text(res$cp, max(data), pos=4, col="red")
  }
  
  plot(d$C_values,type="l", lwd=2, col="blue",
       main="Online CUSUM Statistic",
       xlab="Game", ylab="CUSUM Value")
  
  abline(v = res$cp, col="red", lwd=2, lty=2)   # threshold
  print(paste("This change point was at ",res$cp))
}
  
  
  
#--------------------------------Tuning for online changepoint detection--------------------------------------------------------------------------------------------
tune_online <- function(mu_0 = 0, mu_1 = 5, sigma = 3, h, max_T = 5000) {
  S <- 0
  for (t in 1:max_T) {
    x <- rnorm(1, mean = mu_0, sd = sigma)  # no-change data
    k <- (mu_0 + mu_1) / 2
    S <- max(0, S + (x - k))
    if (S > h) return(t)  # false alarm time
  }
  return(max_T)  # censored if no alarm
}

estimate_ARL0 <- function(h, sims = 1000) {
  rl <- replicate(sims, tune_online(h = h))
  mean(rl)
}
  
hs <- c(5, 8, 10, 12, 15)
ARL0_vals <- sapply(hs, estimate_ARL0)
data.frame(h = hs, ARL0 = ARL0_vals)
# 

#-------------------------------------OFFICIAL ONLINE ANALYSIS CONDUCTING FOR IN GAMES STATS FOR CAVS----------------------------------------------------------------------------
library(stringr)
#We will be conducting online
cavs_a <- cavs %>% filter(yearSeason == 2016)
cavs_ids <- unique(cavs_a$idGame)
#location game and game_id
cavs_b <- cavs_a[,c("idGame","locationGame")]

pbp <- play_by_play(game_ids = cavs_ids)


# cavs_pbp_merged <- pbp %>%
#   left_join(cavs_b, by = "idGame")
# 
# 
# cavs_pbp_merged <- cavs_pbp_merged %>%
#   mutate(
#     points_cavs = if_else(locationGame == "H", scoreHome, scoreAway),
#     points_opp  = if_else(locationGame == "H", scoreAway, scoreHome)
# )
# 
# 
# cavs_pbp2 <- cavs_pbp_merged %>%
#   arrange(idGame, numberPeriod, desc(secondsRemainingQuarter)) %>%
#   mutate(
#     score_diff_cavs = if_else(locationGame == "H", marginScore, -marginScore)
#   ) %>%
#   group_by(idGame) %>%
#   # NEW: fill upwards first to fill beginning-of-game NA
#   fill(score_diff_cavs, .direction = "up") %>%
#   # Then fill downward
#   fill(score_diff_cavs, .direction = "down") %>%
#   ungroup()
# 
# 
# cavs_pbp2 <- cavs_pbp2 %>%
#   mutate(
#     cavs_turnover = case_when(
#       locationGame == "H" & str_detect(descriptionPlayHome, "Turnover") ~ 1,
#       locationGame == "A" & str_detect(descriptionPlayVisitor, "Turnover") ~ 1,
#       TRUE ~ 0
#     ),
#     
#     cavs_3pa = case_when(
#       locationGame == "H" & str_detect(descriptionPlayHome, "3PT") ~ 1,
#       locationGame == "A" & str_detect(descriptionPlayVisitor, "3PT") ~ 1,
#       TRUE ~ 0
#     )
#   )
# 
# cavs_pbp2 <- cavs_pbp2 %>%
#   group_by(idGame) %>%
#   arrange(numberPeriod, desc(secondsRemainingQuarter)) %>%
#   mutate(
#     # Does this event END the Cavs' possession?
#     cavs_poss_end =
#       cavs_turnover == 1 |                  # turnover
#       (cavs_shot_attempt == 1 & cavs_shot_made == 1) |  # made shot
#       (cavs_shot_attempt == 1 & cavs_shot_made == 0 & opp_def_rebound == 1),
#     
#     # Possession starts when the previous row ends a possession
#     possession_start = lag(cavs_poss_end, default = TRUE),
#     
#     possession_index = cumsum(possession_start)
#   ) %>%
#   ungroup()
# 
# cavs_possessions <- cavs_pbp2 %>%
#   group_by(idGame, possession_index) %>%
#   summarise(
#     net_margin_change = last(score_diff_cavs) - first(score_diff_cavs),  # can be -ve, 0 or +ve
#     turnovers         = sum(cavs_turnover),
#     threepa           = sum(cavs_3pa),
#     .groups = "drop"
#   )
# 
# 
# cavs_segments <- cavs_possessions %>%
#   group_by(idGame) %>%
#   mutate(segment_index = ceiling(possession_index / 5)) %>%
#   ungroup() %>%
#   group_by(idGame, segment_index) %>%
#   summarise(
#     seg_net_margin = sum(net_margin_change),  # Gaussian-ish metric
#     seg_turnovers  = sum(turnovers),         # Poisson
#     seg_3pa        = sum(threepa),           # Poisson
#     .groups = "drop"
#   ) %>%
#   arrange(idGame, segment_index)
# 
# 
# cavs_pbp2 %>%
#   select(idGame, numberPeriod, secondsRemainingQuarter,
#          score_diff_cavs, cavs_turnover,
#          possession_start, possession_index) %>%
#   head(20)
# 
# 
# 
# cavs_pbp2 %>%
#   group_by(idGame) %>%
#   filter(possession_start == TRUE) %>%
#   select(idGame, numberPeriod, secondsRemainingQuarter,
#          score_diff_cavs, cavs_turnover, possession_start, possession_index) %>% head(20)


# 
# library(hoopR)
# library(dplyr)
# progressr::with_progress({
#   nba_pbp <- hoopR::load_nba_pbp()
# })
# 
# games_2015 <- hoopR::nba_schedule(season = 2015)
# 
# cavs_games <- games_2015 %>%
#   filter(home_team_tricode == "CLE" | away_team_tricode == "CLE")
# cavs_ids <- cavs_games$game_id
# 
# options(timeout = 1000)
# 
# pbp_list <- list()
# 
# for (g in cavs_ids) {
#   message("Downloading ", g)
#   Sys.sleep(0.5)   # pacing
#   out <- try(nba_pbp(game_id = g), silent = TRUE)
#   if(!inherits(out, "try-error")) {
#     pbp_list[[g]] <- out
#   }
# }
# 
# cavs_pbp <- bind_rows(pbp_list)
# 
# 
# res <- FOCuS(cavs_a$plusminusTeam,18)
# str(res)
# 
# 
# #------------------------------BIVARIATE ANALYSIS----------------------------------------------------------------------------------
# teams <- sort(unique(df$nameTeam))[-c(4,20,21)]
# length(teams)  # should be 30
# 
# compute_team_cusum <- function(team_name, metric = "plusminusTeam",
#                                sigma2_method = var, threshold_c = 50) {
#   
#   y <- df %>%
#     filter(nameTeam == team_name) %>%
#     arrange(idGame) %>%  
#     pull(metric)
#   
#   sigma2 <- sigma2_method(diff(y))    # variance estimate
#   
#   res <- detect_change_cusum(y, sigma2 = sigma2, c = threshold_c)
#   
#   list(
#     y = y,
#     cusum = res$C_values,
#     cp = res$change_point
#   )
# }
# team_stats <- map(teams, ~ compute_team_cusum(.x))
# names(team_stats) <- teams
# pair_list <- combn(teams, 2, simplify = FALSE)
# length(pair_list) 
# # 
# 
# teams10 <- c("Atlanta Hawks","Boston Celtics","Brooklyn Nets","Chicago Bulls","Cleveland Cavaliers",
#              "Dallas Mavericks","Denver Nuggets","Golden State Warriors","Los Angeles Lakers","Miami Heat")
# 
# grid_10 <- expand_grid(
#   row_team = teams10,
#   col_team = teams10
# )
# 
# 
# pair_df_10 <- grid_10 %>%
#   mutate(data = purrr::map2(row_team, col_team, function(r, c) {
#     
#     if (r == c) {
#       # diagonal: empty cell
#       return(tibble(idx = integer(0),
#                     value = numeric(0),
#                     team = character(0)))
#     }
#     
#     cusA <- team_stats[[r]]$cusum
#     cusB <- team_stats[[c]]$cusum
#     
#     if (length(cusA) == 0 || length(cusB) == 0) {
#       return(tibble(idx = integer(0),
#                     value = numeric(0),
#                     team = character(0)))
#     }
#     
#     tibble(
#       idx   = c(seq_along(cusA), seq_along(cusB)),
#       value = c(cusA, cusB),
#       team  = c(rep(r, length(cusA)),
#                 rep(c, length(cusB)))
#     )
#   }))
# 
# 
# pair_df_long <- pair_df_10 %>%
#   unnest(data)
# 
# ggplot(pair_df_long, aes(idx, value, color = team)) +
#   geom_line(linewidth = 0.4, na.rm = TRUE) +
#   facet_grid(row_team ~ col_team, drop = FALSE) +
#   theme_minimal(base_size = 7) +
#   theme(
#     legend.position = "none",
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     strip.text = element_text(size = 7)
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #---------------------BIVARIATE TRACE FOR TEAMS PAIRWISE-----------------------------------------------------------------
# pair_df_10 <- grid_10 %>%
#   mutate(data = purrr::map2(row_team, col_team, function(r, c) {
#     
#     if (r == c) {
#       return(tibble(idx = integer(0),
#                     value = numeric(0)))
#     }
#     
#     yA <- team_stats[[r]]$cusum
#     yB <- team_stats[[c]]$cusum
#     
#     if (length(yA) == 0 || length(yB) == 0) {
#       return(tibble(idx = integer(0),
#                     value = numeric(0)))
#     }
#     
#     n <- min(length(yA), length(yB))
#     Y <- cbind(yA[1:n], yB[1:n])
#     
#     res <- bivariate_online_cusum_trace(Y)
#     
#     tibble(
#       idx   = seq_along(res$M),
#       value = res$M
#     )
#   }))
# pair_df_long <- pair_df_10 %>%
#   unnest(data)
# 
# ggplot(pair_df_long, aes(idx, value)) +
#   geom_line(linewidth = 0.4, na.rm = TRUE) +
#   facet_grid(row_team ~ col_team, drop = FALSE) +
#   theme_minimal(base_size = 7) +
#   theme(
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     strip.text = element_text(size = 7)
#   )
# #--------------------------------------------------------------------------------------------------------------------
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# compute_online_cusum_path <- function(y, sigma2, c) {
#   
#   n <- length(y)
#   out <- vector("list", n)
#   
#   for (t in 2:n) {
#     res <- detect_change_cusum(y[1:t], sigma2, c)
#     
#     out[[t]] <- data.frame(
#       iteration = t,
#       idx = seq_along(res$C_values),
#       C = res$C_values
#     )
#   }
#   
#   do.call(rbind, out)
# }
# 
# 
# online_df <- compute_online_cusum_path(cavs$plusminusTeam, 1, 100)
# step <- 20   # try 10, 20, or 25
# online_df_thin <- online_df[online_df$iteration %% step == 0, ]
# 
# 
# p_bottom <- ggplot(online_df) +
#   geom_segment(aes(
#     x = 1,
#     xend = idx,
#     y = iteration,
#     yend = iteration
#   ),
#   linewidth = 0.4,
#   alpha = 0.6) +
#   scale_y_reverse() +
#   labs(
#     x = "Candidate changepoint location",
#     y = "Iteration (time)"
#   ) +
#   theme_minimal()
# 
# final_res <- detect_change_cusum(cavs$plusminusTeam, 1, 100)
# 
# final_df <- data.frame(
#   idx = seq_along(final_res$C_values),
#   C = final_res$C_values
# )
# 
# layout(matrix(c(1, 2), nrow = 2), heights = c(1, 3))
# par(mar = c(4, 4, 2, 1))
# plot(final_df$idx, final_df$C,
#      type = "l",
#      lwd = 2,
#      col = "black",
#      xlab = "",
#      ylab = "CUSUM",
#      main = "Final (offline) CUSUM")
# plot(
#   range(online_df$idx),
#   range(online_df_thin$iteration),
#   type = "n",
#   xlab = "Candidate changepoint location",
#   ylab = "Iteration (time)",
#   main = "Online CUSUM evolution"
# )
# 
# segments(
#   x0 = 1,
#   y0 = online_df_thin$iteration,
#   x1 = online_df_thin$idx,
#   y1 = online_df_thin$iteration,
#   col = "grey40"
# )
# 
# 
# 
# 
# 
# 
# 
# # for (pair in pair_list[1:6]) {
# #   print(plot_team_pair(pair[[1]], pair[[2]], team_stats))
# # }
# # 
# # 
# # for (i in (1:30)){
# #   print(paste("Loading length for ",teams[i]))
# #   print(nrow(df %>% filter(nameTeam == teams[i]) %>% filter(yearSeason >= 2014 & yearSeason <= 2022)))
# #   
# # }
# # 
# # print(plot_team_pair(pair_list[[1]][1], pair_list[[1]][2], team_stats))
# 
# 
# par(mfrow = c(1, 1))
# par(mar = c(5, 5, 4, 2))
# 
# 
# teamA_plusminus <- df %>% filter(nameTeam == "Boston Celtics") %>% pull(plusminusTeam)
# teamB_plusminus <- df %>% filter(nameTeam == "Cleveland Cavaliers") %>% pull(plusminusTeam)
# 
# library(FOCuS)
# 
# bivariate_online_cusum_trace_fast <- function(Y) {
#   n <- nrow(Y)
#   d <- ncol(Y)
#   S <- apply(Y, 2, cumsum)
#   
#   M <- rep(NA_real_, n)
#   tauhat <- rep(NA_integer_, n)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     k_star <- NA_integer_
#     
#     for (k in 1:(t - 1)) {
#       m1 <- S[k, ] / k
#       m2 <- (S[t, ] - S[k, ]) / (t - k)
#       Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
#       
#       if (Ct > Ct_max) {
#         Ct_max <- Ct
#         k_star <- k
#       }
#     }
#     
#     M[t] <- Ct_max
#     tauhat[t] <- k_star
#   }
#   
#   list(M = M, tauhat = tauhat)
# }
# library(FOCuS)
# 
# bivar_focus_stopping_time_null <- function(
#     n,
#     h,
#     combine = c("sum", "max"),
#     mu0 = NA,
#     K = Inf
# ) {
#   combine <- match.arg(combine)
#   
#   y1 <- rnorm(n)
#   y2 <- rnorm(n)
#   
#   out1 <- FOCuS(y1, thres = h, mu0 = mu0, K = K)
#   out2 <- FOCuS(y2, thres = h, mu0 = mu0, K = K)
#   
#   s1 <- out1$maxs
#   s2 <- out2$maxs
#   
#   M <- if (combine == "sum") s1 + s2 else pmax(s1, s2)
#   
#   T <- which(M > h)[1]
#   if (is.na(T)) n + 1 else T
# }
# 
# estimate_ARL0_surrogate <- function(
#     n,
#     h,
#     sims = 500,
#     combine = "sum",
#     mu0 = NA,
#     K = Inf
# ) {
#   Ts <- replicate(
#     sims,
#     bivar_focus_stopping_time_null(n, h, combine, mu0, K)
#   )
#   mean(Ts)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# calibrate_h_unified <- function(
#     target_arl0 = 500,
#     n_null = 3000,
#     sims = 500,
#     h_grid = seq(10, 120, by = 5),
#     combine = "sum",
#     mu0 = NA,
#     K = Inf
# ) {
#   arls <- sapply(h_grid, function(h)
#     estimate_ARL0_surrogate(
#       n = n_null,
#       h = h,
#       sims = sims,
#       combine = combine,
#       mu0 = mu0,
#       K = K
#     )
#   )
#   
#   h_grid[which.min(abs(arls - target_arl0))]
# }
# h <- calibrate_h_unified(
#   target_arl0 = 500,
#   n_null = 3000,
#   sims = 500,
#   combine = "sum"
# )
# 
# saveRDS(h, "h_ARL500_bivar_focus_sum.rds")
# 
# 
# 
# 
# 
# h <- readRDS("h_ARL500_bivar_focus_sum.rds")
# 
# Y <- scale(cbind(teamA_plusminus, teamB_plusminus))
# y1 <- Y[,1]
# y2 <- Y[,2]
# 
# out1 <- FOCuS(y1, thres = h)
# out2 <- FOCuS(y2, thres = h)
# 
# M <- out1$maxs + out2$maxs
# T <- which(M > h)[1]
# 
# plot(M, type = "l", ylab = "Combined FOCuS statistic")
# abline(h = h, lty = 2)
# if (!is.na(T)) abline(v = T, lwd = 2)
# 
# 
# 
# 
# 
# 
# #----------------------------------plotting ts and cusum for cavs-------------------------------
# 
# y <- cavs$plusminusTeam
# f <- tune_cusum_penalty(length(y),replicates = 20000)
# vv <- detect_change_cusum(y,1,f)
# graphplot(cavs$plusminusTeam,vv)
# 
# #--------------CONTINUING 1ST TASK - CREATING WIN PERCENTAGE FOR THE 2016 SZN------------------------
# 
# wp <- df %>% filter(yearSeason == 2016) %>% select("nameTeam","isWin")
# table(wp)
# win_pct_df <- wp %>%
#   group_by(nameTeam) %>%
#   summarise(
#     games_played = n(),
#     wins = sum(isWin),
#     win_percentage = (wins / games_played)*100,
#     .groups = "drop"
#   )
# 
# win_pct_ranked <- win_pct_df %>%
#   arrange(desc(win_percentage)) %>%
#   mutate(rank = row_number())
# 
# 
# win_pct_ranked <- win_pct_ranked %>%
#   mutate(
#     tier = case_when(
#       win_percentage >= quantile(win_percentage, 0.80) ~ "High-performing",
#       win_percentage <= quantile(win_percentage, 0.20) ~ "Low-performing",
#       TRUE                                             ~ "Mid-performing"
#     )
#   )
# 
# 
# w_data <- as.data.frame(win_pct_ranked)
# 
# #-------------------------------------------------------------------------------------------------
# 
# run_offline_cusum_for_team <- function(y,
#                                        alpha = al[i],
#                                        replicates = 1000) {
#   n <- length(y)
#   
#   # Null model: demeaned series
#   
#   # Tune penalty
#   penalty_std <- tune_cusum_penalty(
#     n = n,
#     alpha = alpha,
#     replicates = replicates
#   )
#   
#   penalty_raw <- penalty_std
#   
#   # Run offline detection
#   res <- detect_change_cusum(
#     y,
#     1,
#     c = penalty_raw
#   )
#   
#   list(
#     n = n,
#     changepoints = res$change_point
#   )
# }
# 
# df1 <- df %>% filter(yearSeason == 2016)
# teams <- w_data$nameTeam
# 
# team_series <- df1 %>%
#   arrange(numberGameTeamSeason) %>%   # ensure correct order
#   group_split(nameTeam)
# 
# 
# x<- list()
# al <-c(0.01,0.05,0.1)
# 
# for (i in c(1:4)){
# 
#   results_list <- lapply(team_series, function(team_df) {
#     y <- team_df$plusminusTeam
#     
#     run_offline_cusum_for_team(y,alpha = al[i])
#   })
#   
#   
#   cusum_summary <- data.frame(
#     nameTeam = teams,
#     n_games = sapply(results_list, function(x) x$n),
#     changepoints = I(lapply(results_list, function(x) x$changepoints))
#   )
#   x[[i]] <- cusum_summary
# }
# 
# #--------------------------------------------------------ONLINE MULTIPLE CP-------------------------------------------
# #-------------------------------------------------
#   
#   
# y <- cavs$plusminusTeam   # ordered chronologically
# win_loss <- ifelse(y > 0, "Won", "Lost")
# t <- seq_along(y)
# season_breaks <- c(82, 164)
# 
# ggplot(data.frame(t, y, win_loss),
#        aes(x = t, y = y, fill = win_loss)) +
#   geom_col(width = 0.9) +
#   scale_fill_manual(values = c("Won" = "steelblue", "Lost" = "firebrick")) +
#   geom_vline(xintercept = season_breaks, linetype = "dashed") +
#   labs(x = "Matches ordered chronologically",
#        y = "Point differential") +
#   theme_minimal()
# 
# 
# 
# results <- lapply(team_series, function(team_df) {
#   y <- scale(team_df$plusminusTeam)[,1]
#   focus(y, model = "meanvar", penalty = beta)$changepoints
# })
# #------------------------------------------------------------------------------------------------------------------
# bivar_cusum_trace <- function(Y) {
#   n <- nrow(Y)
#   d <- ncol(Y)
#   
#   S <- apply(Y, 2, cumsum)
#   
#   M <- rep(NA_real_, n)
#   tauhat <- rep(NA_integer_, n)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     k_star <- NA_integer_
#     
#     for (k in 1:(t - 1)) {
#       m1 <- S[k, ] / k
#       m2 <- (S[t, ] - S[k, ]) / (t - k)
#       
#       Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
#       
#       if (Ct > Ct_max) {
#         Ct_max <- Ct
#         k_star <- k
#       }
#     }
#     
#     M[t] <- Ct_max
#     tauhat[t] <- k_star
#   }
#   
#   list(M = M, tauhat = tauhat)
# }
# 
# Y <- scale(cbind(teamA_plusminus, teamB_plusminus))
# 
# tr <- bivar_cusum_trace(Y)
# 
# plot(tr$M, type = "l",
#      xlab = "Time",
#      ylab = "Max CUSUM statistic",
#      main = "Trace of max bivariate CUSUM over time")
# 
# 
# 
# 
# 
# 
# bivar_online_cusum_stop <- function(Y, h) {
#   n <- nrow(Y)
#   d <- ncol(Y)
#   
#   # cumulative sums
#   S <- apply(Y, 2, cumsum)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     tau_hat <- NA_integer_
#     
#     for (k in 1:(t - 1)) {
#       m1 <- S[k, ] / k
#       m2 <- (S[t, ] - S[k, ]) / (t - k)
#       
#       Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
#       
#       if (Ct > Ct_max) {
#         Ct_max <- Ct
#         tau_hat <- k
#       }
#     }
#     
#     if (Ct_max > h) {
#       return(list(
#         T = t,
#         tau_hat = tau_hat,
#         stat = Ct_max
#       ))
#     }
#   }
#   
#   list(T = NA, tau_hat = NA, stat = NA)
# }
# 
# 
# 
# 
# bivar_online_cusum_stop <- function(Y, h) {
#   n <- nrow(Y)
#   d <- ncol(Y)
#   
#   # cumulative sums
#   S <- apply(Y, 2, cumsum)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     tau_hat <- NA_integer_
#     
#     for (k in 1:(t - 1)) {
#       m1 <- S[k, ] / k
#       m2 <- (S[t, ] - S[k, ]) / (t - k)
#       
#       Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
#       
#       if (Ct > Ct_max) {
#         Ct_max <- Ct
#         tau_hat <- k
#       }
#     }
#     
#     if (Ct_max > h) {
#       return(list(
#         T = t,
#         tau_hat = tau_hat,
#         stat = Ct_max
#       ))
#     }
#   }
#   
#   list(T = NA, tau_hat = NA, stat = NA)
# }
# 
# bivar_online_cusum_with_trace <- function(Y, h = Inf) {
#   n <- nrow(Y)
#   d <- ncol(Y)
#   
#   # cumulative sums
#   S <- apply(Y, 2, cumsum)
#   
#   M <- rep(NA_real_, n)      # online statistic
#   tauhat <- rep(NA_integer_, n)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     k_star <- NA_integer_
#     
#     for (k in 1:(t - 1)) {
#       m1 <- S[k, ] / k
#       m2 <- (S[t, ] - S[k, ]) / (t - k)
#       
#       Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
#       
#       if (Ct > Ct_max) {
#         Ct_max <- Ct
#         k_star <- k
#       }
#     }
#     
#     M[t] <- Ct_max
#     tauhat[t] <- k_star
#     
#     # optional early stopping
#     if (Ct_max > h) break
#   }
#   
#   list(M = M, tauhat = tauhat)
# }
# 
# h <- 20  # example threshold (use calibrated value in practice)
# 
# out <- bivar_online_cusum_with_trace(Y, h)
# 
# plot(out$M[1:T], type = "l",
#      xlab = "Time",
#      ylab = "Online bivariate CUSUM",
#      main = "Sequential online bivariate CUSUM")
# 
# abline(h = h, lty = 2)
# 
# T <- which(out$M > h)[1]
# if (!is.na(T)) {
#   abline(v = T, col = "black", lwd = 2)
#   abline(v = out$tauhat[T], col = "blue", lwd = 2)
# }
# 
# 
# 
# 
# univariate_cusum_trace <- function(y) {
#   n <- length(y)
#   
#   # cumulative sums
#   S <- cumsum(y)
#   
#   M <- rep(NA_real_, n)      # max CUSUM statistic at each time
#   tauhat <- rep(NA_integer_, n)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     k_star <- NA_integer_
#     
#     for (k in 1:(t - 1)) {
#       m1 <- S[k] / k
#       m2 <- (S[t] - S[k]) / (t - k)
#       
#       Ct <- (k * (t - k) / t) * (m1 - m2)^2
#       
#       if (Ct > Ct_max) {
#         Ct_max <- Ct
#         k_star <- k
#       }
#     }
#     
#     M[t] <- Ct_max
#     tauhat[t] <- k_star
#   }
#   
#   list(M = M, tauhat = tauhat)
# }
# univariate_cusum_trace <- function(y) {
#   n <- length(y)
#   
#   # cumulative sums
#   S <- cumsum(y)
#   
#   M <- rep(NA_real_, n)      # max CUSUM statistic at each time
#   tauhat <- rep(NA_integer_, n)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     k_star <- NA_integer_
#     
#     for (k in 1:(t - 1)) {
#       m1 <- S[k] / k
#       m2 <- (S[t] - S[k]) / (t - k)
#       
#       Ct <- (k * (t - k) / t) * (m1 - m2)^2
#       
#       if (Ct > Ct_max) {
#         Ct_max <- Ct
#         k_star <- k
#       }
#     }
#     
#     M[t] <- Ct_max
#     tauhat[t] <- k_star
#   }
#   
#   list(M = M, tauhat = tauhat)
# }
# 
# 
# y <- scale(cavs$plusminusTeam)
# tr <- univariate_cusum_trace(y)
# plot(tr$M, type = "l",
#      xlab = "Time (games)",
#      ylab = "Max CUSUM statistic",
#      main = "Univariate CUSUM trace (change in mean)")
# 
# 
# univariate_online_cusum <- function(y, h) {
#   n <- length(y)
#   S <- cumsum(y)
#   
#   M_online <- rep(NA_real_, n)
#   
#   for (t in 2:n) {
#     Ct_max <- -Inf
#     for (k in 1:(t - 1)) {
#       m1 <- S[k] / k
#       m2 <- (S[t] - S[k]) / (t - k)
#       Ct <- (k * (t - k) / t) * (m1 - m2)^2
#       if (Ct > Ct_max) Ct_max <- Ct
#     }
#     
#     M_online[t] <- Ct_max
#     if (Ct_max > h) break
#   }
#   
#   M_online
# }
# y_std <- as.numeric(scale(cavs$plusminusTeam))
# M_trace <- univariate_cusum_trace(y_std)
# 
# # choose threshold (example; use calibrated value in practice)
# h <- 20
# 
# # BOTTOM: online detection output
# M_online <- univariate_online_cusum(y_std, h)
# 
# # detection time
# T <- which(M_online > h)[1]
# 
# layout(matrix(c(1, 2), nrow = 2),
#        heights = c(1, 4))  # top is small, bottom is large
# 
# par(mar = c(1, 4, 2, 1))
# 
# # ---- TOP: marginal CUSUM trace ----
# plot(
#   M_trace$M,
#   type = "l",
#   xaxt = "n",
#   xlab = "",
#   ylab = "Max CUSUM",
#   main = "CUSUM trace (marginal) and online detection"
# )
# 
# abline(h = h, lty = 2, col = "grey60")
# abline(v = T, lwd = 2)
# 
# par(mar = c(4, 4, 2, 1))
# 
# plot(
#   M_online[1:T],
#   type = "l",
#   xlab = "Time (games)",
#   ylab = "Online CUSUM",
#   main = "Online detection with marginal CUSUM trace"
# )
# 
# abline(h = h, lty = 2, col = "grey50")
# abline(v = T, lwd = 2)
# 
# 
# 
# 
# 
# 
# 
# 
# #-----------------WIN PROBABILITY ANALYSIS--------------------------------------------------------------------------------------------------------
library(tictoc)
tictoc::tic()
progressr::with_progress({
  nba_pbp <- hoopR::load_nba_pbp(season = 2016)
})
tictoc::toc()


data <- nba_pbp %>% filter(away_team_abbrev == "CLE" | home_team_abbrev == "CLE",period <= 4) %>% arrange(game_id,period_number,desc(clock_minutes),desc(clock_seconds))

cavs_a <- cavs %>% filter(yearSeason == 2016)
cavs_ids <- unique(cavs_a$idGame)
cavs_b <- cavs_a[,c("idGame","locationGame")]
pbp <- play_by_play(game_ids = cavs_ids)


msg <- lapply(pbp$descriptionPlayHome, function(x){
  if (is.na(x)) {
    NA}
  else{
    str_split_1(x," ")[1]
    }
  })

num <- pbp$numberEventMessageType

tictoc::toc()
progressr::with_progress({new_data <- as.data.frame(data.table(msgtype = msg,
                       msgnum  = num))})
tictoc::toc()
pbp_reg <- data %>%
  mutate(
    cavs_home = home_team_id == cavs_ids,
    margin = ifelse(
      cavs_home,
      home_score - away_score,
      away_score - home_score
    )
  )
cavs_id <- unique(pbp_reg$home_team_id[pbp_reg$home_team_name == "Cleveland"])

pbp_reg <- pbp_reg %>%
  mutate(
    seconds_remaining = end_game_seconds_remaining
  )
game_outcomes <- pbp_reg %>%
  group_by(game_id) %>%
  summarise(
    cavs_win = as.integer(last(margin) > 0),
    .groups = "drop"
  )
wp_data <- pbp_reg %>%
  left_join(game_outcomes, by = "game_id") %>%
  select(
    cavs_win,
    margin,
    seconds_remaining,
    cavs_home
  )

wp_model <- glm(
  cavs_win ~ margin * seconds_remaining + cavs_home,
  data = wp_data,
  family = binomial(link = "logit")
)
pbp_reg <- pbp_reg %>%
  mutate(
    win_prob = predict(wp_model, newdata = ., type = "response")
  )


game_id_use <- unique(pbp_reg$game_id)[1]

game_wp <- pbp_reg %>%
  filter(game_id == game_id_use) %>%
  arrange(desc(seconds_remaining))
game_wp <- game_wp %>%
  mutate(
    minutes_remaining = seconds_remaining / 60
  )



ggplot(game_wp, aes(x = minutes_remaining, y = win_prob)) +
  geom_hline(yintercept = 0.5, color = "gray50", linetype = "dashed", size = 1) + # Added the 50% line
  geom_line(linewidth = 0.8, color = "steelblue") +
  scale_x_reverse() +
  scale_y_continuous(labels = scales::percent_format()) + # Makes the Y-axis look like 75% instead of 0.75
  labs(
    x = "Minutes remaining (regulation)",
    y = "Win probability",
    title = "Cavaliers win probability over game time"
  ) +
  theme_minimal(base_size = 13)





# #-----------------EDA GRAPHS--------------------------------------------------------------------------------------------
# 
# 
# cavs_data <- cavs[,c("dateGame","plusminusTeam")]
# 
# change_g_cavs3 <- detect_change_cusum(cavs_data$plusminusTeam,1,Inf)
# change_mv_cavs3 <- detect_change_mean_variance(cavs_data$plusminusTeam,Inf)
# 
# c_emp1 <- tune_cusum_penalty(length(cavs_data$plusminusTeam))
# c_emp2 <- tune_threshold_mean_variance(cavs_data$plusminusTeam) #dont length the input variable
# 
# #Comparing plus minus scores for the mean and variance
# 
# par(mfrow = c(3,1))
# 
# plot.ts(detect_change_mean_variance(cavs_data$plusminusTeam,c_emp2$threshold)$LLR_values, type = "l",main="Plus Minus Changepoints Comparison",
#         xlab="Time", ylab="LLR")
# abline(v = detect_change_mean_variance(cavs_data$plusminusTeam,c_emp2$threshold)$change_point, col="blue", lty=2)  # mean-only
# 
# plot.ts(detect_change_cusum(cavs_data$plusminusTeam,1,c_emp1)$C_values, type="l",,xlab="Time",
#         ylab ="CUMSUM")
# 
# abline(v = detect_change_cusum(cavs_data$plusminusTeam,1,c_emp1)$change_point, col="red", lty=2)   # mean-vari
# 
# 
# df %>% filter(idTeam == 1610612739) %>%
#   ggplot(aes(x = dateGame, y = plusminusTeam)) + geom_point()


library(dplyr)
library(hoopR)

# 1. Load the previous season's team box scores (2014-15)
box_15 <- load_nba_team_box(season = 2015)

# 2. Calculate Possessions and Net Rating per team
team_ratings <- box_15 %>%
  # Note: Column names might vary slightly depending on your hoopR version
  # Standardizing to typical ESPN API column names
  mutate(
    possessions = field_goals_attempted - offensive_rebounds + turnovers + (0.44 * free_throws_attempted)
  ) %>%
  group_by(team_abbreviation) %>%
  summarise(
    total_pts = sum(team_score, na.rm = TRUE),
    total_poss = sum(possessions, na.rm = TRUE),
    # Offensive Rating: Points scored per 100 possessions
    off_rtg = (total_pts / total_poss) * 100,
    .groups = "drop"
  )

# To get Defensive Rating, we need to know what opponents scored against them.
# For simplicity in this script, we approximate Net Rating by comparing a team's Off Rtg to the league average.
league_avg_off_rtg <- mean(team_ratings$off_rtg)

team_ratings <- team_ratings %>%
  mutate(net_rating = off_rtg - league_avg_off_rtg)

# 3. Extract Cavs and Bulls Ratings
cavs_rtg <- team_ratings$net_rating[team_ratings$team_abbreviation == "CLE"]
bulls_rtg <- team_ratings$net_rating[team_ratings$team_abbreviation == "CHI"]

# 4. Calculate Net Rating Differential and convert to Pre-Game Probability
# Cavs are the away team, so we subtract the standard 3.0 point Home Court Advantage for the Bulls
net_diff <- cavs_rtg - bulls_rtg - 3.0 

# Logistic conversion: k = 0.15 is the standard NBA scaling factor
pre_game_prob_cavs <- 1 / (1 + exp(-0.15 * net_diff))

pbp_reg <- nba_pbp %>% 
  filter(away_team_abbrev == "CLE" | home_team_abbrev == "CLE") %>%
  arrange(game_id, period, desc(clock_minutes), desc(clock_seconds)) %>%
  fill(home_score, away_score, .direction = "down") 

# 2. Time Engineering & Margin
pbp_reg <- pbp_reg %>%
  mutate(
    period_length_mins = ifelse(period <= 4, 12, 5),
    sec_elapsed_in_period = (period_length_mins * 60) - (clock_minutes * 60 + clock_seconds),
    prior_sec_elapsed = case_when(
      period == 1 ~ 0,
      period == 2 ~ 720,
      period == 3 ~ 1440,
      period == 4 ~ 2160,
      period > 4 ~ 2880 + (period - 5) * 300
    ),
    total_sec_elapsed = prior_sec_elapsed + sec_elapsed_in_period,
    minutes_elapsed = total_sec_elapsed / 60,
    model_seconds_remaining = case_when(
      period <= 4 ~ 2880 - total_sec_elapsed,
      period > 4 ~ (2880 + (period - 4) * 300) - total_sec_elapsed
    ),
    cavs_home = ifelse(home_team_abbrev == "CLE", 1, 0),
    margin = ifelse(cavs_home == 1, home_score - away_score, away_score - home_score),
    
    # The mathematical fix for late-game volatility
    adjusted_margin = margin / sqrt(model_seconds_remaining + 1)
  )

# 3. Game Outcomes
game_outcomes <- pbp_reg %>%
  group_by(game_id) %>%
  summarise(
    final_margin = last(margin[!is.na(margin)]),
    cavs_win = ifelse(final_margin > 0, 1, 0),
    .groups = "drop"
  )

wp_data <- pbp_reg %>%
  left_join(game_outcomes, by = "game_id") %>%
  filter(!is.na(margin)) %>%
  # Convert our pre_game_prob into log-odds so the logistic model can read it
  mutate(baseline_logit = log(pre_game_prob_cavs / (1 - pre_game_prob_cavs)))

# 4. The Anchored Model
# We remove the intercept (- 1) and use our baseline as the anchor (offset)
wp_model <- glm(
  cavs_win ~ offset(baseline_logit) + adjusted_margin - 1,
  data = wp_data,
  family = binomial(link = "logit")
)

wp_data <- wp_data %>%
  mutate(win_prob = predict(wp_model, newdata = ., type = "response"))

# 5. Plotting
game_id_use <- unique(wp_data$game_id)[1] # Opening night game
game_wp <- wp_data %>% filter(game_id == game_id_use)

ggplot(game_wp, aes(x = minutes_elapsed, y = win_prob)) +
  geom_hline(yintercept = 0.5, color = "gray50", linetype = "dashed", linewidth = 1) + 
  geom_line(linewidth = 0.8, color = "steelblue") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) + 
  labs(
    x = "Minutes elapsed",
    y = "Win probability",
    title = "Cavaliers Win Probability vs. Bulls (Opening Night)"
  ) +
  theme_minimal(base_size = 13)


#--------------------------------------------------
# -------------------------------------------------------------------
# 2. MONTE CARLO THRESHOLD CALIBRATION (FIXED FOR SPORTS RUNS)
# -------------------------------------------------------------------
cat("Calibrating Coach's Threshold using constrained Monte Carlo (GLR)...\n")
set.seed(123)
n_sims <- 1000
baseline_window <- 40
window_size <- 40 
min_window <- 10 # THE FIX: Force the algorithm to look for sustained runs, not single shots

max_glr_values <- numeric(n_sims)

for (i in 1:n_sims) {
  null_play_net <- rnorm(total_plays, mean = 0, sd = 0.8)
  
  mu_0 <- mean(null_play_net[1:baseline_window])
  sigma_0 <- sd(null_play_net[1:baseline_window])
  if(sigma_0 == 0) sigma_0 <- 1
  
  max_lr_sim <- 0
  
  for (t in (baseline_window + 1):total_plays) {
    start_idx <- max(baseline_window + 1, t - window_size)
    
    for (k in start_idx:(t-1)) {
      N <- t - k
      if (N < min_window) next # Ignores random 1-play spikes!
      
      mu_1 <- mean(null_play_net[(k+1):t])
      lr <- N * (mu_1 - mu_0)^2 / (2 * sigma_0^2)
      
      if (lr > max_lr_sim) {
        max_lr_sim <- lr
      }
    }
  }
  max_glr_values[i] <- max_lr_sim
}

coach_h <- quantile(max_glr_values, 0.80) 
cat("Coach's Calibrated GLR Threshold (h):", round(coach_h, 3), "\n\n")

# -------------------------------------------------------------------
# 3. THE CONSTRAINED ONLINE GLR ENGINE
# -------------------------------------------------------------------
online_glr <- function(margin_seq, baseline_window = 40, h_threshold, window_size = 40, min_window = 10, lockout_length = 15) {
  n <- length(margin_seq)
  play_net <- c(0, diff(margin_seq)) 
  
  mu_0 <- mean(play_net[1:baseline_window])
  sigma_0 <- sd(play_net[1:baseline_window])
  if(sigma_0 == 0) sigma_0 <- 1
  
  glr_stats <- numeric(n) 
  alarms <- c()
  lockout_timer <- 0 
  
  cat("Running Constrained Online GLR Scan on Live Data...\n")
  for (t in (baseline_window + 1):n) {
    if (lockout_timer > 0) {
      lockout_timer <- lockout_timer - 1
      next
    }
    
    start_idx <- max(baseline_window + 1, t - window_size)
    max_lr <- 0
    
    for (k in start_idx:(t-1)) {
      N <- t - k
      if (N < min_window) next # Hunt for momentum, not moments
      
      mu_1 <- mean(play_net[(k+1):t])
      lr <- N * (mu_1 - mu_0)^2 / (2 * sigma_0^2)
      
      if (lr > max_lr) {
        max_lr <- lr
      }
    }
    
    glr_stats[t] <- max_lr
    
    if (glr_stats[t] > h_threshold) {
      alarms <- c(alarms, t)
      lockout_timer <- lockout_length 
    }
  }
  
  return(list(alarms = alarms, glr_stats = glr_stats))
}

res <- online_glr(pbp_data$margin, h_threshold = coach_h, min_window = min_window)






# =============================================================================
# IN-GAME ONLINE CHANGEPOINT DETECTION
# Early Warning System for Coaching Strategy
#
# Architecture:
#   1. Win probability model (narrative layer)
#   2. Scoring margin rate in fixed time bins (detection layer)
#   3. Online GLR with resets (multiple changepoints per game)
#   4. Contextual output combining detection with WP narrative
#
# The GLR monitors the binned margin rate series. When an alarm fires,
# we report:
#   - When the shift started (estimated changepoint tau_hat)
#   - The magnitude of the scoring rate change
#   - The win probability before and after the shift
#   - An actionable recommendation
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# SECTION 1: WIN PROBABILITY
#
# The WP model is fitted ONCE across all games in the season (outside this
# pipeline), producing wp_data with a win_prob column. This pipeline takes
# wp_data filtered to a single game — no refitting needed.
#
# The user's existing code does:
#   wp_model <- glm(cavs_win ~ offset(baseline_logit) + adjusted_margin - 1,
#                   data = wp_data, family = binomial)
#   wp_data$win_prob <- predict(wp_model, newdata = wp_data, type = "response")
#
# We simply consume that output.
# =============================================================================


# =============================================================================
# SECTION 2: MARGIN RATE SERIES (Detection Layer)
# =============================================================================

#' Bin the game into fixed time intervals and compute scoring margin per bin
#'
#' @param pbp           Play-by-play dataframe (with margin column)
#' @param bin_seconds   Width of each time bin in seconds (default 60 = 1 min)
#' @return Dataframe with one row per bin: bin_index, bin_start, bin_end,
#'         margin_start, margin_end, margin_rate (points scored per bin),
#'         and the win probability at the bin midpoint
bin_margin_rate <- function(pbp, bin_seconds = 60) {
  
  # Total game length in seconds
  max_sec <- max(pbp$total_sec_elapsed, na.rm = TRUE)
  
  # Create bin breaks
  breaks <- seq(0, max_sec + bin_seconds, by = bin_seconds)
  
  pbp_scored <- pbp %>%
    filter(!is.na(margin), !is.na(total_sec_elapsed)) %>%
    mutate(bin = cut(total_sec_elapsed, breaks = breaks,
                     include.lowest = TRUE, labels = FALSE))
  
  bins <- pbp_scored %>%
    group_by(bin) %>%
    summarise(
      bin_start_sec  = min(total_sec_elapsed),
      bin_end_sec    = max(total_sec_elapsed),
      margin_start   = first(margin),
      margin_end     = last(margin),
      margin_change  = last(margin) - first(margin),
      wp_start       = first(win_prob),
      wp_end         = last(win_prob),
      wp_mid         = median(win_prob, na.rm = TRUE),
      n_events       = n(),
      .groups = "drop"
    ) %>%
    arrange(bin) %>%
    mutate(
      bin_index      = row_number(),
      minutes_start  = bin_start_sec / 60,
      minutes_end    = bin_end_sec / 60
    )
  
  return(bins)
}


# =============================================================================
# SECTION 3: ONLINE GLR WITH RESETS
# =============================================================================

#' GLR statistic at a single time t (univariate — same as your existing code)
glr_at_t <- function(y, t, sigma2, k_min = 3) {
  
  if (t < 2 * k_min) return(list(stat = 0, best_k = NA))
  
  cumS <- cumsum(y[1:t])
  S_t  <- cumS[t]
  
  max_lr <- 0
  best_k <- NA
  
  for (k in k_min:(t - k_min)) {
    ybar1 <- cumS[k] / k
    ybar2 <- (S_t - cumS[k]) / (t - k)
    C_t2  <- (k * (t - k) / t) * (ybar1 - ybar2)^2
    lr    <- C_t2 / sigma2
    if (!is.na(lr) && lr > max_lr) {
      max_lr <- lr
      best_k <- k
    }
  }
  
  list(stat = max_lr, best_k = best_k)
}


#' Online GLR with resets — finds multiple changepoints in a single game
#'
#' After each alarm, the detector resets: the post-changepoint data becomes
#' the new baseline, and scanning restarts from the estimated changepoint.
#'
#' @param y        Margin rate series (one value per time bin)
#' @param h        Detection threshold
#' @param k_min    Minimum segment length (in bins)
#' @return List of detected changepoints with alarm times and diagnostics
online_glr_with_resets <- function(y, h, k_min = 3) {
  
  n <- length(y)
  alarms <- list()
  alarm_count <- 0
  
  # Track the full GLR sequence for plotting
  glr_full <- numeric(n)
  
  # Current segment start (shifts after each reset)
  seg_start <- 1
  
  for (t_abs in 1:n) {
    
    # Relative position within current segment
    t_rel <- t_abs - seg_start + 1
    
    if (t_rel < 2 * k_min) {
      glr_full[t_abs] <- 0
      next
    }
    
    # Extract current segment
    y_seg <- y[seg_start:t_abs]
    
    # Estimate sigma from current segment (robust, same as your code)
    if (length(y_seg) > 2) {
      sigma2 <- (mad(diff(y_seg)) / sqrt(2))^2
      if (sigma2 == 0 || is.na(sigma2)) sigma2 <- var(y_seg)
      if (sigma2 == 0 || is.na(sigma2)) sigma2 <- 1  # fallback
    } else {
      sigma2 <- 1
    }
    
    res <- glr_at_t(y_seg, t_rel, sigma2, k_min)
    glr_full[t_abs] <- res$stat
    
    # Check for alarm
    if (!is.na(res$stat) && res$stat > h) {
      alarm_count <- alarm_count + 1
      
      # Estimated changepoint in absolute coordinates
      tau_abs <- seg_start + res$best_k - 1
      
      # Compute mean shift
      pre_mean  <- mean(y_seg[1:res$best_k])
      post_mean <- mean(y_seg[(res$best_k + 1):t_rel])
      shift     <- post_mean - pre_mean
      
      alarms[[alarm_count]] <- list(
        alarm_bin  = t_abs,
        tau_bin    = tau_abs,
        glr_stat   = res$stat,
        pre_mean   = pre_mean,
        post_mean  = post_mean,
        shift      = shift,
        sigma2     = sigma2,
        seg_start  = seg_start,
        seg_length = t_rel
      )
      
      # Reset: new segment starts at the estimated changepoint
      seg_start <- tau_abs + 1
    }
  }
  
  list(
    alarms   = alarms,
    n_alarms = alarm_count,
    glr_full = glr_full,
    h        = h,
    n        = n
  )
}


# =============================================================================
# SECTION 4: THRESHOLD CALIBRATION FOR IN-GAME DETECTION
# =============================================================================

#' Calibrate threshold for the in-game margin rate GLR
#'
#' Simulates null margin rate series (no changepoint) and finds the
#' (1 - alpha) quantile of the max GLR over the game.
#'
#' @param n_bins      Number of time bins per game (e.g., 48 for 1-min bins)
#' @param sigma2      Noise variance for simulation (default 1, standardised)
#' @param alpha       Significance level
#' @param replicates  Number of MC replicates
#' @param k_min       Minimum segment length
#' @return Calibrated threshold h
tune_ingame_threshold <- function(n_bins     = 48,
                                  sigma2     = 1,
                                  alpha      = 0.05,
                                  replicates = 2000,
                                  k_min      = 3,
                                  seed       = 42) {
  set.seed(seed)
  cat(sprintf("Calibrating in-game GLR threshold (n_bins=%d, alpha=%.2f)...\n",
              n_bins, alpha))
  
  max_glr_null <- numeric(replicates)
  
  for (r in 1:replicates) {
    y_null <- rnorm(n_bins, mean = 0, sd = sqrt(sigma2))
    
    max_g <- 0
    for (t in 1:n_bins) {
      if (t < 2 * k_min) next
      res <- glr_at_t(y_null, t, sigma2, k_min)
      if (res$stat > max_g) max_g <- res$stat
    }
    max_glr_null[r] <- max_g
  }
  
  h <- quantile(max_glr_null, probs = 1 - alpha)
  cat(sprintf("-> Threshold h = %.4f\n\n", h))
  return(as.numeric(h))
}


# =============================================================================
# SECTION 5: CONTEXTUAL OUTPUT — ACTIONABLE INSIGHTS
# =============================================================================

#' Generate coaching insights from detected changepoints
#'
#' @param alarms     List of alarms from online_glr_with_resets()
#' @param bins       Binned margin rate dataframe from bin_margin_rate()
#' @param focal_team Abbreviation of the focal team
#' @return Dataframe of insights with timing, magnitude, WP context, and advice
generate_insights <- function(alarms, bins, focal_team) {
  
  if (length(alarms) == 0) {
    cat("No significant momentum shifts detected in this game.\n")
    return(NULL)
  }
  
  insights <- list()
  
  for (i in seq_along(alarms)) {
    a <- alarms[[i]]
    
    # Look up WP context from bins
    tau_bin   <- a$tau_bin
    alarm_bin <- a$alarm_bin
    
    # WP at changepoint and at alarm
    wp_at_cp    <- if (tau_bin <= nrow(bins)) bins$wp_mid[tau_bin] else NA
    wp_at_alarm <- if (alarm_bin <= nrow(bins)) bins$wp_mid[alarm_bin] else NA
    wp_change   <- wp_at_alarm - wp_at_cp
    
    # Timing
    time_at_cp    <- if (tau_bin <= nrow(bins)) bins$minutes_start[tau_bin] else NA
    time_at_alarm <- if (alarm_bin <= nrow(bins)) bins$minutes_start[alarm_bin] else NA
    
    # Determine period from minutes
    period_at_cp <- case_when(
      time_at_cp <= 12 ~ "Q1",
      time_at_cp <= 24 ~ "Q2",
      time_at_cp <= 36 ~ "Q3",
      time_at_cp <= 48 ~ "Q4",
      TRUE             ~ paste0("OT", ceiling((time_at_cp - 48) / 5))
    )
    
    # Direction: positive shift = focal team improving
    direction <- ifelse(a$shift > 0, "positive", "negative")
    
    # Generate recommendation
    if (direction == "negative") {
      # Focal team is losing momentum
      magnitude <- abs(a$shift)
      if (magnitude > 4) {
        advice <- "URGENT: Major scoring run against. Call timeout immediately. Consider full lineup change."
      } else if (magnitude > 2) {
        advice <- "Significant momentum shift against. Consider timeout or defensive substitution."
      } else {
        advice <- "Minor negative trend detected. Monitor closely — consider adjusting defensive scheme."
      }
    } else {
      # Focal team is gaining momentum
      magnitude <- abs(a$shift)
      if (magnitude > 4) {
        advice <- "Strong run in progress. Maintain current lineup. Opponent likely to call timeout."
      } else if (magnitude > 2) {
        advice <- "Positive momentum building. Keep current rotation — avoid disrupting rhythm."
      } else {
        advice <- "Slight positive trend. Maintain current approach."
      }
    }
    
    insights[[i]] <- data.frame(
      shift_num       = i,
      shift_started   = round(time_at_cp, 1),
      alarm_at        = round(time_at_alarm, 1),
      period          = period_at_cp,
      detection_delay = round(time_at_alarm - time_at_cp, 1),
      margin_shift    = round(a$shift, 2),
      direction       = direction,
      wp_before       = round(wp_at_cp, 3),
      wp_after        = round(wp_at_alarm, 3),
      wp_change       = round(wp_change, 3),
      recommendation  = advice,
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, insights)
  
  # Print formatted output
  for (i in 1:nrow(result)) {
    r <- result[i, ]
    cat(sprintf("\n--- Momentum Shift #%d ---\n", r$shift_num))
    cat(sprintf("  Detected at:     %.1f min (%s)\n", r$alarm_at, r$period))
    cat(sprintf("  Shift started:   %.1f min (delay: %.1f min)\n",
                r$shift_started, r$detection_delay))
    cat(sprintf("  Scoring change:  %+.1f pts/min (%s)\n",
                r$margin_shift, r$direction))
    cat(sprintf("  Win probability: %.0f%% -> %.0f%% (%+.0f%%)\n",
                r$wp_before * 100, r$wp_after * 100, r$wp_change * 100))
    cat(sprintf("  >> %s\n", r$recommendation))
  }
  
  return(result)
}


# =============================================================================
# SECTION 6: PLOTTING
# =============================================================================

plot_ingame_detection <- function(bins, detection_result, focal_team, opp_team) {
  
  alarms <- detection_result$alarms
  h      <- detection_result$h
  
  par(mfrow = c(3, 1), mar = c(4, 4.5, 3, 1))
  
  n <- nrow(bins)
  
  # --- Panel 1: Win Probability ---
  plot(bins$minutes_start, bins$wp_mid, type = "l", col = "steelblue", lwd = 1.5,
       main = paste0(focal_team, " Win Probability"),
       xlab = "Minutes Elapsed", ylab = "Win Prob",
       ylim = c(0, 1))
  abline(h = 0.5, col = "grey60", lty = 2)
  
  # Mark detected shifts
  for (a in alarms) {
    if (a$tau_bin <= n) {
      abline(v = bins$minutes_start[a$tau_bin], col = "red", lty = 2, lwd = 1.5)
    }
  }
  
  # Quarter boundaries
  abline(v = c(12, 24, 36, 48), col = "grey80", lty = 3)
  
  # --- Panel 2: Cumulative Margin ---
  plot(bins$minutes_start, bins$margin_end, type = "l", col = "darkorange", lwd = 1.5,
       main = paste0(focal_team, " Scoring Margin"),
       xlab = "Minutes Elapsed", ylab = "Margin")
  abline(h = 0, col = "grey60", lty = 2)
  
  for (a in alarms) {
    if (a$tau_bin <= n) {
      abline(v = bins$minutes_start[a$tau_bin], col = "red", lty = 2, lwd = 1.5)
    }
  }
  abline(v = c(12, 24, 36, 48), col = "grey80", lty = 3)
  
  # --- Panel 3: GLR Statistic ---
  glr <- detection_result$glr_full
  plot(1:length(glr), glr, type = "l", col = "black", lwd = 1.5,
       main = paste0("Online GLR Statistic (h = ", round(h, 2), ")"),
       xlab = "Time Bin", ylab = expression(G[t]))
  abline(h = h, col = "blue", lty = 3, lwd = 2)
  
  for (a in alarms) {
    abline(v = a$alarm_bin, col = "red", lty = 2, lwd = 1.5)
    abline(v = a$tau_bin, col = "purple", lty = 3, lwd = 1.5)
  }
  
  legend("topright",
         legend = c("GLR", "Threshold", "Alarm", "Est. CP"),
         col = c("black", "blue", "red", "purple"),
         lty = c(1, 3, 2, 3), lwd = 2, cex = 0.6, bg = "white")
}


# =============================================================================
# SECTION 7: FULL GAME ANALYSIS PIPELINE
# =============================================================================

#' Run the complete in-game detection pipeline for a single game
#'
#' @param game_wp      wp_data filtered to a single game (must have: margin,
#'                     total_sec_elapsed, minutes_elapsed, win_prob)
#' @param focal_team   Focal team abbreviation
#' @param opp_team     Opponent abbreviation
#' @param bin_seconds  Time bin width in seconds (default 60)
#' @param h            GLR threshold (if NULL, calibrates automatically)
#' @param alpha        Significance level for threshold calibration
#' @param k_min        Minimum segment length in bins
#' @param plot         Whether to produce plots
#' @return List with bins, detection results, and insights
analyse_game <- function(game_wp, focal_team, opp_team,
                         bin_seconds = 60,
                         h = NULL, alpha = 0.05, k_min = 3,
                         plot = TRUE) {
  
  cat(sprintf("\n=== IN-GAME ANALYSIS: %s vs %s ===\n", focal_team, opp_team))
  
  # Check required columns
  required <- c("margin", "total_sec_elapsed", "win_prob")
  missing  <- setdiff(required, names(game_wp))
  if (length(missing) > 0) {
    stop("Missing columns in game_wp: ", paste(missing, collapse = ", "),
         "\nMake sure you pass wp_data filtered to a single game.")
  }
  
  cat(sprintf("  %d play-by-play events\n", nrow(game_wp)))
  
  # Step 1: Bin the margin rate
  cat(sprintf("Binning margin rate (%d-second bins)...\n", bin_seconds))
  bins <- bin_margin_rate(game_wp, bin_seconds = bin_seconds)
  n_bins <- nrow(bins)
  cat(sprintf("  %d time bins created\n", n_bins))
  
  # Step 2: Calibrate threshold if needed
  if (is.null(h)) {
    h <- tune_ingame_threshold(n_bins = n_bins, alpha = alpha, k_min = k_min)
  }
  
  # Step 3: Run online GLR with resets on the margin_change series
  cat("Running online GLR with resets...\n")
  detection <- online_glr_with_resets(bins$margin_change, h = h, k_min = k_min)
  cat(sprintf("  %d momentum shifts detected\n", detection$n_alarms))
  
  # Step 4: Generate insights
  insights <- generate_insights(detection$alarms, bins, focal_team)
  
  # Step 5: Plot
  if (plot) {
    plot_ingame_detection(bins, detection, focal_team, opp_team)
  }
  
  return(list(
    bins       = bins,
    detection  = detection,
    insights   = insights,
    h          = h,
    focal_team = focal_team,
    opp_team   = opp_team
  ))
}


# =============================================================================
# SECTION 8: EXAMPLE USAGE
# =============================================================================
# Prerequisites: run your existing WP pipeline first to produce wp_data
# with win_prob column. Then:
#
# # Pick a game
game_id_use <- unique(wp_data$game_id)[1]  # Opening night
game_wp <- wp_data %>% filter(game_id == game_id_use)

# Determine opponent
opp_team <- ifelse(game_wp$home_team_abbrev[1] == "CLE",
                   game_wp$away_team_abbrev[1],
                   game_wp$home_team_abbrev[1])

# Run the full pipeline
result <- analyse_game(
  game_wp    = game_wp,
  focal_team = "CLE",
  opp_team   = opp_team,
  bin_seconds = 60,    # 1-minute bins
  alpha      = 0.05,
  k_min      = 3,      # minimum 3-bin segment (~3 minutes)
  plot       = TRUE
)

# Access the insights
print(result$insights)

# Try different bin sizes for sensitivity
result_2min <- analyse_game(
  game_wp    = game_wp,
  focal_team = "CLE",
  opp_team   = opp_team,
  bin_seconds = 120,   # 2-minute bins
  alpha = 0.05, k_min = 3, plot = TRUE
)
# 
# library(ggplot2)
# 
# # Assuming 'cavs_pm' is a vector of the Cavaliers' chronological plus-minus values
# # (the same y_1, ..., y_n sequence you use throughout)
# 
# # Base R ACF plot
# acf_result <- acf(cavs$plusminusTeam, lag.max = 20, plot = FALSE)
# 
# # ggplot2 version to match your existing figures
# acf_df <- data.frame(
#   lag = acf_result$lag[-1],    # remove lag 0 (always 1)
#   acf = acf_result$acf[-1]
# )
# 
# n <- length(cavs$plusminusTeam)
# ci_bound <- qnorm(0.975) / sqrt(n)  # 95% confidence band under H0 of independence
# 
# ggplot(acf_df, aes(x = lag, y = acf)) +
#   geom_hline(yintercept = 0, colour = "grey50") +
#   geom_hline(yintercept = c(-ci_bound, ci_bound), 
#              linetype = "dashed", colour = "red") +
#   geom_segment(aes(xend = lag, yend = 0), linewidth = 0.8) +
#   geom_point(size = 2) +
#   labs(
#     title = "Autocorrelation of Game-Level Plus-Minus (Cleveland Cavaliers, 2010–2025)",
#     subtitle = "Dashed lines indicate 95% confidence bounds under independence.",
#     x = "Lag (games)",
#     y = "ACF"
#   ) +
#   theme_minimal()
# 
# # Ljung-Box test (report this one-liner in the text)
# box_test <- Box.test(cavs$plusminusTeam, lag = 10, type = "Ljung-Box")
# print(box_test)
# 
# 
# #---------------------
# # After estimating changepoints, create segment-demeaned residuals
# # Suppose your changepoint locations split the series into segments
# # e.g., cps <- c(0, tau1, tau2, n) for two changepoints
# 
# cps <- c(0,45, 60,n)
# 
# segments <- cut(1:length(cavs$plusminusTeam), breaks = cps, include.lowest = TRUE)
# segment_means <- tapply(cavs$plusminusTeam, segments, mean)
# residuals <- cavs$plusminusTeam - segment_means[segments]
# 
# acf(residuals, lag.max = 20, main = "ACF of Residuals After Segment De-meaning")
# Box.test(residuals, lag = 10, type = "Ljung-Box")
# 

#--------------------------
## ============================================================================
## Offline Univariate Changepoint Detection: All NBA Teams (2010–2025)
## Produces a league-wide summary table mapping detected changepoints
## to game dates and pre/post segment means.


library(dplyr)
library(tidyr)
library(knitr)

# ---- 1. Source your changepoint functions ------------------------------------
# Adjust this path to wherever your functions file lives

# ---- 2. Data loading ---------------------------------------------------------
# Adjust to your actual data pipeline. Expects a dataframe with:
#   slugTeam, dateGame, plusminusTeam, seasonYear
#
# Example:
# all_games <- readRDS("nba_gamelogs_2010_2025.rds")

# ---- 3. Configuration --------------------------------------------------------
REPLICATES     <- 20000
ALPHA          <- 0.05
KMIN           <- 5
RUN_MEAN_VAR   <- TRUE
MV_REPLICATES  <- 500
all_games <- df
# ---- 4. Team list ------------------------------------------------------------
nba_teams <- sort(unique(all_games$slugTeam))
cat("Teams found:", length(nba_teams), "\n")

# ---- 5. Main analysis loop ---------------------------------------------------

results_list <- list()

for (team in nba_teams) {
  
  cat(sprintf("\n=== Processing: %s ===\n", team))
  
  team_data <- all_games %>%
    filter(slugTeam == team) %>%
    arrange(dateGame)
  
  seasons <- sort(unique(team_data$yearSeason))
  
  for (season in seasons) {
    
    season_df <- team_data %>%
      filter(yearSeason == season) %>%
      arrange(dateGame)
    
    y <- season_df$plusminusTeam
    n <- length(y)
    
    if (n < 2 * KMIN + 1) {
      cat(sprintf("  %s %d: skipped (n=%d too short)\n", team, season, n))
      next
    }
    
    # MLE variance — used only for normalising the statistic across teams
    sigma2 <- var(y) * (n - 1) / n
    
    # =====================================================================
    # MEAN MODEL
    # =====================================================================
    
    # tune_cusum_penalty(n, sigma2=1, alpha=0.05, replicates=1000)
    # Your function simulates at unit scale internally, so pass sigma2=1
    threshold_mean <- tune_cusum_penalty(n)
    
    # detect_change_cusum(y, sigma2, c)
    # Call with sigma2=1 to match the threshold scale
    result_mean <- detect_change_cusum(y, 1, threshold_mean)
    
    detected_mean <- !is.null(result_mean$change_point)
    cp_mean       <- result_mean$change_point
    stat_mean     <- result_mean$C_max / sigma2
    
    if (detected_mean) {
      pre_mean  <- mean(y[1:cp_mean])
      post_mean <- mean(y[(cp_mean + 1):n])
      cp_date   <- as.character(season_df$dateGame[cp_mean])
    } else {
      cp_mean   <- NA
      pre_mean  <- NA
      post_mean <- NA
      cp_date   <- NA
    }
    
    row_mean <- data.frame(
      team        = team,
      season      = season,
      n_games     = n,
      model       = "mean",
      reject      = detected_mean,
      changepoint = ifelse(detected_mean, cp_mean, NA),
      cp_date     = cp_date,
      statistic   = round(stat_mean, 4),
      threshold   = round(threshold_mean, 4),
      pre_mean    = round(pre_mean, 2),
      post_mean   = round(post_mean, 2),
      delta       = round(ifelse(detected_mean, post_mean - pre_mean, NA), 2),
      stringsAsFactors = FALSE
    )
    
    results_list[[length(results_list) + 1]] <- row_mean
    
    cat(sprintf("  %s %d (mean):     %s | Game %s | stat=%.3f\n",
                team, season,
                ifelse(detected_mean, "DETECTED", "no change"),
                ifelse(detected_mean, as.character(cp_mean), "-"),
                stat_mean))
    
    # =====================================================================
    # MEAN-VARIANCE MODEL
    # =====================================================================
    if (RUN_MEAN_VAR) {
      
      # tune_threshold_mean_variance(y, M=500, alpha=0.95)
      # alpha is the quantile level, so 0.95 = 95th percentile = 5% significance
      mv_tune      <- tune_threshold_mean_variance(length(y))
      threshold_mv <- mv_tune$threshold
      
      # detect_change_mean_variance(y, c)
      result_mv <- detect_change_mean_variance(y, threshold_mv)
      
      detected_mv <- !is.null(result_mv$change_point)
      cp_mv       <- result_mv$change_point
      stat_mv     <- result_mv$LLR
      
      if (detected_mv) {
        pre_mean_mv  <- mean(y[1:cp_mv])
        post_mean_mv <- mean(y[(cp_mv + 1):n])
        pre_var_mv   <- var(y[1:cp_mv])
        post_var_mv  <- var(y[(cp_mv + 1):n])
        cp_date_mv   <- as.character(season_df$dateGame[cp_mv])
      } else {
        cp_mv        <- NA
        pre_mean_mv  <- NA
        post_mean_mv <- NA
        pre_var_mv   <- NA
        post_var_mv  <- NA
        cp_date_mv   <- NA
      }
      
      row_mv <- data.frame(
        team        = team,
        season      = season,
        n_games     = n,
        model       = "mean-var",
        reject      = detected_mv,
        changepoint = ifelse(detected_mv, cp_mv, NA),
        cp_date     = cp_date_mv,
        statistic   = round(stat_mv, 4),
        threshold   = round(as.numeric(threshold_mv), 4),
        pre_mean    = round(pre_mean_mv, 2),
        post_mean   = round(post_mean_mv, 2),
        delta       = round(ifelse(detected_mv, post_mean_mv - pre_mean_mv, NA), 2),
        stringsAsFactors = FALSE
      )
      
      results_list[[length(results_list) + 1]] <- row_mv
      
      cat(sprintf("  %s %d (mean-var): %s | Game %s | LLR=%.3f\n",
                  team, season,
                  ifelse(detected_mv, "DETECTED", "no change"),
                  ifelse(detected_mv, as.character(cp_mv), "-"),
                  stat_mv))
    }
  }
}

# ---- 6. Assemble results -----------------------------------------------------
results_df <- bind_rows(results_list)

# ---- 7. Output tables ---------------------------------------------------------

# 7a. FULL TABLE (appendix) -- every team-season-model combination
full_table <- results_df %>%
  arrange(team, season, model)

write.csv(full_table, "appendix_full_changepoint_results.csv", row.names = FALSE)
cat("\nFull results saved: appendix_full_changepoint_results.csv\n")

# 7b. DETECTED ONLY -- ranked by statistic strength
detected_only <- results_df %>%
  filter(reject == TRUE) %>%
  arrange(desc(statistic))

write.csv(detected_only, "detected_changepoints.csv", row.names = FALSE)

# 7c. TOP CHANGEPOINT PER TEAM (mean model) -- expanded Table 2
top_by_team <- results_df %>%
  filter(model == "mean", reject == TRUE) %>%
  group_by(team) %>%
  slice_max(statistic, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(statistic))

cat("\n========== Top Changepoints by Team (Mean Model) ==========\n")
print(
  top_by_team %>%
    select(team, season, changepoint, cp_date, statistic, pre_mean, post_mean, delta),
  n = 30
)

# 7d. MODEL COMPARISON -- mean vs mean-var agreement
if (RUN_MEAN_VAR) {
  model_comparison <- results_df %>%
    select(team, season, model, reject, changepoint, statistic) %>%
    pivot_wider(
      names_from  = model,
      values_from = c(reject, changepoint, statistic),
      names_sep   = "_"
    ) %>%
    mutate(
      agreement = case_when(
        `reject_mean` & `reject_mean-var` &
          abs(`changepoint_mean` - `changepoint_mean-var`) <= 10 ~ "agree (within 10 games)",
        `reject_mean` & `reject_mean-var` ~ "both detect, locations differ",
        `reject_mean` & !`reject_mean-var`  ~ "mean only",
        !`reject_mean` & `reject_mean-var`  ~ "mean-var only",
        TRUE                                ~ "neither"
      )
    ) %>%
    arrange(team, season)
  
  write.csv(model_comparison, "model_comparison.csv", row.names = FALSE)
  
  cat("\n========== Model Agreement Summary ==========\n")
  print(table(model_comparison$agreement))
}

# 7e. DETECTION RATE BY SEASON
detection_rate <- results_df %>%
  filter(model == "mean") %>%
  group_by(season) %>%
  summarise(
    n_teams     = n(),
    n_detected  = sum(reject),
    detect_rate = round(n_detected / n_teams, 3),
    avg_stat    = round(mean(statistic), 3),
    .groups     = "drop"
  )

cat("\n========== Detection Rate by Season ==========\n")
print(detection_rate, n = 20)

# 7f. LaTeX TABLE for main text
cat("\n========== LaTeX: Top Changepoints Per Team ==========\n")
top_by_team %>%
  select(Team = team, Season = season, Game = changepoint,
         Date = cp_date, Statistic = statistic,
         Pre = pre_mean, Post = post_mean,
         Delta = delta) %>%
  kable(format = "latex", booktabs = TRUE,
        caption = "Strongest detected offline changepoint per team (mean model, 2010--2025).") %>%
  cat()

cat("\n\nDone.\n")


#--gfhrtyjtjyht-------------------------

glr_at_t <- function(cumS, t, sigma2, k_min = 5) {
  # cumS: cumulative sum of y up to current t (length >= t)
  # We only look at k in [k_min, t - k_min] — same buffer as offline
  if (t < 2 * k_min) return(list(stat = 0, best_k = NA))
  
  S_t    <- cumS[t]
  max_lr <- 0
  best_k <- NA
  
  for (k in k_min:(t - k_min)) {
    ybar1 <- cumS[k] / k
    ybar2 <- (S_t - cumS[k]) / (t - k)
    C_t2  <- (k * (t - k) / t) * (ybar1 - ybar2)^2
    lr    <- C_t2 / sigma2
    if (!is.na(lr) && lr > max_lr) {
      max_lr <- lr
      best_k <- k
    }
  }
  
  list(stat = max_lr, best_k = best_k)
}

# ----------------------------------------------------------------------
# 2b. Full sequential online GLR detector
#     Runs game-by-game and stops at first alarm
# ----------------------------------------------------------------------
online_glr_detect <- function(y, h, k_min = 5) {
  # y      : full season plus-minus vector (in game order)
  # h      : calibrated threshold
  # k_min  : buffer (5 games, same as offline)
  
  y     <- y[!is.na(y)]
  n     <- length(y)
  
  # Estimate sigma2 robustly — IDENTICAL to offline procedure
  sigma2 <- (mad(diff(y)) / sqrt(2))^2
  if (sigma2 == 0) sigma2 <- var(y)
  
  cumS      <- cumsum(y)
  glr_seq   <- numeric(n)   # GLR statistic at each game t
  k_seq     <- rep(NA_integer_, n)  # best k at each game t
  
  T_alarm   <- NA_integer_  # stopping time
  tau_hat   <- NA_integer_  # estimated changepoint at alarm
  
  for (t in 1:n) {
    res          <- glr_at_t(cumS, t, sigma2, k_min)
    glr_seq[t]   <- res$stat
    k_seq[t]     <- res$best_k
    
    # Stop at first crossing — sequential detection
    if (!is.na(res$stat) && res$stat > h && is.na(T_alarm)) {
      T_alarm  <- t
      tau_hat  <- res$best_k
      break   # genuine sequential: stop as soon as alarm fires
    }
  }
  
  list(
    T_alarm   = T_alarm,
    tau_hat   = tau_hat,
    glr_seq   = glr_seq,
    k_seq     = k_seq,
    sigma2    = sigma2,
    h         = h,
    n         = n
  )
}

# ======================================================================
# SECTION 3: MONTE CARLO THRESHOLD CALIBRATION
#
# We calibrate h so that under the null (no changepoint),
# the probability of a false alarm within one season is alpha = 0.05.
#
# This mirrors the offline calibration in tune_cusum_penalty():
#   - Simulate N seasons of length n under N(mu, sigma2) with NO change
#   - For each simulation, compute max_t G_t (the max GLR over the season)
#   - h = (1 - alpha) quantile of that distribution
#
# Critically: we simulate with mean drawn from the EMPIRICAL league
# distribution (not fixed at 0) to avoid the mu_0 = 0 assumption.
# However, because the GLR statistic only depends on DIFFERENCES between
# segment means, any fixed mean cancels out — the null distribution is
# invariant to mu. We simulate with mean = 0 for simplicity, which is
# mathematically equivalent to any other fixed mean for this statistic.
# (This mirrors the same reasoning that validates tune_cusum_penalty.)
# ======================================================================

tune_online_glr_threshold <- function(n          = 82,
                                      sigma2     = 1,
                                      alpha      = 0.05,
                                      replicates = 2000,
                                      k_min      = 5,
                                      seed       = 123) {
  set.seed(seed)
  cat(sprintf("Calibrating online GLR threshold (n=%d, alpha=%.2f, reps=%d)...\n",
              n, alpha, replicates))
  
  max_glr_null <- numeric(replicates)
  
  for (r in 1:replicates) {
    # Null: iid N(0,1) — mean = 0 is valid because the statistic is
    # invariant to location (same justification as tune_cusum_penalty)
    y_null <- rnorm(n, mean = 0, sd = sqrt(sigma2))
    cumS   <- cumsum(y_null)
    
    max_g  <- 0
    for (t in 1:n) {
      res <- glr_at_t(cumS, t, sigma2 = sigma2, k_min = k_min)
      if (res$stat > max_g) max_g <- res$stat
    }
    max_glr_null[r] <- max_g
  }
  
  h <- quantile(max_glr_null, probs = 1 - alpha)
  cat(sprintf("-> Calibrated threshold h = %.4f\n\n", h))
  return(as.numeric(h))
}



#-------------------------------------






## ============================================================================
## Online Univariate GLR Changepoint Detection: All NBA Teams (2010-2025)
## Generalises the study-case online GLR to a full league-wide analysis.
## Produces summary tables comparing online vs offline results.
## ============================================================================

library(dplyr)
library(tidyr)
library(knitr)

# ---- 1. Source your functions ------------------------------------------------
# This file should contain:
#   glr_at_t(), online_glr_detect(), tune_online_glr_threshold()
# AND your offline functions:
#   detect_change_cusum(), tune_cusum_penalty()
 # adjust path as needed

# ---- 2. Data loading ---------------------------------------------------------
# Expects a dataframe `df` (or `all_games`) with:
#   slugTeam, yearSeason (or seasonYear), numberGameTeamSeason, plusminusTeam, dateGame
#
# Example:
# df <- readRDS("nba_gamelogs_2010_2025.rds")

# ---- 3. Configuration --------------------------------------------------------
REPLICATES_ONLINE  <- 2000
REPLICATES_OFFLINE <- 20000
ALPHA              <- 0.05
KMIN               <- 5
AGREEMENT_MARGIN   <- 10   # within 10 games = agreement

# ---- 4. Helper: extract a team-season series ---------------------------------
# Adjust column names to match your dataframe
get_series <- function(df, team, season) {
  season_df <- df %>%
    filter(slugTeam == team, yearSeason == season) %>%
    arrange(numberGameTeamSeason)
  
  y     <- season_df$plusminusTeam
  dates <- season_df$dateGame
  valid <- !is.na(y)
  
  list(y = y[valid], dates = dates[valid])
}

# ---- 5. Calibrate online thresholds -----------------------------------------
# Different season lengths need different thresholds:
#   n = 82 (standard), n = 72 (2020-21 COVID season)
# Pre-calibrate both to avoid recalculating inside the loop

cat("=== Calibrating online GLR thresholds ===\n")

h_82 <- tune_online_glr_threshold(
  n = 82, sigma2 = 1, alpha = ALPHA,
  replicates = REPLICATES_ONLINE, k_min = KMIN, seed = 123
)

h_72 <- tune_online_glr_threshold(
  n = 72, sigma2 = 1, alpha = ALPHA,
  replicates = REPLICATES_ONLINE, k_min = KMIN, seed = 456
)

# ---- 6. Team list and season list --------------------------------------------
nba_teams <- sort(unique(df$slugTeam))
cat("Teams found:", length(nba_teams), "\n\n")

# ---- 7. Main analysis loop ---------------------------------------------------

results_list <- list()

for (team in nba_teams) {
  
  cat(sprintf("\n=== Processing: %s ===\n", team))
  
  team_data <- df %>%
    filter(slugTeam == team)
  
  seasons <- sort(unique(team_data$yearSeason))
  
  for (season in seasons) {
    
    series <- get_series(df, team, season)
    y      <- series$y
    dates  <- series$dates
    n      <- length(y)
    
    if (n < 2 * KMIN + 1) {
      cat(sprintf("  %s %d: skipped (n=%d)\n", team, season, n))
      next
    }
    
    # --- Pick the right pre-calibrated threshold ---
    if (n <= 72) {
      h <- h_72
    } else {
      h <- h_82
    }
    
    # =================================================================
    # ONLINE GLR
    # =================================================================
    res_online <- online_glr_detect(y, h = h, k_min = KMIN)
    
    T_alarm     <- res_online$T_alarm
    tau_online  <- res_online$tau_hat
    detected_on <- !is.na(T_alarm)
    
    if (detected_on) {
      delay       <- T_alarm - tau_online
      cp_date_on  <- as.character(dates[tau_online])
      alarm_date  <- as.character(dates[T_alarm])
      pre_mean_on <- mean(y[1:tau_online])
      post_mean_on <- mean(y[(tau_online + 1):n])
    } else {
      delay        <- NA
      cp_date_on   <- NA
      alarm_date   <- NA
      pre_mean_on  <- NA
      post_mean_on <- NA
    }
    
    # Peak GLR even if no alarm (useful for near-miss analysis)
    peak_glr <- max(res_online$glr_seq, na.rm = TRUE)
    
    # =================================================================
    # OFFLINE CUSUM (for comparison)
    # =================================================================
    sigma2_offline <- var(y) * (n - 1) / n
    threshold_off  <- tune_cusum_penalty(n)
    res_offline    <- detect_change_cusum(y, 1, threshold_off)
    
    detected_off <- !is.null(res_offline$change_point)
    tau_offline  <- res_offline$change_point
    stat_offline <- res_offline$C_max / sigma2_offline
    
    if (detected_off) {
      cp_date_off   <- as.character(dates[tau_offline])
      pre_mean_off  <- mean(y[1:tau_offline])
      post_mean_off <- mean(y[(tau_offline + 1):n])
    } else {
      tau_offline   <- NA
      cp_date_off   <- NA
      pre_mean_off  <- NA
      post_mean_off <- NA
    }
    
    # =================================================================
    # AGREEMENT between online and offline
    # =================================================================
    if (detected_on && detected_off) {
      agree <- abs(tau_online - tau_offline) <= AGREEMENT_MARGIN
    } else {
      agree <- FALSE
    }
    
    agreement_cat <- case_when(
      detected_on & detected_off & agree   ~ "agree (within 10 games)",
      detected_on & detected_off & !agree  ~ "both detect, locations differ",
      detected_on & !detected_off          ~ "online only",
      !detected_on & detected_off          ~ "offline only",
      TRUE                                 ~ "neither"
    )
    
    # =================================================================
    # Store row
    # =================================================================
    row <- data.frame(
      team             = team,
      season           = season,
      n_games          = n,
      # Online results
      online_detected  = detected_on,
      online_tau       = ifelse(detected_on, tau_online, NA),
      online_cp_date   = cp_date_on,
      alarm_time       = ifelse(detected_on, T_alarm, NA),
      alarm_date       = alarm_date,
      detection_delay  = ifelse(detected_on, delay, NA),
      peak_glr         = round(peak_glr, 4),
      threshold        = round(h, 4),
      glr_ratio        = round(peak_glr / h, 4),
      online_pre_mean  = round(pre_mean_on, 2),
      online_post_mean = round(post_mean_on, 2),
      # Offline results
      offline_detected = detected_off,
      offline_tau      = ifelse(detected_off, tau_offline, NA),
      offline_cp_date  = cp_date_off,
      offline_stat     = round(stat_offline, 4),
      offline_pre_mean = round(pre_mean_off, 2),
      offline_post_mean = round(post_mean_off, 2),
      # Comparison
      agreement        = agreement_cat,
      stringsAsFactors = FALSE
    )
    
    results_list[[length(results_list) + 1]] <- row
    
    cat(sprintf("  %s %d: online=%s (T=%s, tau=%s, delay=%s) | offline=%s (tau=%s) | %s\n",
                team, season,
                ifelse(detected_on,  "ALARM", "---"),
                ifelse(detected_on,  as.character(T_alarm), "-"),
                ifelse(detected_on,  as.character(tau_online), "-"),
                ifelse(detected_on,  as.character(delay), "-"),
                ifelse(detected_off, "DETECTED", "---"),
                ifelse(detected_off, as.character(tau_offline), "-"),
                agreement_cat))
  }
}

# ---- 8. Assemble results -----------------------------------------------------
results_df <- bind_rows(results_list)

# ---- 9. Output tables ---------------------------------------------------------

# 9a. FULL TABLE (appendix)
write.csv(results_df, "appendix_online_vs_offline_full.csv", row.names = FALSE)
cat("\nFull results saved: appendix_online_vs_offline_full.csv\n")

# 9b. ONLINE DETECTIONS ONLY — ranked by GLR strength
online_detected <- results_df %>%
  filter(online_detected == TRUE) %>%
  arrange(desc(peak_glr))

write.csv(online_detected, "online_detected_changepoints.csv", row.names = FALSE)

cat(sprintf("\nOnline detections: %d out of %d team-seasons (%.1f%%)\n",
            nrow(online_detected), nrow(results_df),
            100 * nrow(online_detected) / nrow(results_df)))

# 9c. AGREEMENT SUMMARY
cat("\n========== Online vs Offline Agreement ==========\n")
print(table(results_df$agreement))

# 9d. NEAR-MISSES — offline detected but online didn't fire,
#     ranked by how close the GLR got to the threshold
near_misses <- results_df %>%
  filter(offline_detected == TRUE, online_detected == FALSE) %>%
  arrange(desc(glr_ratio)) %>%
  select(team, season, peak_glr, threshold, glr_ratio,
         offline_tau, offline_cp_date, offline_stat)

cat("\n========== Near Misses (offline detected, online missed) ==========\n")
print(near_misses, n = 20)

write.csv(near_misses, "online_near_misses.csv", row.names = FALSE)

# 9e. DETECTION DELAY SUMMARY (for cases where online fired)
if (nrow(online_detected) > 0) {
  cat("\n========== Detection Delay Summary ==========\n")
  cat(sprintf("  Median delay: %.1f games\n", median(online_detected$detection_delay, na.rm = TRUE)))
  cat(sprintf("  Mean delay:   %.1f games\n", mean(online_detected$detection_delay, na.rm = TRUE)))
  cat(sprintf("  Min delay:    %d games\n",   min(online_detected$detection_delay, na.rm = TRUE)))
  cat(sprintf("  Max delay:    %d games\n",   max(online_detected$detection_delay, na.rm = TRUE)))
}

# 9f. TOP ONLINE CHANGEPOINTS PER TEAM
top_online_by_team <- results_df %>%
  filter(online_detected == TRUE) %>%
  group_by(team) %>%
  slice_max(peak_glr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(peak_glr))

cat("\n========== Strongest Online Detection Per Team ==========\n")
print(
  top_online_by_team %>%
    select(team, season, online_tau, online_cp_date,
           alarm_time, detection_delay, peak_glr,
           online_pre_mean, online_post_mean),
  n = 30
)

# 9g. DETECTION RATE BY SEASON
det_rate <- results_df %>%
  group_by(season) %>%
  summarise(
    n_teams         = n(),
    n_online        = sum(online_detected),
    n_offline       = sum(offline_detected),
    n_both          = sum(online_detected & offline_detected),
    online_rate     = round(n_online / n_teams, 3),
    offline_rate    = round(n_offline / n_teams, 3),
    .groups = "drop"
  )

cat("\n========== Detection Rate by Season ==========\n")
print(det_rate, n = 20)

# 9h. LaTeX TABLE — comparative results (main text Table 3 replacement)
cat("\n========== LaTeX: Online vs Offline Comparison ==========\n")
results_df %>%
  filter(online_detected | offline_detected) %>%
  select(Team = team, Season = season,
         `Offline CP` = offline_tau,
         `Online CP`  = online_tau,
         `Stop Time`  = alarm_time,
         Delay        = detection_delay,
         Agreement    = agreement) %>%
  kable(format = "latex", booktabs = TRUE,
        caption = "Comparative analysis of retrospective vs.\\ sequential changepoints across all teams (2010--2025).") %>%
  cat()

cat("\n\nDone.\n")