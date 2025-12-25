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

# simpsons_episodes <- read_csv("https://www.lancaster.ac.uk/~romano/teaching/2425MATH337/datasets/simpsons_episodes.csv")
# simpsons_episodes <- simpsons_episodes |> 
#   mutate(Episode = id + 1, Season = as.factor(season), Rating = tmdb_rating)
# simpsons_episodes <- simpsons_episodes[-nrow(simpsons_episodes), ]
# simp <- simpsons_ratings[,c("Episode","Rating")]

df <- read.csv("NBA_2010_2020.csv") #Reading NBA data
df[df == "LA Clippers"] <- "Los Angeles Clippers"
df %>% filter(idTeam == 1610612739) %>%
  ggplot(aes(x = dateGame, y = pctFG3Team)) + geom_point()

cavs <- filter(df, idTeam == "1610612739")

cavs_data <- cavs[,c("dateGame","plusminusTeam")]

change_g_cavs3 <- detect_change_cusum(cavs_data$plusminusTeam,1,Inf)
change_mv_cavs3 <- detect_change_mean_variance(cavs_data$plusminusTeam,Inf)

c_emp1 <- tune_cusum_penalty(length(cavs_data$plusminusTeam))
c_emp2 <- tune_threshold_mean_variance(cavs_data$plusminusTeam) #dont length the input variable

#Comparing plus minus scores for the mean and variance

par(mfrow = c(3,1))

plot.ts(detect_change_mean_variance(cavs_data$plusminusTeam,c_emp2$threshold)$LLR_values, type = "l",main="Plus Minus Changepoints Comparison",
        xlab="Time", ylab="LLR")
abline(v = detect_change_mean_variance(cavs_data$plusminusTeam,c_emp2$threshold)$change_point, col="blue", lty=2)  # mean-only

plot.ts(detect_change_cusum(cavs_data$plusminusTeam,1,c_emp1)$C_values, type="l",,xlab="Time",
        ylab ="CUMSUM")

abline(v = detect_change_cusum(cavs_data$plusminusTeam,1,c_emp1)$change_point, col="red", lty=2)   # mean-vari


# 


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




#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





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



library(hoopR)
library(dplyr)
progressr::with_progress({
  nba_pbp <- hoopR::load_nba_pbp()
})

games_2015 <- hoopR::nba_schedule(season = 2015)

cavs_games <- games_2015 %>%
  filter(home_team_tricode == "CLE" | away_team_tricode == "CLE")
cavs_ids <- cavs_games$game_id

options(timeout = 1000)

pbp_list <- list()

for (g in cavs_ids) {
  message("Downloading ", g)
  Sys.sleep(0.5)   # pacing
  out <- try(nba_pbp(game_id = g), silent = TRUE)
  if(!inherits(out, "try-error")) {
    pbp_list[[g]] <- out
  }
}

cavs_pbp <- bind_rows(pbp_list)


res <- FOCuS(cavs_a$plusminusTeam,18)
str(res)


#------------------------------BIVARIATE ANALYSIS----------------------------------------------------------------------------------
teams <- sort(unique(df$nameTeam))[-c(4,20,21)]
length(teams)  # should be 30

compute_team_cusum <- function(team_name, metric = "plusminusTeam",
                               sigma2_method = var, threshold_c = 50) {
  
  y <- df %>%
    filter(nameTeam == team_name) %>%
    arrange(idGame) %>%  
    pull(metric)
  
  sigma2 <- sigma2_method(diff(y))    # variance estimate
  
  res <- detect_change_cusum(y, sigma2 = sigma2, c = threshold_c)
  
  list(
    y = y,
    cusum = res$C_values,
    cp = res$change_point
  )
}
team_stats <- map(teams, ~ compute_team_cusum(.x))
names(team_stats) <- teams
pair_list <- combn(teams, 2, simplify = FALSE)
length(pair_list) 
# 

teams10 <- c("Atlanta Hawks","Boston Celtics","Brooklyn Nets","Chicago Bulls","Cleveland Cavaliers",
             "Dallas Mavericks","Denver Nuggets","Golden State Warriors","Los Angeles Lakers","Miami Heat")

grid_10 <- expand_grid(
  row_team = teams10,
  col_team = teams10
)


pair_df_10 <- grid_10 %>%
  mutate(data = purrr::map2(row_team, col_team, function(r, c) {
    
    if (r == c) {
      # diagonal: empty cell
      return(tibble(idx = integer(0),
                    value = numeric(0),
                    team = character(0)))
    }
    
    cusA <- team_stats[[r]]$cusum
    cusB <- team_stats[[c]]$cusum
    
    if (length(cusA) == 0 || length(cusB) == 0) {
      return(tibble(idx = integer(0),
                    value = numeric(0),
                    team = character(0)))
    }
    
    tibble(
      idx   = c(seq_along(cusA), seq_along(cusB)),
      value = c(cusA, cusB),
      team  = c(rep(r, length(cusA)),
                rep(c, length(cusB)))
    )
  }))


pair_df_long <- pair_df_10 %>%
  unnest(data)

ggplot(pair_df_long, aes(idx, value, color = team)) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  facet_grid(row_team ~ col_team, drop = FALSE) +
  theme_minimal(base_size = 7) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 7)
  )









#---------------------BIVARIATE TRACE FOR TEAMS PAIRWISE-----------------------------------------------------------------
pair_df_10 <- grid_10 %>%
  mutate(data = purrr::map2(row_team, col_team, function(r, c) {
    
    if (r == c) {
      return(tibble(idx = integer(0),
                    value = numeric(0)))
    }
    
    yA <- team_stats[[r]]$cusum
    yB <- team_stats[[c]]$cusum
    
    if (length(yA) == 0 || length(yB) == 0) {
      return(tibble(idx = integer(0),
                    value = numeric(0)))
    }
    
    n <- min(length(yA), length(yB))
    Y <- cbind(yA[1:n], yB[1:n])
    
    res <- bivariate_online_cusum_trace(Y)
    
    tibble(
      idx   = seq_along(res$M),
      value = res$M
    )
  }))
pair_df_long <- pair_df_10 %>%
  unnest(data)

ggplot(pair_df_long, aes(idx, value)) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  facet_grid(row_team ~ col_team, drop = FALSE) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 7)
  )
#--------------------------------------------------------------------------------------------------------------------



























compute_online_cusum_path <- function(y, sigma2, c) {
  
  n <- length(y)
  out <- vector("list", n)
  
  for (t in 2:n) {
    res <- detect_change_cusum(y[1:t], sigma2, c)
    
    out[[t]] <- data.frame(
      iteration = t,
      idx = seq_along(res$C_values),
      C = res$C_values
    )
  }
  
  do.call(rbind, out)
}


online_df <- compute_online_cusum_path(cavs$plusminusTeam, 1, 100)
step <- 20   # try 10, 20, or 25
online_df_thin <- online_df[online_df$iteration %% step == 0, ]


p_bottom <- ggplot(online_df) +
  geom_segment(aes(
    x = 1,
    xend = idx,
    y = iteration,
    yend = iteration
  ),
  linewidth = 0.4,
  alpha = 0.6) +
  scale_y_reverse() +
  labs(
    x = "Candidate changepoint location",
    y = "Iteration (time)"
  ) +
  theme_minimal()

final_res <- detect_change_cusum(cavs$plusminusTeam, 1, 100)

final_df <- data.frame(
  idx = seq_along(final_res$C_values),
  C = final_res$C_values
)

layout(matrix(c(1, 2), nrow = 2), heights = c(1, 3))
par(mar = c(4, 4, 2, 1))
plot(final_df$idx, final_df$C,
     type = "l",
     lwd = 2,
     col = "black",
     xlab = "",
     ylab = "CUSUM",
     main = "Final (offline) CUSUM")
plot(
  range(online_df$idx),
  range(online_df_thin$iteration),
  type = "n",
  xlab = "Candidate changepoint location",
  ylab = "Iteration (time)",
  main = "Online CUSUM evolution"
)

segments(
  x0 = 1,
  y0 = online_df_thin$iteration,
  x1 = online_df_thin$idx,
  y1 = online_df_thin$iteration,
  col = "grey40"
)







# for (pair in pair_list[1:6]) {
#   print(plot_team_pair(pair[[1]], pair[[2]], team_stats))
# }
# 
# 
# for (i in (1:30)){
#   print(paste("Loading length for ",teams[i]))
#   print(nrow(df %>% filter(nameTeam == teams[i]) %>% filter(yearSeason >= 2014 & yearSeason <= 2022)))
#   
# }
# 
# print(plot_team_pair(pair_list[[1]][1], pair_list[[1]][2], team_stats))


par(mfrow = c(1, 1))
par(mar = c(5, 5, 4, 2))


teamA_plusminus <- df %>% filter(nameTeam == "Boston Celtics") %>% pull(plusminusTeam)
teamB_plusminus <- df %>% filter(nameTeam == "Cleveland Cavaliers") %>% pull(plusminusTeam)

library(FOCuS)

bivariate_online_cusum_trace_fast <- function(Y) {
  n <- nrow(Y)
  d <- ncol(Y)
  S <- apply(Y, 2, cumsum)
  
  M <- rep(NA_real_, n)
  tauhat <- rep(NA_integer_, n)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    k_star <- NA_integer_
    
    for (k in 1:(t - 1)) {
      m1 <- S[k, ] / k
      m2 <- (S[t, ] - S[k, ]) / (t - k)
      Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
      
      if (Ct > Ct_max) {
        Ct_max <- Ct
        k_star <- k
      }
    }
    
    M[t] <- Ct_max
    tauhat[t] <- k_star
  }
  
  list(M = M, tauhat = tauhat)
}
library(FOCuS)

bivar_focus_stopping_time_null <- function(
    n,
    h,
    combine = c("sum", "max"),
    mu0 = NA,
    K = Inf
) {
  combine <- match.arg(combine)
  
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  
  out1 <- FOCuS(y1, thres = h, mu0 = mu0, K = K)
  out2 <- FOCuS(y2, thres = h, mu0 = mu0, K = K)
  
  s1 <- out1$maxs
  s2 <- out2$maxs
  
  M <- if (combine == "sum") s1 + s2 else pmax(s1, s2)
  
  T <- which(M > h)[1]
  if (is.na(T)) n + 1 else T
}

estimate_ARL0_surrogate <- function(
    n,
    h,
    sims = 500,
    combine = "sum",
    mu0 = NA,
    K = Inf
) {
  Ts <- replicate(
    sims,
    bivar_focus_stopping_time_null(n, h, combine, mu0, K)
  )
  mean(Ts)
}








calibrate_h_unified <- function(
    target_arl0 = 500,
    n_null = 3000,
    sims = 500,
    h_grid = seq(10, 120, by = 5),
    combine = "sum",
    mu0 = NA,
    K = Inf
) {
  arls <- sapply(h_grid, function(h)
    estimate_ARL0_surrogate(
      n = n_null,
      h = h,
      sims = sims,
      combine = combine,
      mu0 = mu0,
      K = K
    )
  )
  
  h_grid[which.min(abs(arls - target_arl0))]
}
h <- calibrate_h_unified(
  target_arl0 = 500,
  n_null = 3000,
  sims = 500,
  combine = "sum"
)

saveRDS(h, "h_ARL500_bivar_focus_sum.rds")





h <- readRDS("h_ARL500_bivar_focus_sum.rds")

Y <- scale(cbind(teamA_plusminus, teamB_plusminus))
y1 <- Y[,1]
y2 <- Y[,2]

out1 <- FOCuS(y1, thres = h)
out2 <- FOCuS(y2, thres = h)

M <- out1$maxs + out2$maxs
T <- which(M > h)[1]

plot(M, type = "l", ylab = "Combined FOCuS statistic")
abline(h = h, lty = 2)
if (!is.na(T)) abline(v = T, lwd = 2)































#----------------------------------plotting ts and cusum for cavs-------------------------------

y <- cavs$plusminusTeam
f <- tune_cusum_penalty(length(y),replicates = 20000)
vv <- detect_change_cusum(y,1,f)
graphplot(cavs$plusminusTeam,vv)

#--------------CONTINUING 1ST TASK - CREATING WIN PERCENTAGE FOR THE 2016 SZN------------------------

wp <- df %>% filter(yearSeason == 2016) %>% select("nameTeam","isWin")
table(wp)
win_pct_df <- wp %>%
  group_by(nameTeam) %>%
  summarise(
    games_played = n(),
    wins = sum(isWin),
    win_percentage = (wins / games_played)*100,
    .groups = "drop"
  )

win_pct_ranked <- win_pct_df %>%
  arrange(desc(win_percentage)) %>%
  mutate(rank = row_number())


win_pct_ranked <- win_pct_ranked %>%
  mutate(
    tier = case_when(
      win_percentage >= quantile(win_percentage, 0.80) ~ "High-performing",
      win_percentage <= quantile(win_percentage, 0.20) ~ "Low-performing",
      TRUE                                             ~ "Mid-performing"
    )
  )


w_data <- as.data.frame(win_pct_ranked)

#-------------------------------------------------------------------------------------------------

run_offline_cusum_for_team <- function(y,
                                       alpha = al[i],
                                       replicates = 1000) {
  n <- length(y)
  
  # Null model: demeaned series
  
  # Tune penalty
  penalty_std <- tune_cusum_penalty(
    n = n,
    alpha = alpha,
    replicates = replicates
  )
  
  penalty_raw <- penalty_std
  
  # Run offline detection
  res <- detect_change_cusum(
    y,
    1,
    c = penalty_raw
  )
  
  list(
    n = n,
    changepoints = res$change_point
  )
}

df1 <- df %>% filter(yearSeason == 2016)
teams <- w_data$nameTeam

team_series <- df1 %>%
  arrange(numberGameTeamSeason) %>%   # ensure correct order
  group_split(nameTeam)


x<- list()
al <-c(0.01,0.05,0.1)

for (i in c(1:4)){

  results_list <- lapply(team_series, function(team_df) {
    y <- team_df$plusminusTeam
    
    run_offline_cusum_for_team(y,alpha = al[i])
  })
  
  
  cusum_summary <- data.frame(
    nameTeam = teams,
    n_games = sapply(results_list, function(x) x$n),
    changepoints = I(lapply(results_list, function(x) x$changepoints))
  )
  x[[i]] <- cusum_summary
}

#--------------------------------------------------------ONLINE MULTIPLE CP-------------------------------------------
#-------------------------------------------------
  
  
y <- cavs$plusminusTeam   # ordered chronologically
win_loss <- ifelse(y > 0, "Won", "Lost")
t <- seq_along(y)
season_breaks <- c(82, 164)

ggplot(data.frame(t, y, win_loss),
       aes(x = t, y = y, fill = win_loss)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = c("Won" = "steelblue", "Lost" = "firebrick")) +
  geom_vline(xintercept = season_breaks, linetype = "dashed") +
  labs(x = "Matches ordered chronologically",
       y = "Point differential") +
  theme_minimal()



results <- lapply(team_series, function(team_df) {
  y <- scale(team_df$plusminusTeam)[,1]
  focus(y, model = "meanvar", penalty = beta)$changepoints
})
#------------------------------------------------------------------------------------------------------------------
bivar_cusum_trace <- function(Y) {
  n <- nrow(Y)
  d <- ncol(Y)
  
  S <- apply(Y, 2, cumsum)
  
  M <- rep(NA_real_, n)
  tauhat <- rep(NA_integer_, n)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    k_star <- NA_integer_
    
    for (k in 1:(t - 1)) {
      m1 <- S[k, ] / k
      m2 <- (S[t, ] - S[k, ]) / (t - k)
      
      Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
      
      if (Ct > Ct_max) {
        Ct_max <- Ct
        k_star <- k
      }
    }
    
    M[t] <- Ct_max
    tauhat[t] <- k_star
  }
  
  list(M = M, tauhat = tauhat)
}

Y <- scale(cbind(teamA_plusminus, teamB_plusminus))

tr <- bivar_cusum_trace(Y)

plot(tr$M, type = "l",
     xlab = "Time",
     ylab = "Max CUSUM statistic",
     main = "Trace of max bivariate CUSUM over time")






bivar_online_cusum_stop <- function(Y, h) {
  n <- nrow(Y)
  d <- ncol(Y)
  
  # cumulative sums
  S <- apply(Y, 2, cumsum)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    tau_hat <- NA_integer_
    
    for (k in 1:(t - 1)) {
      m1 <- S[k, ] / k
      m2 <- (S[t, ] - S[k, ]) / (t - k)
      
      Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
      
      if (Ct > Ct_max) {
        Ct_max <- Ct
        tau_hat <- k
      }
    }
    
    if (Ct_max > h) {
      return(list(
        T = t,
        tau_hat = tau_hat,
        stat = Ct_max
      ))
    }
  }
  
  list(T = NA, tau_hat = NA, stat = NA)
}




bivar_online_cusum_stop <- function(Y, h) {
  n <- nrow(Y)
  d <- ncol(Y)
  
  # cumulative sums
  S <- apply(Y, 2, cumsum)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    tau_hat <- NA_integer_
    
    for (k in 1:(t - 1)) {
      m1 <- S[k, ] / k
      m2 <- (S[t, ] - S[k, ]) / (t - k)
      
      Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
      
      if (Ct > Ct_max) {
        Ct_max <- Ct
        tau_hat <- k
      }
    }
    
    if (Ct_max > h) {
      return(list(
        T = t,
        tau_hat = tau_hat,
        stat = Ct_max
      ))
    }
  }
  
  list(T = NA, tau_hat = NA, stat = NA)
}

bivar_online_cusum_with_trace <- function(Y, h = Inf) {
  n <- nrow(Y)
  d <- ncol(Y)
  
  # cumulative sums
  S <- apply(Y, 2, cumsum)
  
  M <- rep(NA_real_, n)      # online statistic
  tauhat <- rep(NA_integer_, n)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    k_star <- NA_integer_
    
    for (k in 1:(t - 1)) {
      m1 <- S[k, ] / k
      m2 <- (S[t, ] - S[k, ]) / (t - k)
      
      Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
      
      if (Ct > Ct_max) {
        Ct_max <- Ct
        k_star <- k
      }
    }
    
    M[t] <- Ct_max
    tauhat[t] <- k_star
    
    # optional early stopping
    if (Ct_max > h) break
  }
  
  list(M = M, tauhat = tauhat)
}

h <- 20  # example threshold (use calibrated value in practice)

out <- bivar_online_cusum_with_trace(Y, h)

plot(out$M[1:T], type = "l",
     xlab = "Time",
     ylab = "Online bivariate CUSUM",
     main = "Sequential online bivariate CUSUM")

abline(h = h, lty = 2)

T <- which(out$M > h)[1]
if (!is.na(T)) {
  abline(v = T, col = "black", lwd = 2)
  abline(v = out$tauhat[T], col = "blue", lwd = 2)
}




univariate_cusum_trace <- function(y) {
  n <- length(y)
  
  # cumulative sums
  S <- cumsum(y)
  
  M <- rep(NA_real_, n)      # max CUSUM statistic at each time
  tauhat <- rep(NA_integer_, n)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    k_star <- NA_integer_
    
    for (k in 1:(t - 1)) {
      m1 <- S[k] / k
      m2 <- (S[t] - S[k]) / (t - k)
      
      Ct <- (k * (t - k) / t) * (m1 - m2)^2
      
      if (Ct > Ct_max) {
        Ct_max <- Ct
        k_star <- k
      }
    }
    
    M[t] <- Ct_max
    tauhat[t] <- k_star
  }
  
  list(M = M, tauhat = tauhat)
}
univariate_cusum_trace <- function(y) {
  n <- length(y)
  
  # cumulative sums
  S <- cumsum(y)
  
  M <- rep(NA_real_, n)      # max CUSUM statistic at each time
  tauhat <- rep(NA_integer_, n)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    k_star <- NA_integer_
    
    for (k in 1:(t - 1)) {
      m1 <- S[k] / k
      m2 <- (S[t] - S[k]) / (t - k)
      
      Ct <- (k * (t - k) / t) * (m1 - m2)^2
      
      if (Ct > Ct_max) {
        Ct_max <- Ct
        k_star <- k
      }
    }
    
    M[t] <- Ct_max
    tauhat[t] <- k_star
  }
  
  list(M = M, tauhat = tauhat)
}


y <- scale(cavs$plusminusTeam)
tr <- univariate_cusum_trace(y)
plot(tr$M, type = "l",
     xlab = "Time (games)",
     ylab = "Max CUSUM statistic",
     main = "Univariate CUSUM trace (change in mean)")


univariate_online_cusum <- function(y, h) {
  n <- length(y)
  S <- cumsum(y)
  
  M_online <- rep(NA_real_, n)
  
  for (t in 2:n) {
    Ct_max <- -Inf
    for (k in 1:(t - 1)) {
      m1 <- S[k] / k
      m2 <- (S[t] - S[k]) / (t - k)
      Ct <- (k * (t - k) / t) * (m1 - m2)^2
      if (Ct > Ct_max) Ct_max <- Ct
    }
    
    M_online[t] <- Ct_max
    if (Ct_max > h) break
  }
  
  M_online
}
y_std <- as.numeric(scale(cavs$plusminusTeam))
M_trace <- univariate_cusum_trace(y_std)

# choose threshold (example; use calibrated value in practice)
h <- 20

# BOTTOM: online detection output
M_online <- univariate_online_cusum(y_std, h)

# detection time
T <- which(M_online > h)[1]

layout(matrix(c(1, 2), nrow = 2),
       heights = c(1, 4))  # top is small, bottom is large

par(mar = c(1, 4, 2, 1))

# ---- TOP: marginal CUSUM trace ----
plot(
  M_trace$M,
  type = "l",
  xaxt = "n",
  xlab = "",
  ylab = "Max CUSUM",
  main = "CUSUM trace (marginal) and online detection"
)

abline(h = h, lty = 2, col = "grey60")
abline(v = T, lwd = 2)

par(mar = c(4, 4, 2, 1))

plot(
  M_online[1:T],
  type = "l",
  xlab = "Time (games)",
  ylab = "Online CUSUM",
  main = "Online detection with marginal CUSUM trace"
)

abline(h = h, lty = 2, col = "grey50")
abline(v = T, lwd = 2)







#------------------------START OF DISSO ANALYSIS---------------------------------------------------------------------

#-----------------EDA--------------------------------------------------------------------------------------------------------
library(tictoc)
tictoc::tic()
progressr::with_progress({
  nba_pbp <- hoopR::load_nba_pbp(season = 2016)
})
tictoc::toc()


df2 <- nba_pbp %>% filter(season ==  2022)