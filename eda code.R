source("C:/Users/joten/Documents/nbacode.R")


# #Declaring NBA data and the cavs data
# df <- read.csv("NBA_2010_2020.csv") #Reading NBA data
# df[df == "LA Clippers"] <- "Los Angeles Clippers"
# df[df == "Charlotte Bobcats"] <- "Charlotte Hornets"
# df[df == "New Orleans Hornets"] <- "New Orleans Pelicans"
# df <- df %>% filter(idTeam != "1610612751")
# 
# cavs <- filter(df, idTeam == "1610612739")
# 
# cavs_data <- cavs[,c("dateGame","plusminusTeam")]
# 
# 
# tictoc::tic()
# progressr::with_progress({
#   nba_pbp <- hoopR::load_nba_pbp(season = 2016)
# })
# tictoc::toc()
# 
# 
# data <- nba_pbp %>% filter(away_team_abbrev == "CLE" | home_team_abbrev == "CLE",period <= 4) %>% arrange(game_id,period_number,desc(clock_minutes),desc(clock_seconds))
# 
# 

# 
bivariate_cusum_single <- function(Y) {
  # Y is an n x 2 matrix
  n <- nrow(Y)
  Cvals <- numeric(n - 1)

  for (k in 1:(n - 1)) {
    mu1 <- colMeans(Y[1:k, , drop = FALSE])
    mu2 <- colMeans(Y[(k + 1):n, , drop = FALSE])
    diff <- mu1 - mu2
    Cvals[k] <- (k * (n - k) / n) * sum(diff^2)
  }

  list(
    tau_hat = which.max(Cvals),
    C_max   = max(Cvals),
    C_vals  = Cvals
  )
}


tune_bivariate_cusum_threshold <- function(n, alpha = 0.05,replicates = 2000){
  Cmax_vals <- numeric(replicates)

  for (r in seq_len(replicates)) {
    # simulate bivariate Gaussian noise
    Y <- matrix(rnorm(2 * n), ncol = 2)

    # compute offline bivariate CUSUM
    res <- bivariate_cusum_single(Y)

    # store maximum statistic
    Cmax_vals[r] <- res$C_max
  }

  # empirical (1 - alpha) quantile
  quantile(Cmax_vals, probs = 1 - alpha, names = FALSE)
}

# 
alpha <- 0.05
n_season <- 82

cutoffs <- list(
  mean = tune_cusum_penalty(
    n = n_season,
    alpha = alpha,
    replicates = 2000
  ),
  meanvar = tune_threshold_mean_variance(
    n = n_season,
    alpha = alpha,
    replicates = 2000
  ),
  bivar = tune_bivariate_cusum_threshold(
    n = n_season,
    alpha = alpha,
    replicates = 2000
  )
)

saveRDS(cutoffs, "global_cutoffs.rds")




bivariate_cusum_test <- function(Y, alpha = 0.05, replicates = 2000) {
  Y <- scale(Y)
  n <- nrow(Y)
  
  res <- bivariate_cusum_single(Y)
  
  c_alpha <- cutoffs$bivar
  
  list(
    tau_hat = res$tau_hat,
    C_max   = res$C_max,
    cutoff  = c_alpha,
    reject  = res$C_max > c_alpha
  )
}























#Creating variables for gaussian change and mean-variance change
change_g_cavs3 <- detect_change_cusum(cavs_data$plusminusTeam,1,Inf)
change_mv_cavs3 <- detect_change_mean_variance(cavs_data$plusminusTeam,Inf)
c_emp1 <- tune_cusum_penalty(length(cavs_data$plusminusTeam))
c_emp2 <- tune_threshold_mean_variance(cavs_data$plusminusTeam) #dont length the input variable



#Comparing plus minus scores for the mean and variance (CAVS) with log likelihood values and cumsum scores
plt1 <- function (){
  par(mfrow = c(2,1))

  plot.ts(detect_change_mean_variance(cavs_data$plusminusTeam,c_emp2$threshold)$LLR_values, type = "l",main="Plus Minus Changepoints Comparison",
          xlab="Time", ylab="LLR")
  abline(v = detect_change_mean_variance(cavs_data$plusminusTeam,c_emp2$threshold)$change_point, col="blue", lty=2)  # mean-only
  
  plot.ts(detect_change_cusum(cavs_data$plusminusTeam,1,c_emp1)$C_values, type="l",,xlab="Time",
          ylab ="CUMSUM")
  
  abline(v = detect_change_cusum(cavs_data$plusminusTeam,1,c_emp1)$change_point, col="red", lty=2)   # mean-vari
}




#time series of the Cavs 
plt2 <- function(df){
  df %>% filter(idTeam == 1610612739) %>%
  ggplot(aes(x = dateGame, y = plusminusTeam)) + geom_point()
}


#Theoretical values of Monte Carlo tuning algorithm
plt3 <- function(){

  theoretical_threshold <- function(n, alpha, scale_factor = 1) {
    a_n <- 1 / sqrt(2 * log(log(n)))
    b_n <- 2 * log(log(n)) + 0.5 * log(log(log(n))) - 0.5 * log(pi)
    # Gumbel quantile for upper tail 1 - alpha
    q_alpha <- -log(-0.5 * log(1 - alpha))
    # Scale to match empirical statistic
    return((b_n + a_n * q_alpha)^2)
  }
  
  n_values <- c(100, 500, 1000, 10000)
  alpha_values <- c(0.1, 0.05, 0.01)
  sigma2 <- 1
  replicates <- 1000
  
  results <- expand.grid(n = n_values, alpha = alpha_values)
  results$empirical_c <- NA
  results$theoretical_c <- NA
  
  for (i in 1:nrow(results)) {
  n <- results$n[i]
  alpha <- results$alpha[i]
  
  cat("Running simulation for n =", n, "alpha =", alpha, "\n")
  
  # Empirical threshold
  results$empirical_c[i] <- tune_cusum_penalty(n, sigma2, alpha, replicates)
  
  # Theoretical threshold (choose one form, e.g. 2 * log(log(n)))
  results$theoretical_c[i] <- theoretical_threshold(n, alpha)
  }
  
  
  ggplot(results, aes(x = n)) +
    geom_line(aes(y = empirical_c, color = as.factor(alpha), linetype = "Empirical"), linewidth = 1) +
    geom_line(aes(y = theoretical_c, color = as.factor(alpha), linetype = "Theoretical"), linewidth = 1, alpha = 0.7) +
    scale_x_log10() +
    labs(
      title = "CUSUM Thresholds: Monte Carlo vs Theoretical",
      x = "Signal length (n)",
      y = "CUSUM threshold (c)",
      color = "Alpha",
      linetype = "Type"
    ) +
    theme_minimal(base_size = 14)
}











#------------------------------BIVARIATE ANALYSIS OF PLUSMINUS SCORES ON GRID---------------------------------------------------------------------------------
plt4 <- function(){
  teams <- sort(unique(df$nameTeam))
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
  names(team_stats) <- unique(df$slugTeam)
  pair_list <- combn(teams, 2, simplify = FALSE)
  length(pair_list)
  #
  
  teams10 <- c("ATL","BOS","BKN","CHI","CLE",
               "DAL","DEN","GSW","LAL","MIA")
  
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
}
  
  
  
#Time series and CHANGEPOINT PLOT FOR CAVALIERS----------------------------------------
plt111 <- function(){
  y <- cavs$plusminusTeam
  f <- tune_cusum_penalty(length(y),replicates = 20000)
  vv <- detect_change_cusum(y,1,f)
  par(mfrow=c(2,1))
  
  plot.ts(y,
          main = paste("Plus Minus Over Time"),
          ylab = "Plus Minus Score",
          xlab = "t",
          lwd = 2)
  
  abline(v = vv$change_point, col = "green", lty = 2,lwd = 3)
  #328 ,  656
  
  plot.ts(vv$C_values,
          main = paste("CUSUM Statistics for","Plus Minus Score"),
          ylab = expression(C[t]^2),
          xlab = "t",
          lwd = 2)
  
  
  abline(v = vv$change_point, col = "green", lty = c(1,2,2),lwd = 3)

} 
  
  
#-----------------POINT DIFFERENTIAL WIN LOSS FOR CAVS----------------------------------------------------
plt11<- function(){
  y <- cavs$plusminusTeam
  win_loss <- ifelse(y > 0, "Won", "Lost")
  t <- seq_along(y)
  season_breaks <- c(328,656)
  
  ggplot(data.frame(t, y, win_loss),
         aes(x = t, y = y, fill = win_loss)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = c("Won" = "steelblue", "Lost" = "firebrick")) +
    geom_vline(xintercept = season_breaks, linetype = "dashed") +
    labs(x = "Matches ordered chronologically",
         y = "Point differential") +
    theme_minimal()
}








#---------------------------------------UNIVARIATE CUSUM MAX TRACE-----------------------------------------
plt5 <- function(){
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
}



# #-------------------------------------TRACE OF MAX BIVARIATE CUSUM----------------------------------------------------------------------------
# plt6 <- function(){
#   #trace of max bivariate cusum
#   par(mfrow = c(1, 1))
#   par(mar = c(5, 5, 4, 2))
#    
#   bivar_cusum_trace <- function(Y) {
#     n <- nrow(Y)
#     d <- ncol(Y)
#   
#     S <- apply(Y, 2, cumsum)
#   
#     M <- rep(NA_real_, n)
#     tauhat <- rep(NA_integer_, n)
#   
#     for (t in 2:n) {
#       Ct_max <- -Inf
#       k_star <- NA_integer_
#   
#       for (k in 1:(t - 1)) {
#         m1 <- S[k, ] / k
#         m2 <- (S[t, ] - S[k, ]) / (t - k)
#   
#         Ct <- (k * (t - k) / t) * sum((m1 - m2)^2)
#   
#         if (Ct > Ct_max) {
#           Ct_max <- Ct
#           k_star <- k
#         }
#       }
#   
#       M[t] <- Ct_max
#       tauhat[t] <- k_star
#     }
#   
#     list(M = M, tauhat = tauhat)
#   }
#   teamA_plusminus <- df %>% filter(nameTeam == "Boston Celtics") %>% pull(plusminusTeam)
#   teamB_plusminus <- df %>% filter(nameTeam == "Cleveland Cavaliers") %>% pull(plusminusTeam)
#   Y <- scale(cbind(teamA_plusminus, teamB_plusminus))
#   
#   tr <- bivar_cusum_trace(Y)
#   
#   plot(tr$M, type = "l",
#        xlab = "Time",
#        ylab = "Max CUSUM statistic",
#        main = "Trace of max bivariate CUSUM over time")
#   
# }
# 
# 
# 
# 
# #---------------------------SEQUENTIAL ONLINE CUMSUM-----------------------------------------------
# plt7 <- function(){
# 
#   out <- bivar_online_cusum_with_trace(Y, h)
#   
#   plot(out$M[1:T], type = "l",
#        xlab = "Time",
#        ylab = "Online bivariate CUSUM",
#        main = "Sequential online bivariate CUSUM")
#   
#   abline(h = h, lty = 2)
#   
#   T <- which(out$M > h)[1]
#   if (!is.na(T)) {
#     abline(v = T, col = "black", lwd = 2)
#     abline(v = out$tauhat[T], col = "blue", lwd = 2)
#   }
# }
# 
# 
# 
# 
# #---------------------WIN PROBABILITY GRAPH FOR CAVS-----------------------------------------
# 
# plt8 <- function(){
#   data <- nba_pbp %>% filter(away_team_abbrev == "CLE" | home_team_abbrev == "CLE",period <= 4) %>% arrange(game_id,period_number,desc(clock_minutes),desc(clock_seconds))
#   
#   
#   cavs_id <- unique(data$home_team_id[data$home_team_name == "Cleveland"])
#   pbp_reg <- data %>%
#     mutate(
#       cavs_home = home_team_id == cavs_id,
#       margin = ifelse(
#         cavs_home,
#         home_score - away_score,
#         away_score - home_score
#       )
#     )
#   pbp_reg <- pbp_reg %>%
#     mutate(
#       seconds_remaining = end_game_seconds_remaining
#     )
#   game_outcomes <- pbp_reg %>%
#     group_by(game_id) %>%
#     summarise(
#       cavs_win = as.integer(last(margin) > 0),
#       .groups = "drop"
#     )
#   wp_data <- pbp_reg %>%
#     left_join(game_outcomes, by = "game_id") %>%
#     select(
#       cavs_win,
#       margin,
#       seconds_remaining,
#       cavs_home
#     )
#   
#   wp_model <- glm(
#     cavs_win ~ margin + seconds_remaining + cavs_home,
#     data = wp_data,
#     family = binomial(link = "logit")
#   )
#   pbp_reg <- pbp_reg %>%
#     mutate(
#       win_prob = predict(wp_model, newdata = ., type = "response")
#     )
#   
#   
#   game_id_use <- unique(pbp_reg$game_id)[1]
#   
#   game_wp <- pbp_reg %>%
#     filter(game_id == game_id_use) %>%
#     arrange(desc(seconds_remaining))
#   game_wp <- game_wp %>%
#     mutate(
#       minutes_remaining = seconds_remaining / 60
#     )
#   
#   ggplot(game_wp, aes(x = minutes_remaining, y = win_prob)) +
#     geom_line(linewidth = 0.8, color = "steelblue") +
#     scale_x_reverse() +
#     labs(
#       x = "Minutes remaining (regulation)",
#       y = "Win probability",
#       title = "CAVS win probability over game time vs CHI"
#     ) +
#     theme_minimal(base_size = 13)
# }


#----------------------ROLLING MEAN AND VARIANCE------------------------------------------------------------
plt9 <- function(){
  
  cavs_df <- df %>%
    filter(slugTeam == "CLE") %>%
    arrange(numberGameTeamSeason) %>%   # or date, but must be ordered
    mutate(idGame = row_number()) %>%
    select(idGame, plusminusTeam)
  
  
  window <- 10
  
  cavs_df <- cavs_df %>%
    mutate(
      roll_mean = rollapply(
        plusminusTeam,
        width = window,
        FUN = mean,
        align = "right",
        fill = NA
      ),
      roll_var = rollapply(
        plusminusTeam,
        width = window,
        FUN = var,
        align = "right",
        fill = NA
      )
    )
  cavs_long <- cavs_df %>%
    select(idGame, roll_mean, roll_var) %>%
    pivot_longer(-idGame)
  
  ggplot(cavs_long, aes(idGame, value)) +
    geom_line() +
    facet_wrap(~ name, scales = "free_y", ncol = 1) +
    labs(
      title = "Cavaliers Rolling Mean and Variance Diagnostics",
      x = "Game",
      y = NULL
    ) +
    theme_minimal()
  #look at these graph, apparent missing v
}
  
#-----TEAM SUMMARY TABLE FOR PM ETC--------------------------------------
team_season_summary <- df %>%
filter(!is.na(plusminusTeam)) %>%
arrange(nameTeam, yearSeason, numberGameTeamSeason) %>%
group_by(nameTeam, yearSeason) %>%
summarise(
  n_games = n(),
  
  mean_pm = mean(plusminusTeam),
  var_pm  = var(plusminusTeam),
  sd_pm   = sd(plusminusTeam),
  
  max_pm  = max(plusminusTeam),
  min_pm  = min(plusminusTeam),
  max_abs_pm = max(abs(plusminusTeam)),
  
  # rolling variance (window = 5)
  max_roll_var_5 = max(
    rollapply(
      plusminusTeam,
      width = 5,
      FUN = var,
      fill = NA,
      align = "right"
    ),
    na.rm = TRUE
  ),
  
  acf1 = acf(plusminusTeam, plot = FALSE)$acf[2],
  
  .groups = "drop"
)


#---------------BOXPLOTS FOR PLUSMINUS SCORES VARIANCE-------------------------------------
plt10 <- function(){
  
  ggplot(team_season_summary,
         aes(x = factor(yearSeason), y = var_pm)) +
    geom_boxplot(outlier.alpha = 0.3) +
    labs(
      x = "Season",
      y = "Variance of plus–minus",
      title = "Distribution of plus–minus variance across teams"
    ) +
    theme_minimal()
}
 #----------------------HEATMAP FOR PLUSMINUS SCORESS-------------------------------
plt12 <- function(){
  
  ggplot(team_season_summary,
         aes(x = var_pm)) +
    geom_histogram(bins = 30, fill = "grey80") +
    geom_vline(
      data = subset(team_season_summary, nameTeam == "Cavaliers"),
      aes(xintercept = var_pm),
      color = "red",
      linewidth = 1
    ) +
    labs(
      x = "Variance of plus–minus",
      title = "Cavaliers volatility in league context"
    ) +
    theme_minimal()
  
  ggplot(team_season_summary,
         aes(x = factor(yearSeason),
             y = nameTeam,
             fill = var_pm)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(
      x = "Season",
      y = "Team",
      fill = "Var(+/-)",
      title = "Heatmap of plus–minus volatility"
    ) +
    theme_minimal()
}



# ---------------PRE,PRESENT AND POS LEBRON--------------------------------------
# plt13 <- function(){
#   cavs_df <- df %>%
#     filter(slugTeam == "CLE") %>%
#     mutate(
#       period = case_when(
#         yearSeason >= 2010 & yearSeason <= 2014 ~ "Pre-LeBron",
#         yearSeason >= 2015 & yearSeason <= 2018 ~ "LeBron era",
#         yearSeason > 2018                  ~ "Post-LeBron",
#         TRUE                            ~ NA_character_),
# 
#       period = factor(
#         period,
#         levels = c(
#           "Pre-LeBron",
#           "LeBron era",
#           "Post-LeBron"
#         )
#       )
#   )
# 
# 
# 
# 
#   ggplot(cavs_df, aes(x = period, y =plusminusTeam)) +
#     geom_boxplot(outlier.alpha = 0.3) +
#     labs(
#       title = "Cavaliers Plus/Minus: Pre vs Post LeBron",
#       x = NULL,
#       y = "Game-level Plus/Minus"
#     ) +
#     theme_minimal()
# }


#---------------OFFFLINE BIVARIATE CUMSUM-------------------------
prepare_bivariate_series <- function(df, teamA, teamB, season_year) {

  A <- df %>%
    filter(slugTeam == teamA, yearSeason %in% season_year) %>%
    arrange(numberGameTeamSeason) %>%
    mutate(t = row_number()) %>%
    select(t, yA = plusminusTeam)

  B <- df %>%
    filter(slugTeam == teamB, yearSeason %in% season_year) %>%
    arrange(numberGameTeamSeason) %>%
    mutate(t = row_number()) %>%
    select(t, yB = plusminusTeam)

  inner_join(A, B, by = "t")
}


Y_df <- prepare_bivariate_series(df, "CLE", "GSW", 2010:2018)



Y_mat <- as.matrix(Y_df[, c("yA", "yB")])
Y_mat <- scale(Y_mat)

res <- bivariate_cusum_test(Y_mat)

test_res <- bivariate_cusum_test(Y_mat)

test_res$reject    # TRUE / FALSE
test_res$tau_hat





#League wide bivariate analysis for 2018-------------------

teams <- unique(df$slugTeam)
pairs <- combn(teams, 2, simplify = FALSE)

results <- lapply(pairs, function(p) {
  Y_df <- prepare_bivariate_series(df, p[1], p[2], 2018)
  if (nrow(Y_df) < 60) return(NULL)

  res <- bivariate_cusum_test(scale(as.matrix(Y_df[, c("yA","yB")])))
  tibble(
    teamA = p[1],
    teamB = p[2],
    tau_hat = res$tau_hat,
    C_max = res$C_max
  )
})

# pair_results <- bind_rows(results)
# c_alpha <- cutoffs$bivar
# 
# pair_results <- pair_results %>%
#   mutate(reject = C_max > c_alpha)
# 
# sig_pairs <-  pair_results %>% filter(reject)
# 
# 
# 









#Online bivariate detection algorithm for plotting


bivariate_online_trace <- function(Y) {
  Y <- scale(Y)
  n <- nrow(Y)
  S <- apply(Y, 2, cumsum)

  M <- numeric(n)
  tauhat <- integer(n)

  for (t in 2:n) {
    Ct <- numeric(t - 1)
    for (k in 1:(t - 1)) {
      m1 <- S[k, ] / k
      m2 <- (S[t, ] - S[k, ]) / (t - k)
      Ct[k] <- (k * (t - k) / t) * sum((m1 - m2)^2)
    }
    tauhat[t] <- which.max(Ct)
    M[t] <- max(Ct)
  }

  list(M = M, tauhat = tauhat)
}







plot_bivariate_pair <- function(Y_df, res, teamA, teamB) {

  Y <- scale(as.matrix(Y_df[, c("yA", "yB")]))
  tau <- res$tau_hat

  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

  # Panel A: Team A
  plot(Y[,1], type = "l", lwd = 2,
       ylab = paste(teamA, "plus/minus (std)"),
       xlab = "Game",
       main = paste("Bivariate changepoint:", teamA, "&", teamB))
  abline(v = tau, col = "red", lty = 2, lwd = 2)

  # Panel B: Team B
  plot(Y[,2], type = "l", lwd = 2,
       ylab = paste(teamB, "plus/minus (std)"),
       xlab = "Game")
  abline(v = tau, col = "red", lty = 2, lwd = 2)

  # Panel C: Bivariate CUSUM trace
  trace <- bivariate_online_trace(Y)
  plot(trace$M, type = "l", lwd = 2,
       ylab = "Bivariate CUSUM",
       xlab = "Game")
  abline(v = tau, col = "red", lty = 2, lwd = 2)
}








# 
# plt23 <- function(){
#   Y_df <- prepare_bivariate_series(df, "CLE", "GSW", 2018)
#   Y_mat <- scale(as.matrix(Y_df[, c("yA", "yB")]))
#   res <- bivariate_cusum_single(Y_mat)
# 
#   plot_bivariate_pair(Y_df, res, "CLE", "GSW")
# }

