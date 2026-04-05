source("C:/Users/joten/Documents/Jevon Oteng Lancaster/Changepoint-Detection/eda code.R")

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

detect_change_cusum3 <- function(y, c, sigma2 = 1) {
  n <- length(y)
  # Enforce Buffer: We need at least 10 games to reasonably split
  if (n < 10) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  S <- cumsum(y)
  sigma <- sqrt(sigma2)
  
  # Calculate CUSUM statistic
  C_stat <- (max(S) - min(S)) / (sigma * sqrt(n))
  
  # Find Changepoint (Argmax)
  # We restrict the search to be away from the very edges (Buffer of 5)
  # This prevents identifying "Game 1" or "Game 82" as structural breaks
  candidate_tau <- which.max(abs(S))
  
  # Logic: If the candidate is within the buffer zone, we invalidate it
  # (Optional strictness, but good for consistency with Mean-Var)
  if (candidate_tau < 5 || candidate_tau > (n - 5)) {
    # If the max deviation is at the edge, treat as no change or search inward
    # For simplicity here, we stick to the raw statistic but note the location
  }
  
  list(
    reject = (C_stat > c),
    statistic = C_stat,
    change_point = if (C_stat > c) candidate_tau else NA
  )
}







# ==============================================================================
# 1. MEAN DETECTOR (CUSUM) - With Safety Buffer
# ==============================================================================
detect_change_cusum1 <- function(y, c, sigma2 = 1) {
  n <- length(y)
  if (n < 10) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  # FIX 1: CENTER the data first (critical for proper CUSUM)
  y_centered <- y - mean(y)
  sigma <- sqrt(sigma2)
  
  # FIX 2: Use proper two-sample mean test statistic
  get_cusum_stat <- function(tau) {
    mean1 <- mean(y_centered[1:tau])
    mean2 <- mean(y_centered[(tau+1):n])
    
    # Standardized difference between segments
    stat <- abs(mean1 - mean2) * sqrt(tau * (n - tau) / n) / sigma
    return(stat)
  }
  
  # Search over valid taus (buffer of 3 games)
  taus <- 3:(n-3)
  if (length(taus) == 0) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  cusum_stats <- sapply(taus, get_cusum_stat)
  C_stat <- max(cusum_stats)
  best_tau <- taus[which.max(cusum_stats)]
  
  list(
    reject = (C_stat > c),
    statistic = C_stat,
    change_point = if (C_stat > c) best_tau else NA
  )
}

# ==============================================================================
# 2. MEAN-VARIANCE DETECTOR - STRICT BUFFER ENFORCED
# ==============================================================================
detect_change_mean_variance0 <- function(y, c, sigma2 = 1) {
  n <- length(y)
  
  # Safety: If season is too short to have a buffer, skip it
  if (n < 12) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  # 1. Precompute sums for fast calculation
  cs <- cumsum(y)
  cs2 <- cumsum(y^2)
  
  # 2. Likelihood Ratio Function
  get_lr <- function(tau) {
    # Segment 1 (1 to tau)
    mean1 <- cs[tau] / tau
    var1  <- (cs2[tau] - tau * mean1^2) / tau
    var1  <- max(var1, 1e-12) # Floor variance
    
    # Segment 2 (tau+1 to n)
    len2  <- n - tau
    mean2 <- (cs[n] - cs[tau]) / len2
    var2  <- (cs2[n] - cs2[tau] - len2 * mean2^2) / len2
    var2  <- max(var2, 1e-12)
    
    # Global (Null Hypothesis)
    mean0 <- cs[n] / n
    var0  <- (cs2[n] - n * mean0^2) / n
    var0  <- max(var0, 1e-12)
    
    # Log-Likelihood Ratio
    llr <- n * log(var0) - (tau * log(var1) + (n - tau) * log(var2))
    return(llr)
  }
  
  # 3. CRITICAL UPDATE: The Safety Buffer
  # We only calculate LLR for taus between 5 and n-5.
  # This makes it IMPOSSIBLE to detect a change at Game 2 or Game 81.
  taus <- 5:(n-5)
  
  if (length(taus) == 0) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  # Calculate statistic for valid taus only
  lr_values <- sapply(taus, get_lr)
  
  # Find max
  max_lr <- max(lr_values)
  best_index <- which.max(lr_values)
  best_tau <- taus[best_index] # Map back to actual game number
  
  list(
    reject = (max_lr > c),
    statistic = max_lr,
    change_point = if (max_lr > c) best_tau else NA
  )
}

tune_penalties_global <- function(n_sim = 2000, n_games = 82) {
  print("--- TUNING PENALTIES (Monte Carlo) ---")
  set.seed(123)
  stats_mean <- numeric(n_sim)
  stats_var  <- numeric(n_sim)
  
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  for (i in 1:n_sim) {
    # Simulate random season
    y <- rnorm(n_games, 0, 1)
    
    # Run detectors with c=0 just to capture the statistic
    stats_mean[i] <- detect_change_cusum1(y,1, c = 0)$statistic
    stats_var[i]  <- detect_change_mean_variance0(y, c = 0)$statistic
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Save thresholds (95th percentile) to Global Environment
  C_MEAN <<- quantile(stats_mean, 0.90, na.rm=TRUE)
  C_VAR  <<- quantile(stats_var, 0.90, na.rm=TRUE)
  
  print(" ")
  print(paste("Global C_MEAN set to:", round(C_MEAN, 3)))
  print(paste("Global C_VAR set to:", round(C_VAR, 3)))
}



#--------------------------COMPARING ALPHA VALUES THE GAUSSIAN AND MEAN-VARIANCE MODELS-----------------
safe_cp <- function(x) {
  if (is.null(x)) NA_integer_ else x
}

tune_penalties_global()
# ==============================================================================
# 3. ANALYZER FUNCTION (The Driver)
# ==============================================================================
analyze_cavs_season <- function(df, season) {
  
  # --- A. Data Prep ---
  y_raw <- df %>%
    dplyr::filter(slugTeam == "CLE", yearSeason == season) %>%
    dplyr::arrange(numberGameTeamSeason) %>%
    dplyr::pull(plusminusTeam)
  
  n <- length(y_raw)
  if (n < 20) return(NULL) # Skip very short seasons
  
  # --- B. Robust Scaling ---
  sigma_robust <- mad(diff(y_raw)) / sqrt(2)
  if (sigma_robust == 0) sigma_robust <- sd(y_raw)
  
  # Standardize (Variance forced to 1)
  y_std <- (y_raw - mean(y_raw)) / sigma_robust
  
  # --- C. Run Detectors ---
  if (!exists("C_MEAN") || !exists("C_VAR")) {
    stop("Global penalties (C_MEAN, C_VAR) are missing. Run tuning first.")
  }
  
  # Run both algorithms on standardized data
  res_mean <- detect_change_cusum1(y_std, c = C_MEAN, sigma2 = 1)
  res_mvar <- detect_change_mean_variance0(y_std, c = C_VAR, sigma2 = 1)
  
  # --- D. Safety Wrapper (Prevents "differing rows" error) ---
  safe_val <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA)
    return(x)
  }
  
  # --- E. Build Result ---
  data.frame(
    season = season,
    model = c("mean", "mean-variance"),
    
    reject = c(
      safe_val(res_mean$reject),
      safe_val(res_mvar$reject)
    ),
    
    changepoint = c(
      safe_val(res_mean$change_point),
      safe_val(res_mvar$change_point)
    ),
    
    statistic = c(
      safe_val(res_mean$statistic),
      safe_val(res_mvar$statistic)
    ),
    
    cutoff = c(C_MEAN, C_VAR),
    sigma_scaling_factor = c(sigma_robust, sigma_robust)
  )
}

# HOW TO TEST:
# 1. Run this whole block.
# 2. Run: analyze_cavs_season(df, 2022)
# 3. Check if changepoint is now 5 or NA (instead of 2).

# cavs_seasons <- sort(unique(df$yearSeason))
# 
# cavs_results <- cavs_seasons %>%
#   map_dfr(~ analyze_cavs_season(df, season = .x))
# View(cavs_results)
# # saveRDS(cavs_results, "cavs_results.rds")
# 

#----------------------------COMPARING DIFF METRICS ACROSS MODELS---------------------------------------

library(dplyr)
library(purrr)

run_offline_all_metrics_<- function(df_net) {
  
  keys <- df_net %>% distinct(yearSeason, slugTeam)
  
  # --- HELPER: The Stability Scan (Simulated MBIC) ---
  # This function runs your detector multiple times with strict penalties.
  # It only returns TRUE if the changepoint is stable across >50% of scans.
  scan_for_change <- function(y, model_type, penalty_range) {
    votes <- c()
    
    for (p in penalty_range) {
      if (model_type == "mean") {
        # Using YOUR function
        res <- detect_change_cusum1(y, sigma2 = 1, c = p)
        cp <- res$change_point
      } else if (model_type == "meanvar") {
        # Using YOUR function
        res <- detect_change_mean_variance0(y, c = p)
        cp <- res$change_point
      } else if (model_type == "poisson") {
        # Using YOUR function
        res <- detect_changepoint_poisson(y, c = p)
        cp <- res$change_point
      }
      
      # Collect valid changepoints (ignore NA or empty)
      if (!is.null(cp) && !is.na(cp) && length(cp) > 0) {
        votes <- c(votes, cp[1])
      }
    }
    
    # Stability Check
    if (length(votes) == 0) return(list(reject = FALSE, tau = NA))
    
    freq <- table(votes)
    best_cp <- as.numeric(names(which.max(freq)))
    count <- max(freq)
    
    # Must appear in >50% of the strict scans to be "Real"
    if (count >= (length(penalty_range) / 2)) {
      return(list(reject = TRUE, tau = best_cp))
    } else {
      return(list(reject = FALSE, tau = NA))
    }
  }
  
  # --- MAIN LOOP ---
  map_dfr(seq_len(nrow(keys)), function(i) {
    
    season <- keys$yearSeason[i]
    team   <- keys$slugTeam[i]
    
    df_ts <- df_net %>%
      filter(yearSeason == season, slugTeam == team) %>%
      arrange(numberGameTeamSeason)
    
    n <- nrow(df_ts)
    if (n < 40) return(NULL) # Skip short seasons
    
    out <- list()
    
    ## ---------- Continuous Metrics (PlusMinus, NetRating) ----------
    for (metric in c("plusminusTeam", "net_rating")) {
      
      # Robust Scaling
      val <- df_ts[[metric]]
      sigma <- mad(diff(val)) / sqrt(2)
      if (sigma == 0) sigma <- sd(val)
      y <- (val - mean(val)) / sigma
      
      # 1. MEAN MODEL (Scan 4.0 - 12.0)
      # These are strict CUSUM thresholds simulating MBIC
      res_m <- scan_for_change(y, "mean", seq(4.0, 12.0, by = 1.0))
      
      out[[length(out)+1]] <- tibble(
        yearSeason = season, team = team, metric = metric, model = "mean",
        reject = res_m$reject, tau_hat = res_m$tau
      )
      
      # 2. MEAN-VARIANCE MODEL (Scan 4.0 - 12.0)
      res_mv <- scan_for_change(y, "meanvar", seq(10.0, 30.0, by = 2.0))
      
      out[[length(out)+1]] <- tibble(
        yearSeason = season, team = team, metric = metric, model = "mean-variance",
        reject = res_mv$reject, tau_hat = res_mv$tau
      )
    }
    
    ## ---------- Discrete Metrics (Turnovers, Rebounds) ----------
    # Note: Poisson thresholds are scale-dependent, but for normalized CUSUM
    # or Likelihood Ratios, strict penalties (e.g., 4-10) often work similarly.
    # Assuming your detect_changepoint_poisson uses a Likelihood Ratio stat:
    
    for (metric in c("tovTeam", "trebTeam")) {
      y <- df_ts[[metric]]
      y <- y[!is.na(y)]
      if (length(y) < 40) next
      
      # We use a similar strict scan to avoid "Game 2" noise in turnovers
      # Range: 3.0 to 10.0 (slightly lower base than continuous but still strict)
      res_p <- scan_for_change(y, "poisson", seq(3.0, 10.0, by = 1.0))
      
      out[[length(out)+1]] <- tibble(
        yearSeason = season, team = team, metric = metric, model = "poisson",
        reject = res_p$reject, tau_hat = res_p$tau
      )
    }
    
    bind_rows(out)
  })
}

offline_results <- run_offline_all_metrics_(df_joined)
saveRDS(offline_results, "offline_results.rds")



















# ==============================================================================
# TUNING FUNCTION - Sets Global Penalties
# ==============================================================================
tune_penalties_global <- function(n_sim = 2000, n_games = 82) {
  print("--- TUNING PENALTIES (Monte Carlo) ---")
  
  stats_mean <- numeric(n_sim)
  stats_var  <- numeric(n_sim)
  
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  for (i in 1:n_sim) {
    # Simulate random season
    y <- rnorm(n_games, 0, 1)
    
    # Run detectors with c=0 just to capture the statistic
    stats_mean[i] <- detect_change_cusum1(y, c = 0)$statistic
    stats_var[i]  <- detect_change_mean_variance0(y, c = 0)$statistic
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Save thresholds (90th percentile for less conservative detection)
  C_MEAN <<- quantile(stats_mean, 0.90, na.rm=TRUE)
  C_VAR  <<- quantile(stats_var, 0.90, na.rm=TRUE)
  
  print(" ")
  print(paste("Global C_MEAN set to:", round(C_MEAN, 3)))
  print(paste("Global C_VAR set to:", round(C_VAR, 3)))
  
  # DIAGNOSTIC: Show distribution
  print(paste("CUSUM stats - Min:", round(min(stats_mean, na.rm=TRUE), 3),
              "Max:", round(max(stats_mean, na.rm=TRUE), 3),
              "Median:", round(median(stats_mean, na.rm=TRUE), 3)))
  print(paste("MeanVar stats - Min:", round(min(stats_var, na.rm=TRUE), 3),
              "Max:", round(max(stats_var, na.rm=TRUE), 3),
              "Median:", round(median(stats_var, na.rm=TRUE), 3)))
}

# ==============================================================================
# 1. MEAN DETECTOR (CUSUM) - FIXED VERSION
# ==============================================================================
detect_change_cusum1 <- function(y, c, sigma2 = 1) {
  n <- length(y)
  if (n < 10) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  # FIX 1: CENTER the data first (critical for proper CUSUM)
  y_centered <- y - mean(y)
  sigma <- sqrt(sigma2)
  
  # FIX 2: Use proper two-sample mean test statistic
  get_cusum_stat <- function(tau) {
    mean1 <- mean(y_centered[1:tau])
    mean2 <- mean(y_centered[(tau+1):n])
    
    # Standardized difference between segments
    stat <- abs(mean1 - mean2) * sqrt(tau * (n - tau) / n) / sigma
    return(stat)
  }
  
  # Search over valid taus (buffer of 3 games)
  taus <- 3:(n-3)
  if (length(taus) == 0) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  cusum_stats <- sapply(taus, get_cusum_stat)
  C_stat <- max(cusum_stats)
  best_tau <- taus[which.max(cusum_stats)]
  
  list(
    reject = (C_stat > c),
    statistic = C_stat,
    change_point = if (C_stat > c) best_tau else NA
  )
}

# ==============================================================================
# 2. MEAN-VARIANCE DETECTOR - STRICT BUFFER ENFORCED
# ==============================================================================
detect_change_mean_variance1 <- function(y, c, sigma2 = 1) {
  n <- length(y)
  
  # Safety: If season is too short to have a buffer, skip it
  if (n < 10) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  # 1. Precompute sums for fast calculation
  cs <- cumsum(y)
  cs2 <- cumsum(y^2)
  
  # 2. Likelihood Ratio Function
  get_lr <- function(tau) {
    # Segment 1 (1 to tau)
    mean1 <- cs[tau] / tau
    var1  <- (cs2[tau] - tau * mean1^2) / tau
    var1  <- max(var1, 1e-12) # Floor variance
    
    # Segment 2 (tau+1 to n)
    len2  <- n - tau
    mean2 <- (cs[n] - cs[tau]) / len2
    var2  <- (cs2[n] - cs2[tau] - len2 * mean2^2) / len2
    var2  <- max(var2, 1e-12)
    
    # Global (Null Hypothesis)
    mean0 <- cs[n] / n
    var0  <- (cs2[n] - n * mean0^2) / n
    var0  <- max(var0, 1e-12)
    
    # Log-Likelihood Ratio
    llr <- n * log(var0) - (tau * log(var1) + (n - tau) * log(var2))
    return(llr)
  }
  
  # 3. The Safety Buffer - Reduced from 5 to 3
  taus <- 3:(n-3)
  
  if (length(taus) == 0) return(list(reject = FALSE, statistic = 0, change_point = NA))
  
  # Calculate statistic for valid taus only
  lr_values <- sapply(taus, get_lr)
  
  # Find max
  max_lr <- max(lr_values)
  best_index <- which.max(lr_values)
  best_tau <- taus[best_index] # Map back to actual game number
  
  list(
    reject = (max_lr > c),
    statistic = max_lr,
    change_point = if (max_lr > c) best_tau else NA
  )
}

# ==============================================================================
# 3. ANALYZER FUNCTION WITH FULL DIAGNOSTICS
# ==============================================================================
analyze_cavs_season <- function(df, season, team = "ATL", verbose = TRUE) {
  
  # --- A. Data Prep ---
  y_raw <- df %>%
    dplyr::filter(slugTeam == team, yearSeason == season) %>%
    dplyr::arrange(numberGameTeamSeason) %>%
    dplyr::pull(plusminusTeam)
  
  n <- length(y_raw)
  
  # DIAGNOSTIC OUTPUT
  if (verbose) {
    cat(sprintf("\n========== Season %d (%s) ==========\n", season, team))
    cat(sprintf("Games: %d\n", n))
  }
  
  if (n < 20) {
    if (verbose) cat("  -> SKIPPED (too short)\n")
    return(NULL)
  }
  
  # --- B. Robust Scaling ---
  sigma_robust <- mad(diff(y_raw)) / sqrt(2)
  if (sigma_robust == 0) sigma_robust <- sd(y_raw)
  
  if (verbose) {
    cat(sprintf("Raw data - Mean: %.2f, SD: %.2f, Sigma_robust: %.3f\n", 
                mean(y_raw), sd(y_raw), sigma_robust))
  }
  
  # FIX 3: RE-CENTER the data (this is crucial!)
  y_std <- (y_raw - mean(y_raw)) / sigma_robust
  
  if (verbose) {
    cat(sprintf("Standardized - Mean: %.3f, SD: %.3f\n", 
                mean(y_std), sd(y_std)))
  }
  
  # --- C. Run Detectors ---
  if (!exists("C_MEAN") || !exists("C_VAR")) {
    stop("Global penalties (C_MEAN, C_VAR) are missing. Run tune_penalties_global() first.")
  }
  
  # Run both algorithms on standardized data
  res_mean <- detect_change_cusum1(y_std, c = C_MEAN, sigma2 = 1)
  res_mvar <- detect_change_mean_variance1(y_std, c = C_VAR, sigma2 = 1)
  
  # DIAGNOSTIC OUTPUT
  if (verbose) {
    cat(sprintf("CUSUM:\n"))
    cat(sprintf("  Statistic: %.3f | Threshold: %.3f | Reject: %s | CP: %s\n",
                res_mean$statistic, C_MEAN, res_mean$reject,
                ifelse(is.na(res_mean$change_point), "None", res_mean$change_point)))
    
    cat(sprintf("Mean-Variance:\n"))
    cat(sprintf("  Statistic: %.3f | Threshold: %.3f | Reject: %s | CP: %s\n",
                res_mvar$statistic, C_VAR, res_mvar$reject,
                ifelse(is.na(res_mvar$change_point), "None", res_mvar$change_point)))
  }
  
  # --- D. Safety Wrapper ---
  safe_val <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA)
    return(x)
  }
  
  # --- E. Build Result ---
  data.frame(
    season = season,
    team = team,
    n_games = n,
    model = c("mean", "mean-variance"),
    
    reject = c(
      safe_val(res_mean$reject),
      safe_val(res_mvar$reject)
    ),
    
    changepoint = c(
      safe_val(res_mean$change_point),
      safe_val(res_mvar$change_point)
    ),
    
    statistic = c(
      safe_val(res_mean$statistic),
      safe_val(res_mvar$statistic)
    ),
    
    cutoff = c(C_MEAN, C_VAR),
    sigma_scaling_factor = c(sigma_robust, sigma_robust)
  )
}

# ==============================================================================
# 4. SYNTHETIC DATA VALIDATION TEST
# ==============================================================================
test_synthetic_changepoint <- function() {
  cat("\n========== TESTING ON SYNTHETIC DATA ==========\n")
  
  # Create fake season with clear break at game 41
  set.seed(123)
  y_fake <- c(rnorm(41, mean = 5, sd = 10),   # Good first half
              rnorm(41, mean = -5, sd = 10))  # Bad second half
  
  # Standardize
  y_fake_std <- (y_fake - mean(y_fake)) / sd(y_fake)
  
  # Test both methods
  cat("\nTrue changepoint: Game 41\n")
  
  res_cusum <- detect_change_cusum1(y_fake_std, c = C_MEAN, sigma2 = 1)
  cat(sprintf("CUSUM detected: Game %s (statistic: %.3f)\n", 
              ifelse(is.na(res_cusum$change_point), "None", res_cusum$change_point),
              res_cusum$statistic))
  
  res_mvar <- detect_change_mean_variance1(y_fake_std, c = C_VAR, sigma2 = 1)
  cat(sprintf("Mean-Var detected: Game %s (statistic: %.3f)\n",
              ifelse(is.na(res_mvar$change_point), "None", res_mvar$change_point),
              res_mvar$statistic))
  
  cat("\nIf both methods detect around Game 41, the code is working correctly!\n")
}

# ==============================================================================
# 5. RUN FULL ANALYSIS ACROSS MULTIPLE SEASONS
# ==============================================================================
run_full_analysis <- function(df, team = "ATL", seasons = 2010:2025, verbose = TRUE, test_synthetic = TRUE) {
  
  # First, tune the penalties
  tune_penalties_global(n_sim = 2000, n_games = 82)
  
  # Optional: Test on synthetic data
  if (test_synthetic) {
    test_synthetic_changepoint()
  }
  
  # Then analyze each season
  results_list <- lapply(seasons, function(s) {
    analyze_cavs_season(df, season = s, team = team, verbose = verbose)
  })
  
  # Combine all results
  results <- dplyr::bind_rows(results_list)
  
  # Summary
  if (verbose) {
    cat("\n\n========== SUMMARY ==========\n")
    cat(sprintf("Total seasons analyzed: %d\n", length(unique(results$season))))
    cat(sprintf("CUSUM detections: %d\n", sum(results$reject[results$model == "mean"], na.rm = TRUE)))
    cat(sprintf("Mean-Var detections: %d\n", sum(results$reject[results$model == "mean-variance"], na.rm = TRUE)))
    
    # Show distribution of changepoints
    cp_cusum <- results$changepoint[results$model == "mean" & !is.na(results$changepoint)]
    if (length(cp_cusum) > 0) {
      cat(sprintf("\nCUSUM changepoint distribution:\n"))
      cat(sprintf("  Min: %d, Max: %d, Median: %.0f\n", 
                  min(cp_cusum), max(cp_cusum), median(cp_cusum)))
    }
  }
  
  return(results)
}

# ==============================================================================
# USAGE EXAMPLES:
# ==============================================================================
# Run with full diagnostics and synthetic test
# results <- run_full_analysis(df, team = "CLE", seasons = 2010:2025, verbose = TRUE, test_synthetic = TRUE)
# 
# # Run quietly without synthetic test
# results <- run_full_analysis(df, team = "ATL", seasons = 2010:2025, verbose = FALSE, test_synthetic = FALSE)
# 
# # View results
# View(results)
# 
# # Filter to see only detections
# detections <- results %>% filter(reject == TRUE)
# View(detections)
# ==============================================================================
# 4. RUN ANALYSIS ACROSS MULTIPLE SEASONS
# ==============================================================================
run_full_analysis <- function(df, team = "CLE", seasons = 2010:2025, verbose = TRUE) {
  
  # First, tune the penalties
  tune_penalties_global(n_sim = 2000, n_games = 82)
  
  # Then analyze each season
  results_list <- lapply(seasons, function(s) {
    analyze_cavs_season(df, season = s, team = team, verbose = verbose)
  })
  
  # Combine all results
  results <- dplyr::bind_rows(results_list)
  
  # Summary
  if (verbose) {
    cat("\n\n========== SUMMARY ==========\n")
    cat(sprintf("Total seasons analyzed: %d\n", length(unique(results$season))))
    cat(sprintf("CUSUM detections: %d\n", sum(results$reject[results$model == "mean"], na.rm = TRUE)))
    cat(sprintf("Mean-Var detections: %d\n", sum(results$reject[results$model == "mean-variance"], na.rm = TRUE)))
  }
  
  return(results)
}

# ==============================================================================
# USAGE EXAMPLE:
# ==============================================================================
# results <- run_full_analysis(df, team = "ATL", seasons = 2010:2025, verbose = TRUE)
# 
# # To analyze quietly:
# results <- run_full_analysis(df, team = "CLE", seasons = 2010:2025, verbose = TRUE)
# 
# # View results
# View(results)
# 
# 



# ==============================================================================
# STEP 3: DETERMINE BEST METHOD BASED ON VALIDATION
# ==============================================================================
evaluate_methods <- function(results_df) {
  
  cat("\n\n========== METHOD COMPARISON ACROSS TOP TEAMS ==========\n")
  
  # Count detections by method
  cusum_detections <- results_df %>%
    dplyr::filter(model == "mean", reject == TRUE) %>%
    dplyr::group_by(team) %>%
    dplyr::summarise(cusum_count = n())
  
  meanvar_detections <- results_df %>%
    dplyr::filter(model == "mean-variance", reject == TRUE) %>%
    dplyr::group_by(team) %>%
    dplyr::summarise(meanvar_count = n())
  
  comparison <- dplyr::full_join(cusum_detections, meanvar_detections, by = "team") %>%
    dplyr::mutate(
      cusum_count = ifelse(is.na(cusum_count), 0, cusum_count),
      meanvar_count = ifelse(is.na(meanvar_count), 0, meanvar_count)
    )
  
  print(comparison)
  
  # Recommendation
  total_cusum <- sum(comparison$cusum_count)
  total_meanvar <- sum(comparison$meanvar_count)
  
  cat(sprintf("\nTotal CUSUM detections: %d\n", total_cusum))
  cat(sprintf("Total Mean-Variance detections: %d\n", total_meanvar))
  
  cat("\nRECOMMENDATION:\n")
  cat("Use CUSUM for initial screening (higher sensitivity to mean shifts).\n")
  cat("Use Mean-Variance for confirmation (detects fundamental regime changes).\n")
  cat("Prioritize investigating changepoints confirmed by BOTH methods.\n")
  
  return(comparison)
}

# ==============================================================================
# STEP 4: VISUALIZATION FUNCTION
# ==============================================================================
plot_top_teams_changepoints <- function(df, results_df, seasons = 2010:2025) {
  
  library(ggplot2)
  
  top_teams <- unique(results_df$team)
  
  for (team in top_teams) {
    # Get all plus-minus data
    team_data <- df %>%
      dplyr::filter(slugTeam == team, yearSeason %in% seasons) %>%
      dplyr::arrange(yearSeason, numberGameTeamSeason) %>%
      dplyr::mutate(game_index = row_number())
    
    # Get detected changepoints
    changepoints_cusum <- results_df %>%
      dplyr::filter(team == !!team, model == "mean", reject == TRUE) %>%
      dplyr::left_join(
        df %>%
          dplyr::filter(slugTeam == !!team) %>%
          dplyr::group_by(yearSeason) %>%
          dplyr::summarise(season_start = min(numberGameTeamSeason)),
        by = c("season" = "yearSeason")
      ) %>%
      dplyr::mutate(global_cp = cumsum(c(0, diff(season) != 0)) * 82 + changepoint)
    
    # Plot
    p <- ggplot(team_data, aes(x = game_index, y = plusminusTeam)) +
      geom_line(alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_vline(data = changepoints_cusum, 
                 aes(xintercept = global_cp), 
                 color = "red", linewidth = 1, alpha = 0.7) +
      labs(
        title = sprintf("%s: Plus-Minus Over Time (2010-2025)", team),
        subtitle = "Red lines indicate CUSUM-detected changepoints",
        x = "Game Number (Across All Seasons)",
        y = "Plus-Minus Score"
      ) +
      theme_minimal()
    
    print(p)
  }
}





# ==============================================================================
# DETAILED INVESTIGATION OF TOP 3 TEAMS
# ==============================================================================
investigate_top_teams_detailed <- function(df, seasons = 2010:2025) {
  
  top_teams <- c("HOU", "POR", "MIA")
  
  detailed_results <- list()
  
  for (team in top_teams) {
    cat(sprintf("\n\n########################################\n"))
    cat(sprintf("### DETAILED ANALYSIS: %s ###\n", team))
    cat(sprintf("########################################\n"))
    
    team_changes <- data.frame()
    
    for (season in seasons) {
      # Get data
      y_raw <- df %>%
        dplyr::filter(slugTeam == team, yearSeason == season) %>%
        dplyr::arrange(numberGameTeamSeason) %>%
        dplyr::pull(plusminusTeam)
      
      n <- length(y_raw)
      if (n < 20) next
      
      # Standardize
      sigma_robust <- mad(diff(y_raw)) / sqrt(2)
      if (sigma_robust == 0) sigma_robust <- sd(y_raw)
      y_std <- (y_raw - mean(y_raw)) / sigma_robust
      
      # Run both methods
      res_cusum <- detect_change_cusum1(y_std, c = C_MEAN, sigma2 = 1)
      res_mvar <- detect_change_mean_variance1(y_std, c = C_VAR, sigma2 = 1)
      
      # If either detected, analyze the segments
      if (res_cusum$reject || res_mvar$reject) {
        
        cp_cusum <- res_cusum$change_point
        cp_mvar <- res_mvar$change_point
        
        cat(sprintf("\n=== SEASON %d ===\n", season))
        cat(sprintf("CUSUM: %s (stat=%.3f) | Mean-Var: %s (stat=%.3f)\n",
                    ifelse(is.na(cp_cusum), "No change", paste("Game", cp_cusum)),
                    res_cusum$statistic,
                    ifelse(is.na(cp_mvar), "No change", paste("Game", cp_mvar)),
                    res_mvar$statistic))
        
        # Use CUSUM changepoint if available, otherwise mean-variance
        cp <- ifelse(!is.na(cp_cusum), cp_cusum, cp_mvar)
        
        if (!is.na(cp)) {
          # Calculate before/after statistics
          before <- y_raw[1:cp]
          after <- y_raw[(cp+1):n]
          
          cat(sprintf("\nBEFORE Changepoint (Games 1-%d):\n", cp))
          cat(sprintf("  Mean: %.2f | SD: %.2f | N: %d\n", 
                      mean(before), sd(before), length(before)))
          
          cat(sprintf("AFTER Changepoint (Games %d-%d):\n", cp+1, n))
          cat(sprintf("  Mean: %.2f | SD: %.2f | N: %d\n", 
                      mean(after), sd(after), length(after)))
          
          cat(sprintf("CHANGE: Mean Δ = %.2f | SD Δ = %.2f\n", 
                      mean(after) - mean(before),
                      sd(after) - sd(before)))
          
          # Store result
          team_changes <- rbind(team_changes, data.frame(
            team = team,
            season = season,
            changepoint_game = cp,
            cusum_detected = !is.na(cp_cusum),
            meanvar_detected = !is.na(cp_mvar),
            both_confirmed = !is.na(cp_cusum) & !is.na(cp_mvar),
            mean_before = mean(before),
            mean_after = mean(after),
            sd_before = sd(before),
            sd_after = sd(after),
            mean_delta = mean(after) - mean(before),
            sd_delta = sd(after) - sd(before)
          ))
        }
      }
    }
    
    detailed_results[[team]] <- team_changes
    
    # Summary
    cat(sprintf("\n\n--- %s SUMMARY ---\n", team))
    if (nrow(team_changes) > 0) {
      confirmed <- sum(team_changes$both_confirmed)
      cat(sprintf("Total changepoints: %d\n", nrow(team_changes)))
      cat(sprintf("Confirmed by both methods: %d\n", confirmed))
      cat(sprintf("Average mean delta: %.2f\n", mean(team_changes$mean_delta)))
      cat(sprintf("Average SD delta: %.2f\n", mean(team_changes$sd_delta)))
    }
  }
  
  return(dplyr::bind_rows(detailed_results))
}

# ==============================================================================
# VISUALIZATION: BEFORE/AFTER COMPARISON
# ==============================================================================
plot_changepoint_effects <- function(detailed_results) {
  library(ggplot2)
  
  # Plot 1: Mean changes
  p1 <- ggplot(detailed_results, aes(x = factor(season), y = mean_delta, fill = both_confirmed)) +
    geom_col() +
    facet_wrap(~team, scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "orange"),
                      labels = c("TRUE" = "Both methods", "FALSE" = "One method")) +
    labs(title = "Mean Performance Change at Detected Changepoints",
         subtitle = "Negative = Performance declined after changepoint",
         x = "Season",
         y = "Change in Mean Plus-Minus",
         fill = "Confirmation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  
  # Plot 2: Variance changes
  p2 <- ggplot(detailed_results, aes(x = factor(season), y = sd_delta, fill = both_confirmed)) +
    geom_col() +
    facet_wrap(~team, scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "orange"),
                      labels = c("TRUE" = "Both methods", "FALSE" = "One method")) +
    labs(title = "Variance Change at Detected Changepoints",
         subtitle = "Negative = More consistent after changepoint",
         x = "Season",
         y = "Change in SD",
         fill = "Confirmation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p2)
}

# ==============================================================================
# RUN DETAILED INVESTIGATION
# ==============================================================================
# detailed_analysis <- investigate_top_teams_detailed(df, seasons = 2010:2025)
# 
# # View detailed results
# View(detailed_analysis)
# 
# # Create visualizations
# plot_changepoint_effects(detailed_analysis)

# # Export for further analysis
# write.csv(detailed_analysis, "top3_teams_changepoint_details.csv", row.names = FALSE)
# 
# 


library(dplyr)
library(purrr)

# ==============================================================================
# 1. HELPER: SCANNER FOR MAX CUSUM (Ranking Logic)
# ==============================================================================
# This function calculates the max CUSUM statistic for a single team across all seasons
# without printing verbose output. It is used to rank the teams.

get_team_max_cusum <- function(df, team_slug, seasons = 2010:2025) {
  
  max_stat <- 0
  best_season <- NA
  
  for (s in seasons) {
    # 1. Data Prep (Reusing the robust logic from your main function)
    y_raw <- df %>%
      dplyr::filter(slugTeam == team_slug, yearSeason == s) %>%
      dplyr::arrange(numberGameTeamSeason) %>%
      dplyr::pull(plusminusTeam)
    
    n <- length(y_raw)
    if (n < 20) next # Skip short seasons
    
    # 2. Robust Scaling & Centering (Critical for valid comparison)
    sigma_robust <- mad(diff(y_raw)) / sqrt(2)
    if (sigma_robust == 0) sigma_robust <- sd(y_raw)
    y_std <- (y_raw - mean(y_raw)) / sigma_robust
    
    # 3. Get CUSUM Statistic (We don't need the cutoff here, just the raw stat)
    # We pass c=0 just to get the statistic returned
    res <- detect_change_cusum1(y_std, c = 0, sigma2 = 1)
    
    if (!is.na(res$statistic) && res$statistic > max_stat) {
      max_stat <- res$statistic
      best_season <- s
    }
  }
  
  return(data.frame(team = team_slug, max_cusum = max_stat, peak_season = best_season))
}

# ==============================================================================
# 2. RANKING FUNCTION
# ==============================================================================
rank_teams_by_max_cusum <- function(df, seasons = 2010:2025) {
  cat("\n--- SCANNING LEAGUE FOR MAX STRUCTURAL BREAKS ---\n")
  
  all_teams <- unique(df$slugTeam)
  # Remove "null" teams if any exist in data
  all_teams <- all_teams[!is.na(all_teams)]
  
  # Progress bar since this scans 30 teams * 15 seasons
  pb <- txtProgressBar(min = 0, max = length(all_teams), style = 3)
  
  rankings <- map_dfr(seq_along(all_teams), function(i) {
    t <- all_teams[i]
    res <- get_team_max_cusum(df, team_slug = t, seasons = seasons)
    setTxtProgressBar(pb, i)
    return(res)
  })
  close(pb)
  
  # Sort by highest CUSUM statistic descending
  ranked <- rankings %>%
    dplyr::arrange(desc(max_cusum))
  
  return(ranked)
}

# ==============================================================================
# 3. ANALYSIS FUNCTION (Modified to use Rankings)
# ==============================================================================
analyze_top_teams <- function(df, rankings_df, n_top = 3, seasons = 2010:2025) {
  
  # Select the top N teams from the ranking dataframe
  top_n_teams <- head(rankings_df$team, n_top)
  
  cat(sprintf("\n\n========== ANALYZING TOP %d TEAMS WITH HIGHEST MAX CUSUM ==========\n", n_top))
  print(head(rankings_df, n_top))
  
  all_results <- list()
  
  for (team in top_n_teams) {
    cat(sprintf("\n\n##################################################\n"))
    cat(sprintf("########## %s (Ranked via Max CUSUM) ##########\n", team))
    cat(sprintf("##################################################\n"))
    
    # Run full analysis with both methods
    # NOTE: ensure analyze_cavs_season is defined in your environment!
    team_results <- lapply(seasons, function(s) {
      # Use tryCatch to prevent loop stopping on data errors
      tryCatch({
        analyze_cavs_season(df, season = s, team = team, verbose = TRUE)
      }, error = function(e) return(NULL))
    })
    
    # Combine results
    team_df <- dplyr::bind_rows(team_results)
    
    if (nrow(team_df) > 0) {
      all_results[[team]] <- team_df
      
      # --- Summary ---
      cat(sprintf("\n--- %s SUMMARY ---\n", team))
      cat(sprintf("Total seasons analyzed: %d\n", length(unique(team_df$season))))
      
      # Count detections
      n_mean <- sum(team_df$reject[team_df$model == "mean"], na.rm = TRUE)
      n_mvar <- sum(team_df$reject[team_df$model == "mean-variance"], na.rm = TRUE)
      
      cat(sprintf("CUSUM detections: %d\n", n_mean))
      cat(sprintf("Mean-Var detections: %d\n", n_mvar))
      
      # Find Concordance (Where BOTH models detected a change)
      confirmed <- team_df %>%
        dplyr::filter(reject == TRUE) %>%
        dplyr::group_by(season) %>%
        dplyr::filter(n_distinct(model) == 2) %>% # Must have both models present
        dplyr::summarise(
          mean_cp = changepoint[model=="mean"],
          mvar_cp = changepoint[model=="mean-variance"],
          mean_stat = statistic[model=="mean"],
          mvar_stat = statistic[model=="mean-variance"]
        )
      
      if (nrow(confirmed) > 0) {
        cat("\nHIGH-CONFIDENCE SEASONS (Both models triggered):\n")
        print(confirmed)
      } else {
        cat("\nNo seasons confirmed by both methods simultaneously.\n")
      }
    } else {
      cat("\nNo valid data found for this team.\n")
    }
  }
  
  return(dplyr::bind_rows(all_results, .id = "team"))
}

# ==============================================================================
# 4. MASTER WORKFLOW (Run This)
# ==============================================================================
run_league_analysis_workflow <- function(df) {
  
  # Step 1: Tune Penalties (Critical first step)
  # Uses the calibration function you already have
  tune_penalties_global(n_sim = 2000, n_games = 82)
  
  # Step 2: Rank Teams
  rankings <- rank_teams_by_max_cusum(df)
  
  # Step 3: Analyze the Top 3
  final_results <- analyze_top_teams(df, rankings, n_top = 3)
  
  return(list(rankings = rankings, results = final_results))
}

# Assuming 'df' is your main NBA dataframe
output <- run_league_analysis_workflow(df)

# # View the full rankings to see who missed the cut
# print(output$rankings)
# 
# # View the detailed analysis of the winners
# print(output$results)
# 
# 
# 
# 












library(ggplot2)
library(dplyr)
library(ggrepel)

# ==============================================================================
# 1. HELPER: FIND GLOBAL CHANGEPOINT
# ==============================================================================
find_global_changepoint <- function(y) {
  # Simple CUSUM on the full stitched dataset
  n <- length(y)
  if (n < 50) return(NA)
  
  # Center the data globally to find deviation from the 15-year average
  y_centered <- y - mean(y)
  
  # CUSUM Statistic Vector
  # |Mean(1:k) - Mean(k+1:n)| * weighting
  get_stat <- function(k) {
    if (k < 50 || k > (n-50)) return(0) # Buffer for 15-year data
    mean1 <- mean(y_centered[1:k])
    mean2 <- mean(y_centered[(k+1):n])
    abs(mean1 - mean2) * sqrt(k * (n - k) / n)
  }
  
  stats <- sapply(1:(n-1), get_stat)
  best_k <- which.max(stats)
  
  return(best_k)
}

# ==============================================================================
# 2. MAIN PLOTTING FUNCTION
# ==============================================================================
plot_franchise_trajectory <- function(df, team_slug) {
  
  # --- A. Data Prep ---
  # Get all games for the team chronologically
  long_data <- df %>%
    dplyr::filter(slugTeam == team_slug) %>%
    dplyr::arrange(yearSeason, numberGameTeamSeason) %>%
    dplyr::mutate(
      GameIndex = row_number(),
      GlobalMean = mean(plusminusTeam),
      # The CUSUM Path: Cumulative sum of how much they beat/missed the 15-year average
      CUSUM = cumsum(plusminusTeam - GlobalMean)
    )
  
  if (nrow(long_data) == 0) return(NULL)
  
  # --- B. Find the Structural Break ---
  y_vec <- long_data$plusminusTeam
  cp_index <- find_global_changepoint(y_vec)
  
  # Get metadata for the label (Season and Game Number)
  cp_info <- long_data %>% dplyr::filter(GameIndex == cp_index)
  label_text <- paste0("Global Break\n", cp_info$yearSeason, " (Game ", cp_info$numberGameTeamSeason, ")")
  
  # --- C. Find Season Boundaries (for dotted grid lines) ---
  season_starts <- long_data %>%
    dplyr::group_by(yearSeason) %>%
    dplyr::summarise(StartIndex = min(GameIndex)) %>%
    dplyr::pull(StartIndex)
  
  # --- D. Plotting ---
  p <- ggplot(long_data, aes(x = GameIndex, y = CUSUM)) +
    
    # 1. The Trajectory Line
    geom_line(color = "#2c3e50", size = 1) +
    
    # 2. Season Grid Lines (Light Grey)
    geom_vline(xintercept = season_starts, color = "grey85", linetype = "dotted") +
    
    # 3. The Structural Break Line (Red)
    geom_vline(xintercept = cp_index, color = "#e74c3c", size = 1.2, linetype = "longdash") +
    
    # 4. The Label (Smartly placed)
    annotate("label", x = cp_index, y = max(long_data$CUSUM) * 0.9, 
             label = label_text, fill = "#e74c3c", color = "white", fontface = "bold", alpha = 0.9) +
    
    # 5. Aesthetics
    theme_minimal() +
    labs(
      title = paste0(team_slug, ": 15-Year Performance Trajectory (2010–2025)"),
      subtitle = "Cumulative deviations from the 15-year average. Kinks indicate structural breaks.",
      x = "Cumulative Game Number",
      y = "Cumulative Performance (CUSUM)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold")
    )
  
  print(p)
}

# ==============================================================================
# 3. EXECUTE FOR TARGET TEAMS
# # ==============================================================================
# 
# # HOUSTON: Expect a massive 'V' shape or sharp turn around the Harden trade
# plot_franchise_trajectory(df, "HOU")
# 
# # MIAMI: Expect the 'Big 3' volatility and the post-LeBron adjustment
# plot_franchise_trajectory(df, "MIA")
# 
# # PORTLAND: Expect the steady Lillard climb
# plot_franchise_trajectory(df, "POR")
# 
# 










library(dplyr)
library(ggplot2)
library(gridExtra)

# ==============================================================================
# 1. THE ONLINE CUSUM ALGORITHM (Page's CUSUM)
# ==============================================================================
run_online_cusum2 <- function(y, k = 0.5, h = 5) {
  n <- length(y)
  
  # We track two statistics simultaneously:
  # S_pos: Accumulates evidence of IMPROVEMENT (Positive Shift)
  # S_neg: Accumulates evidence of DECLINE (Negative Shift)
  s_pos <- numeric(n)
  s_neg <- numeric(n)
  alarms <- integer(0)
  
  for (t in 1:n) {
    # Get previous values (start at 0 for Game 1)
    prev_pos <- if (t == 1) 0 else s_pos[t-1]
    prev_neg <- if (t == 1) 0 else s_neg[t-1]
    
    # Update Page's CUSUM Statistic
    # Formula: Max(0, Previous + (Observation - Target) - Slack)
    # We assume standardized data, so Target Mean = 0
    s_pos[t] <- max(0, prev_pos + y[t] - k)
    s_neg[t] <- max(0, prev_neg - y[t] - k) # Note: we subtract y[t] for negative detection
    
    # Trigger Alarm if threshold is breached
    if (s_pos[t] > h || s_neg[t] > h) {
      alarms <- c(alarms, t)
    }
  }
  
  return(list(pos = s_pos, neg = s_neg, alarms = alarms))
}

# ==============================================================================
# 2. CALIBRATION (Finding 'h' for ARL = 82)
# ==============================================================================
# We tune 'h' so a stable team only triggers a false alarm once per season on average.
calibrate_online_threshold <- function(target_arl = 82, k = 0.5, n_sim = 2000) {
  message(paste("--- Calibrating Online Threshold (h) for ARL =", target_arl, "---"))
  
  h_candidates <- seq(3, 8, by = 0.1) # Search range
  best_h <- NA
  best_diff <- Inf
  
  pb <- txtProgressBar(min = 0, max = length(h_candidates), style = 3)
  
  for (idx in seq_along(h_candidates)) {
    h <- h_candidates[idx]
    run_lengths <- numeric(n_sim)
    
    for (i in 1:n_sim) {
      # Simulate stable random season (Null Hypothesis)
      y_sim <- rnorm(target_arl * 4, 0, 1) 
      res <- run_online_cusum2(y_sim, k = k, h = h)
      
      if (length(res$alarms) > 0) {
        run_lengths[i] <- res$alarms[1]
      } else {
        run_lengths[i] <- target_arl * 4 # Censored
      }
    }
    
    # Compare simulated ARL to target (82)
    avg_run_length <- mean(run_lengths)
    diff <- abs(avg_run_length - target_arl)
    
    if (diff < best_diff) {
      best_diff <- diff
      best_h <- h
    }
    setTxtProgressBar(pb, idx)
  }
  close(pb)
  
  message(paste("\nOptimal Threshold h:", best_h))
  return(best_h)
}

# ==============================================================================
# 3. REPLAY EXPERIMENT (The 4 Case Studies)
# ==============================================================================
replay_season_online <- function(df, team_slug, season, offline_cp, h_threshold) {
  
  # A. Data Prep
  y_raw <- df %>%
    dplyr::filter(slugTeam == team_slug, yearSeason == season) %>%
    dplyr::arrange(numberGameTeamSeason) %>%
    dplyr::pull(plusminusTeam)
  
  if(length(y_raw) < 20) return(NULL)
  
  # B. Robust Scaling (Simulating known variance)
  sigma_robust <- mad(diff(y_raw)) / sqrt(2)
  if(sigma_robust == 0) sigma_robust <- sd(y_raw)
  y_std <- (y_raw - mean(y_raw)) / sigma_robust
  
  # C. Run Online Detector
  res <- run_online_cusum2(y_std, k = 0.5, h = h_threshold)
  
  # D. Analyze Results
  first_alarm <- if(length(res$alarms) > 0) res$alarms[1] else NA
  delay <- if(!is.na(first_alarm)) (first_alarm - offline_cp) else NA
  
  # E. Create Plot
  plot_data <- data.frame(
    Game = 1:length(y_raw),
    Improvement = res$pos,
    Decline = res$neg
  )
  
  p <- ggplot(plot_data, aes(x = Game)) +
    # Threshold Line
    geom_hline(yintercept = h_threshold, color = "grey50", linetype = "dashed") +
    annotate("text", x = 5, y = h_threshold + 0.2, label = paste("Threshold h =", h_threshold), color = "grey50") +
    
    # Trends
    geom_line(aes(y = Improvement, color = "Improvement"), size = 1) +
    geom_line(aes(y = Decline, color = "Decline"), size = 1) +
    
    # Events
    geom_vline(xintercept = offline_cp, color = "blue", linetype = "dotted", size = 1) +
    # Alarm Point
    geom_point(data = data.frame(x = first_alarm, y = h_threshold), 
               aes(x = x, y = y), color = "red", size = 4, shape = 18, na.rm = TRUE) +
    
    scale_color_manual(values = c("Improvement" = "#2ecc71", "Decline" = "#e74c3c")) +
    labs(title = paste(team_slug, season),
         subtitle = paste("Detection Delay:", ifelse(is.na(delay), "NO ALARM", paste(delay, "Games"))),
         y = "Online CUSUM Score", x = "Game") +
    theme_minimal() +
    theme(legend.position = "top", legend.title = element_blank())
  
  return(list(team = team_slug, season = season, alarm = first_alarm, delay = delay, plot = p))
}

# ==============================================================================
# 4. EXECUTION BLOCK
# ==============================================================================

# Step 1: Calibrate (Run once to set global threshold)
# Note: For speed, you can manually set ONLINE_H <- 4.8 if calibration takes too long
ONLINE_H <- calibrate_online_threshold(target_arl = 82)

# Step 2: Run the 4 Case Studies
# 1. Cleveland 2022 (Game 33 - Variance Shift?)
res_cle <- replay_season_online(df, "CLE", 2022, offline_cp = 33, h_threshold = ONLINE_H)

# 2. Houston 2020 (Game 59 - Collapse)
res_hou <- replay_season_online(df, "HOU", 2020, offline_cp = 59, h_threshold = ONLINE_H)

# 3. Miami 2014 (Game 76 - Fatigue)
res_mia <- replay_season_online(df, "MIA", 2014, offline_cp = 76, h_threshold = ONLINE_H)

# 4. Portland 2022 (Game 18 - Injury)
res_por <- replay_season_online(df, "POR", 2022, offline_cp = 18, h_threshold = ONLINE_H)

# Step 3: Print Summary
summary_table <- data.frame(
  Team = c("CLE", "HOU", "MIA", "POR"),
  Season = c(2022, 2020, 2014, 2022),
  Type = c("Variance (Defensive)", "Variance/Mean (Collapse)", "Mean (Fatigue)", "Mean (Injury)"),
  Offline_CP = c(33, 59, 76, 18),
  Online_Alarm = c(res_cle$alarm, res_hou$alarm, res_mia$alarm, res_por$alarm),
  Delay = c(res_cle$delay, res_hou$delay, res_mia$delay, res_por$delay)
)

f <- print(summary_table)

# Step 4: Show Plots
x <- grid.arrange(res_cle$plot, res_hou$plot, res_mia$plot, res_por$plot, ncol=2)



###--------------------------------------------------------------------------------
# =============================================================================
# ONLINE BIVARIATE GLR Changepoint Detection of NBA Player Transitions
#
# This is a direct bivariate extension of the univariate online GLR detector.
# The univariate online GLR at time t is:
#
#   G_t = max_{k} [k(t-k)/t] * (ybar_1:k - ybar_{k+1:t})^2 / sigma^2
#
# The bivariate extension replaces the scalar squared difference with
# the bivariate quadratic form:
#
#   G_t = max_{k} [k(t-k)/t] * (ybar_1:k - ybar_{k+1:t})^T Sigma^{-1} (ybar_1:k - ybar_{k+1:t})
#
# Same code structure, same cumulative sums, just two dimensions.
# =============================================================================

library(dplyr)

# =============================================================================
# SECTION 1: DATA PREPARATION
# =============================================================================

prepare_bivariate_series <- function(df, team_a, team_b, season) {
  
  extract_team <- function(slug) {
    df %>%
      filter(slugTeam == slug, yearSeason == season) %>%
      arrange(numberGameTeamSeason) %>%
      select(game_num = numberGameTeamSeason, dateGame, pm = plusminusTeam)
  }
  
  series_a <- extract_team(team_a)
  series_b <- extract_team(team_b)
  
  merged <- inner_join(series_a, series_b, by = "game_num", suffix = c("_a", "_b"))
  
  data.frame(
    game_index = merged$game_num,
    dateGame_a = merged$dateGame_a,
    dateGame_b = merged$dateGame_b,
    pm_a       = merged$pm_a,
    pm_b       = merged$pm_b
  )
}


# =============================================================================
# SECTION 2: ONLINE BIVARIATE GLR
# =============================================================================

# ----------------------------------------------------------------------
# 2a. Bivariate GLR statistic at a single time t
#     Direct extension of your glr_at_t() to two dimensions.
#
#     Instead of:
#       C_t2 = [k(t-k)/t] * (ybar1 - ybar2)^2
#       lr   = C_t2 / sigma2
#
#     We compute:
#       delta  = ybar1 - ybar2          (2-vector)
#       C_t2   = [k(t-k)/t] * delta^T Sigma_inv delta
#       (no division by sigma2 needed — already in the quadratic form)
# ----------------------------------------------------------------------

bivariate_glr_at_t <- function(cumS_A, cumS_B, t, sigma_A, sigma_B, k_min = 5) {
  # cumS_A, cumS_B: cumulative sums for each team up to current t
  # sigma_A, sigma_B: robust SD estimates for each component
  # k_min: minimum segment length (same buffer as offline)
  
  if (t < 2 * k_min) return(list(stat = 0, best_k = NA))
  
  S_A_t <- cumS_A[t]
  S_B_t <- cumS_B[t]
  
  max_lr <- 0
  best_k <- NA
  
  for (k in k_min:(t - k_min)) {
    # Segment means for team A
    ybar1_A <- cumS_A[k] / k
    ybar2_A <- (S_A_t - cumS_A[k]) / (t - k)
    
    # Segment means for team B
    ybar1_B <- cumS_B[k] / k
    ybar2_B <- (S_B_t - cumS_B[k]) / (t - k)
    
    # Standardised differences (using identity covariance after standardising)
    delta_A <- (ybar1_A - ybar2_A) / sigma_A
    delta_B <- (ybar1_B - ybar2_B) / sigma_B
    
    # Bivariate quadratic form with Sigma = I (after standardisation)
    # This is the direct analogue of C_t2 / sigma2 in your univariate code
    weight <- (k * (t - k)) / t
    lr     <- weight * (delta_A^2 + delta_B^2)
    
    if (!is.na(lr) && lr > max_lr) {
      max_lr <- lr
      best_k <- k
    }
  }
  
  list(stat = max_lr, best_k = best_k)
}


# ----------------------------------------------------------------------
# 2b. Full sequential online bivariate GLR detector
#     Runs game-by-game and stops at first alarm.
#     Direct extension of your online_glr_detect().
# ----------------------------------------------------------------------

online_bivariate_glr_detect <- function(Y, h, k_min = 5) {
  # Y      : n x 2 matrix (col 1 = team A PM, col 2 = team B PM)
  # h      : calibrated threshold
  # k_min  : buffer (5 games, same as offline)
  
  n <- nrow(Y)
  
  # Estimate sigma robustly — IDENTICAL to your offline procedure
  # Applied independently to each component
  sigma_A <- mad(diff(Y[, 1])) / sqrt(2)
  if (sigma_A == 0) sigma_A <- sd(Y[, 1])
  
  sigma_B <- mad(diff(Y[, 2])) / sqrt(2)
  if (sigma_B == 0) sigma_B <- sd(Y[, 2])
  
  # Cumulative sums for each component
  cumS_A <- cumsum(Y[, 1])
  cumS_B <- cumsum(Y[, 2])
  
  glr_seq <- numeric(n)        # GLR statistic at each game t
  k_seq   <- rep(NA_integer_, n)  # best k at each game t
  
  T_alarm <- NA_integer_       # stopping time
  tau_hat <- NA_integer_       # estimated changepoint at alarm
  
  for (t in 1:n) {
    res        <- bivariate_glr_at_t(cumS_A, cumS_B, t, sigma_A, sigma_B, k_min)
    glr_seq[t] <- res$stat
    k_seq[t]   <- res$best_k
    
    # Stop at first crossing — genuine sequential detection
    if (!is.na(res$stat) && res$stat > h && is.na(T_alarm)) {
      T_alarm <- t
      tau_hat <- res$best_k
      break
    }
  }
  
  list(
    T_alarm  = T_alarm,
    tau_hat  = tau_hat,
    glr_seq  = glr_seq,
    k_seq    = k_seq,
    sigma_A  = sigma_A,
    sigma_B  = sigma_B,
    h        = h,
    n        = n
  )
}


# Also provide a version that runs through the FULL season (no early stopping)
# for plotting purposes — so we can see the full S_t trajectory
online_bivariate_glr_full <- function(Y, h, k_min = 5) {
  
  n <- nrow(Y)
  
  sigma_A <- mad(diff(Y[, 1])) / sqrt(2)
  if (sigma_A == 0) sigma_A <- sd(Y[, 1])
  
  sigma_B <- mad(diff(Y[, 2])) / sqrt(2)
  if (sigma_B == 0) sigma_B <- sd(Y[, 2])
  
  cumS_A <- cumsum(Y[, 1])
  cumS_B <- cumsum(Y[, 2])
  
  glr_seq <- numeric(n)
  k_seq   <- rep(NA_integer_, n)
  
  T_alarm <- NA_integer_
  tau_hat <- NA_integer_
  
  for (t in 1:n) {
    res        <- bivariate_glr_at_t(cumS_A, cumS_B, t, sigma_A, sigma_B, k_min)
    glr_seq[t] <- res$stat
    k_seq[t]   <- res$best_k
    
    # Record first alarm but DON'T stop — continue for plotting
    if (!is.na(res$stat) && res$stat > h && is.na(T_alarm)) {
      T_alarm <- t
      tau_hat <- res$best_k
    }
  }
  
  list(
    T_alarm  = T_alarm,
    tau_hat  = tau_hat,
    glr_seq  = glr_seq,
    k_seq    = k_seq,
    sigma_A  = sigma_A,
    sigma_B  = sigma_B,
    h        = h,
    n        = n
  )
}


# =============================================================================
# SECTION 3: MONTE CARLO THRESHOLD CALIBRATION
#
# Direct extension of your tune_online_glr_threshold() to 2 dimensions.
#
# Under H0 (no changepoint), simulate bivariate iid N(0, I_2) data,
# compute max_t G_t over the season, and take the (1 - alpha) quantile.
#
# Same justification as your code: the GLR statistic depends only on
# DIFFERENCES between segment means, so it is invariant to the true
# mean — we simulate with mean = 0 for simplicity.
# =============================================================================

tune_online_bivariate_glr_threshold <- function(n          = 82,
                                                alpha      = 0.05,
                                                replicates = 2000,
                                                k_min      = 5,
                                                seed       = 123) {
  set.seed(seed)
  cat(sprintf("Calibrating online bivariate GLR threshold (n=%d, alpha=%.2f, reps=%d)...\n",
              n, alpha, replicates))
  
  max_glr_null <- numeric(replicates)
  
  for (r in 1:replicates) {
    # Null: bivariate iid N(0, I_2)
    Y_null <- matrix(rnorm(n * 2), ncol = 2)
    
    cumS_A <- cumsum(Y_null[, 1])
    cumS_B <- cumsum(Y_null[, 2])
    
    # sigma = 1 under the null (since we simulated N(0,1))
    max_g <- 0
    for (t in 1:n) {
      res <- bivariate_glr_at_t(cumS_A, cumS_B, t,
                                sigma_A = 1, sigma_B = 1, k_min = k_min)
      if (res$stat > max_g) max_g <- res$stat
    }
    max_glr_null[r] <- max_g
  }
  
  h <- quantile(max_glr_null, probs = 1 - alpha)
  cat(sprintf("-> Calibrated threshold h = %.4f\n\n", h))
  return(as.numeric(h))
}


# =============================================================================
# SECTION 4: PLOTTING
# =============================================================================

plot_online_bivariate_glr <- function(biv_data, result, team_a, team_b,
                                      trade_game_idx = NULL,
                                      player_name = NULL) {
  
  par(mfrow = c(3, 1), mar = c(4, 4.5, 3, 1))
  
  n     <- nrow(biv_data)
  alarm <- result$T_alarm
  tau   <- result$tau_hat
  
  # Panel 1: Sending team raw plus-minus
  plot(1:n, biv_data$pm_a, type = "l", col = "steelblue", lwd = 1.5,
       main = paste0(team_a, " Plus-Minus [Sending]"),
       xlab = "Game Number", ylab = "Plus-Minus")
  abline(h = 0, col = "grey70", lty = 3)
  if (!is.null(trade_game_idx)) abline(v = trade_game_idx, col = "darkgreen", lwd = 2)
  if (!is.na(alarm)) abline(v = alarm, col = "red", lty = 2, lwd = 2)
  if (!is.na(tau))   abline(v = tau, col = "purple", lty = 3, lwd = 1.5)
  
  # Panel 2: Receiving team raw plus-minus
  plot(1:n, biv_data$pm_b, type = "l", col = "darkorange", lwd = 1.5,
       main = paste0(team_b, " Plus-Minus [Receiving]"),
       xlab = "Game Number", ylab = "Plus-Minus")
  abline(h = 0, col = "grey70", lty = 3)
  if (!is.null(trade_game_idx)) abline(v = trade_game_idx, col = "darkgreen", lwd = 2)
  if (!is.na(alarm)) abline(v = alarm, col = "red", lty = 2, lwd = 2)
  if (!is.na(tau))   abline(v = tau, col = "purple", lty = 3, lwd = 1.5)
  
  # Panel 3: Running GLR statistic
  glr <- result$glr_seq
  plot(1:length(glr), glr, type = "l", col = "black", lwd = 1.5,
       main = bquote(G[t] ~ " | h = " ~ .(round(result$h, 2))),
       xlab = "Game Number", ylab = expression(G[t]))
  abline(h = result$h, col = "blue", lty = 3, lwd = 2)
  if (!is.null(trade_game_idx)) abline(v = trade_game_idx, col = "darkgreen", lwd = 2)
  if (!is.na(alarm)) abline(v = alarm, col = "red", lty = 2, lwd = 2)
  if (!is.na(tau))   abline(v = tau, col = "purple", lty = 3, lwd = 1.5)
  
  # Legend
  legend_labels <- c(
    if (!is.na(alarm)) paste0("Alarm t=", alarm) else "No alarm",
    paste0("h = ", round(result$h, 1))
  )
  legend_cols <- c("red", "blue")
  legend_ltys <- c(2, 3)
  if (!is.na(tau)) {
    legend_labels <- c(legend_labels, paste0("Est. CP = ", tau))
    legend_cols <- c(legend_cols, "purple")
    legend_ltys <- c(legend_ltys, 3)
  }
  if (!is.null(player_name)) {
    legend_labels <- c(legend_labels, paste0("Trade: ", player_name))
    legend_cols <- c(legend_cols, "darkgreen")
    legend_ltys <- c(legend_ltys, 1)
  }
  legend("topright", legend = legend_labels, col = legend_cols,
         lty = legend_ltys, lwd = 2, cex = 0.6, bg = "white")
}


# =============================================================================
# SECTION 5: ANALYSIS PIPELINE
# =============================================================================

analyse_online_glr <- function(df, team_from, team_to, season,
                               trade_date, player_name,
                               h = NULL, alpha = 0.05,
                               mc_reps = 2000, k_min = 5,
                               hit_window = 5,
                               plot = TRUE) {
  
  cat("\n==============================\n")
  cat(player_name, "(", team_from, "->", team_to, ", season", season, ")\n")
  cat("==============================\n")
  
  # Prepare data
  biv_data <- prepare_bivariate_series(df, team_from, team_to, season)
  n <- nrow(biv_data)
  cat("Series length:", n, "games\n")
  
  Y <- cbind(biv_data$pm_a, biv_data$pm_b)
  
  # Calibrate threshold if not provided
  if (is.null(h)) {
    h <- tune_online_bivariate_glr_threshold(
      n = n, alpha = alpha, replicates = mc_reps, k_min = k_min
    )
  }
  
  # Run FULL season (for plotting) — records first alarm but doesn't stop
  result <- online_bivariate_glr_full(Y, h = h, k_min = k_min)
  
  # Trade game index
  dates_a   <- as.Date(biv_data$dateGame_a)
  trade_idx <- which.min(abs(dates_a - as.Date(trade_date)))
  
  # Evaluate
  alarm   <- result$T_alarm
  tau_hat <- result$tau_hat
  
  if (!is.na(alarm)) {
    delay <- alarm - trade_idx
    hit   <- abs(alarm - trade_idx) <= hit_window
    
    cat(sprintf("\n*** ALARM at game %d ***\n", alarm))
    cat(sprintf("Estimated changepoint (tau_hat): game %d\n", tau_hat))
    cat(sprintf("Known trade date: game %d\n", trade_idx))
    cat(sprintf("Detection delay (alarm - trade): %+d games\n", delay))
    cat(sprintf("Changepoint accuracy (tau_hat - trade): %+d games\n",
                tau_hat - trade_idx))
  } else {
    delay <- NA
    hit   <- FALSE
    cat("\nNo alarm raised during the season.\n")
    cat(sprintf("Known trade date: game %d\n", trade_idx))
    cat(sprintf("Max GLR reached: %.2f (threshold: %.2f, ratio: %.3f)\n",
                max(result$glr_seq, na.rm = TRUE), h,
                max(result$glr_seq, na.rm = TRUE) / h))
  }
  
  # Mean shifts at trade point
  if (trade_idx > 1 && trade_idx < n) {
    cat(sprintf("\nMean shifts at trade date:\n"))
    cat(sprintf("  %s: before=%.1f, after=%.1f (shift=%+.1f)\n",
                team_from,
                mean(biv_data$pm_a[1:trade_idx]),
                mean(biv_data$pm_a[(trade_idx + 1):n]),
                mean(biv_data$pm_a[(trade_idx + 1):n]) - mean(biv_data$pm_a[1:trade_idx])))
    cat(sprintf("  %s: before=%.1f, after=%.1f (shift=%+.1f)\n",
                team_to,
                mean(biv_data$pm_b[1:trade_idx]),
                mean(biv_data$pm_b[(trade_idx + 1):n]),
                mean(biv_data$pm_b[(trade_idx + 1):n]) - mean(biv_data$pm_b[1:trade_idx])))
  }
  
  # Plot
  if (plot) {
    plot_online_bivariate_glr(biv_data, result, team_from, team_to,
                              trade_idx, player_name)
  }
  
  return(list(
    player    = player_name,
    team_from = team_from,
    team_to   = team_to,
    season    = season,
    n         = n,
    h         = h,
    T_alarm   = alarm,
    tau_hat   = tau_hat,
    trade_idx = trade_idx,
    delay     = delay,
    hit       = hit,
    detected  = !is.na(alarm),
    result    = result,
    biv_data  = biv_data
  ))
}


# =============================================================================
# SECTION 6: FOUR-CASE TRADE REGISTER
# =============================================================================

four_cases <- data.frame(
  player     = c("Pascal Siakam", "DeMarcus Cousins",
                 "Nikola Vucevic", "Kyrie Irving"),
  from       = c("TOR", "SAC", "ORL", "BKN"),
  to         = c("IND", "NOP", "CHI", "DAL"),
  date       = as.Date(c("2024-01-16", "2017-02-20",
                         "2021-03-25", "2023-02-06")),
  season     = c(2024, 2017, 2021, 2023),
  est_impact = c(7, 8, 5, 5),
  stringsAsFactors = FALSE
)


# =============================================================================
# SECTION 7: RUN THE ANALYSIS
# =============================================================================
# Uncomment to run.
# NOTE: Check team slugs — NOP may be "NOH" in older data.
# =============================================================================

# df <- read.csv("nba_game_logs.csv")
#
#
# # ---- 7a. Calibrate thresholds at both alpha levels -----------------------
h_82_05 <- tune_online_bivariate_glr_threshold(n = 82, alpha = 0.05,
                                                 replicates = 3000, k_min = 5)
h_82_10 <- tune_online_bivariate_glr_threshold(n = 82, alpha = 0.10,
                                                 replicates = 3000, k_min = 5)
# COVID-shortened season
h_72_05 <- tune_online_bivariate_glr_threshold(n = 72, alpha = 0.05,
                                                 replicates = 3000, k_min = 5)
h_72_10 <- tune_online_bivariate_glr_threshold(n = 72, alpha = 0.10,
                                                 replicates = 3000, k_min = 5)


# ---- 7b. Run all four cases at both alpha levels -------------------------
all_results <- list()

for (alpha_level in c(0.05, 0.10)) {
  for (i in 1:nrow(four_cases)) {
    tr <- four_cases[i, ]

    # Pick threshold for season length and alpha
    if (tr$season == 2021) {
      h_use <- if (alpha_level == 0.05) h_72_05 else h_72_10
    } else {
      h_use <- if (alpha_level == 0.05) h_82_05 else h_82_10
    }

    label <- paste0(tr$player, "_a", gsub("0\\.", "", as.character(alpha_level)))

    out <- analyse_online_glr(
      df, team_from = tr$from, team_to = tr$to, season = tr$season,
      trade_date = tr$date, player_name = tr$player,
      h = h_use, k_min = 5, hit_window = 5, plot = FALSE
    )
    out$alpha <- alpha_level

    all_results[[label]] <- out
  }
}


# ---- 7c. Summary table ---------------------------------------------------
summary_table <- do.call(rbind, lapply(all_results, function(r) {
  data.frame(
    Player     = r$player,
    Season     = r$season,
    Alpha      = r$alpha,
    n          = r$n,
    h          = round(r$h, 2),
    Max_GLR    = round(max(r$result$glr_seq, na.rm = TRUE), 2),
    Ratio      = round(max(r$result$glr_seq, na.rm = TRUE) / r$h, 3),
    Alarm      = ifelse(is.na(r$T_alarm), "—", as.character(r$T_alarm)),
    tau_hat    = ifelse(is.na(r$tau_hat), "—", as.character(r$tau_hat)),
    Trade_At   = r$trade_idx,
    Delay      = ifelse(is.na(r$delay), "—", sprintf("%+d", r$delay)),
    CP_error   = ifelse(is.na(r$tau_hat), "—",
                        sprintf("%+d", r$tau_hat - r$trade_idx)),
    Detected   = r$detected,
    stringsAsFactors = FALSE
  )
}))
rownames(summary_table) <- NULL

cat("\n\n=== FULL RESULTS (both alpha levels) ===\n\n")
print(summary_table, row.names = FALSE)


# ---- 7d. Publication plots (use alpha = 0.10 for main figures) -----------
for (i in 1:nrow(four_cases)) {
  tr <- four_cases[i, ]

  if (tr$season == 2021) {
    h_use <- h_72_10
  } else {
    h_use <- h_82_10
  }

  cat("\n--- Plot:", tr$player, "(alpha = 0.10) ---\n")
  analyse_online_glr(
    df, team_from = tr$from, team_to = tr$to, season = tr$season,
    trade_date = tr$date, player_name = tr$player,
    h = h_use, k_min = 5, hit_window = 5, plot = TRUE
  )
}