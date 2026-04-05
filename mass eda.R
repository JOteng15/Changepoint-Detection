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

source("C:/Users/joten/Documents/Jevon Oteng Lancaster/Changepoint-Detection/load_data.R")

# Narrative: The NBA is a high-variance environment. 
# This plot compares the distribution of Net Ratings across all teams.
library(ggplot2)
library(dplyr)


df %>%
  ggplot(aes(x = reorder(slugTeam, plusminusTeam, FUN = median), y = plusminusTeam, fill = (slugTeam == "CLE"))) +
  geom_boxplot(outlier.alpha = 0.1) +
  scale_fill_manual(values = c("TRUE" = "#860038", "FALSE" = "grey80")) +
  coord_flip() +
  labs(title = "NBA Performance Volatility by Team (2010-2025)",
       subtitle = "The Cavaliers exhibit some of the widest performance spreads in the league.",
       x = "Team", y = "Game Plus/Minus") +
  theme_minimal() + theme(legend.position = "none")



df_cavs <- df %>% 
  filter(slugTeam == "CLE") %>% 
  arrange(yearSeason, numberGameTeamSeason) %>%
  mutate(cum_pm = cumsum(plusminusTeam),
         game_id = row_number())

ggplot(df_cavs, aes(x = game_id, y = cum_pm)) +
  geom_line(color = "#860038", linewidth = 1) +
  geom_vline(xintercept = c(328, 656), linetype = "dashed", alpha = 0.5) + # Approximated LeBron entry/exit
  annotate("text", x = 150, y = 500, label = "The 'Rebuild' Era", color = "grey40") +
  annotate("text", x = 480, y = 1500, label = "The 'Contention' Era", color = "grey40") +
  labs(title = "Cavaliers Cumulative Plus/Minus (2010-2025)",
       subtitle = "Shifts in the 'slope' suggest structural changes in team quality.",
       x = "Total Games Played", y = "Cumulative Plus/Minus") +
  theme_classic()


df %>%
  mutate(is_cavs = ifelse(slugTeam == "CLE", "Cavaliers", "Rest of League")) %>%
  ggplot(aes(x = tovTeam, y = ptsTeam, color = isWin)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~is_cavs) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("W" = "#00471B", "L" = "#EEE1C6")) +
  labs(title = "Turnovers vs. Points Scored",
       subtitle = "Comparing fundamental basketball relationships across groups.") +
  theme_minimal()

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
      title = "Cavs Rolling Mean and Variance on Plus Minus score (10 game window)",
      x = "Game",
      y = NULL
    ) +
    theme_minimal()
  #look at these graph, apparent missing v
}

plt9()





#--NEW TIMELINE FOR THE PAPER - EDA--------------------------------------------------------------

library(ggplot2)
library(dplyr)

# Highlight CLE in a league-wide boxplot
ggplot(df, aes(x = reorder(slugTeam, plusminusTeam, median), y = plusminusTeam, fill = (slugTeam == "CLE"))) +
  geom_boxplot(outlier.alpha = 0.1, color = "grey30") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#860038", "FALSE" = "grey80")) +
  labs(title = "NBA Performance Distribution (2010-2025)",
       subtitle = "Cavaliers show one of the widest spreads in the league.",
       x = "Team", y = "Game Plus/Minus") +
  theme_minimal() + theme(legend.position = "none")


team_stats <- df %>%
  group_by(slugTeam) %>%
  summarise(
    mean_pm = mean(plusminusTeam),
    sd_pm = sd(plusminusTeam)
  )

ggplot(team_stats, aes(x = mean_pm, y = sd_pm)) +
  geom_point(aes(color = (slugTeam == "CLE")), size = 4) +
  geom_text(aes(label = slugTeam), vjust = -1, size = 3, check_overlap = TRUE) +
  scale_color_manual(values = c("TRUE" = "#860038", "FALSE" = "grey70")) +
  labs(title = "Team Success vs. Performance Volatility",
       subtitle = "High SD suggests structural shifts (Potential Changepoints).",
       x = "Average Plus/Minus", y = "Standard Deviation of Performance") +
  theme_light() + theme(legend.position = "none")


ggplot(df, aes(x = fgaTeam, y = ptsTeam)) +
  geom_bin2d(bins = 30) +
  geom_smooth(data = filter(df, slugTeam == "CLE"), aes(color = "Cavaliers"), method = "lm", se = F) +
  scale_fill_gradient(low = "grey90", high = "grey40") +
  scale_color_manual(values = c("Cavaliers" = "#FDBB30")) +
  labs(title = "Scoring Efficiency: League vs. Cavaliers",
       x = "Field Goal Attempts", y = "Points") +
  theme_minimal()


library(zoo)

df_cavs <- df %>% 
  filter(slugTeam == "CLE") %>% 
  arrange(dateGame) %>%
  mutate(rolling_avg = rollmean(plusminusTeam, k = 41, fill = NA)) # 41 games = half season

ggplot(df_cavs, aes(x = seq_along(dateGame), y = plusminusTeam)) +
  geom_point(alpha = 0.1, color = "#860038") +
  geom_line(aes(y = rolling_avg), color = "#FDBB30", size = 1.2) +
  annotate("text", x = 200, y = 30, label = "The First Rebuild", fontface = "italic") +
  annotate("text", x = 500, y = 30, label = "The Championship Window", fontface = "italic") +
  labs(title = "Visualizing the 'Eras' of Cavaliers Basketball",
       subtitle = "The 41-game rolling average reveals clear structural shifts.",
       x = "Games Played (Cumulative)", y = "Plus/Minus") +
  theme_classic()





# GLOBAL PLOT 1: Testing the Normality Assumption
# Narrative: Verifying if NBA point differentials follow the Gaussian distribution 
# assumed by standard CUSUM models.

ggplot(df, aes(x = plusminusTeam)) +
  # 1. The Real Data (Histogram)
  geom_histogram(aes(y = ..density..), bins = 60, fill = "steelblue", alpha = 0.4, color = "white") +
  
  # 2. The Theoretical Normal Distribution (Red Line)
  stat_function(fun = dnorm, 
                args = list(mean = mean(df$plusminusTeam), sd = sd(df$plusminusTeam)), 
                color = "#860038", size = 1.2, linetype = "dashed") +
  
  labs(title = "Distributional Properties of NBA Margins",
       subtitle = "Data approximates Normality but exhibits 'heavy tails' (extreme blowouts).",
       x = "Game Point Differential", y = "Density") +
  theme_minimal()


# GLOBAL PLOT 2: The Home Court Advantage
# Narrative: Quantifying the structural bias introduced by game location.

ggplot(df, aes(x = locationGame, y = plusminusTeam, fill = locationGame)) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.1) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "white") + # Mark the Mean
  scale_fill_manual(values = c("A" = "#E41A1C", "H" = "#377EB8"), 
                    labels = c("Away", "Home")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Structural Bias: The Home Court Advantage",
       subtitle = "Home teams (Blue) enjoy a significant mean performance shift vs. Away teams (Red).",
       x = "Game Location", y = "Point Differential") +
  theme_bw()


# 1.1 League Evolution: The Changing Baseline
# Narrative: The NBA is not a static environment. 
# We observe a structural drift in Offensive Rating and 3-Point Volume across the league.

league_evolution <- df %>%
  # Calculate per-game metrics first
  group_by(yearSeason, idGame) %>%
  summarise(
    game_pts = mean(ptsTeam),
    game_3pa = mean(fg3aTeam),
    game_pace = mean(fgaTeam + 0.44*ftaTeam + tovTeam),
    .groups = 'drop'
  ) %>%
  # Now average by season
  group_by(yearSeason) %>%
  summarise(
    avg_pts = mean(game_pts),
    avg_3pa = mean(game_3pa),
    avg_pace = mean(game_pace)
  ) %>%
  pivot_longer(cols = starts_with("avg"), names_to = "metric", values_to = "value")

# Create Faceted Plot
metric_labels <- c(
  "avg_3pa" = "3-Point Attempts",
  "avg_pace" = "Pace (Possessions)",
  "avg_pts" = "Points per Game"
)

ggplot(league_evolution, aes(x = yearSeason, y = value)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_y", labeller = as_labeller(metric_labels)) +
  labs(title = "Non-Stationarity in League Averages (2010-2025)",
       subtitle = "Global drifts in scoring and strategy complicate static baselines.",
       x = "Season", y = "League Average") +
  theme_bw()


# 1.2 Drivers of Performance: Correlation Matrix
# Narrative: Identifying which metrics strongly correlate with Plus/Minus (Winning).

library(ggcorrplot)

# Select numeric columns relevant to performance
cor_data <- df %>%
  select(plusminusTeam, ptsTeam, fg3mTeam, fg3aTeam, 
         astTeam, tovTeam, orebTeam, drebTeam, pfTeam) %>%
  cor(use = "complete.obs")

# Plot Heatmap
ggcorrplot(cor_data, 
           method = "square", 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           colors = c("#B2182B", "white", "#2166AC"),
           title = "Correlation of Game Metrics with Performance")


# METRIC 4: The Efficiency Frontier (3P Volume vs Accuracy)
# Narrative: The Cavs broke the 'law of diminishing returns' during the Championship era.

# 1. Prepare Data
shooting_trend <- df %>%
  group_by(yearSeason, slugTeam) %>%
  summarise(
    vol_3pa = mean(fg3aTeam),
    acc_3p = mean(pctFG3Team),
    .groups = 'drop'
  )

# 2. Isolate Cavs and League Trend
cavs_shooting <- shooting_trend %>% filter(slugTeam == "CLE") %>% arrange(yearSeason)
league_trend <- shooting_trend %>% group_by(yearSeason) %>% summarise(vol_3pa = mean(vol_3pa), acc_3p = mean(acc_3p))

# 3. Plot
ggplot(shooting_trend, aes(x = vol_3pa, y = acc_3p)) +
  # Background: The "Cloud" of all teams
  geom_point(color = "grey85", alpha = 0.4) +
  
  # The Cavs Path
  geom_path(data = cavs_shooting, color = "#860038", size = 1.2, 
            arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  geom_point(data = cavs_shooting, color = "#FDBB30", size = 3) +
  
  # Labels for Context
  geom_text(data = filter(cavs_shooting, yearSeason %in% c(2011, 2015, 2017, 2019, 2022)), 
            aes(label = yearSeason), vjust = -1.5, color = "#860038", fontface = "bold") +
  
  # Annotations for "Regimes"
  annotate("text", x = 32, y = 0.385, label = "The 'Super-Team' Peak\n(High Vol, High Acc)", 
           color = "#860038", size = 3, fontface = "italic") +
  annotate("text", x = 15, y = 0.33, label = "The 'Dark Ages'\n(Low Vol, Low Acc)", 
           color = "grey50", size = 3, fontface = "italic") +
  
  labs(title = "Breaking the Trade-off: 3-Point Volume vs. Accuracy",
       subtitle = "Most teams sacrifice accuracy for volume. The 2016-17 Cavs maximized BOTH.",
       x = "Volume (3PA per Game)", 
       y = "Accuracy (3P%)") +
  theme_minimal()



# METRIC 1: 3-Point Evolution (Strategy)
# Narrative: Tracking the 'LeBron Effect' on shot selection.

# Calculate League Average vs Cavs per Season
trend_3pt <- df %>%
  group_by(yearSeason, is_cavs = (slugTeam == "CLE")) %>%
  summarise(avg_3pa = mean(fg3aTeam), .groups = 'drop')

ggplot(trend_3pt, aes(x = yearSeason, y = avg_3pa, color = is_cavs)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#860038"), 
                     labels = c("League Average", "Cavaliers")) +
  geom_vline(xintercept = c(2014.5, 2018.5), linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 2016.5, y = 33, label = "The 'Space' Era", color = "#860038", fontface="italic") +
  labs(title = "Tactical Evolution: 3-Point Attempts per Game",
       subtitle = "The Cavaliers' adoption of the 3-point shot highlights a strategic structural break.",
       x = "Season", y = "3PA per Game", color = "Team") +
  theme_minimal()



# METRIC 3: The Identity Path (Offensive vs Defensive Efficiency)
# Narrative: Tracing the migration of the Cavaliers through performance quadrants.

# 1. Calculate Efficiency Metrics
team_efficiency <- df %>%
  mutate(possessions = fgaTeam + 0.44*ftaTeam + tovTeam,
         off_rtg = 100 * (ptsTeam / possessions),
         # Def Rtg = Points Allowed / Possessions
         pts_allowed = ptsTeam - plusminusTeam, 
         def_rtg = 100 * (pts_allowed / possessions)) %>%
  group_by(slugTeam, yearSeason) %>%
  summarise(off_rtg = mean(off_rtg), 
            def_rtg = mean(def_rtg), .groups = 'drop')

# 2. Isolate Cavs for the "Path"
cavs_path <- team_efficiency %>% 
  filter(slugTeam == "CLE") %>%
  arrange(yearSeason)

# 3. Plot
ggplot(team_efficiency, aes(x = off_rtg, y = def_rtg)) +
  # Background: All teams as faint grey dots
  geom_point(color = "grey85", alpha = 0.5) +
  
  # The Cavs Path: Connected line with dots
  geom_path(data = cavs_path, color = "#860038", size = 1, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_point(data = cavs_path, color = "#FDBB30", size = 3) +
  
  # Label key years to show the "Journey"
  geom_text(data = filter(cavs_path, yearSeason %in% c(2011, 2016, 2018, 2019, 2022)), 
            aes(label = yearSeason), vjust = -1, fontface = "bold", color = "#860038") +
  
  # Invert Y axis because Lower Def Rating is BETTER
  scale_y_reverse() + 
  
  labs(title = "The Cavaliers' Identity Path (2010-2025)",
       subtitle = "Tracing the structural migration from Rebuild (Bottom Left) to Contender (Top Right).",
       x = "Offensive Rating (Pts/100 Poss)", 
       y = "Defensive Rating (Allowed/100 Poss) [Inverted]") +
  theme_bw()



team_stability <- df %>%
  # Step A: Calculate the average Plus/Minus for EACH SEASON first
  group_by(slugTeam, yearSeason) %>%
  summarise(
    season_mean = mean(plusminusTeam, na.rm = TRUE), 
    .groups = 'drop'
  ) %>%
  
  # Step B: Calculate the SD of those season averages (The Structural Volatility)
  group_by(slugTeam) %>%
  summarise(
    mean_overall = mean(season_mean, na.rm = TRUE),      # The X-axis
    sd_of_seasons = sd(season_mean, na.rm = TRUE)        # The Y-axis
  )