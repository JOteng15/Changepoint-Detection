# Sequential Optimization & Online GLR Changepoint Detection

## Overview
This repository implements an statistical monitoring suite in R, designed to detect structural breaks in sequential data streams. Moving beyond retrospective (offline) analysis, the core of this project focuses on **Online Sequential Detection**—monitoring live data streams to identify distributional shifts with minimal detection delay.


## Core Architecture & Algorithms

* **Online Bivariate GLR ($G_t$):** Implements a real-time sequential detector using a bivariate quadratic form. It monitors two simultaneous data streams, dynamically updating the covariance matrix to trigger alarms when the joint distribution structurally shifts.
* **Page's CUSUM with Resets:** A bounded-memory sequential algorithm that tracks positive and negative shifts in data streams. Upon breaching a calibrated threshold ($h$), the system registers the alarm, isolates the estimated changepoint ($\tau$), and resets its baseline to continue monitoring.
* **Monte Carlo Threshold Calibration:** Replaces arbitrary hyperparameters with statistically rigorous tuning. The code uses Monte Carlo simulations of Gaussian noise to dynamically calibrate the detection threshold ($h$), explicitly balancing the Average Run Length ($ARL_0$) against the False Alarm Rate ($\alpha$).
