# Trans-Ridge
This repository contains the codes for our simulations and second real data example. All codes are written in R. A quick navigation is given as

- Simulations
  - Robustness Git.R : Contains the codes to perform simulations when there is a shift in the testing data distributions. The codes re-produce figure 6. in the manuscript.
  - Simulation Git.R : Contains the codes to perform comparisons of the theoretical risks obtained from Theorem 2.2 and Theorem 5.2, and risks based on synthetic data. The codes reproduce Figure 1. and Figure 4. in the manuscript.
- Real Data
  - TLridge Git.R: It contains the codes for estimation population parameters including $\rho_{kk'}, \sigma_k, \alpha_k$ as described in the supplement to the manuscript. In addition, it contains the codes to estimate the optimal prediction weights and the optimal estimation weight. At last, the leave one out AUC performance for both weights are reported.
  - Naive Git.R: It contains the codes to report the leave one out AUC performance for ridge regressions fitted with only target population data and with all data pooled together.
  - DF.rda: Raw data file of our second real data example.
