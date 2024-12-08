library(tidyverse)
library(foreach)

rho_sim <- 0.5
alpha_sq_sim <- 1
sigma_sq_sim <- 1
n_sim <- 100
K_sim <- 5
p_sim <- 150
n_vec <- seq(150, 50, length.out = K_sim+1)

##########################################
################mean shift################
##########################################
lambda_grid_nonI_r <- expand.grid(lambda_sim = seq(0.5, 3, 0.2),
                                  ite = 1:50)
X_cov_list_I <- map(c(150), ~{diag(..1)})
X_cov_foreach_I <- X_cov_list_I[map_dbl(lambda_grid_nonI_r$p_sim, ~{which(unique(lambda_grid_nonI_r$p_sim) == .)})]

X_cov_list <- map(c(p_sim), exp_decay_cov, width_actual = 101, c = 0.07)
X_mean_list <- map(1:length(X_cov_list_I[[1]]), ~{
  1
})
set.seed(111)

pred_empirical_mean <- foreach(lambda_sim = lambda_grid_nonI_r$lambda_sim,
                               X_cov = map(1:nrow(lambda_grid_nonI_r), 
                                           ~{X_cov_list[[1]]}),
                               mean_sh = X_mean_list) %do%{
                                 print(p_sim)
                                 gamma_vec <- p_sim / n_vec
                                 rho_mat <- matrix(alpha_sq_sim * sigma_sq_sim * rho_sim/p_sim, nrow = K_sim+1, ncol = K_sim+1)
                                 diag(rho_mat) <- alpha_sq_sim * sigma_sq_sim / p_sim
                                 beta_mat <- MASS::mvrnorm(p_sim, rep(0, K_sim+1), Sigma = rho_mat)
                                 
                                 X_list <- map2(1:(K_sim+1), n_vec, ~{MASS::mvrnorm(.y, rep(0, p_sim), Sigma = X_cov)})
                                 
                                 
                                 Y_list <- map2(lapply(seq_len(ncol(beta_mat)), function(i) beta_mat[,i]), 
                                                X_list, 
                                                ~{..2 %*% as.matrix(..1) + rnorm(nrow(..2))})
                                 
                                 
                                 #######################  #######################  #######################  #######################
                                 beta_hat_list <- map2(X_list, Y_list, ~{
                                   solve(t(..1) %*% ..1 + nrow(..1) / (nrow(..1) / n_vec[K_sim+1]) * lambda_sim * diag(p_sim)) %*% t(..1) %*% ..2
                                 })
                                 #######################  #######################  #######################  #######################
                                 
                                 lambda_use <- map_dbl(n_vec, ~{lambda_sim / (.x / n_vec[K_sim+1])})
                                 
                                 Sigma_hat <- map2(X_list, n_vec, ~{t(.x) %*% .x / .y})
                                 resolvent <- map2(Sigma_hat, lambda_use, ~{solve(.x + 
                                                                                    diag(.y, p_sim))})
                                 m <- map_dbl(resolvent, ~{sum(diag(.x)/p_sim)}) 
                                 m_per <- map_dbl(resolvent, ~{sum(diag(.x %*% .x)/p_sim)})
                                 v <- gamma_vec * (m - 1 / lambda_use) + (1 / lambda_use)
                                 v_per <- gamma_vec * (m_per - 1 / lambda_use^2) + (1 / lambda_use^2)
                                 
                                 order_0 <- 1
                                 order_1 <- (1/gamma_vec) *  (1 /(v*lambda_use) - 1)
                                 order_2 <- (1/gamma_vec) * (v - lambda_use * v_per) / (lambda_use * v)^2
                                 # order_11 <- (m - lambda_use * m_per) / 
                                 #   (1 - gamma_vec + gamma_vec * lambda_use^2 * m_per)
                                 
                                 cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                 diag(cov_mat_beta) <- 1
                                 
                                 C_ij <- matrix(1, K_sim + 1, K_sim+1)
                                 for (x in 1:(K_sim + 1)) {
                                   for (y in 1:(K_sim + 1)) {
                                     C_ij[x, y] <- order_0 - lambda_use[x] * order_1[x] - lambda_use[y] * order_1[y] +
                                       (lambda_use[x] * m[x] - lambda_use[y] * m[y]) / ((m[y] - m[x]))
                                   }
                                 }
                                 
                                 #return(list(order_0, order_1, order_2, lambda_use))
                                 
                                 D_i <- matrix(order_0 - lambda_use * order_1, K_sim + 1, 1)
                                 F_i <- diag(order_1 - lambda_use* order_2)
                                 diag(C_ij) <- (order_0 - 2 * lambda_use * order_1 + lambda_use^2 * order_2)
                                 
                                 D_f <- D_i * cov_mat_beta[, (K_sim+1)]* alpha_sq_sim * sigma_sq_sim
                                 F_f <- gamma_vec * sigma_sq_sim * F_i
                                 C_f <- C_ij * cov_mat_beta * alpha_sq_sim * sigma_sq_sim
                                 
                                 weight_temp <- as.numeric(solve(C_f + F_f) %*% D_f)
                                 
                                 
                                 beta_pool <- map2(beta_hat_list, weight_temp, ~{..1 * ..2}) %>% 
                                   Reduce( '+', .)
                                 
                                 #####
                                 n_vec_try <- n_vec
                                 p_sim_try <- p_sim
                                 resolvent_list <-  map2(X_list, Y_list, ~{
                                   solve(t(..1) %*% ..1 / nrow(..1) + 1 / (nrow(..1) / n_vec_try[K_sim+1]) * lambda_sim * diag(p_sim_try))
                                 })
                                 
                                 m_vec <- map_dbl(resolvent_list, ~{sum(diag(.x)/p_sim_try)})
                                 m_per_vec <- map_dbl(resolvent_list, ~{sum(diag(.x %*% .x)/p_sim_try)})
                                 gamma_vec <- p_sim_try/(n_vec_try)
                                 lambda_vec <- lambda_sim / (n_vec_try / n_vec_try[K_sim+1])
                                 
                                 cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                 diag(cov_mat_beta) <- 1
                                 
                                 cross_t <- outer(1:(K_sim+1), 1:(K_sim+1), function(x, y){
                                   1  + (lambda_vec[x] - lambda_vec[y]) *
                                     (m_vec[x] * m_vec[y] / (m_vec[x] - m_vec[y]))
                                   
                                   #(1 - lambda_vec[x] * m_vec[x]) * (1 - lambda_vec[y] * m_vec[y])
                                 })
                                 diag(cross_t) <- 0
                                 
                                 V <- sigma_sq_sim * cov_mat_beta[, K_sim+1] * alpha_sq_sim * (1 - lambda_vec * m_vec)
                                 A_ij <- cov_mat_beta * sigma_sq_sim * alpha_sq_sim * cross_t
                                 A_ii <- sigma_sq_sim * alpha_sq_sim * (1 - 2 * lambda_vec * m_vec + lambda_vec^2 * m_per_vec)
                                 R_ii <- sigma_sq_sim * (gamma_vec * m_vec - gamma_vec * lambda_vec * m_per_vec)
                                 diag(A_ij) <- A_ii
                                 R <- diag(R_ii)
                                 A <- A_ij
                                 res <- solve(A + R) %*% as.matrix(V)
                                 
                                 beta_pool_r <- map2(beta_hat_list, res, ~{..1 * ..2}) %>%
                                   Reduce( '+', .)
                                 
                                 X_test <- MASS::mvrnorm(1000, rep(0, p_sim), Sigma = X_cov)
                                 X_test_t <- MASS::mvrnorm(1000, rep(mean_sh, p_sim), Sigma = X_cov)
                                 
                                 
                                 y_test <- X_test %*% beta_mat[,K_sim+1] + rnorm(1000)
                                 y_test_r <- X_test_t %*% beta_mat[,K_sim+1] + rnorm(1000)
                                 
                                 y_hat_TL <- X_test %*% beta_pool
                                 y_hat_TL_rob <- X_test_t %*% beta_pool_r
                                 
                                 y_hat_naive <- X_test %*% beta_hat_list[[K_sim+1]]
                                 y_hat_naive_r <- X_test_t %*% beta_hat_list[[K_sim+1]]
                                 
                                 
                                 data.frame(TL_risk= mean((y_hat_TL - y_test)^2),
                                            TL_risk_r = mean((y_hat_TL - y_test_r)^2),
                                            TL_risk_rob = mean((y_hat_TL_rob - y_test_r)^2),
                                            naive_risk = mean((as.numeric(y_hat_naive - y_test)^2)),
                                            naive_risk_r = mean((as.numeric(y_hat_naive_r - y_test)^2))
                                 )
                                 
                               }

pred_empirical_mean %>%
  reduce(rbind) %>%
  cbind(lambda_grid_nonI_r) %>%
  pivot_longer(1:5) %>%
  filter(name != "naive_risk") %>%
  mutate(name = case_when(name == "naive_risk_r" ~ "Naive Ridge (Shifted)",
                          name == "TL_risk" ~ "Optimal Prediction Weight",
                          name == "TL_risk_r" ~ "Optimal Prediction Weight (Shifted)",
                          name == "TL_risk_rob" ~ "Optimal Estimation Weight (Shifted)")) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(lambda_sim), y = value,  fill = name), lwd = 0.3, outlier.size = 0.3) +
  labs(x = expression("lambda ("~lambda~")"), y = "Prediction Risk", fill = "") +
  ylim(0, 3)  +
  theme_bw(base_size = 10)+
  theme(legend.position="bottom")+ 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


################################################
################Covariance shift################
################################################

lambda_grid_nonI_r <- expand.grid(lambda_sim = seq(0.5, 3, 0.2),
                                  ite = 1:50)

X_cov_list_shift <- map(c(p_sim), exp_decay_cov, width_actual = 21, c = 0.7)
X_cov_list <- map(c(p_sim), exp_decay_cov, width_actual = 101, c = 0.07)
set.seed(111)

pred_empirical_shift_cov <- foreach(lambda_sim = lambda_grid_nonI_r$lambda_sim,
                                    X_cov = map(1:nrow(lambda_grid_nonI_r), 
                                                ~{X_cov_list[[1]]}),
                                    X_cov_t = map(1:nrow(lambda_grid_nonI_r), 
                                                  ~{X_cov_list_shift[[1]]})) %do%{
                                                    print(p_sim)
                                                    gamma_vec <- p_sim / n_vec
                                                    rho_mat <- matrix(alpha_sq_sim * sigma_sq_sim * rho_sim/p_sim, nrow = K_sim+1, ncol = K_sim+1)
                                                    diag(rho_mat) <- alpha_sq_sim * sigma_sq_sim / p_sim
                                                    beta_mat <- MASS::mvrnorm(p_sim, rep(0, K_sim+1), Sigma = rho_mat)
                                                    
                                                    X_list <- map2(1:(K_sim+1), n_vec, ~{MASS::mvrnorm(.y, rep(0, p_sim), Sigma = X_cov)})
                                                    
                                                    
                                                    Y_list <- map2(lapply(seq_len(ncol(beta_mat)), function(i) beta_mat[,i]), 
                                                                   X_list, 
                                                                   ~{..2 %*% as.matrix(..1) + rnorm(nrow(..2))})
                                                    
                                                    
                                                    #######################  #######################  #######################  #######################
                                                    beta_hat_list <- map2(X_list, Y_list, ~{
                                                      solve(t(..1) %*% ..1 + nrow(..1) / (nrow(..1) / n_vec[K_sim+1]) * lambda_sim * diag(p_sim)) %*% t(..1) %*% ..2
                                                    })
                                                    #######################  #######################  #######################  #######################
                                                    
                                                    lambda_use <- map_dbl(n_vec, ~{lambda_sim / (.x / n_vec[K_sim+1])})
                                                    
                                                    Sigma_hat <- map2(X_list, n_vec, ~{t(.x) %*% .x / .y})
                                                    resolvent <- map2(Sigma_hat, lambda_use, ~{solve(.x + 
                                                                                                       diag(.y, p_sim))})
                                                    m <- map_dbl(resolvent, ~{sum(diag(.x)/p_sim)}) 
                                                    m_per <- map_dbl(resolvent, ~{sum(diag(.x %*% .x)/p_sim)})
                                                    v <- gamma_vec * (m - 1 / lambda_use) + (1 / lambda_use)
                                                    v_per <- gamma_vec * (m_per - 1 / lambda_use^2) + (1 / lambda_use^2)
                                                    
                                                    order_0 <- 1
                                                    order_1 <- (1/gamma_vec) *  (1 /(v*lambda_use) - 1)
                                                    order_2 <- (1/gamma_vec) * (v - lambda_use * v_per) / (lambda_use * v)^2
                                                    # order_11 <- (m - lambda_use * m_per) / 
                                                    #   (1 - gamma_vec + gamma_vec * lambda_use^2 * m_per)
                                                    
                                                    cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                                    diag(cov_mat_beta) <- 1
                                                    
                                                    C_ij <- matrix(1, K_sim + 1, K_sim+1)
                                                    for (x in 1:(K_sim + 1)) {
                                                      for (y in 1:(K_sim + 1)) {
                                                        C_ij[x, y] <- order_0 - lambda_use[x] * order_1[x] - lambda_use[y] * order_1[y] +
                                                          (lambda_use[x] * m[x] - lambda_use[y] * m[y]) / ((m[y] - m[x]))
                                                      }
                                                    }
                                                    
                                                    #return(list(order_0, order_1, order_2, lambda_use))
                                                    
                                                    D_i <- matrix(order_0 - lambda_use * order_1, K_sim + 1, 1)
                                                    F_i <- diag(order_1 - lambda_use* order_2)
                                                    diag(C_ij) <- (order_0 - 2 * lambda_use * order_1 + lambda_use^2 * order_2)
                                                    
                                                    D_f <- D_i * cov_mat_beta[, (K_sim+1)]* alpha_sq_sim * sigma_sq_sim
                                                    F_f <- gamma_vec * sigma_sq_sim * F_i
                                                    C_f <- C_ij * cov_mat_beta * alpha_sq_sim * sigma_sq_sim
                                                    
                                                    weight_temp <- as.numeric(solve(C_f + F_f) %*% D_f)
                                                    
                                                    
                                                    beta_pool <- map2(beta_hat_list, weight_temp, ~{..1 * ..2}) %>% 
                                                      Reduce( '+', .)
                                                    
                                                    #####
                                                    n_vec_try <- n_vec
                                                    p_sim_try <- p_sim
                                                    resolvent_list <-  map2(X_list, Y_list, ~{
                                                      solve(t(..1) %*% ..1 / nrow(..1) + 1 / (nrow(..1) / n_vec_try[K_sim+1]) * lambda_sim * diag(p_sim_try))
                                                    })
                                                    
                                                    m_vec <- map_dbl(resolvent_list, ~{sum(diag(.x)/p_sim_try)})
                                                    m_per_vec <- map_dbl(resolvent_list, ~{sum(diag(.x %*% .x)/p_sim_try)})
                                                    gamma_vec <- p_sim_try/(n_vec_try)
                                                    lambda_vec <- lambda_sim / (n_vec_try / n_vec_try[K_sim+1])
                                                    
                                                    cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                                    diag(cov_mat_beta) <- 1
                                                    
                                                    cross_t <- outer(1:(K_sim+1), 1:(K_sim+1), function(x, y){
                                                      1  + (lambda_vec[x] - lambda_vec[y]) *
                                                        (m_vec[x] * m_vec[y] / (m_vec[x] - m_vec[y]))
                                                      
                                                      #(1 - lambda_vec[x] * m_vec[x]) * (1 - lambda_vec[y] * m_vec[y])
                                                    })
                                                    diag(cross_t) <- 0
                                                    
                                                    V <- sigma_sq_sim * cov_mat_beta[, K_sim+1] * alpha_sq_sim * (1 - lambda_vec * m_vec)
                                                    A_ij <- cov_mat_beta * sigma_sq_sim * alpha_sq_sim * cross_t
                                                    A_ii <- sigma_sq_sim * alpha_sq_sim * (1 - 2 * lambda_vec * m_vec + lambda_vec^2 * m_per_vec)
                                                    R_ii <- sigma_sq_sim * (gamma_vec * m_vec - gamma_vec * lambda_vec * m_per_vec)
                                                    diag(A_ij) <- A_ii
                                                    R <- diag(R_ii)
                                                    A <- A_ij
                                                    res <- solve(A + R) %*% as.matrix(V)
                                                    
                                                    beta_pool_r <- map2(beta_hat_list, res, ~{..1 * ..2}) %>%
                                                      Reduce( '+', .)
                                                    
                                                    X_test <- MASS::mvrnorm(1000, rep(0, p_sim), Sigma = X_cov)
                                                    X_test_t <- MASS::mvrnorm(1000, rep(0, p_sim), Sigma = X_cov_t)
                                                    
                                                    
                                                    y_test <- X_test %*% beta_mat[,K_sim+1] + rnorm(1000)
                                                    y_test_r <- X_test_t %*% beta_mat[,K_sim+1] + rnorm(1000)
                                                    
                                                    y_hat_TL <- X_test %*% beta_pool
                                                    y_hat_TL_rob <- X_test_t %*% beta_pool_r
                                                    
                                                    y_hat_naive <- X_test %*% beta_hat_list[[K_sim+1]]
                                                    y_hat_naive_r <- X_test_t %*% beta_hat_list[[K_sim+1]]
                                                    
                                                    
                                                    data.frame(TL_risk= mean((y_hat_TL - y_test)^2),
                                                               TL_risk_r = mean((y_hat_TL - y_test_r)^2),
                                                               TL_risk_rob = mean((y_hat_TL_rob - y_test_r)^2),
                                                               naive_risk = mean((as.numeric(y_hat_naive - y_test)^2)),
                                                               naive_risk_r = mean((as.numeric(y_hat_naive_r - y_test)^2))
                                                    )
                                                    
                                                  }

pred_empirical_shift_cov %>%
  reduce(rbind) %>%
  cbind(lambda_grid_nonI_r) %>%
  pivot_longer(1:5) %>%
  filter(name != "naive_risk") %>%
  mutate(name = case_when(name == "naive_risk_r" ~ "Naive Ridge (Shifted)",
                          name == "TL_risk" ~ "Optimal Prediction Weight",
                          name == "TL_risk_r" ~ "Optimal Prediction Weight (Shifted)",
                          name == "TL_risk_rob" ~ "Optimal Estimation Weight (Shifted)")) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(lambda_sim), y = value,  fill = name), lwd=0.3, outlier.size = 0.3) +
  labs(x = expression("lambda ("~lambda~")"), y = "Prediction Risk", fill = "") +
  ylim(0, 3)  +
  theme_bw(base_size = 10)+
  theme(legend.position="bottom")+ 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


