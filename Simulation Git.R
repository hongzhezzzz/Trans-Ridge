t_w_any <- function(cov_mat){
  t <- eigen(cov_mat, symmetric = T)$values
  p_sim <- nrow(cov_mat)
  w <- rep(1/p_sim, p_sim)
  return(list(t = t, 
              w = w))
}
compute_ST <- function(w, t, gamma, grid_size = 1e5){
  
  v <- seq(1/grid_size, 1e3, length.out = grid_size)
  z <- rep(1, grid_size)
  
  for (i in 1:grid_size) {
    z[i] = -1/v[i] + gamma * sum(w * t/(1 + v[i] * t))
  }
  
  v <- v[z<0]
  z <- z[z<0]
  
  lambda_seq <- -z
  lambda_ind <- (lambda_seq<50) & (lambda_seq>1e-2)
  lambda_seq <- lambda_seq[lambda_ind]
  
  v <- v[lambda_ind]
  z <- z[lambda_ind]
  m <- v/gamma + (1/gamma-1)/z
  
  L <- length(lambda_seq)
  v_prime <- rep(1,L)
  
  for (i in 1:L) {
    v_prime[i] = 1/(1/v[i]^2 - gamma * sum( w*t^2 / (1 + t*v[i])^2))
  }
  
  m_prime <- v_prime/gamma - (1/gamma-1)/z^2
  
  return(list(lambda_seq = lambda_seq, 
              m = m, 
              v = v, 
              m_prime = m_prime, 
              v_prime = v_prime))
}
library(tidyverse)
library(foreach)
lambda_sim <- 1
rho_sim <- 0.5
alpha_sq_sim <- 1
sigma_sq_sim <- 1
K_sim <- 5
################################################################################
######################    Estimation Risk       ################################
################################################################################

########################################
#Same number of Observations Covariance#
########################################
set.seed(111)
n_sim <- 500
lambda_grid <- expand.grid(lambda_sim = seq(0.5, 11, 1),
                           p_sim = c(250, 500, 750),
                           ite = 1:50)
X_cov_list_I <- map(c(250, 500, 750), ~{diag(..1)}) 
X_cov_foreach_I <- X_cov_list_I[map_dbl(lambda_grid$p_sim, 
                                        ~{which(unique(lambda_grid$p_sim) == .)})]


gamma_empirical_oldV <- foreach(lambda_sim = lambda_grid$lambda_sim,
                                p_sim = lambda_grid$p_sim,
                                X_cov = X_cov_foreach_I) %do%{
                                  print(p_sim)
                                  
                                  rho_mat <- matrix(alpha_sq_sim * sigma_sq_sim * rho_sim/p_sim, nrow = K_sim+1, ncol = K_sim+1)
                                  diag(rho_mat) <- alpha_sq_sim * sigma_sq_sim / p_sim
                                  beta_mat <- MASS::mvrnorm(p_sim, rep(0, K_sim+1), Sigma = rho_mat)
                                  
                                  X_list <- map(1:(K_sim+1), ~{MASS::mvrnorm(n_sim, rep(0, p_sim), Sigma = X_cov)})

                                  
                                  Sigma_hat_pool <- t(reduce(X_list, rbind)) %*% reduce(X_list, rbind) / (n_sim * K_sim + n_sim)
                                  Sigma_hat_naive <- t(X_list[[K_sim+1]]) %*% X_list[[K_sim+1]] / (n_sim)
                                  
                                  Y_list <- map2(lapply(seq_len(ncol(beta_mat)), function(i) beta_mat[,i]), 
                                                 X_list, 
                                                 ~{..2 %*% as.matrix(..1) + rnorm(n_sim)})
                                  
                                  #######################  #######################  #######################  #######################
                                  beta_hat_list <- map2(X_list, Y_list, ~{
                                    solve(t(..1) %*% ..1 + n_sim * lambda_sim * diag(p_sim)) %*% t(..1) %*% ..2
                                  })
                                  #######################  #######################  #######################  #######################
                                  
                                  resolvent <- solve(Sigma_hat_pool + diag(lambda_sim, p_sim))
                                  m <- sum(diag(resolvent)/p_sim)
                                  m_per <- sum(diag(resolvent %*% resolvent))/ p_sim
                                  
                                  gamma_ <- p_sim/(n_sim)
                                  

                                  A_ij <- rho_sim * sigma_sq_sim * alpha_sq_sim * (1 - lambda_sim * m)^2
                                  A_ii <- sigma_sq_sim * alpha_sq_sim * (1 - 2 * lambda_sim * m + lambda_sim^2 * m_per)
                                  V_i <- sigma_sq_sim * rho_sim * alpha_sq_sim * (1 - lambda_sim * m)
                                  V_K <- sigma_sq_sim * alpha_sq_sim * (1 - lambda_sim * m)
                                  R_ii <- sigma_sq_sim * (gamma_ * m - gamma_ * lambda_sim * m_per)
                                  
                                  V <- c(rep(V_i, K_sim), V_K)
                                  R <- diag(R_ii, K_sim+1)
                                  A <- matrix(A_ij, nrow = K_sim+1, ncol = K_sim+1)
                                  diag(A) <- A_ii
                                  res <- solve(A + R) %*% as.matrix(V)
                                  
                                  beta_pool <- map2(beta_hat_list, res, ~{..1 * ..2}) %>% 
                                    Reduce( '+', .)
                                  
                                  data.frame(TL_risk= sum((as.numeric(beta_pool) - beta_mat[,K_sim+1])^2),
                                             naive_risk = sum((as.numeric(beta_hat_list[[K_sim+1]]) - beta_mat[, K_sim+1])^2),
                                             base_risk = sum(beta_mat[,K_sim+1]^2))
                                  
                                }

v_t_list_I <- map(X_cov_list_I, t_w_any)[map_dbl(lambda_grid$p_sim, 
                                   ~{which(unique(lambda_grid$p_sim) == .)})]

p_grid_asymp <- expand.grid(p_sim = c(250, 500, 750))
gamma_asymp_oldV <- foreach(p_sim = p_grid_asymp$p_sim,
                            alpha_sq_sim = rep(alpha_sq_sim, nrow(p_grid_asymp)),
                            n_sim = rep(n_sim, nrow(p_grid_asymp)),
                            rho_sim = rep(rho_sim, nrow(p_grid_asymp)),
                            t_w = v_t_list_I) %do% {
                              
                              gamma_ <- p_sim /n_sim
                              v_m_res <- compute_ST(w = t_w$w, t = t_w$t, gamma = gamma_)
                              
                              res <- pmap_dbl(list(v_m_res$m,
                                                   v_m_res$m_prime,
                                                   v_m_res$lambda_seq), function(m, m_per, lambda_sim){
                                                     
                                                     A_ij <- rho_sim * sigma_sq_sim * alpha_sq_sim * (1 - lambda_sim * m)^2
                                                     A_ii <- sigma_sq_sim * alpha_sq_sim * (1 - 2 * lambda_sim * m + lambda_sim^2 * m_per)
                                                     V_i <- sigma_sq_sim * rho_sim * alpha_sq_sim * (1 - lambda_sim * m)
                                                     V_K <- sigma_sq_sim * alpha_sq_sim * (1 - lambda_sim * m)
                                                     R_ii <- sigma_sq_sim * (gamma_ * m - gamma_ * lambda_sim * m_per)
                                                     
                                                     V <- c(rep(V_i, K_sim), V_K)
                                                     R <- diag(R_ii, K_sim+1)
                                                     A <- matrix(A_ij, nrow = K_sim+1, ncol = K_sim+1)
                                                     diag(A) <- A_ii
                                                     res2 <- (alpha_sq_sim * sigma_sq_sim - t(as.matrix(V)) %*% solve(A + R) %*% as.matrix(V))
                                                     
                                                     return(res2)
                                                   })
                              
                              return(data.frame(lambda = v_m_res$lambda_seq, risk = res))
                            }

gamma_empirical_oldV %>%
  reduce(rbind) %>%
  cbind(lambda_grid) %>%
  mutate(gamma = p_sim / n_sim) %>%
  left_join(bind_rows(map2(p_grid_asymp$p_sim, gamma_asymp_oldV, ~{
    mutate(..2, p_sim = ..1) %>%
      filter(lambda < 30) %>%
      filter(lambda > 0.5)
  }))) %>%
  mutate(gamma = factor(gamma, 
                        levels = c(0.5, 1, 1.5),
                        ordered = TRUE, 
                        labels=c(expression(gamma ~ " = 0.5"), 
                                 expression(gamma ~ " = 1"), 
                                 expression(gamma ~ " = 1.5")))) %>%
  ggplot() +
  geom_boxplot(aes(x = (lambda_sim), y = TL_risk, group = lambda_sim), lwd=0.3, outlier.size = 0.3) +
  geom_line(aes(x = lambda, y = risk), lwd = 0.3, color = "red") +
  facet_grid(~ gamma, labeller = label_parsed) +
  labs(x = expression("lambda ("~lambda~")") , y = "Estimation Risk", color = "") +
  theme_bw(base_size = 10) +
  coord_cartesian(xlim = c(0, 10), expand = FALSE)+ 
  scale_x_continuous(breaks=seq(0, 10, 2.5),labels=c("0", "2.5", "5", "7.5", "10"))  +
  theme(panel.spacing = unit(1, "lines"))

############################
#Identity Covariance Matrix#
############################
set.seed(111)
n_vec <- seq(750, 250, length.out = K_sim+1)
lambda_grid_uneven <- expand.grid(lambda_sim = seq(0.5, 11, 1),
                                  ite = 1:50)

X_cov_list_I_uneven <- map(c(750), ~{diag(..1)}) 
X_cov_150 <- map(1:nrow(lambda_grid_uneven), ~{X_cov_list_I_uneven[[1]]})

gamma_empirical_uneven <- foreach(lambda_sim = lambda_grid_uneven$lambda_sim,
                                  X_cov = X_cov_150) %do%{
                                    print(p_sim)
                                    
                                    rho_mat <- matrix(alpha_sq_sim * sigma_sq_sim * rho_sim/p_sim, nrow = K_sim+1, ncol = K_sim+1)
                                    diag(rho_mat) <- alpha_sq_sim * sigma_sq_sim / p_sim
                                    beta_mat <- MASS::mvrnorm(p_sim, rep(0, K_sim+1), Sigma = rho_mat)
                                    
                                    X_list <- map2(1:(K_sim+1), n_vec, ~{MASS::mvrnorm(.y, rep(0, p_sim), Sigma = X_cov)})
                                    
                                    Y_list <- map2(lapply(seq_len(ncol(beta_mat)), function(i) beta_mat[,i]), 
                                                   X_list, 
                                                   ~{..2 %*% as.matrix(..1) + rnorm(nrow(..2))})
                                    
                                    #######################  #######################  #######################  #######################
                                    beta_hat_list <- map2(X_list, Y_list, ~{
                                      solve(t(..1) %*% ..1 + nrow(..1)/ (nrow(..1) / n_vec[K_sim+1]) * lambda_sim * diag(p_sim)) %*% t(..1) %*% ..2
                                    })
                                    #######################  #######################  #######################  #######################
                                    
                                    resolvent_list <-  map2(X_list, Y_list, ~{
                                      solve(t(..1) %*% ..1 / nrow(..1) + 1 / (nrow(..1) / n_vec[K_sim+1]) * lambda_sim * diag(p_sim))
                                    })
                                    
                                    m_vec <- map_dbl(resolvent_list, ~{sum(diag(.x)/p_sim)})
                                    m_per_vec <- map_dbl(resolvent_list, ~{sum(diag(.x %*% .x)/p_sim)})
                                    gamma_vec <- p_sim/(n_vec)
                                    lambda_vec <- lambda_sim / (n_vec / n_vec[K_sim+1])
                                    
                                    cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                    diag(cov_mat_beta) <- 1
                                    
                                    cross_t <- outer(1:(K_sim+1), 1:(K_sim+1), function(x, y){
                                      # 1  + (lambda_vec[x] - lambda_vec[y]) *
                                      #   (m_vec[x] * m_vec[y] / (m_vec[x] - m_vec[y]))
                                      (1 - lambda_vec[x] * m_vec[x]) * (1 - lambda_vec[y] * m_vec[y])
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
                                    
                                    beta_pool <- map2(beta_hat_list, res, ~{..1 * ..2}) %>% 
                                      Reduce( '+', .)
                                    
                                    data.frame(TL_risk= sum((as.numeric(beta_pool) - beta_mat[,K_sim+1])^2),
                                               naive_risk = sum((as.numeric(beta_hat_list[[K_sim+1]]) - beta_mat[, K_sim+1])^2),
                                               base_risk = sum(beta_mat[,K_sim+1]^2))
                                    
                                  }

X_cov_list_I_uneven <- map(c(750), ~{diag(..1)}) 
v_t_list_I_aaa <- map(X_cov_list_I_uneven, t_w_any)
v_m_res_matched_list <- map(v_t_list_I_aaa, ~{
  t_w <- .x
  
  lam_ex <- map(p_sim/n_vec, ~{
    compute_ST(w = t_w$w, t = t_w$t, gamma = .x, grid_size = 1e5)
  }) %>%
    purrr::transpose()
  
  close_indx <- map2(lam_ex$lambda_seq[-(K_sim + 1)],
                     n_vec[1:K_sim] / n_vec[K_sim+1],
                     function(x1, x2){
                       map_dbl(lam_ex$lambda_seq[[(K_sim + 1)]], ~{
                         which.min(abs(x1 - (.)/x2))
                       })
                     })
  
  close_indx[[(K_sim + 1)]] <- 1:length(lam_ex$lambda_seq[[(K_sim + 1)]])
  
  map(lam_ex, function(l1){
    map2(l1, close_indx, ~{.x[.y]})
  }) %>%
    map(purrr::transpose)
  
})

gamma_asymp_uneven <- foreach(alpha_sq_sim = rep(alpha_sq_sim, length(v_m_res_matched_list)),
                              rho_sim = rep(rho_sim, length(v_m_res_matched_list)),
                              v_m_res_matched = v_m_res_matched_list) %do% {
                                
                                gamma_vec <- p_sim /n_vec
                                
                                res <- pmap_dbl(list(v_m_res_matched$m,
                                                     v_m_res_matched$m_prime,
                                                     v_m_res_matched$lambda_seq), function(m, m_per, lambda_sim){
                                                       
                                                       m_vec <- unlist(m)
                                                       m_per_vec <- unlist(m_per)
                                                       lambda_vec <- unlist(lambda_sim)
                                                       
                                                       cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                                       diag(cov_mat_beta) <- 1
                                                       
                                                       cross_t <- outer(1:(K_sim+1), 1:(K_sim+1), function(x, y){
                                                         #   1 + (lambda_vec[x] - lambda_vec[y]) *
                                                         #     (m_vec[x] * m_vec[y] / (m_vec[x] - m_vec[y]))
                                                         (1 - lambda_vec[x] * m_vec[x]) * (1 - lambda_vec[y] * m_vec[y])
                                                       })
                                                       diag(cross_t) <- 0
                                                       
                                                       V <- alpha_sq_sim * sigma_sq_sim * cov_mat_beta[, K_sim+1] *  (1 - lambda_vec * m_vec)
                                                       A_ij <- cov_mat_beta * sigma_sq_sim * alpha_sq_sim * cross_t
                                                       A_ii <- sigma_sq_sim * alpha_sq_sim * (1 - 2 * lambda_vec * m_vec + lambda_vec^2 * m_per_vec)
                                                       R_ii <- sigma_sq_sim * (gamma_vec * m_vec - gamma_vec * lambda_vec * m_per_vec)
                                                       diag(A_ij) <- A_ii
                                                       R <- diag(R_ii)
                                                       A <- A_ij
                                                       
                                                       res2 <- alpha_sq_sim * sigma_sq_sim - t(as.matrix(V)) %*% solve(A + R) %*% as.matrix(V)
                                                       
                                                       return(res2)
                                                     })
                                
                                return(data.frame(lambda = map_dbl(v_m_res_matched$lambda_seq, ~{.[[K_sim + 1]]}), 
                                                  risk = res))
                              }


gamma_empirical_uneven %>%
  reduce(rbind) %>%
  cbind(lambda_grid_uneven) %>%
  ggplot() +
  geom_boxplot(aes(x = lambda_sim, y = TL_risk, group = lambda_sim), lwd=0.3, outlier.size = 0.3) +
  geom_line(data = gamma_asymp_uneven[[1]] %>%
              filter(lambda <= 40) %>%
              filter(lambda >= 0.5)
            , aes(x = lambda, y = risk), lwd=0.3, color = "red") +
  labs(x = expression("lambda ("~lambda~")") , y = "Estimation Risk", color = "") +
  theme_bw(base_size = 10)+
  coord_cartesian(xlim = c(0, 10), expand = FALSE) + 
  scale_x_continuous(breaks = seq(0, 10, 2.5),labels=c("0", "2.5", "5", "7.5", "10"))

###########################
#General Covariance Matrix#
###########################

set.seed(111)
p_sim_try <- 750
X_cov_exp <- exp_decay_cov(p_sim_try, width_actual = 101, c = 0.07)
#X_cov_exp <- exp_decay_cov(p_sim_try)
X_cov_150_exp <- map(1:nrow(lambda_grid_uneven), ~{X_cov_exp})
n_vec_try <- seq(750, 250, length.out = K_sim+1)

gamma_empirical_uneven_exp <- foreach(lambda_sim = lambda_grid_uneven$lambda_sim,
                                      X_cov = X_cov_150_exp) %do%{
                                        print(p_sim_try)
                                        
                                        rho_mat <- matrix(alpha_sq_sim * sigma_sq_sim * rho_sim/p_sim_try, nrow = K_sim+1, ncol = K_sim+1)
                                        diag(rho_mat) <- alpha_sq_sim * sigma_sq_sim / p_sim_try
                                        beta_mat <- MASS::mvrnorm(p_sim_try, rep(0, K_sim+1), Sigma = rho_mat)
                                        
                                        X_list <- map2(1:(K_sim+1), n_vec_try, ~{MASS::mvrnorm(.y, rep(0, p_sim_try), Sigma = X_cov)})
                                        
                                        Y_list <- map2(lapply(seq_len(ncol(beta_mat)), function(i) beta_mat[,i]), 
                                                       X_list, 
                                                       ~{..2 %*% as.matrix(..1) + rnorm(nrow(..2))})
                                        
                                        #######################  #######################  #######################  #######################
                                        beta_hat_list <- map2(X_list, Y_list, ~{
                                          solve(t(..1) %*% ..1 + nrow(..1)/ (nrow(..1) / n_vec_try[K_sim+1]) * lambda_sim * diag(p_sim_try)) %*% t(..1) %*% ..2
                                        })
                                        #######################  #######################  #######################  #######################
                                        
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
                                        
                                        beta_pool <- map2(beta_hat_list, res, ~{..1 * ..2}) %>% 
                                          Reduce( '+', .)
                                        
                                        data.frame(TL_risk= sum((as.numeric(beta_pool) - beta_mat[,K_sim+1])^2),
                                                   naive_risk = sum((as.numeric(beta_hat_list[[K_sim+1]]) - beta_mat[, K_sim+1])^2),
                                                   base_risk = sum(beta_mat[,K_sim+1]^2))
                                        
                                      }

v_t_list_I_exp <- map(list(X_cov_exp), t_w_any)
v_m_res_matched_list_exp_try <- map(v_t_list_I_exp, ~{
  t_w <- .x
  
  lam_ex <- map(p_sim_try/n_vec_try, ~{
    compute_ST(w = t_w$w, t = t_w$t, gamma = .x, grid_size = 1e5)
  }) %>%
    purrr::transpose()
  
  close_indx <- map2(lam_ex$lambda_seq[-(K_sim + 1)],
                     n_vec_try[1:K_sim] / n_vec_try[K_sim+1],
                     function(x1, x2){
                       map_dbl(lam_ex$lambda_seq[[(K_sim + 1)]], ~{
                         which.min(abs(x1 - (.)/x2))
                       })
                     })
  
  close_indx[[(K_sim + 1)]] <- 1:length(lam_ex$lambda_seq[[(K_sim + 1)]])
  
  map(lam_ex, function(l1){
    map2(l1, close_indx, ~{.x[.y]})
  }) %>%
    map(purrr::transpose)
  
})

gamma_asymp_uneven_exp <- foreach(alpha_sq_sim = rep(alpha_sq_sim, length(v_m_res_matched_list_exp_try)),
                                  rho_sim = rep(rho_sim, length(v_m_res_matched_list_exp_try)),
                                  v_m_res_matched = v_m_res_matched_list_exp_try) %do% {
                                    
                                    gamma_vec <- p_sim_try /n_vec_try
                                    
                                    res <- pmap_dbl(list(v_m_res_matched$m,
                                                         v_m_res_matched$m_prime,
                                                         v_m_res_matched$lambda_seq,
                                                         v_m_res_matched$v), function(m, m_per, lambda_sim, v){
                                                           
                                                           m_vec <- unlist(m)
                                                           m_per_vec <- unlist(m_per)
                                                           lambda_vec <- unlist(lambda_sim)
                                                           
                                                           cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                                           diag(cov_mat_beta) <- 1
                                                           
                                                           cross_t <- outer(1:(K_sim+1), 1:(K_sim+1), function(x, y){

                                                             1 +  (lambda_vec[x] - lambda_vec[y]) *
                                                               (m_vec[x] * m_vec[y] / (m_vec[x] - m_vec[y]))
                                                          
                                                           })
                                                           diag(cross_t) <- 0
                                                           
                                                           V <- alpha_sq_sim * sigma_sq_sim * cov_mat_beta[, K_sim+1] *  (1 - lambda_vec * m_vec)
                                                           A_ij <- cov_mat_beta * cross_t * sigma_sq_sim * alpha_sq_sim 
                                                           A_ii <- sigma_sq_sim * alpha_sq_sim * (1 - 2 * lambda_vec * m_vec + lambda_vec^2 * m_per_vec)
                                                           R_ii <- sigma_sq_sim * (gamma_vec * m_vec - gamma_vec * lambda_vec * m_per_vec)
                                                           diag(A_ij) <- A_ii
                                                           R <- diag(R_ii)
                                                           A <- A_ij
                                                           
                                                           res2 <- alpha_sq_sim * sigma_sq_sim - t(as.matrix(V)) %*% solve(A + R) %*% as.matrix(V)
                                                           
                                                           return(res2)
                                                         })
                                    
                                    return(data.frame(lambda = map_dbl(v_m_res_matched$lambda_seq, ~{.[[K_sim + 1]]}), 
                                                      risk = res))
                                  }



gamma_empirical_uneven_exp %>%
  reduce(rbind) %>%
  cbind(lambda_grid_uneven) %>%
  ggplot() +
  geom_boxplot(aes(x = lambda_sim, y = TL_risk, group = lambda_sim), lwd=0.3, outlier.size = 0.3) +
  geom_line(data = gamma_asymp_uneven_exp[[1]] %>%
              filter(lambda <= 30) %>%
              filter(lambda >= 0.5)
            , aes(x = lambda, y = risk), lwd=0.3, color = "red") +
  labs(x = expression("lambda ("~lambda~")") , y = "Estimation Risk", color = "") +
  theme_bw(base_size = 10) +
  coord_cartesian(xlim = c(0, 10), expand = FALSE)+ 
  scale_x_continuous(breaks = seq(0, 10, 2.5),labels=c("0", "2.5", "5", "7.5", "10"))

################################################################################
######################    Prediction Risk       ################################
################################################################################

########################################
#Same number of Observations Covariance#
########################################
set.seed(111)
lambda_grid <- expand.grid(lambda_sim = seq(0.5, 11, 1),
                           p_sim = c(250, 500, 750),
                           ite = 1:50)

X_cov_list_I_uneven <- map(c(250, 500, 750), ~{diag(..1)}) 
v_t_list_I_aaa <- map(X_cov_list_I_uneven, t_w_any)
pred_empirical_new <- foreach(lambda_sim = lambda_grid$lambda_sim,
                              p_sim = lambda_grid$p_sim,
                              X_cov = X_cov_foreach_I) %do%{
                                print(p_sim)
                                
                                gamma_ <- p_sim / n_sim
                                rho_mat <- matrix(alpha_sq_sim * sigma_sq_sim * rho_sim/p_sim, nrow = K_sim+1, ncol = K_sim+1)
                                diag(rho_mat) <- alpha_sq_sim * sigma_sq_sim / p_sim
                                beta_mat <- MASS::mvrnorm(p_sim, rep(0, K_sim+1), Sigma = rho_mat)
                                
                                X_list <- map(1:(K_sim+1), ~{MASS::mvrnorm(n_sim, rep(0, p_sim), Sigma = X_cov)})
                                
                                # X_list <- map(1:(K_sim+1), ~{
                                #   matrix(rnorm(n_sim * p_sim), n_sim, p_sim) %*% expm::sqrtm(X_cov)
                                #   })
                                
                                Y_list <- map2(lapply(seq_len(ncol(beta_mat)), function(i) beta_mat[,i]), 
                                               X_list, 
                                               ~{..2 %*% as.matrix(..1) + rnorm(n_sim)})
                                
                                #######################  #######################  #######################  #######################
                                beta_hat_list <- map2(X_list, Y_list, ~{
                                  solve(t(..1) %*% ..1 + n_sim * lambda_sim * diag(p_sim)) %*% t(..1) %*% ..2
                                })
                                #######################  #######################  #######################  #######################
                                
                                #beta_naive_hat <- solve(t(X_list[[K_sim+1]]) %*% X_list[[K_sim+1]] + n_sim * lambda_sim * diag(p_sim)) %*% t(X_list[[K_sim+1]]) %*% Y_list[[K_sim+1]]
                                
                                
                                Sigma_hat <- t(X_list[[K_sim+1]]) %*% X_list[[K_sim+1]] / n_sim
                                resolvent <- solve(Sigma_hat + diag(lambda_sim, p_sim))
                                m <- sum(diag(resolvent)/p_sim)
                                m_per <- sum(diag(resolvent %*% resolvent))/ p_sim
                                v <- gamma_ * (m - 1 / lambda_sim) + (1 / lambda_sim)
                                v_per <- gamma_ * (m_per - 1 / lambda_sim^2) + (1 / lambda_sim^2)
                                
                                order_0 <- 1
                                order_1 <- (1/gamma_) *  (1 /(v*lambda_sim) - 1)
                                order_2 <- (1/gamma_) * (v - lambda_sim * v_per) / (lambda_sim * v)^2
                                order_11 <- (m - lambda_sim * m_per) / 
                                  (1 - gamma_ + gamma_ * lambda_sim^2 * m_per)
                                
                                cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                diag(cov_mat_beta) <- 1
                                
                                C_ij <- matrix(order_0 - 2 * lambda_sim * order_1 + lambda_sim^2 * order_11, K_sim + 1, K_sim + 1)
                                diag(C_ij) <- (order_0 - 2 * lambda_sim * order_1 + lambda_sim^2 * order_2)
                                D_i <- matrix(order_0 - lambda_sim * order_1, K_sim + 1, 1)
                                F_i <- diag((order_1 - lambda_sim* order_2),  K_sim + 1)
                                
                                D_f <- D_i * cov_mat_beta[, (K_sim+1)]* alpha_sq_sim * sigma_sq_sim
                                F_f <- gamma_ * sigma_sq_sim * F_i
                                C_f <- C_ij * cov_mat_beta * alpha_sq_sim * sigma_sq_sim
                                
                                weight_temp <- as.numeric(solve(C_f + F_f) %*% D_f)
                                
                                beta_pool <- map2(beta_hat_list, weight_temp, ~{..1 * ..2}) %>% 
                                  Reduce( '+', .)
                                
                                X_test <- MASS::mvrnorm(1000, rep(0, p_sim), Sigma = X_cov)
                                y_test <- X_test %*% beta_mat[,K_sim+1] + rnorm(1000)
                                y_hat_TL <- X_test %*% beta_pool
                                y_hat_naive <- X_test %*% beta_hat_list[[K_sim+1]]
                                
                                
                                data.frame(TL_risk= mean((y_hat_TL - y_test)^2),
                                           naive_risk = mean((as.numeric(y_hat_naive - y_test)^2)),
                                           base_risk = mean(y_test^2))
                                
                              }

p_grid_asymp <- expand.grid(p_sim = c(250, 500, 750))

pred_asymp_opt <- foreach(p_sim = p_grid_asymp$p_sim,
                          alpha_sq_sim = rep(alpha_sq_sim, nrow(p_grid_asymp)),
                          n_sim = rep(n_sim, nrow(p_grid_asymp)),
                          rho_sim = rep(rho_sim, nrow(p_grid_asymp)),
                          t_w = v_t_list_I_aaa) %do% {
                            
                            gamma_ <- p_sim /n_sim
                            v_m_res <- compute_ST(w = t_w$w, t = t_w$t, gamma = gamma_)
                            
                            res <- pmap_dbl(list(v_m_res$m,
                                                 v_m_res$m_prime,
                                                 v_m_res$v,
                                                 v_m_res$v_prime,
                                                 v_m_res$lambda_seq), function(m, m_per, v, v_per, lambda_sim){
                                                   
                                                   order_0 <- 1
                                                   order_1 <- (1/gamma_) *  (1 /(v*lambda_sim) - 1)
                                                   order_2 <- (1/gamma_) * (v - lambda_sim * v_per) / (lambda_sim * v)^2
                                                   order_11 <- (m - lambda_sim * m_per) / 
                                                     (1 - gamma_ + gamma_ * lambda_sim^2 * m_per)
                                                   
                                                   cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                                   diag(cov_mat_beta) <- 1
                                                   
                                                   C_ij <- matrix(order_0 - 2 * lambda_sim * order_1 + lambda_sim^2 * order_11, K_sim + 1, K_sim + 1)
                                                   diag(C_ij) <- (order_0 - 2 * lambda_sim * order_1 + lambda_sim^2 * order_2)
                                                   D_i <- matrix(order_0 - lambda_sim * order_1, K_sim + 1, 1)
                                                   F_i <- diag((order_1 - lambda_sim* order_2),  K_sim + 1)
                                                   
                                                   D_f <- D_i * cov_mat_beta[, (K_sim+1)]* alpha_sq_sim * sigma_sq_sim
                                                   F_f <- gamma_ * sigma_sq_sim * F_i
                                                   C_f <- C_ij * cov_mat_beta * alpha_sq_sim * sigma_sq_sim
                                                   
                                                   res1 <- alpha_sq_sim * sigma_sq_sim - t(D_f) %*% solve(F_f + C_f) %*% D_f + order_0
                                                   return(res1)
                                                 })
                            
                            return(data.frame(lambda = v_m_res$lambda_seq, risk = res))
                          }


pred_empirical_new %>%
  reduce(rbind) %>%
  cbind(lambda_grid) %>%
  filter(lambda_sim <= 20) %>%
  mutate(gamma = p_sim / n_sim) %>%
  left_join(bind_rows(map2(p_grid_asymp$p_sim, pred_asymp_opt, ~{
    mutate(..2, p_sim = ..1) %>%
      filter(lambda < 15) %>%
      filter(lambda > 0.5)
  }))) %>%
  mutate(gamma = factor(gamma, 
                        levels = c(0.5, 1, 1.5),
                        ordered = TRUE, 
                        labels=c(expression(gamma ~ " = 0.5"), 
                                 expression(gamma ~ " = 1"), 
                                 expression(gamma ~ " = 1.5")))) %>%
  ggplot() +
  geom_boxplot(aes(x = (lambda_sim), y = TL_risk, group = lambda_sim), lwd=0.3, outlier.size = 0.3) +
  geom_line(aes(x = lambda, y = risk), lwd=0.3, color = "red") +
  facet_grid(~ gamma, labeller = label_parsed) +
  labs(x = expression("lambda ("~lambda~")") , y = "Prediction Risk", color = "") +
  theme_bw(base_size = 10) +
  ylim(1.1, 2)+
  coord_cartesian(xlim = c(0, 10), expand = FALSE) + 
  scale_x_continuous(breaks=seq(0, 10, 2.5),labels=c("0", "2.5", "5", "7.5", "10")) +
  theme(panel.spacing = unit(1, "lines"))

############################
#Identity Covariance Matrix#
############################
set.seed(111)
p_sim <- 750
lambda_grid_uneven <- expand.grid(lambda_sim = seq(0.5, 11, 1),
                                  ite = 1:50)
X_cov_list_I_uneven <- map(c(p_sim), ~{diag(..1)}) 
X_cov_150 <- map(1:nrow(lambda_grid_uneven), ~{X_cov_list_I_uneven[[1]]})

pred_empirical_I_opt <- foreach(lambda_sim = lambda_grid_uneven$lambda_sim,
                                X_cov = X_cov_150) %do%{
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
                                        lambda_use[x]* m[x]* m[y] * lambda_use[y]
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
                                  
                                  X_test <- MASS::mvrnorm(1000, rep(0, p_sim), Sigma = X_cov)
                                  
                                  y_test <- X_test %*% beta_mat[,K_sim+1] + rnorm(1000)
                                  y_hat_TL <- X_test %*% beta_pool
                                  y_hat_naive <- X_test %*% beta_hat_list[[K_sim+1]]
                                  
                                  
                                  data.frame(TL_risk= mean((y_hat_TL - y_test)^2),
                                             naive_risk = mean((as.numeric(y_hat_naive - y_test)^2)),
                                             base_risk = mean(y_test^2))
                                  
                                }


v_m_res_matched_fixed <- map(list(t_w_any(diag(p_sim))), ~{
  t_w <- .x
  
  lam_ex <- map(p_sim/n_vec, ~{
    compute_ST(w = t_w$w, t = t_w$t, gamma = .x, grid_size = 1e5)
  }) %>%
    purrr::transpose()
  
  close_indx <- map2(lam_ex$lambda_seq[-(K_sim + 1)],
                     n_vec[1:K_sim] / n_vec[K_sim+1],
                     function(x1, x2){
                       map_dbl(lam_ex$lambda_seq[[(K_sim + 1)]], ~{
                         which.min(abs(x1 - (.)/x2))
                       })
                     })
  
  close_indx[[(K_sim + 1)]] <- 1:length(lam_ex$lambda_seq[[(K_sim + 1)]])
  
  map(lam_ex, function(l1){
    map2(l1, close_indx, ~{.x[.y]})
  }) %>%
    map(purrr::transpose)
  
})

pred_asymp_I_opt <- foreach(alpha_sq_sim = 1,
                            rho_sim = rho_sim,
                            v_m_res_matched = v_m_res_matched_fixed) %do% {
                              
                              gamma_vec <- p_sim / n_vec
                              
                              res <- pmap_dbl(list(v_m_res_matched$m,
                                                   v_m_res_matched$m_prime,
                                                   v_m_res_matched$v,
                                                   v_m_res_matched$v_prime,
                                                   v_m_res_matched$lambda_seq), function(m, m_per, v, v_per, lambda_sim){
                                                     
                                                     m <- unlist(m)
                                                     m_per <- unlist(m_per)
                                                     v <- unlist(v)
                                                     v_per <- unlist(v_per)
                                                     lambda_sim <- unlist(lambda_sim)
                                                     
                                                     order_0 <- 1
                                                     order_1 <- (1/gamma_vec) *  (1 /(v*lambda_sim) - 1)
                                                     order_2 <- (1/gamma_vec) * (v - lambda_sim * v_per) / (lambda_sim * v)^2
                                                     order_11 <- outer(1:(K_sim+1), 1:(K_sim+1), function(x, y){
                                                       lambda_sim[x] * m[x] * lambda_sim[y] * m[y]
                                                     })
                                                     
                                                     diag(order_11) <- 0
                                                     
                                                     cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                                     diag(cov_mat_beta) <- 1
                                                     
                                                     C_ij <- matrix(1, K_sim + 1, K_sim+1)
                                                     
                                                     for (x in 1:(K_sim + 1)) {
                                                       for (y in 1:(K_sim + 1)) {
                                                         C_ij[x, y] <- order_0 - lambda_sim[x] * order_1[x] - lambda_sim[y] * order_1[y] + 
                                                           lambda_sim[x]* m[x]* m[y] * lambda_sim[y]
                                                       }
                                                     }
                                                     
                                                     diag(C_ij) <- (order_0 - 2 * lambda_sim * order_1 + lambda_sim^2 * order_2)
                                                     
                                                     D_i <- matrix(order_0 - lambda_sim * order_1, K_sim + 1, 1)
                                                     F_i <- diag(order_1 - lambda_sim* order_2)
                                                     
                                                     D_f <- D_i * cov_mat_beta[, (K_sim+1)]* alpha_sq_sim * sigma_sq_sim
                                                     F_f <- gamma_vec * sigma_sq_sim * F_i
                                                     C_f <- C_ij * cov_mat_beta * alpha_sq_sim * sigma_sq_sim
                                                     
                                                     res1 <- order_0 * sigma_sq_sim * alpha_sq_sim -
                                                       t(D_f) %*% 
                                                       solve(C_f + F_f) %*% 
                                                       D_f +
                                                       sigma_sq_sim
                                                     
                                                     return(res1)
                                                   })
                              
                              return(data.frame(lambda = map_dbl(v_m_res_matched$lambda_seq, ~{.[[K_sim + 1]]}),
                                                risk = res))
                            }


pred_empirical_I_opt %>%
  reduce(rbind) %>%
  cbind(lambda_grid_uneven) %>%
  ggplot() +
  geom_boxplot(aes(x = lambda_sim, y = TL_risk, group = lambda_sim), lwd=0.3, outlier.size = 0.3) +
  geom_line(data = pred_asymp_I_opt[[1]] %>%
              filter(lambda <= 30) %>%
              filter(lambda >= 0.4)
            , aes(x = lambda, y = risk), lwd=0.3, color = "red") +
  #this is from uneven n Sigma
  labs(x = expression("lambda ("~lambda~")") , y = "Prediction Risk", color = "") +
  theme_bw(base_size = 10) +
  ylim(1.4, 2.0) +
  coord_cartesian(xlim = c(0, 10), expand = FALSE) + 
  scale_x_continuous(breaks=seq(0, 10, 2.5),labels=c("0", "2.5", "5", "7.5", "10")) 


###########################
#General Covariance Matrix#
###########################
lambda_grid_nonI <- expand.grid(lambda_sim = seq(0.5, 11, 1),
                                ite = 1:50)
p_sim <- 750
n_vec <- seq(750, 250, length.out = K_sim+1)
X_cov_list <- map(c(p_sim), exp_decay_cov, width_actual = 701, c = 0.07)
set.seed(111)
pred_empirical_opt_nonI <- foreach(lambda_sim = lambda_grid_nonI$lambda_sim,
                                   X_cov = map(1:nrow(lambda_grid_nonI), 
                                               ~{X_cov_list[[1]]})) %do%{
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
                                                     # C_ij[x, y] <- order_0 - lambda_use[x] * order_1[x] - lambda_use[y] * order_1[y] + 
                                                     #   lambda_use[x]* m[x]* m[y] * lambda_use[y]
                                                     # C_ij[x, y] <- order_0 - lambda_use[x] * order_1[x] - lambda_use[y] * order_1[y] +
                                                     #   (lambda_use[x] * m[x] - lambda_use[y] * m[y]) / (lambda_use[x] * lambda_use[y] * (m[y] - m[x]))
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
                                                 
                                                 X_test <- MASS::mvrnorm(1000, rep(0, p_sim), Sigma = X_cov)
                                                 
                                                 y_test <- X_test %*% beta_mat[,K_sim+1] + rnorm(1000)
                                                 y_hat_TL <- X_test %*% beta_pool
                                                 y_hat_naive <- X_test %*% beta_hat_list[[K_sim+1]]
                                                 
                                                 
                                                 data.frame(TL_risk= mean((y_hat_TL - y_test)^2),
                                                            naive_risk = mean((as.numeric(y_hat_naive - y_test)^2)),
                                                            base_risk = mean(y_test^2))
                                                 
                                               }

v_t_list_I_exp <- map(X_cov_list, t_w_any)
v_m_res_matched_exp <- map(v_t_list_I_exp, ~{
  t_w <- .x
  lam_ex <- map(p_sim/n_vec, ~{
    compute_ST(w = t_w$w, t = t_w$t, gamma = .x, grid_size = 1e5)
  }) %>%
    purrr::transpose()
  
  close_indx <- map2(lam_ex$lambda_seq[-(K_sim + 1)],
                     n_vec[1:K_sim] / n_vec[K_sim+1],
                     function(x1, x2){
                       map_dbl(lam_ex$lambda_seq[[(K_sim + 1)]], ~{
                         which.min(abs(x1 - (.)/x2))
                       })
                     })
  
  close_indx[[(K_sim + 1)]] <- 1:length(lam_ex$lambda_seq[[(K_sim + 1)]])
  
  map(lam_ex, function(l1){
    map2(l1, close_indx, ~{.x[.y]})
  }) %>%
    map(purrr::transpose)
  
})


pred_asymp_opt_nonI <- foreach(alpha_sq_sim = rep(alpha_sq_sim, length(v_m_res_matched_exp)),
                               rho_sim = rep(rho_sim, length(v_m_res_matched_exp)),
                               v_m_res_matched = v_m_res_matched_exp) %do% {
                                 
                                 gamma_vec <- p_sim / n_vec
                                 
                                 res <- pmap_dbl(list(v_m_res_matched$m,
                                                      v_m_res_matched$m_prime,
                                                      v_m_res_matched$v,
                                                      v_m_res_matched$v_prime,
                                                      v_m_res_matched$lambda_seq), function(m, m_per, v, v_per, lambda_sim){
                                                        
                                                        m <- unlist(m)
                                                        m_per <- unlist(m_per)
                                                        v <- unlist(v)
                                                        v_per <- unlist(v_per)
                                                        lambda_sim <- unlist(lambda_sim)
                                                        
                                                        order_0 <- 1
                                                        order_1 <- (1/gamma_vec) *  (1 /(v*lambda_sim) - 1)
                                                        order_2 <- (1/gamma_vec) * (v - lambda_sim * v_per) / (lambda_sim * v)^2
                                                        
                                                        cov_mat_beta <- matrix(rho_sim, K_sim + 1, K_sim + 1)
                                                        diag(cov_mat_beta) <- 1
                                                        
                                                        C_ij <- matrix(1, K_sim + 1, K_sim+1)
                                                        
                                                        for (x in 1:(K_sim + 1)) {
                                                          for (y in 1:(K_sim + 1)) {
                                                            C_ij[x, y] <- order_0 - lambda_sim[x] * order_1[x] - lambda_sim[y] * order_1[y] +
                                                              (lambda_sim[x] * m[x] - lambda_sim[y] * m[y]) / (m[y] - m[x])
 
                                                          }
                                                        }
                                                        
                                                        diag(C_ij) <- (order_0 - 2 * lambda_sim * order_1 + lambda_sim^2 * order_2)
                                                        
                                                        D_i <- matrix(order_0 - lambda_sim * order_1, K_sim + 1, 1)
                                                        F_i <- diag(order_1 - lambda_sim * order_2)
                                                        D_f <- D_i * cov_mat_beta[, (K_sim+1)]
                                                        F_f <- gamma_vec * F_i
                                                        C_f <- C_ij * cov_mat_beta
                                                        
                                                        res1 <- order_0 -
                                                          t(D_f) %*%
                                                          solve(C_f + F_f) %*%
                                                          D_f +
                                                          sigma_sq_sim
                                                        
                                                        return(res1)
                                                      })
                                 
                                 return(data.frame(lambda = map_dbl(v_m_res_matched$lambda_seq, ~{.[[K_sim + 1]]}),
                                                   risk = res))
                               }


pred_empirical_opt_nonI %>%
  reduce(rbind) %>%
  cbind(lambda_grid_nonI) %>%
  ggplot() +
  geom_boxplot(aes(x = lambda_sim, y = TL_risk, group = lambda_sim), lwd=0.3, outlier.size = 0.3) +
  #this is from uneven n Sigma
  geom_line(data = pred_asymp_opt_nonI[[1]] %>%
              filter(lambda <= 30) %>%
              filter(lambda >= 0.5)
            , aes(x = lambda, y = risk), lwd=0.3, color = "red") +
  labs(x = expression("lambda ("~lambda~")") , y = "Prediction Risk", color = "") +
  theme_bw(base_size = 10) +
  coord_cartesian(xlim = c(0, 10), expand = FALSE) +
  ylim(1, 1.7)+
  scale_x_continuous(breaks=seq(0, 10, 2.5),labels=c("0", "2.5", "5", "7.5", "10")) 


