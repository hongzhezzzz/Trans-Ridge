############################################################################################
#A few notes for running the real data codes:
#
#- Please set current working directory to the folder where this code file locate
#
#- Please make sure the raw data in downloaded into the same folder
#
#- To switch the target population, you can simply switch the "target_name". The rest
# will be automatically taken care of.
#
#- To speed up computation, you only need to run the Bayesian estimation codes once.
#
#- Please make sure you have the Bayesian estimation results before running the "Naive
#Git. R" file
############################################################################################

library(glmnet)
library(foreach)
library(tidyverse)
library(rstan)
library(pROC)

set.seed(111)

ind.set <- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re, 1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}
# Define the Stan model
stan_model_code <- "
data {
  int<lower=0> K;  // number of models
  int<lower=0> n[K];  // number of observations for each model
  int<lower=0> p;  // number of predictors
  matrix[n[1], p] X[K];  // design matrices
  vector[n[1]] y[K];  // response variables
}
parameters {
  matrix[p, K] beta;  // Coefficients for each model
  vector<lower=0, upper=0.5>[K] alpha;  // Coefficient standard deviations for each model
  vector<lower=0, upper=1>[K] sigma;  // Error standard deviations for each model
  vector<lower=-1, upper=1>[K * (K - 1) / 2] rho_raw;  // Off-diagonal elements of the correlation matrix
}
transformed parameters {
  corr_matrix[K] Rho;
  cov_matrix[K] Sigma_out;
  {
    int idx;
    idx = 1;
    for (k in 1:K) {  // Loop from 1 to K-1 for proper indexing
      for (j in 1:K) {
      if(j > k){
        Rho[k, j] = rho_raw[idx];
        Rho[j, k] = rho_raw[idx];
        idx = idx + 1;
        }
      }
    }
    for (k in 1:K) {
      Rho[k, k] = 1;
    }
    
    for (k in 1:K) {
      for (j in 1:K) {
        Sigma_out[k, j] = Rho[k, j] * alpha[k] * alpha[j] * sigma[k] * sigma[j] / p;
      }
    }
  }
}
model {
  vector[n[1]] mu[K];
  for (k in 1:K) {
    mu[k] = X[k] * beta[, k];  // Mean of y[k]
    y[k] ~ normal(mu[k], sigma[k]);  // Likelihood
  }

  // Priors
  for (j in 1:p) {
    beta[j, ] ~ multi_normal(rep_vector(0, K), Sigma_out);
  }
  
  alpha ~ normal(0, 1);  // Half-normal prior for non-negative values
  sigma ~ normal(0, 1);  // Half-normal prior for non-negative values
  for (i in 1:(K * (K - 1) / 2)) {
    target += log_mix(0.5,
                      normal_lpdf(rho_raw[i] | -0.5, 0.2),
                      normal_lpdf(rho_raw[i] | 0.5, 0.2));
  }
}
"
# Compile the Stan model
stan_model_ <- stan_model(model_code = stan_model_code)

target_name <- c("Zackular", "Zeller", "Baxter")[1]

if(target_name == "Zackular"){
  name_vec <- c("Zeller", "Baxter", "Zackular")
  set_indx <- c(210-83, 698-210, 83)
  load("DF.rda")
  DF.wk <- rbind(DF.wk[83:210, ],
                 DF.wk[211:698, ],
                 DF.wk[1:83, ])
}else if(target_name == "Zeller"){
  name_vec <- c("Baxter", "Zackular", "Zeller")
  set_indx <- c(698-210, 83, 210-83)
  load("DF.rda")
  DF.wk <- rbind(DF.wk[211:698, ],
                 DF.wk[1:83, ],
                 DF.wk[83:210, ])
}else{
  load("DF.rda")
  name_vec <- c("Zackular", "Zeller", "Baxter")
  set_indx <- c(83, 210-83, 698-210)
  
}

DF.wk <- map(1:3, ~{
  ind.set(set_indx, .x) %>%
    sample(length(.), replace = FALSE)
}) %>%
  map(~{DF.wk[.x, ]}) %>%
  bind_rows()

wk_tb <- as_tibble(DF.wk)
group_index <- rep(1:3, set_indx)
wk_tb$DiseaseState <- as.numeric(scale(ifelse(wk_tb$DiseaseState == "CRC", 1, -1)))

#####
target_ind <- 3
K <- 3
ex <- ind.set(set_indx, target_ind)
train_xy <- wk_tb

naive_set <- train_xy %>%
  slice(ex) %>%
  select(-DiseaseState, -gender, -age, -BMI) %>%
  apply(2, function(x) sum(x == 0) / length(x)) %>%
  {names(.)[which(. <= 0.9)]}

naive_set2 <- train_xy %>%
  slice(ind.set(set_indx, 2)) %>%
  select(-DiseaseState, -gender, -age, -BMI) %>%
  apply(2, function(x) sum(x == 0) / length(x)) %>%
  {names(.)[which(. <= 0.9)]}

naive_set1 <- train_xy %>%
  slice(ind.set(set_indx, 1)) %>%
  select(-DiseaseState, -gender, -age, -BMI) %>%
  apply(2, function(x) sum(x == 0) / length(x)) %>%
  {names(.)[which(. <= 0.9)]}

naive_set_f <-intersect(naive_set, naive_set2) %>%
  intersect(naive_set1)

train_x <- select(train_xy, gender:BMI, all_of(naive_set_f))

target_scale <- map(2: ncol(train_x), ~{
  to_b <- pull(train_x, .x)
  list(m = mean(to_b), s = sd(to_b))
}) %>% 
  {purrr::transpose(.)}

for (i in 2:ncol(train_x)) {
  train_x[, i] <- (train_x[, i] - target_scale$m[[i-1]]) / target_scale$s[[i-1]]
}

train_x <- model.matrix(~., data=train_x)
train_y <- train_xy$DiseaseState

##########################################################################
#################Estimating Variance and Correlation######################
##########################################################################
X_list <- map(1:3, ~{
  ind.set(set_indx, .x)
})

x_train_list <- map(X_list, ~{train_x[.x, ]})
y_train_list <- map(X_list, ~{train_y[.x]})
n_vec <- map_dbl(x_train_list, nrow)

pad_vec <- max(n_vec) - n_vec + 1
X_padded <- map2(x_train_list, pad_vec, ~{scale(rbind(.x, matrix(0, .y, ncol(.x))))})
y_padded <- map2(y_train_list, pad_vec, ~{as.numeric(scale(c(.x, rep(0, .y))))})

micro_list <- list(K = 3,
                   n = rep(max(n_vec) + 1, K),
                   p = ncol(x_train_list[[1]]),
                   X = X_padded,
                   y = y_padded)

start_ <- Sys.time()
rho_fit <- sampling(stan_model_, data = micro_list, iter = 12000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99, max_treedepth = 16))
end_ <- Sys.time()

rho_mat <- get_posterior_mean(rho_fit,par=c("Rho")) %>%
  apply(1, mean) %>%
  matrix(3,3)

colnames(rho_mat) <- name_vec
rownames(rho_mat) <- name_vec

alpha_sigma <- get_posterior_mean(rho_fit, par=c("alpha")) %>%
  apply(1, mean) *
  get_posterior_mean(rho_fit,par=c("sigma")) %>%
  apply(1, mean)

names(alpha_sigma) <- name_vec

alpha_ <- get_posterior_mean(rho_fit, par=c("alpha")) %>%
  apply(1, mean)

sigma_ <- get_posterior_mean(rho_fit,par=c("sigma")) %>%
  apply(1, mean)

names(alpha_) <- name_vec
names(sigma_) <- name_vec

write.table(alpha_sigma, file="alpha_sigma.Rdata")
write.table(rho_mat, file="rho_mat.Rdata")
write.table(alpha_, file="alpha_.Rdata")
write.table(sigma_, file="sigma_.Rdata")
##########################################################
#############################end##########################
##########################################################
#prediction weight
train_cv_list <- ind.set(set_indx, 3)
lambda_factor_vector <- 10^seq(-4, 4)
score_list <- list()

for (index in train_cv_list) {
  X_list <- map(1:3, ~{
    ind.set(set_indx, .x)
  })
  
  X_list[[3]] <- X_list[[3]][-which(index == X_list[[3]])]
  
  x_train_list <- map(X_list, ~{train_x[.x, ]})
  y_train_list <- map(X_list, ~{train_y[.x]})
  
  x_test_i <- train_x[index,]
  y_test_i <- train_y[index]
  
  n_vec <- map_dbl(x_train_list, nrow)
  
  sigma_ <-  read.table("sigma_.Rdata") %>% 
    {setNames(pull(., 1), rownames(.))}
  
  sigma_ <- sigma_[name_vec]
  
  alpha_ <-  read.table("alpha_.Rdata") %>% 
    {setNames(pull(., 1), rownames(.))}
  
  alpha_ <- alpha_[name_vec]
  
  alpha_sigma <- read.table("alpha_sigma.Rdata") %>% 
    {setNames(pull(., 1), rownames(.))}
  
  alpha_sigma <- alpha_sigma[name_vec]
  
  rho_mat <- read.table("rho_mat.Rdata") %>% as.matrix() %>% 
    .[name_vec, name_vec]
  
  gamma_ <- n_vec / ncol(x_train_list[[1]])
  
  train_XY <- map2(x_train_list, y_train_list, ~{
    list(.x, .y)
  })
  
  score_vec <- c()
  acc_vec <- c()
  #model_list <- list()
  #for (lambda_factor in c(0.01)) {
  for (lambda_factor in lambda_factor_vector) {
    lambda_snips <- gamma_ / alpha_^2
    lambda_snips <- lambda_factor * lambda_snips
    
    
    inverse_list_p <- foreach(X_Y = train_XY,
                              lambda = lambda_snips) %do% {
                                train_x_ <- X_Y[[1]]
                                
                                temp_dgC <- (t(train_x_) %*% train_x_/nrow(train_x_) +  diag(lambda, ncol(train_x_))) %>%
                                  solve()
                                
                                
                                return(temp_dgC)
                              }
    
    m_K_vec <- map_dbl(inverse_list_p, ~{mean(diag(.x))})
    m_K_prime_vec <- map_dbl(inverse_list_p, ~{mean(diag(.x)^2)})
    v_K_vec <- gamma_ * (m_K_vec - 1 / lambda_snips) + (1 / lambda_snips)
    v_K_prime_vec <- gamma_ * (m_K_prime_vec - 1 / lambda_snips^2) + (1 / lambda_snips^2)
    
    order_0 <- 1
    order_1 <- (1/gamma_) *  (1 /(v_K_vec*lambda_snips) - 1)
    order_2 <- (1/gamma_) * (v_K_vec - lambda_snips * v_K_prime_vec) / (lambda_snips * v_K_vec)^2
    
    C_ij <- matrix(1, K, K)
    for (x in 1:(K)) {
      for (y in 1:(K)) {
        C_ij[x, y] <- order_0 - lambda_snips[x] * order_1[x] -  lambda_snips[y] * order_1[y] +
          (lambda_snips[x] * m_K_vec[x] - lambda_snips[y] * m_K_vec[y]) / ((m_K_vec[y] - m_K_vec[x]))
        
      }
    }
    
    diag(C_ij) <- (order_0 - 2 * lambda_snips * order_1 + lambda_snips^2 * order_2)
    
    D_i <- matrix(order_0 - lambda_snips * order_1, K, 1)
    F_i <- diag((order_1 - lambda_snips* order_2),  K)
    
    D_f <- D_i * rho_mat[, K]* alpha_sigma^2
    F_f <- gamma_ * alpha_^2 * F_i
    C_f <- C_ij * rho_mat * matrix(alpha_sigma) %*% alpha_sigma
    
    weight_ <- as.numeric(abs(solve(C_f + F_f) %*% D_f))
    
    beta_list <- foreach(X_Y = train_XY,
                         lambda = lambda_snips) %do% {
                           train_x_ <- X_Y[[1]]
                           train_y_ <- X_Y[[2]]
                           
                           net_rid <- glmnet(train_x_, 
                                             train_y_, alpha = 0, lambda = lambda, standardize = TRUE)
                           
                           return(list(inter = coef(net_rid)[1],
                                       beta = coef(net_rid)[-1]))
                         }
    
    inter_tl <- mean(unlist(purrr::transpose(beta_list)[[1]]))
    beta_tl <- (purrr::transpose(beta_list)[[2]] %>% reduce(cbind)) %*% weight_
    score_vec <- c(score_vec, x_test_i %*% beta_tl + inter_tl)
    acc_vec <- c(acc_vec, mean(((x_test_i %*% beta_tl + inter_tl) < 0) == (y_test_i < 0)))
  }
  
  score_list[[length(score_list) + 1]] <- score_vec
}

X_list <- map(1:3, ~{
  ind.set(set_indx, .x)
})

data.table::transpose(score_list) %>%
  map_dbl(~{
    roc(train_y[X_list[[3]]], .x) %>%
      {as.numeric(.$auc)}
  }) %>% max

#estimation
est_score_list <- list()
for (index in train_cv_list) {
  X_list <- map(1:3, ~{
    ind.set(set_indx, .x)
  })
  
  X_list[[3]] <- X_list[[3]][-which(index == X_list[[3]])]
  
  x_train_list <- map(X_list, ~{train_x[.x, ]})
  y_train_list <- map(X_list, ~{train_y[.x]})
  
  x_test_i <- train_x[index,]
  y_test_i <- train_y[index]
  
  n_vec <- map_dbl(x_train_list, nrow)
  
  sigma_ <-  read.table("sigma_.Rdata") %>% 
    {setNames(pull(., 1), rownames(.))}
  
  sigma_ <- sigma_[name_vec]
  
  alpha_ <-  read.table("alpha_.Rdata") %>% 
    {setNames(pull(., 1), rownames(.))}
  
  alpha_ <- alpha_[name_vec]
  
  alpha_sigma <- read.table("alpha_sigma.Rdata") %>% 
    {setNames(pull(., 1), rownames(.))}
  
  alpha_sigma <- alpha_sigma[name_vec]
  
  rho_mat <- read.table("rho_mat.Rdata") %>% as.matrix() %>% 
    .[name_vec, name_vec]
  
  gamma_ <- n_vec / ncol(x_train_list[[1]])
  
  train_XY <- map2(x_train_list, y_train_list, ~{
    list(.x,
         .y)
  })
  
  score_vec <- c()

  for (lambda_factor in lambda_factor_vector) {
    lambda_snips <- gamma_ / alpha_^2
    lambda_snips <- lambda_factor * lambda_snips
    
    inverse_list_p <- foreach(X_Y = train_XY,
                              lambda = lambda_snips) %do% {
                                train_x_ <- X_Y[[1]]
                                
                                temp_dgC <- (t(train_x_) %*% train_x_/nrow(train_x_) +  
                                               diag(lambda, ncol(train_x_))) %>%
                                  solve()
                                
                                return(temp_dgC)
                              }
    
    m_K_vec <- map_dbl(inverse_list_p, ~{mean(diag(.x))})
    m_K_prime_vec <- map_dbl(inverse_list_p, ~{mean(diag(.x)^2)})
    
    
    pk_1 <- rho_mat[,3]
    v_i <- alpha_sigma^2 * pk_1 * (1 - lambda_snips * m_K_vec)
    
    A_ij <- matrix(1, K, K)
    for (x in 1:(K)) {
      for (y in 1:(K)) {
        A_ij[x, y] <- 1 + (lambda_snips[x] - lambda_snips[y]) * (m_K_vec[x] * m_K_vec[y]) / (m_K_vec[x] - m_K_vec[y])
      }
    }
    A_ij <- A_ij * rho_mat * matrix(alpha_sigma) %*% alpha_sigma
    
    diag(A_ij) <- alpha_sigma^2 * 
      (1 - 2 * lambda_snips * m_K_vec + lambda_snips^2 * m_K_prime_vec)
    
    R_ii <- diag(alpha_^2 * gamma_ * (m_K_vec - lambda_snips * m_K_prime_vec), K)
    weight_ <- solve(R_ii + A_ij) %*% matrix(v_i)
    
    beta_list <- foreach(X_Y = train_XY,
                         lambda = lambda_snips) %do% {
                           train_x_ <- X_Y[[1]]
                           train_y_ <- X_Y[[2]]
                           
                           net_rid <- glmnet(train_x_, 
                                             train_y_, alpha = 0, lambda = lambda, standardize = TRUE)
                           
                           return(list(inter = coef(net_rid)[1],
                                       beta = coef(net_rid)[-1]))
                         }
    
    inter_tl <- mean(unlist(purrr::transpose(beta_list)[[1]]))
    beta_tl <- (purrr::transpose(beta_list)[[2]] %>% reduce(cbind)) %*% weight_
    
    score_vec <- c(score_vec, (x_test_i %*% beta_tl+ inter_tl))

  }
  est_score_list[[length(est_score_list) + 1]] <- score_vec
}

X_list <- map(1:3, ~{
  ind.set(set_indx, .x)
})
data.table::transpose(est_score_list) %>%
  map_dbl(~{
    roc(train_y[X_list[[3]]], .x) %>%
      {as.numeric(.$auc)}
  }) %>% max












