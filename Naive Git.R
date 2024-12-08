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
library(glmnet)
library(foreach)
library(tidyverse)
library(pROC)

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

set.seed(111)
DF.wk <- map(1:3, ~{
  ind.set(set_indx, .x) %>%
    sample(length(.), replace = FALSE)
}) %>%
  map(~{DF.wk[.x, ]}) %>%
  bind_rows()


wk_tb <- as_tibble(DF.wk)
group_index <- rep(1:3, set_indx)
wk_tb$DiseaseState <- as.numeric(scale(ifelse(wk_tb$DiseaseState == "CRC", 1, -1)))
target_ind <- 3
ex <- ind.set(set_indx, target_ind)
testing_ind <- split(ex, cut_number(ex, 10))[[10]]
new_id <- set_indx
new_id[3] <- new_id[3] - length(testing_ind)
training_ind <- map(1:3, ~{
  ind.set(new_id, .x)
})
train_xy <- wk_tb
naive_set <- train_xy %>%
  select(-DiseaseState, -gender, -age, -BMI) %>%
  apply(2, function(x) sum(x == 0) / length(x)) %>%
  {names(.)[which(. <= 0.9)]}

train_x <- select(train_xy, gender:BMI, all_of(naive_set))
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
##########################################################
train_cv_list <- ind.set(set_indx, 3)
lambda_factor_vector <- 10^seq(-4, 4)

#one should make sure the covariance matrix is aligned with current outcome
alpha_sigma <- read.table("alpha_sigma.Rdata") %>% 
  {setNames(pull(., 1), rownames(.))}

score_list <- list()
for (index in train_cv_list) {
  
  train_index_i <- (1:nrow(train_x)) == index
  
  train_x_i <- train_x[!train_index_i, ]
  test_x_i <- train_x[train_index_i, ]
  train_y_i <- train_y[!train_index_i]
  test_y_i <- train_y[train_index_i]
  score_vec <- c()
  for (lambda_factor in lambda_factor_vector) {
    gamma_ <- ncol(train_x_i) / nrow(train_x_i)
    lambda_snips <- gamma_ / mean(alpha_sigma)
    lambda_snips <- lambda_factor * lambda_snips
    
    ridge_m <- glmnet(train_x_i, 
                      train_y_i, alpha = 0, lambda = lambda_snips, standardize = TRUE)
    score_vec <- c(score_vec, as.numeric(test_x_i %*% ridge_m$beta))
    acc_vec <- c(acc_vec, mean((test_x_i %*% ridge_m$beta < 0 ) == (test_y_i < 0)))
    
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
  })%>% max


##########################################################
#####################Target data Only#####################
##########################################################

set.seed(111)

target_name <- c("Zackular", "Zeller", "Baxter")[2]
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
target_ind <- 3
train_xy <- wk_tb[group_index == target_ind, ]

naive_set <- train_xy %>%
  select(-DiseaseState, -gender, -age, -BMI) %>%
  apply(2, function(x) sum(x == 0) / length(x)) %>%
  {names(.)[which(. <= 0.9)]}

train_x <- select(train_xy, gender:BMI, all_of(naive_set))

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

train_cv_list <- 1:nrow(train_x)
lambda_factor_vector <- 10^seq(-4, 4)

alpha_sigma <- read.table("alpha_sigma.Rdata") %>% 
  {setNames(pull(., 1), rownames(.))}
#one should make sure the covariance matrix is aligned with current outcome

naive_score_list <- list()

for (index in train_cv_list) {
  
  train_index_i <- (1:nrow(train_x)) == index
  
  train_x_i <- train_x[!train_index_i, ]
  test_x_i <- train_x[train_index_i, ]
  train_y_i <- train_y[!train_index_i]
  test_y_i <- train_y[train_index_i]
  
  score_vec <- c()
  acc_vec <- c()
  #model_list <- list()
  for (lambda_factor in lambda_factor_vector) {
    gamma_ <- ncol(train_x_i) / nrow(train_x_i)
    lambda_snips <- gamma_ / mean(alpha_sigma)
    lambda_snips <- lambda_factor * lambda_snips
    
    ridge_m <- glmnet(train_x_i, 
                      train_y_i, alpha = 0, lambda = lambda_snips, standardize = TRUE)
    score_vec <- c(score_vec, as.numeric(test_x_i %*% ridge_m$beta))
  }
  
  naive_score_list[[length(naive_score_list) + 1]] <- score_vec
}

data.table::transpose(naive_score_list) %>%
  map_dbl(~{
    roc(train_y, .x) %>%
      {as.numeric(.$auc)}
  }) %>% max



