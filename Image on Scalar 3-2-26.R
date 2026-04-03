library(plot.matrix)
library(Matrix)
library(fastBayesReg)
library(BayesGPfit)
library(reshape2)
library(gridExtra)
library(locfdr)
library(Rcpp)
library(profvis)
library(foreach)
library(RNifti)
library(oro.nifti)
library(neurobase)

setwd("~/Dissertation/Code")

sourceCpp('bracket.quadratic.cpp')

#evaluate quadratic function
quadratic.f = function(coef, x) {
  return(coef[1]*x^2+coef[2]*x+coef[3])
}

set.seed(1)

voxels_in_regions <- function(region_list){
  unlist(lapply(region_list, `[[`, 1))
}

voxel_array_index <- function(voxels){
  voxel_d1 <- ifelse(voxels %% d1 != 0, voxels %% d1, d1)
  voxel_d2 <- ceiling(((voxels - 0.1) %% (d1 * d2)) / d1)
  voxel_d3 <- ceiling(voxels / (d1 * d2))
  voxel_matrix <- matrix(c(voxel_d1, voxel_d2, voxel_d3), ncol = 3)
  return(voxel_matrix)
}

region_list_to_betahat <- function(region_list){
  betahat <- matrix(0, nrow = 1, ncol = V)
  if(length(region_list) == 0){
    return(betahat)
  }
  else{
    for(i in 1:length(region_list)){
      betahat[region_list[[i]][[1]]] <- region_list[[i]][[2]]
    }
    return(betahat)
  }
}

get_neighbors <- function(region, region_list = NULL, unused_only, mask = NULL){
  if(is.null(region_list)){
    used_voxels <- region
  }
  else{
    used_voxels <- voxels_in_regions(region_list)
  }
  neighbor_list <- numeric()
  for(i in 1:length(region)){
    voxel <- region[i]
    voxel_index <- voxel_array_index(voxel)
    voxel_d1 <- voxel_index[1]
    voxel_d2 <- voxel_index[2]
    voxel_d3 <- voxel_index[3]
    potential_neighbors <- matrix(0, ncol = 3, nrow = 6)
    potential_neighbors[1, ] <- c(voxel_d1 + 1, voxel_d2, voxel_d3)
    potential_neighbors[2, ] <- c(voxel_d1 - 1, voxel_d2, voxel_d3)
    potential_neighbors[3, ] <- c(voxel_d1, voxel_d2 + 1, voxel_d3)
    potential_neighbors[4, ] <- c(voxel_d1, voxel_d2 - 1, voxel_d3)
    potential_neighbors[5, ] <- c(voxel_d1, voxel_d2, voxel_d3 + 1)
    potential_neighbors[6, ] <- c(voxel_d1, voxel_d2, voxel_d3 - 1)
    remove_index <- numeric()
    for(i in 1:6){
      if(0 %in% potential_neighbors[i, ] | potential_neighbors[i, ][1] > d1 | potential_neighbors[i, ][2] > d2 | potential_neighbors[i, ][3] > d3){
        remove_index <- c(remove_index, i)
      }
      else if(unused_only){
        if((potential_neighbors[i, 1] + d1 * (potential_neighbors[i , 2] - 1) + d1 * d2 * (potential_neighbors[i, 3] - 1)) %in% used_voxels){
          remove_index <- c(remove_index, i)
        }
      }
      if(!is.null(mask)){
        if(mask[potential_neighbors[i, 1] + d1 * (potential_neighbors[i , 2] - 1) + d1 * d2 * (potential_neighbors[i, 3] - 1)] == 0){
          remove_index <- c(remove_index, i)
        }
      }
    }
    if(length(remove_index) > 0){
      potential_neighbors <- potential_neighbors[-remove_index, , drop = F]
    }
    neighbors <- potential_neighbors[, 1] + d1 * (potential_neighbors[ , 2] - 1) + d1 * d2 * (potential_neighbors[ , 3] - 1)
    neighbors <- neighbors[!(neighbors %in% region)]
    neighbor_list <- c(neighbor_list, neighbors)
  }
  neighbor_list <- neighbor_list[!duplicated(neighbor_list)]
  return(neighbor_list)
}

# Computes norm(X %*% betahat_full_new, type = "F") ^ 2 - norm(X %*% betahat_full, type = "F") ^ 2
matrix_squared_norm_difference <- function(which_beta, c1, v1, c2 = NULL, v2 = NULL) { 
  XTX_jj <- XTX[which_beta, which_beta]
  cross_term1 <- 2 * c1 * sum(XTXB[which_beta, v1])
  sq_norm_diff <- cross_term1 + length(v1) * (c1 ^ 2) * XTX_jj
  if(!is.null(c2)){
    cross_term2 <- 2 * c2 * sum(XTXB[which_beta, v2])
    sq_norm_diff <- sq_norm_diff + cross_term2 + length(v2) * (c2 ^ 2) * XTX_jj
  }
  return(sq_norm_diff)
}

compute_sad_diff <- function(betahat, c1, v1, c2 = NULL, v2 = NULL) {
  # Get all modified voxels
  mod_voxels <- if (is.null(v2)) v1 else c(v1, v2)
  
  # Find all unique edge indices to process
  edge_indices <- unique(unlist(
    lapply(mod_voxels, function(v) {
      edge_structure$index[[as.character(v)]]
    }),
    use.names = FALSE
  ))
  
  # Initialize difference
  delta_sad <- 0
  
  # Process each affected edge only once
  for (idx in edge_indices) {
    i <- edge_structure$edges[idx, 1]
    j <- edge_structure$edges[idx, 2]
    
    # Original difference
    orig_diff <- abs(betahat[i] - betahat[j])
    
    # Calculate modifications
    mod_i <- betahat[i] + 
      (if (i %in% v1) c1 else 0) + 
      (if (!is.null(v2) && (i %in% v2)) c2 else 0)
    
    mod_j <- betahat[j] + 
      (if (j %in% v1) c1 else 0) + 
      (if (!is.null(v2) && (j %in% v2)) c2 else 0)
    
    # Modified difference
    mod_diff <- abs(mod_i - mod_j)
    
    delta_sad <- delta_sad + (mod_diff - orig_diff)
  }
  
  return(delta_sad)
}

abs_sum_diff <- function(betahat, c1, v1, c2 = NULL, v2 = NULL) {
  # Calculate sum difference for v1 voxels
  sum_diff <- sum(abs(betahat[v1] + c1)) - sum(abs(betahat[v1]))
  
  # Add sum difference for v2 voxels if provided
  if (!is.null(v2)) {
    sum_diff <- sum_diff + sum(abs(betahat[v2] + c2)) - sum(abs(betahat[v2]))
  }
  
  return(sum_diff)
}

loss_difference <- function(betahat, which_beta, c1, v1, c2 = NULL, v2 = NULL){
  changing_voxels <- v1
  change_amount <- rep(c1, length(v1))
  if(!is.null(v2)){
    changing_voxels <- c(changing_voxels, v2)
    change_amount <- c(change_amount, rep(c2, length(v2)))
  }
  term1 <- matrix_squared_norm_difference(which_beta, c1, v1, c2, v2)
  term2 <- XTY_tilde[which_beta, changing_voxels] %*% change_amount
  ms_loss_diff <- (term1 - 2 * term2) / n
  
  neighbor_loss_diff <- compute_sad_diff(betahat, c1, v1, c2, v2)
  
  L1_loss_diff <- abs_sum_diff(betahat, c1, v1, c2, v2)
  loss_diff <- ms_loss_diff + lambda * neighbor_loss_diff + mu * L1_loss_diff
  return(loss_diff)
}


expand_region <- function(region_list, region_index, betahat, which_beta){
  region_voxels <- region_list[[region_index]][[1]]
  region_beta   <- region_list[[region_index]][[2]]
  all_voxels_in_regions <- voxels_in_regions(region_list)
  
  # Fix: use unlist() instead of c() in a loop (O(n) vs O(n^2))
  region_neighbors <- unique(unlist(neighbor_list[region_voxels]))
  region_neighbors <- region_neighbors[!(region_neighbors %in% all_voxels_in_regions)]
  
  if(length(region_neighbors) == 0){
    return(c(Inf, NA))
  }
  else{
    expand_loss_diffs <- numeric(length(region_neighbors))
    for(i in 1:length(region_neighbors)){
      expand_loss_diffs[i] <- loss_difference(betahat, which_beta, region_beta, region_neighbors[i])
    }
    min_expand_loss_diff <- min(expand_loss_diffs)
    min_index           <- which.min(expand_loss_diffs)
    voxel_to_add        <- region_neighbors[min_index]
    return(c(min_expand_loss_diff, voxel_to_add))
  }
}

shrink_region <- function(region_list, region_index, betahat, which_beta){
  region_voxels <- region_list[[region_index]][[1]]
  region_beta   <- region_list[[region_index]][[2]]
  if(length(region_voxels) == 1){
    shrink_loss_diff <- loss_difference(betahat, which_beta, -region_beta, region_voxels)
    return(c(shrink_loss_diff, region_voxels))
  }
  else{
    shrink_loss_diffs <- numeric(length(region_voxels))
    for(i in 1:length(region_voxels)){
      shrink_loss_diffs[i] <- loss_difference(betahat, which_beta, -region_beta, region_voxels[i])
    }
    min_shrink_loss_diff <- min(shrink_loss_diffs)
    min_index            <- which.min(shrink_loss_diffs)
    voxel_to_remove      <- region_voxels[min_index]
    return(c(min_shrink_loss_diff, voxel_to_remove))
  }
}

bracket_quadratic_setup <- function(region_voxels, region_size, region_beta, region_neighbors, betahat, which_beta, XTXB){
  neighbors_size <- length(region_neighbors)
  
  coefs <- matrix(0, nrow = 2 * neighbors_size + 3, ncol = 3)
  coefs[1:neighbors_size, 2] <- -lambda
  coefs[1:neighbors_size, 3] <- lambda * (betahat[region_neighbors] - region_beta)
  coefs[(neighbors_size + 1):(2 * neighbors_size), 2] <- lambda
  coefs[(neighbors_size + 1):(2 * neighbors_size), 3] <- -lambda * (betahat[region_neighbors] - region_beta)
  coefs[2 * neighbors_size + 1, 2] <- -mu * region_size
  coefs[2 * neighbors_size + 1, 3] <- -mu * region_size * region_beta
  coefs[2 * neighbors_size + 2, 2] <- mu * region_size
  coefs[2 * neighbors_size + 2, 3] <- mu * region_size * region_beta
  coefs[2 * neighbors_size + 3, 1] <- region_size * XTX[which_beta, which_beta] / n
  coefs[2 * neighbors_size + 3, 2] <- 2 / n * sum(XTXB[which_beta, region_voxels] - XTY_tilde[which_beta, region_voxels])
  
  b <- matrix(0, nrow = 2 * neighbors_size + 3, ncol = 2)
  b[1:neighbors_size, 1] <- -Inf
  b[1:neighbors_size, 2] <- betahat[region_neighbors] - region_beta
  b[(neighbors_size + 1):(2 * neighbors_size), 1] <- betahat[region_neighbors] - region_beta
  b[(neighbors_size + 1):(2 * neighbors_size), 2] <- Inf
  b[2 * neighbors_size + 1, 1] <- -Inf
  b[2 * neighbors_size + 1, 2] <- region_beta
  b[2 * neighbors_size + 2, 1] <- region_beta
  b[2 * neighbors_size + 2, 2] <- Inf
  b[2 * neighbors_size + 3, 1] <- -Inf
  b[2 * neighbors_size + 3, 2] <- Inf
  return(list(coefs, b))
}

revalue_region <- function(region_list, region_index, betahat, which_beta){ 
  region_voxels <- region_list[[region_index]][[1]]
  region_beta <- region_list[[region_index]][[2]]
  region_size <- length(region_voxels)
  region_neighbors <- numeric()
  for(i in 1:region_size){
    region_neighbors <- c(region_neighbors, neighbor_list[[region_voxels[i]]])
  }
  region_neighbors <- region_neighbors[!duplicated(region_neighbors)]
  region_neighbors <- region_neighbors[!(region_neighbors %in% region_voxels)]
  bracket_quadratic_inputs <- bracket_quadratic_setup(region_voxels, region_size, region_beta, region_neighbors, betahat, which_beta, XTXB)
  beta_optimal_diff <- bracket_quadratic(coef = bracket_quadratic_inputs[[1]], a = 0, b = bracket_quadratic_inputs[[2]])$par
  beta_optimal <- region_beta + beta_optimal_diff
  if(abs(beta_optimal) < 1e-4){
    beta_optimal_diff <- -region_beta
  }
  revalue_loss_diff <- loss_difference(betahat, which_beta, beta_optimal_diff, region_voxels)
  if(abs(revalue_loss_diff) < 1e-6){
    revalue_loss_diff <- 0
  }
  return(c(revalue_loss_diff, beta_optimal))
}


merge_regions <- function(region_list, betahat, num_regions, which_beta, betahat_full){
  if(num_regions == 1){
    return(list(Inf, NA, NA))
  }
  
  merge_loss_diff <- numeric(num_regions - 1)
  beta_optimal    <- numeric(num_regions - 1)
  region_betas    <- numeric(num_regions)
  
  for(i in 1:num_regions){
    region_betas[i] <- region_list[[i]][[2]]
  }
  
  region_ordering        <- order(region_betas)
  regions_ordered_voxels <- vector("list", num_regions)
  regions_ordered_betas  <- numeric(num_regions)
  
  for(i in 1:num_regions){
    regions_ordered_voxels[[i]] <- region_list[[region_ordering[i]]][[1]]
    regions_ordered_betas[i]    <- region_list[[region_ordering[i]]][[2]]
  }
  
  for(i in 1:(num_regions - 1)){
    voxels_to_merge <- c(regions_ordered_voxels[[i]], regions_ordered_voxels[[i + 1]])
    region_size     <- length(voxels_to_merge)
    
    # Fix: use unlist() instead of c() in a loop
    region_neighbors <- unique(unlist(neighbor_list[voxels_to_merge]))
    region_neighbors <- region_neighbors[!(region_neighbors %in% voxels_to_merge)]
    
    mean_beta       <- mean(regions_ordered_betas[i:(i + 1)])
    temp_region_list <- region_list[-c(region_ordering[i], region_ordering[i + 1])]
    temp_region_list[[num_regions - 1]] <- list(voxels_to_merge, mean_beta)
    temp_betahat    <- region_list_to_betahat(temp_region_list)
    
    # Fix: use passed-in betahat_full instead of accessing alphahat_full directly
    temp_betahat_full <- betahat_full
    temp_betahat_full[which_beta, ] <- temp_betahat
    temp_XTXB <- XTX %*% temp_betahat_full
    
    bracket_quadratic_inputs <- bracket_quadratic_setup(
      voxels_to_merge, region_size, mean_beta, region_neighbors,
      temp_betahat, which_beta, temp_XTXB
    )
    beta_optimal[i] <- bracket_quadratic(
      coef = bracket_quadratic_inputs[[1]], a = 0, b = bracket_quadratic_inputs[[2]]
    )$par + mean_beta
    
    if(abs(beta_optimal[i]) < 1e-4){
      beta_optimal[i] <- 0
    }
    
    merge_loss_diff[i] <- loss_difference(
      betahat, which_beta,
      beta_optimal[i] - regions_ordered_betas[i],     regions_ordered_voxels[[i]],
      beta_optimal[i] - regions_ordered_betas[i + 1], regions_ordered_voxels[[i + 1]]
    )
    if(abs(merge_loss_diff[i]) < 1e-6){
      merge_loss_diff[i] <- 0
    }
  }
  
  min_index <- which.min(merge_loss_diff)
  return(list(
    merge_loss_diff[min_index],
    beta_optimal[min_index],
    c(region_ordering[min_index], region_ordering[min_index + 1])
  ))
}

gradient <- function(betahat, voxel, which_beta){
  total <- 2 / n * (XTXB[which_beta, voxel] - XTY_tilde[which_beta, voxel])
  total <- total + lambda * sum(sign(betahat[voxel] - betahat[neighbor_list[[voxel]]]))
  total <- total + mu * sign(betahat[voxel])
  return(total)
}


split_region <- function(region_list, region_index, betahat, which_beta, betahat_full){
  region_voxels <- region_list[[region_index]][[1]]
  region_size   <- length(region_voxels)
  region_beta   <- region_list[[region_index]][[2]]
  
  if(region_size == 1){
    return(list(Inf, NA, NA, NA))
  }
  
  # Split voxels into two halves based on gradient ranking
  gradient_values <- numeric(region_size)
  for(i in 1:region_size){
    gradient_values[i] <- gradient(betahat, region_voxels[i], which_beta)
  }
  halfway             <- floor(region_size / 2)
  first_region_index  <- which(rank(gradient_values) <= halfway)
  second_region_index <- which(rank(gradient_values) > halfway)
  region_1      <- region_voxels[first_region_index]
  region_2      <- region_voxels[second_region_index]
  region_1_size <- length(region_1)
  region_2_size <- length(region_2)
  
  # Fix: use unlist() instead of c() in a loop
  region_1_neighbors <- unique(unlist(neighbor_list[region_1]))
  region_1_neighbors <- region_1_neighbors[!(region_1_neighbors %in% region_1)]
  region_2_neighbors <- unique(unlist(neighbor_list[region_2]))
  region_2_neighbors <- region_2_neighbors[!(region_2_neighbors %in% region_2)]
  
  starting_betahat  <- betahat
  # Fix: use passed-in betahat_full instead of accessing alphahat_full directly
  temp_betahat_full <- betahat_full
  region_1_beta     <- region_beta
  region_2_beta     <- region_beta
  split_loss_diff   <- Inf
  split_iter        <- 0
  
  # Fix: cap iterations to prevent infinite loop
  while(abs(split_loss_diff) > 1e-4 && split_iter <= 100){
    temp_betahat <- starting_betahat
    temp_betahat_full[which_beta, ] <- temp_betahat
    temp_XTXB <- XTX %*% temp_betahat_full
    
    bracket_quadratic_inputs_1 <- bracket_quadratic_setup(
      region_1, region_1_size, region_1_beta,
      region_1_neighbors, temp_betahat, which_beta, temp_XTXB
    )
    region_1_beta_diff <- bracket_quadratic(
      coef = bracket_quadratic_inputs_1[[1]], a = 0, b = bracket_quadratic_inputs_1[[2]]
    )$par
    region_1_beta <- region_1_beta + region_1_beta_diff
    temp_betahat[region_1] <- region_1_beta
    temp_betahat_full[which_beta, ] <- temp_betahat
    temp_XTXB <- XTX %*% temp_betahat_full
    
    bracket_quadratic_inputs_2 <- bracket_quadratic_setup(
      region_2, region_2_size, region_2_beta,
      region_2_neighbors, temp_betahat, which_beta, temp_XTXB
    )
    region_2_beta_diff <- bracket_quadratic(
      coef = bracket_quadratic_inputs_2[[1]], a = 0, b = bracket_quadratic_inputs_2[[2]]
    )$par
    region_2_beta <- region_2_beta + region_2_beta_diff
    temp_betahat[region_2] <- region_2_beta
    
    split_loss_diff  <- loss_difference(starting_betahat, which_beta,
                                        region_1_beta_diff, region_1,
                                        region_2_beta_diff, region_2)
    starting_betahat <- temp_betahat
    split_iter       <- split_iter + 1
  }
  if(split_iter > 100){
    warning(sprintf("split_region failed to converge for region %d", region_index))
  }
  
  if(abs(region_1_beta) < 1e-4){ region_1_beta <- 0 }
  if(abs(region_2_beta) < 1e-4){ region_2_beta <- 0 }
  
  full_split_loss_diff <- loss_difference(betahat, which_beta,
                                          region_1_beta - region_beta, region_1,
                                          region_2_beta - region_beta, region_2)
  return(list(full_split_loss_diff, c(region_1_beta, region_2_beta), region_1, region_2))
}

update_region_list <- function(region_list, which_beta){
  num_regions <- length(region_list)
  if(num_regions == 0 | is.null(region_list[[1]][[1]])){
    terminate_due_to_null_regions <- T
    return_info <- list(T, NULL, NULL, NULL, NULL)
    print(paste("FINISHED UPDATING BETA ", p1_index[which_beta], sep = " "))
    return(list(region_list, return_info))
  }
  num_voxels_in_regions <- length(voxels_in_regions(region_list))
  betahat <- region_list_to_betahat(region_list)
  betahat_full <- alphahat_full[p1_index, , drop = F]
  current_min_loss_diff <- 0
  update <- 0
  terminate_due_to_convergence <- F
  terminate_due_to_null_regions <- F
  if(num_voxels_in_regions < V * 0.2){
    for(i in 1:num_regions){
      expand_loss_diff_region <- expand_region(region_list, i, betahat, which_beta)
      if(expand_loss_diff_region[1] < current_min_loss_diff){
        current_min_loss_diff <- expand_loss_diff_region[1]
        region_to_expand <- i
        voxel_to_add <- expand_loss_diff_region[2]
        update <- 1
      }
    }
  }
  for(i in 1:num_regions){
    shrink_loss_diff_region <- shrink_region(region_list, i, betahat, which_beta)
    if(shrink_loss_diff_region[1] < current_min_loss_diff){
      current_min_loss_diff <- shrink_loss_diff_region[1]
      region_to_shrink <- i
      voxel_to_remove <- shrink_loss_diff_region[2]
      update <- 2
    }
  }
  for(i in 1:num_regions){
    revalue_loss_diff <- revalue_region(region_list, i, betahat, which_beta)
    if(revalue_loss_diff[1] < current_min_loss_diff){
      current_min_loss_diff <- revalue_loss_diff[1]
      region_to_revalue <- i
      optimal_beta_value <- revalue_loss_diff[2]
      update <- 3
    }
  }
  merge_optimal <- merge_regions(region_list, betahat, num_regions, which_beta, betahat_full)
  merge_loss_diff <- merge_optimal[[1]]
  if(merge_loss_diff < current_min_loss_diff){
    current_min_loss_diff <- merge_loss_diff
    regions_to_merge <- merge_optimal[[3]]
    voxels_to_merge <- c(region_list[[regions_to_merge[1]]][[1]], region_list[[regions_to_merge[2]]][[1]])
    merged_betahat <- merge_optimal[[2]]
    update <- 4
  }
  if(num_regions < 15){
    for(i in 1:num_regions){
      split_optimal <- split_region(region_list, i, betahat, which_beta, betahat_full)
      split_loss_diff <- split_optimal[[1]]
      if(split_loss_diff < current_min_loss_diff){
        current_min_loss_diff <- split_loss_diff
        region_to_split <- i
        new_region_1 <- split_optimal[[3]]
        new_region_2 <- split_optimal[[4]]
        new_beta_1 <- split_optimal[[2]][1]
        new_beta_2 <- split_optimal[[2]][2]
        update <- 5
      }
    }
  }
  
  new_region_list <- region_list
  if(update == 1){
    return_info <- list(F, new_region_list[[region_to_expand]][[2]], voxel_to_add, NULL, NULL)
    new_region_list[[region_to_expand]][[1]] <- c(new_region_list[[region_to_expand]][[1]], voxel_to_add)
    print(paste("EXPAND REGION", region_to_expand, "ADD VOXEL", voxel_to_add, sep = " "))
  }
  else if(update == 2){
    return_info <- list(F, -new_region_list[[region_to_shrink]][[2]], voxel_to_remove, NULL, NULL)
    if(length(region_list[[region_to_shrink]][[1]]) == 1){
      new_region_list <- new_region_list[-region_to_shrink]
    }
    else{
      new_region_list[[region_to_shrink]][[1]] <- setdiff(new_region_list[[region_to_shrink]][[1]], voxel_to_remove)
    }
    print(paste("SHRINK REGION", region_to_shrink, "REMOVE VOXEL", voxel_to_remove, sep = " "))
  }
  else if(update == 3){
    return_info <- list(F, optimal_beta_value - new_region_list[[region_to_revalue]][[2]], new_region_list[[region_to_revalue]][[1]], NULL, NULL)
    new_region_list[[region_to_revalue]][[2]] <- optimal_beta_value
    print(paste("REVALUE REGION", region_to_revalue, "TO", round(optimal_beta_value, 4), sep = " "))
  }
  else if(update == 4){
    return_info <- list(F, merged_betahat - region_list[[regions_to_merge[1]]][[2]], region_list[[regions_to_merge[1]]][[1]],
                        merged_betahat - region_list[[regions_to_merge[2]]][[2]], region_list[[regions_to_merge[2]]][[1]])
    new_region_list <- new_region_list[-regions_to_merge]
    new_region_list[[num_regions - 1]] <- list(voxels_to_merge, merged_betahat)
    print(paste("MERGE REGIONS", regions_to_merge[1], "AND", regions_to_merge[2], sep = " "))
  }
  else if(update == 5){
    return_info <- list(F, new_beta_1 - new_region_list[[region_to_split]][[2]], new_region_1,
                        new_beta_2 - new_region_list[[region_to_split]][[2]], new_region_2)
    new_region_list <- new_region_list[-region_to_split]
    new_region_list[[num_regions]] <- list(new_region_1, new_beta_1)
    new_region_list[[num_regions + 1]] <- list(new_region_2, new_beta_2)
    print(paste("SPLIT REGION", region_to_split, sep = " "))
  }
  else{
    terminate_due_to_convergence <- T
    return_info <- list(T, NULL, NULL, NULL, NULL)
    print(paste("FINISHED UPDATING BETA ", p1_index[which_beta], sep = " "))
  }
  if(length(new_region_list) == 0){
    terminate_due_to_null_regions <- T
    return_info <- list(T, NULL, NULL, NULL, NULL)
    print(paste("FINISHED UPDATING BETA ", p1_index[which_beta], sep = " "))
  }
  # Eliminate regions with beta = 0
  else{
    zero_regions <- numeric()
    for(i in 1:length(new_region_list)){
      if(new_region_list[[i]][[2]] == 0){
        zero_regions <- c(zero_regions, i)
      }
    }
    if(length(zero_regions) > 0){
      new_region_list <- new_region_list[-zero_regions]
    }
  }
  if(length(new_region_list) == 0){
    terminate_due_to_null_regions <- T
    return_info <- list(T, NULL, NULL, NULL, NULL)
    print(paste("FINISHED UPDATING BETA ", p1_index[which_beta], sep = " "))
  }
  # arrange regions by beta value
  # merge regions with near-identical betas
  
  if((inner_iter %% 5 == 0 | terminate_due_to_convergence) & !terminate_due_to_null_regions){
    new_region_betas <- numeric(length(new_region_list))
    for(i in 1:length(new_region_list)){
      new_region_betas[i] <- new_region_list[[i]][[2]]
    }
    new_region_betas_ordering <- order(new_region_betas)
    new_region_list <- new_region_list[new_region_betas_ordering]
    new_region_betas_ordered <- new_region_betas[new_region_betas_ordering]
    new_region_betas_duplicated <- duplicated(round(new_region_betas_ordered, 3))
    if(length(new_region_list) > 1){
      for(i in length(new_region_list):1){
        if(new_region_betas_duplicated[i]){
          new_region_list[[i - 1]][[1]] <- c(new_region_list[[i - 1]][[1]], new_region_list[[i]][[1]])
          new_region_list <- new_region_list[-i]
          print(paste("MERGED SIMILAR REGIONS", i - 1, "AND", i, sep = " "))
        }
      }
    }
  }
  
  return(list(new_region_list, return_info))
}

naive_solution <- function(Y, Q){
  mlm <- lm(Y ~ Q)
  mlm_coef <- coef(summary(mlm))
  mua_solution <- matrix(nrow = p, ncol = V)
  for(j in 1:p){
    mua_solution[j , ] <- unname(unlist(lapply(mlm_coef, function(f) f[j, 1]))) * 
      (p.adjust(unname(unlist(lapply(mlm_coef, function(f) f[j, 4]))), "BH") < 0.05)
  }
  return(mua_solution)
}

smooth_vec <- function(vec, times){
  smoothed_vec <- vec
  for(i in 1:times){
    smoothed_vec <- smoothingMatrix %*% smoothed_vec
  }
  return(as.numeric(smoothed_vec))
}


initialize_regions_improved <- function(Y, Q, t_values = NULL, which_alpha, use_locfdr = F, breaks = 120, df = 7, nulltype = 1){
  initial_voxels <- numeric()
  naive_betahat <- least_squares_solution[which_alpha, ]
  if(is.null(t_values)){
    t_values <- unname(unlist(lapply(mlm_coef, function(f) f[which_alpha, 3])))
    t_values[is.na(t_values)] <- 0
  }
  smoothed_t <- smooth_vec(t_values, 2)
  if(use_locfdr){
    smoothed_t_locfdr <- locfdr(smoothed_t, bre = breaks, df = df, nulltype = nulltype, plot = 0)$fdr
    initial_voxels <- which(smoothed_t_locfdr <= 0.05)
  }
  if(length(initial_voxels) == 0){
    order_t <- order(smoothed_t)
    initial_voxels <- c(head(order_t, floor(V / 50)), tail(order_t, floor(V / 50)))
  }
  initial_voxels_positive <- initial_voxels[which(smoothed_t[initial_voxels] > 0)]
  initial_voxels_negative <- initial_voxels[which(smoothed_t[initial_voxels] < 0)]
  region_list_initial <- list()
  if(length(initial_voxels_positive) == 0 & length(initial_voxels_negative) == 0){
    region_list_initial[[1]] <- list(NULL, 0)
  }else if(length(initial_voxels_positive) == 0){
    region_list_initial[[1]] <- list(initial_voxels_negative, mean(naive_betahat[initial_voxels_negative]))
  }else if(length(initial_voxels_negative) == 0){
    region_list_initial[[1]] <- list(initial_voxels_positive, mean(naive_betahat[initial_voxels_positive]))
  }else{
    region_list_initial[[1]] <- list(initial_voxels_negative, mean(naive_betahat[initial_voxels_negative]))
    region_list_initial[[2]] <- list(initial_voxels_positive, mean(naive_betahat[initial_voxels_positive]))
  }
  return(region_list_initial)
}

num_voxels <- function(region_list_full){
  if(all(null_tracker)){
    return(0)
  }
  else{
    total_num_voxels <- 0
    for(i in 1:length(p1_index)){
      region_list <- region_list_full[[p1_index[i]]]
      for(j in 1:length(region_list)){
        total_num_voxels <- total_num_voxels + length(region_list[[j]][[1]]) 
      }
    }
    return(total_num_voxels)
  }
}

num_boundary_voxels <- function(region_list_full){
  if(all(null_tracker)){
    return(0)
  }
  else{
    total_boundary_voxels <- 0
    for(i in 1:length(p1_index)){
      region_list <- region_list_full[[p1_index[i]]]
      for(j in 1:length(region_list)){
        voxels_in_regions <- voxels_in_regions(region_list)
        region_voxels <- region_list[[j]][[1]]
        for(v in region_voxels){
          v_neighbors <- neighbor_list[[v]]
          # if(any(!(v_neighbors %in% region_voxels))){
          #   total_boundary_voxels <- total_boundary_voxels + 1
          # }
          if(any(!(v_neighbors %in% voxels_in_regions))){
            total_boundary_voxels <- total_boundary_voxels + 1
          }
        }
      }
    }
    return(total_boundary_voxels)
  }
}

compute_AIC_BIC <- function(Y, Q, alphahat_full, region_list_full){
  total_num_regions <- 0
  num_voxels <- num_voxels(region_list_full)
  num_boundary_voxels <- num_boundary_voxels(region_list_full)
  for(i in p1_index){
    total_num_regions <- total_num_regions + length(region_list_full[[i]])
  }
  MSE <- mean((Y - Q %*% alphahat_full) ^ 2)
  AIC <- 2 * total_num_regions + V * n * log(MSE)
  BIC_regions <- log(V * n) * total_num_regions + V * n * log(MSE)
  BIC_boundary_voxels <- log(V * n) * num_boundary_voxels + V * n * log(MSE)
  BIC_voxels <- log(V * n) * num_voxels + V * n * log(MSE)
  output <- c(AIC, BIC_regions, BIC_boundary_voxels, BIC_voxels, MSE, total_num_regions, num_voxels, num_boundary_voxels)
  names(output) <- c("AIC", "BIC (Regions)", "BIC (Boundary Voxels)", "BIC (Voxels)", "MSE", "Number of Regions", "Number of Voxels", "Number of Boundary Voxels")
  return(output)
}


region_update_algorithm <- function(initial_region_list, initial_betahat, which_beta){
  terminate <- F
  inner_iter <<- 0
  new_region_list <- initial_region_list
  new_betahat <- initial_betahat
  while(inner_iter <= 500){
    new_info <- update_region_list(new_region_list, which_beta)
    new_region_list <- new_info[[1]]
    new_betahat <- region_list_to_betahat(new_region_list)
    other_info <- new_info[[2]]
    terminate <- other_info[[1]]
    if(terminate){
      return(new_region_list)
    }
    else{
      c1 <- other_info[[2]]
      v1 <- other_info[[3]]
      c2 <- other_info[[4]]
      v2 <- other_info[[5]]
      changing_voxels <- v1
      change_amount <- rep(c1, length(v1))
      if(!is.null(v2)){
        changing_voxels <- c(changing_voxels, v2)
        change_amount <- c(change_amount, rep(c2, length(v2)))
      }
      delta_betahat <- numeric(V)
      delta_betahat[v1] <- c1
      if(!is.null(c2)){
        delta_betahat[v2] <- c2
      }
      XTXB[ , changing_voxels] <<- XTXB[ , changing_voxels] + outer(XTX[ , which_beta], delta_betahat[changing_voxels])
      inner_iter <<- inner_iter + 1
    }
  }
  warning(sprintf("region_update_algorithm hit max iterations (500) for beta index %d without converging.", which_beta))
  return(new_region_list)
}

update_gammahat <- function(alphahat_full, thetahat){
  return(ZTZ1ZTY - ZTZ1ZTX %*% alphahat_full[p1_index , , drop = F] - (ZTZ1ZT %*% thetahat) %*% Psi_star)
}

update_thetahat <- function(X, Z, alphahat_full){
  thetahat <- (Y_Psi_starT - X %*% tcrossprod(alphahat_full[p1_index, , drop = F], Psi_star) - Z %*% tcrossprod(alphahat_full[p2_index, , drop = F], Psi_star)) %*% scaling_matrix_2
  return(thetahat)
}

update_XTY_tilde <- function(gammahat, thetahat){
  XTZ_gammahat <- XTZ %*% gammahat
  XT_thetahat <- crossprod(X, thetahat)
  XT_thetahat_Psi_star <- XT_thetahat %*% Psi_star
  XTY_tilde <- XTY - XTZ_gammahat - XT_thetahat_Psi_star
  return(XTY_tilde)
}


compute_full_loss <- function(alphahat_full, thetahat){
  ms_loss <- sum((Y - X %*% alphahat_full[p1_index, , drop = F] - 
                    Z %*% alphahat_full[p2_index, , drop = F] - 
                    thetahat %*% Psi_star) ^ 2) / n
  beta_diffs <- alphahat_full[p1_index, edge_list[, 1], drop = FALSE] -
    alphahat_full[p1_index, edge_list[, 2], drop = FALSE]
  neighbor_loss <- sum(abs(beta_diffs))
  l1_loss <- sum(abs(alphahat_full[p1_index, ]))
  full_loss <- ms_loss + lambda * neighbor_loss + mu * l1_loss
  return(c(ms_loss, neighbor_loss, l1_loss, full_loss))
}

full_algorithm <- function(max_iter = 50, initial_full_region_list){
  # debug code
  # print(paste0("STARTING LOSS: ", compute_full_loss(alphahat_full, thetahat)[4]))
  #
  outer_iter <<- 0
  null_tracker <<- (1:p) %in% p2_index
  print(paste("ITERATION", outer_iter, sep = " "))
  for(j in 1:p){
    if(j %in% p1_index){
      print(paste("UPDATING BETA", j, sep = " "))
      which_beta <- which(p1_index == j)
      initial_region_list <- initial_full_region_list[[j]]
      initial_betahat <- region_list_to_betahat(initial_region_list)
      region_list_full[[j]] <<- region_update_algorithm(initial_region_list, initial_betahat, which_beta)
      alphahat_full[j , ] <<- region_list_to_betahat(region_list_full[[j]])
      if(length(region_list_full[[j]]) == 0){
        null_tracker[j] <<- T
        region_list_full[[j]] <<- list(list(NULL, 0))
      }
      else if(is.null(region_list_full[[j]][[1]][[1]])){
        null_tracker[j] <<- T
      }
    }
    else{
      region_list_full[[j]] <<- list(list(NULL, 0))
    }
  }
  gammahat <<- update_gammahat(alphahat_full, thetahat)
  alphahat_full[p2_index , ] <<- gammahat
  # debug code
  # print(paste0("AFTER UPDATING GAMMAHAT LOSS: ", compute_full_loss(alphahat_full, thetahat)[4]))
  #
  thetahat <<- update_thetahat(X, Z, alphahat_full)
  # debug code
  # print(paste0("AFTER UPDATING THETAHAT LOSS: ", compute_full_loss(alphahat_full, thetahat)[4]))
  #
  if(all(null_tracker)){
    return(list(region_list_full, alphahat_full, thetahat))
  }
  XTY_tilde <<- update_XTY_tilde(gammahat, thetahat) 
  iterate <- T
  outer_iter <<- outer_iter + 1
  while(iterate & outer_iter <= max_iter & !all(null_tracker)){
    print(paste("ITERATION", outer_iter, sep = " "))
    previous_alphahat_full <- alphahat_full
    for(j in p1_index){
      print(paste("UPDATING BETA", j, sep = " "))
      which_beta <- which(p1_index == j)
      initial_betahat <- region_list_to_betahat(region_list_full[[j]])
      region_list_full[[j]] <<- region_update_algorithm(region_list_full[[j]], initial_betahat, which_beta)
      alphahat_full[j , ] <<- region_list_to_betahat(region_list_full[[j]])
      if(length(region_list_full[[j]]) == 0){
        null_tracker[j] <<- T
        region_list_full[[j]] <<- list(list(NULL, 0))
      }
      else if(is.null(region_list_full[[j]][[1]][[1]])){
        null_tracker[j] <<- T
      }
    }
    gammahat <<- update_gammahat(alphahat_full, thetahat)
    alphahat_full[p2_index , ] <<- gammahat
    # debug code
    # print(paste0("AFTER UPDATING GAMMAHAT LOSS: ", compute_full_loss(alphahat_full, thetahat)[4]))
    #
    thetahat <<- update_thetahat(X, Z, alphahat_full)
    # debug code
    # print(paste0("AFTER UPDATING THETAHAT LOSS: ", compute_full_loss(alphahat_full, thetahat)[4]))
    #
    XTY_tilde <<- update_XTY_tilde(gammahat, thetahat)
    alphahat_norm_diff <- norm(alphahat_full - previous_alphahat_full, type = "F")
    print(paste("NORM DIFF: ", alphahat_norm_diff))
    if(alphahat_norm_diff < 0.01){
      iterate <- F
    }
    outer_iter <<- outer_iter + 1
  }
  return(list(region_list_full, alphahat_full, thetahat))
}
