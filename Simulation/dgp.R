# Lee, Shi and Gao (2020) On Lasso for predictive regression
# Data generating process (DGP) for Monte Carlo experiments

# ------------------- Pure unit roots --------------
generate_data_pure_unit_root <-
  function(n, coef_0, sigma_w_mat, phi = NULL, ur_index = c(1, 2, 3), n_burn_in = 1000) {
    # n_burn_in should be larger than 4
    
    a_0 <- coef_0[1]
    b_0 <- coef_0[-1]
    b_0[ur_index] <- b_0[ur_index] / sqrt(n)
    p <- length(b_0)
    
    # Generate innovations, denoted as w
    w_mat <- mvtnorm::rmvnorm(n = n + n_burn_in, sigma = sigma_w_mat)
    w_mat <- w_mat[-(1:(n_burn_in - 3)), ]
    
    x_mat <- apply(w_mat[, 1:p], 2, cumsum)[-c(1, nrow(w_mat)), ]
    u <- w_mat[-(1:2), p + 1]
    y <- a_0 + x_mat %*% b_0 + u
    
    return(list(
      y_est = y[1:n],
      x_est = x_mat[1:n, ],
      y_pred = y[length(y)],
      x_pred = x_mat[nrow(x_mat), ]
    ))
  }

# ------------------- Generate innovation ------------------
generate_innovation <- function(n, sigma_w_mat, phi, lag_depth = 0, n_burn_in = 1000) {
  w_mat <- mvtnorm::rmvnorm(n = n + n_burn_in, sigma = sigma_w_mat)
  ar_coef <- diag(phi)
  xi_mat <- w_mat
  for (i in 2:(n + n_burn_in)) {
    xi_mat[i, ] <- ar_coef %*% xi_mat[i - 1, ] + w_mat[i, ]
  }
  xi_mat <- xi_mat[-(1:(n_burn_in - (3 + lag_depth))), ]
  
  xi_mat
}

# ------------------- Mixed roots -------------------
generate_data_mixed_roots <-
  function(n,
           coef_0,
           sigma_w_mat,
           phi,
           ur_index = c(8, 9),
           p_z = 3,
           p_c = 4,
           p_x = 5,
           n_burn_in = 1000) {
    
    a_0 <- coef_0[1]
    b_0 <- coef_0[-1]
    b_0[ur_index] <- b_0[ur_index] / sqrt(n)
    p <- p_z + p_c + p_x
    
    # Generate z_i, v_i, e_i, u_i
    xi_mat <-
      generate_innovation(n, sigma_w_mat, phi, n_burn_in = n_burn_in)
    z_mat <- xi_mat[, 1:p_z, drop = FALSE]
    v_mat <- xi_mat[, (p_z + 1):(p_z + p_c), drop = FALSE]
    e_mat <- xi_mat[, (p_z + p_c + 1):(p_z + p_c + p_x), drop = FALSE]
    u <- xi_mat[, p + 1]
    
    # Generate x_1_i, x_c_i
    x_1_mat <- apply(e_mat, 2, cumsum)
    x_c_mat <- v_mat
    x_c_mat[, c(1, 3)] <- apply(v_mat[, c(1, 3)], 2, cumsum)
    x_c_mat[, c(2, 4)] <- x_c_mat[, c(1, 3)] - v_mat[, c(2, 4)]
    
    # Generate x_i, y_i
    x_mat <- cbind(z_mat, x_c_mat, x_1_mat)[-c(1, nrow(xi_mat)),]
    u <- u[-(1:2)]
    y <- a_0 + x_mat %*% b_0 + u
    
    list(
      y_est = y[1:n],
      x_est = x_mat[1:n,],
      y_pred = y[length(y)],
      x_pred = x_mat[nrow(x_mat),]
    )
  }

# ------------------- Local_to_unity -------------
generate_data_lur <-
  function(n,
           coef_0,
           sigma_w_mat,
           phi,
           ur_index = c(8, 11),
           phi_x = c(1 - 7.493427 / n, 1 - 5.915274 / n, 1 - 8.927722 / n, 1, 1),
           phi_v = c(1 - 4.838619 / n, 1),
           p_z = 3,
           p_c = 4,
           p_x = 5,
           n_burn_in = 1000) {
    
    a_0 <- coef_0[1]
    b_0 <- coef_0[-1]
    b_0[ur_index] <- b_0[ur_index] / sqrt(n)
    p <- p_z + p_c + p_x
    
    # Generate z_i, v_i, e_i, u_i
    xi_mat <-
      generate_innovation(n, sigma_w_mat, phi, n_burn_in = n_burn_in)
    z_mat <- xi_mat[, 1:p_z, drop = FALSE]
    v_mat <- xi_mat[, (p_z + 1):(p_z + p_c), drop = FALSE]
    e_mat <- xi_mat[, (p_z + p_c + 1):(p_z + p_c + p_x), drop = FALSE]
    u <- xi_mat[, p + 1]
    
    # Generate x_1_i, x_c_i
    x_1_mat <- e_mat
    x_c_mat <- v_mat
    for (i in 2:nrow(x_1_mat)) {
      x_1_mat[i, ] <- c(diag(phi_x) %*% x_1_mat[i - 1, ]) + e_mat[i,]
      x_c_mat[i, c(1, 3)] <- c(diag(phi_v) %*% x_c_mat[i - 1, c(1, 3)]) + v_mat[i, c(1, 3)]
    }
    x_c_mat[, c(2, 4)] <- x_c_mat[, c(1, 3)] - v_mat[, c(2, 4)]
    
    # Generate x_i, y_i
    x_mat <- cbind(z_mat, x_c_mat, x_1_mat)[-c(1, nrow(xi_mat)),]
    u <- u[-(1:2)]
    y <- a_0 + x_mat %*% b_0 + u
    
    list(
      y_est = y[1:n],
      x_est = x_mat[1:n,],
      y_pred = y[length(y)],
      x_pred = x_mat[nrow(x_mat),]
    )
    
  }

# ---------- stationary auto-regression with mixed roots ------
generate_data_ardl <-
  function(n,
           coef_0,
           sigma_w_mat,
           phi,
           ur_index = c(10, 13),
           phi_x = c(1 - 7.493427 / n, 1 - 5.915274 / n, 1 - 8.927722 / n, 1, 1),
           phi_v = c(1 - 4.838619 / n, 1),
           p_z = 3,
           p_c = 4,
           p_x = 5,
           n_burn_in = 1000) {
    
    a_0 <- coef_0[1]
    b_0 <- coef_0[-1]
    b_0[ur_index] <- b_0[ur_index] / sqrt(n)
    p <- p_z + p_c + p_x
    
    # Generate z_i, v_i, e_i, u_i
    xi_mat <-
      generate_innovation(n, sigma_w_mat, phi, lag_depth = 1, n_burn_in = n_burn_in)
    z_mat <- xi_mat[, 1:p_z, drop = FALSE]
    v_mat <- xi_mat[, (p_z + 1):(p_z + p_c), drop = FALSE]
    e_mat <- xi_mat[, (p_z + p_c + 1):(p_z + p_c + p_x), drop = FALSE]
    u <- xi_mat[, p + 1]
    
    # Generate x_1_i, x_c_i
    x_1_mat <- e_mat
    x_c_mat <- v_mat
    for (i in 2:nrow(x_1_mat)) {
      x_1_mat[i, ] <- c(diag(phi_x) %*% x_1_mat[i - 1, ]) + e_mat[i,]
      x_c_mat[i, c(1, 3)] <- c(diag(phi_v) %*% x_c_mat[i - 1, c(1, 3)]) + v_mat[i, c(1, 3)]
    }
    x_c_mat[, c(2, 4)] <- x_c_mat[, c(1, 3)] - v_mat[, c(2, 4)]
    
    # Generate x_i, y_i
    m <- nrow(x_1_mat)
    x_mat <- cbind(x_c_mat[-c(1, m), ],
                   x_c_mat[-c(m - 1, m), ],
                   x_1_mat[-c(1, m), ],
                   x_1_mat[-c(m - 1, m), ],
                   z_mat[-c(1, m), ],
                   z_mat[-c(m - 1, m), ])
    u <- u[-c(1,2)]
    
    y <- rep(0, length(u))
    y[1] <- a_0 + sum(x_mat[1, ] * b_0[-1]) + u[1]
    for (i in 2:length(y)) {
      y[i] <- a_0 + b_0[1] * y[i - 1] + sum(x_mat[i, ] * b_0[-1]) + u[i]
    }
    
    x_mat <- cbind(y[-length(y)], x_mat[-1, ])
    y <- y[-1]
    
    list(
      y_est = y[1:n],
      x_est = x_mat[1:n,],
      y_pred = y[length(y)],
      x_pred = x_mat[nrow(x_mat),]
    )
  }

