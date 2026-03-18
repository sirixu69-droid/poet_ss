library(matrixStats) 
library(parallel)
library(MASS)
library("SPCAvRP")
library("mvtnorm")
library("mnormt")
library(ICSNP)
library(ICS)
library("SpatialNP")
library("MNM")
library(flare)
#library(glasso)
library(glassoFast)
library(stats)
norm2=function(x){sqrt(sum(x^2))}


# Project a matrix to the PSD cone
eigen_projection <- function(A, epsilon = 1e-5) {
  eig <- eigen(A, symmetric = TRUE)
  values <- pmax(eig$values, epsilon)
  return(eig$vectors %*% diag(values) %*% t(eig$vectors))
}


# Select gamma
select_gamma <- function(A, mu) {
  A_abs_vec <- as.vector(abs(A)) / mu
  A_sorted  <- sort(A_abs_vec, decreasing = TRUE)
  d2 <- length(A_sorted) 
  cums <- cumsum(A_sorted)           # cums[u] = sum_{j=1}^u A_sorted[j] 
  u_idx <- 1:(d2 - 1)
  sum_u      <- cums[u_idx]        - u_idx     * A_sorted[u_idx]
  sum_u_plus <- cums[u_idx + 1L]   - (u_idx+1) * A_sorted[u_idx + 1L]  
  idx <- which(sum_u < 1 & sum_u_plus >= 1)
  if (length(idx) > 0) {
    u <- idx[1]
    gamma <- (cums[u] - 1) / u
    return(list(gamma = gamma, u = u, threshold = mu * gamma))
  } 
  warning("gamma not found, using minimal element.")
  gamma <- min(A_sorted)
  return(list(gamma = gamma, u = d2, threshold = mu * gamma))
}
# smoothed l_inf norm approximation
smooth_linf <- function(A, mu) {
  gamma_obj <- select_gamma(A, mu)
  gamma <- gamma_obj$gamma
  U_hat <- sign(A) * pmax(abs(A) / mu - gamma, 0)
  return(sum(U_hat * A) - mu * norm(U_hat, "F")^2 / 2)
}

Q_value <- function(mu, S_hat, S, M, G, eta, theta) {
  Delta <- S - M
  inner_prod <- sum(G * Delta)
  norm_sq <- norm(Delta, "F")^2
  return(smooth_linf(S_hat - M, mu) + inner_prod + (1 / (2*eta * theta)) * norm_sq)
}

accelerated_proximal_gradient <- function(S_hat, mu, max_iter = 300, tol = 1e-4) {
  p <- nrow(S_hat)
  S0 <- eigen_projection(S_hat)
  S_prev <- S0
  W_prev <- S0     
  for (t in 1:max_iter) {
    theta_t <- 2 / (1 + t)   
    M_t <- (1 - theta_t) * S_prev + theta_t * W_prev   
    Delta <- S_hat - M_t
    gamma <- select_gamma(Delta, mu)$gamma
    G_t <- -sign(Delta) * pmax(abs(Delta / mu) - gamma, 0)
    eta <- mu
    step_size <- eta / theta_t
    W_temp <- eigen_projection(W_prev - step_size * G_t)
    S_temp <- (1 - theta_t) * S_prev + theta_t * W_temp 
    W_prev <- W_temp
    S_curr <- S_temp    
    gap <- abs(smooth_linf(S_hat - S_curr, mu) - smooth_linf(S_hat - S_prev, mu))
    if (gap <= tol * mu) {
      break
    }    
    S_prev <- S_curr
  }
  return(S_curr)
}



adaptive_threshold <- function(S, tau) {
  diag_S <- diag(S)
  diag_pos <- pmax(diag_S, 1e-8)  
  
  tau_mat <- outer(diag_pos, diag_pos, function(x, y) tau * sqrt(x * y))
  S_offdiag <- S
  diag(S_offdiag) <- 0
  
  S_thresh <- sign(S_offdiag) * pmax(abs(S_offdiag) - tau_mat, 0)
  diag(S_thresh) <- diag_S  
  return(S_thresh)
}

robust_variance_diag <- function(Y, delta0 = 1e-4) {
  n <- nrow(Y)
  p <- ncol(Y)
  
  h_func <- function(x) {
    sign(x) * log(1 + abs(x) + x^2 / 2)
  }
  
  
  catoni_m_estimator <- function(z, v, eps = 1 / (max(n, p)^2)) {
    alpha <- sqrt(2 * log(1 / eps) / (n * (v + 2 * v * log(1 / eps) / (n - 2 * log(1 / eps)))))
    h_sum <- function(mu) {
      sum(h_func(alpha * (z - mu)))
    }
    
    mu_hat <- tryCatch({
      uniroot(h_sum, interval = c(min(z) - 10 * sd(z), max(z) + 10 * sd(z)))$root
    }, error = function(e) {
      mean(z)  
    })
    
    return(mu_hat)
  }
  
  
  mu_hat <- numeric(p)
  eta_hat <- numeric(p)
  sigma2_hat <- numeric(p)
  
  for (j in 1:p) {
    yj <- Y[, j]
    vj <- 3 * var(yj)
    
    mu_hat[j] <- catoni_m_estimator(yj, vj)
    yj_sq <- yj^2
    vj2 <- 3 * var(yj^3)
    eta_hat[j] <- catoni_m_estimator(yj_sq, vj2)
    sigma2_hat[j] <- max(eta_hat[j] - mu_hat[j]^2, delta0)
  }
  
  D_hat <- diag(sqrt(sigma2_hat))
  return(D_hat)
}
estimate_Lambda <- function(Y, D_hat, m) {
  n <- nrow(Y)
  p <- ncol(Y)
  upper_idx <- upper.tri(matrix(0, n, n))
  s_matrix <- matrix(NA_real_, nrow = sum(upper_idx), ncol = p)
  for (j in seq_len(p)) {
    Yj <- Y[, j]
    diff_mat <- matrix(Yj, n, n) - matrix(Yj, n, n, byrow = TRUE)
    s_matrix[, j] <- sign(diff_mat)[upper_idx]
  }
  
  coeff <- 2 / (n * (n - 1))
  T_hat <- coeff * crossprod(s_matrix)
  diag(T_hat) <- 1
  R_hat <- sin(pi / 2 * T_hat)
  Sigma1_hat <- D_hat %*% R_hat %*% D_hat
  eig_vals <- eigen(Sigma1_hat, symmetric = TRUE, only.values = TRUE)$values
  lambda_topm <- eig_vals[1:m]
  Lambda_hat <- diag(lambda_topm, nrow = m, ncol = m)
  return(list(
    Lambda_hat = Lambda_hat,
    Sigma1_hat = Sigma1_hat
  ))
}
estimate_mu_huber_cv <- function(Y, c_candidates = c(0.5, 1.0, 1.5), nfold = 3) {
  n <- nrow(Y)  
  p <- ncol(Y)  
  logN <- log(p)
  
  mu_hat <- numeric(p)
  best_c <- numeric(p)
  for (i in 1:p) {
    y <- Y[, i]
    scale_full <- mad(y)
    fold_id <- sample(rep(1:nfold, length.out = n))
    cv_loss <- numeric(length(c_candidates))
    
    for (j in seq_along(c_candidates)) {
      c_val <- c_candidates[j]
      H <- c_val * scale_full * sqrt(n / logN) 
      
      fold_errors <- numeric(nfold)
      for (k in 1:nfold) {
        y_train <- y[fold_id != k]
        y_test  <- y[fold_id == k]
        scale_train <- mad(y_train)
        mu_train <- MASS::huber(y_train, k = H / scale_train)$mu
        fold_errors[k] <- mean((y_test - mu_train)^2)
      }
      cv_loss[j] <- mean(fold_errors)
    }
    j_best <- which.min(cv_loss)
    best_c[i] <- c_candidates[j_best]
    H_best <- best_c[i] * scale_full * sqrt(n / logN)
    mu_hat[i] <- MASS::huber(y, k = H_best / scale_full)$mu
  }
  
  return(list(mu = mu_hat, c_opt = best_c))
}


estimate_Exi2 <- function(Y, 
                          mu_hat, 
                          Sigma_shape = NULL,
                          epsilon = 0.345, 
                          c_candidates = seq(0.3, 1.0, length.out = 10), 
                          nfold = 3) {
  n <- nrow(Y)
  p <- ncol(Y)
  Y_centered <- sweep(Y, 2, mu_hat, FUN = "-")
  if (is.null(Sigma_shape)) {
    Sigma_shape <- diag(p)
  }
  if (isTRUE(all.equal(Sigma_shape, diag(p)))) {
    xi_sq <- rowSums(Y_centered^2) / p
  } else {
    eig_shape <- eigen(Sigma_shape, symmetric = TRUE)
    U   <- eig_shape$vectors
    lam <- eig_shape$values
    
    Z <- Y_centered %*% U     # n × p
    xi_sq <- rowSums(Z^2 * rep(lam, each = n)) / p
  }
  scale_full <- stats::mad(xi_sq)  
  fold_id <- sample(rep(1:nfold, length.out = n))
  cv_loss <- numeric(length(c_candidates))
  for (j in seq_along(c_candidates)) {
    c_val <- c_candidates[j]
    H <- c_val * n^min(1 / (1 + epsilon / 2), 1/2)  
    
    fold_errors <- numeric(nfold)
    for (k in 1:nfold) {
      train_xi <- xi_sq[fold_id != k]
      test_xi  <- xi_sq[fold_id == k]
      scale_train <- stats::mad(train_xi)
      huber_fit <- MASS::huber(train_xi, k = H / scale_train)
      mu_train  <- huber_fit$mu
      fold_errors[k] <- mean((test_xi - mu_train)^2)
    }
    cv_loss[j] <- mean(fold_errors)
  }
  
  j_best <- which.min(cv_loss)
  c_opt  <- c_candidates[j_best]
  H_best <- c_opt * n^min(1 / (1 + epsilon / 2), 1/2)
  huber_fit_final <- MASS::huber(xi_sq, k = H_best / scale_full)
  
  list(
    Exi2      = huber_fit_final$mu, 
    c_opt     = c_opt,
    H_best    = H_best,
    xi_sq     = xi_sq,
    cv_loss   = cv_loss,
    Sigma_shape_used = Sigma_shape
  )
}
regTME_svd <- function(Y, alpha = NULL, tol = 1e-4, maxit = 100, eps = 1e-12) {
  n <- nrow(Y)
  p <- ncol(Y)
  gamma <- p / n
  S_sample <- cov(Y)               
  S_sample_norm <- (p * S_sample) / sum(diag(S_sample))  
  s_max <- norm(S_sample_norm, type = "2")
  if (is.null(alpha)) {
    alpha <- max(0.1, 1.1 * (gamma - 1 + s_max * (1 + sqrt(gamma))^2))
  }
  
  if (n < p) {
    svd_Y <- svd(Y)            
    U <- svd_Y$u               
    D <- diag(svd_Y$d)         
    V <- svd_Y$v               
    Z <- U %*% D               
    proj_mat <- V             
    dim_Z <- n
    use_proj <- TRUE
  } else {
    Z <- Y
    proj_mat <- diag(p)
    dim_Z <- p
    use_proj <- FALSE
  }
  Sigma <- (alpha / (1 + alpha)) * diag(dim_Z)
  for (iter in 1:maxit) {
    Sigma_inv <- solve(Sigma)
    quad_forms <- rowSums((Z %*% Sigma_inv) * Z) + eps   
    W <- 1 / quad_forms                                   
    Zw <- Z * sqrt(W)                                     
    updated <- crossprod(Zw)                              
    Sigma_new <- (dim_Z / n) * updated
    Sigma_new <- Sigma_new / (1 + alpha) + (alpha / (1 + alpha)) * diag(dim_Z)
    if (norm(Sigma_new - Sigma, type = "F") < tol) {
      Sigma <- Sigma_new
      break
    }
    Sigma <- Sigma_new
  }
  A_low <- Sigma - (alpha / (1 + alpha)) * diag(dim_Z)
  if (use_proj) {
    A_full <- proj_mat %*% A_low %*% t(proj_mat)  # p x p
  } else {
    A_full <- A_low
  }
  A_full <- (A_full + t(A_full)) / 2
  den <- sum(diag(A_full))
  if (abs(den) < eps) {
    warning("trace(A_full) 非常小，加入 eps 以避免除零")
    den <- den + sign(den + eps) * eps
  }
  Sigma_out <- p * A_full / den   
  return(Sigma_out)
}

generate_data_four_cases <- function(p, n, m, mi, 
                                     s = c(1, 0.75^2, 0.5^2)) {
  
  # Sigma_u = (0.9^{|i-j|})
  indices <- matrix(1:p, nrow = p, ncol = p)
  Sigma_u <- 0.9^abs(indices - t(indices))
  
  # loading matrix B
  B <- matrix(0, nrow = p, ncol = m)
  for (k in 1:m) {
    B[, k] <- rnorm(p, mean = 0, sd = sqrt(s[k]))
  }
  
  # raw covariance/scatter of y
  Sigma <- B %*% t(B) + Sigma_u
  
  # normalization constant c = p / tr(BB' + Sigma_u)
  c0 <- p / sum(diag(Sigma))
  
  # true covariance of (f,u): Sigma_fu = c * diag(I_m, Sigma_u)
  Sigma_fu <- matrix(0, nrow = m + p, ncol = m + p)
  Sigma_fu[1:m, 1:m] <- c0 * diag(m)
  Sigma_fu[(m + 1):(m + p), (m + 1):(m + p)] <- c0 * Sigma_u
  
  # true covariance of y and u
  Sigma0 <- c0 * Sigma
  Sigma_u0 <- c0 * Sigma_u
  
  # generate (f, u)
  if (mi == 1) {
    # Scenario I: multivariate normal
    joint_sample <- mvtnorm::rmvnorm(n, mean = rep(0, m + p), sigma = Sigma_fu)
    
  } else if (mi == 2) {
    # Scenario II: multivariate t, df = 4
    df <- 4
    joint_sample <- mvtnorm::rmvt(n, sigma = Sigma_fu * (df - 2) / df, df = df)
    
  } else if (mi == 3) {
    # Scenario III: multivariate t, df = 2.2
    df <- 2.2
    joint_sample <- mvtnorm::rmvt(n, sigma = Sigma_fu * (df - 2) / df, df = df)
    
  } else if (mi == 4) {
    # Scenario IV: 0.8 N(0, Sigma_fu) + 0.2 N(0, 10 Sigma_fu)
    z <- rbinom(n, size = 1, prob = 0.2)
    joint_sample <- matrix(0, nrow = n, ncol = m + p)
    
    n1 <- sum(z == 0)
    n2 <- sum(z == 1)
    
    if (n1 > 0) {
      joint_sample[z == 0, ] <- mvtnorm::rmvnorm(
        n1, mean = rep(0, m + p), sigma = Sigma_fu
      )
    }
    if (n2 > 0) {
      joint_sample[z == 1, ] <- mvtnorm::rmvnorm(
        n2, mean = rep(0, m + p), sigma = 10 * Sigma_fu
      )
      
    }
    joint_sample <- joint_sample / sqrt(2.8)
  } else {
    stop("mi must be one of 1, 2, 3, 4.")
  }
  
  # construct y = Bf + u
  f_t <- joint_sample[, 1:m, drop = FALSE]
  u_t <- joint_sample[, (m + 1):(m + p), drop = FALSE]
  Y <- f_t %*% t(B) + u_t
  
  list(
    Y = Y,
    Sigma0 = Sigma0,
    Sigma_u0 = Sigma_u0
  )
}
compute_2_part_nonsparse <- function(Y, m, method = c("SAMPLE", "FLW","TylerSPCA","IPSN")) {
  method <- match.arg(method)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SAMPLE") {
    S_hat <- cov(Y)
    eig <- eigen(S_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- S_hat - low_rank_part
    mu <- sqrt(log(p) / n)
    tau <- 1 * (mu +1/sqrt(p)) 
    R_hat <- adaptive_threshold(R_tilde, tau = tau)
    Sigma_hat <- low_rank_part+R_hat
    Sigma_u_hat <- R_hat
    Sigma_hat  <- (Sigma_hat  + t(Sigma_hat )) / 2
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
  }
  
  if (method == "FLW") {
    D_hat <- robust_variance_diag(Y)
    estimate <- estimate_Lambda(Y, D_hat, m)
    Lambda_hat <- estimate$Lambda_hat
    S_hat <- estimate$Sigma1_hat
    HK<-RCov(Y)
    eig <- eigen(HK, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- S_hat - low_rank_part
    R_tilde1 <- eigen_projection(R_tilde,epsilon = 1e-5)
    max_diff <- max(abs(R_tilde - R_tilde1))
    mu <- sqrt(log(p) / n)
    tau <- 1 * (mu +1/sqrt(p))  
    if (max_diff < 2*mu) {
      R_tilde2 <- R_tilde1
    } else {
      R_tilde2 <- accelerated_proximal_gradient(R_tilde, mu, max_iter = 300, tol = 5e-1)
    }
    R_hat <- adaptive_threshold(R_tilde2, tau = tau)
    Sigma_hat <- low_rank_part+R_hat
    Sigma_u_hat <- R_hat
    Sigma_hat  <- (Sigma_hat  + t(Sigma_hat )) / 2
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
  }
  if (method == "TylerSPCA") {
    S_hat <- SCov(Y)*p
    eig <- eigen(S_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- S_hat - low_rank_part
    mu <- sqrt(log(p) / n)
    tau <- 0.8*(mu +sqrt(log(n)/n))
    R_hat <- adaptive_threshold(R_tilde, tau = tau)
    R_hat1 <- eigen_projection(R_hat,epsilon = 1e-5)
    Sigma_hat <- low_rank_part+R_hat
    mu_hat <- spatial.median(Y)
    Y_centered <- sweep(Y, 2, mu_hat, FUN = "-")  
    Sigma_S_inv <- solve(Sigma_hat)
    quad_forms <- rowSums((Y_centered %*% Sigma_S_inv) * Y_centered) 
    weights <- 1 / quad_forms  
    Sigma_T <- (p / n) * t(Y_centered) %*% diag(weights) %*% Y_centered
    Sigma_T <- (Sigma_T + t(Sigma_T)) / 2  
    eig_Sigma_T <- eigen(Sigma_T, symmetric = TRUE)
    Gamma_hat <- eig_Sigma_T$vectors[, 1:m, drop = FALSE]
    if(length(eig_Sigma_T$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig_Sigma_T$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig_Sigma_T$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde_T <- Sigma_T - low_rank_part
    tau <- 0.8*(mu +sqrt(log(n)/n))
    R_hat_T <- adaptive_threshold(R_tilde_T, tau = tau)
    Sigma_hat <- low_rank_part + R_hat_T
    Sigma_u_hat <- R_hat_T
    Sigma_hat  <- (Sigma_hat  + t(Sigma_hat )) / 2
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Exi2 <- estimate_Exi2(Y, mu_hat, Sigma_shape = Sigma_S_inv)$Exi2
    Sigma_hat <- Exi2*Sigma_hat
    Lambda_hat <- Exi2*Lambda_hat
    Sigma_u_hat <- Exi2*Sigma_u_hat
  }
  if (method == "IPSN") {
    mu_hat <- estimate_mu_huber_cv(Y)$mu
    HK<-RCov(Y)
    eig <- eigen(HK, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE]
    P_I_basis <- qr.Q(qr(Gamma_hat), complete = TRUE)[, (ncol(Gamma_hat) + 1):nrow(Gamma_hat)]
    X_hat <- matrix(0, nrow = n, ncol = p) 
    for (t in 1:n) {
      y_centered <- Y[t, ] - mu_hat  
      proj <- crossprod(P_I_basis, y_centered)  
      norm_proj <- sqrt(sum(proj^2))
      X_hat[t, ] <- sqrt(p) * y_centered / norm_proj
    }
    Sigma_0_raw <- matrix(0, nrow = p, ncol = p)
    for (t in 1:n) {
      x_t <- X_hat[t, , drop = FALSE]  
      Sigma_0_raw <- Sigma_0_raw + crossprod(x_t)  
    }
    Sigma_0_raw <- Sigma_0_raw / n
    eta_hat <- p / sum(diag(Sigma_0_raw))
    Sigma_0_hat <- Sigma_0_raw * eta_hat 
    eig_Sigma0 <- eigen(Sigma_0_hat, symmetric = TRUE)
    Gamma_hat <- eig_Sigma0$vectors[, 1:m] 
    if(length(eig_Sigma0$values[1:m]) == 1) {
      Lambda_hat0 <- matrix(eig_Sigma0$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat0 <- diag(eig_Sigma0$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)
    mu <-  sqrt(log(p)/n)
    R_tilde <- Sigma_0_hat - low_rank_part
    tau <- 1 *(mu +sqrt(log(p)/p)) 
    Sigma_u_hat0 <- adaptive_threshold(R_tilde, tau = tau)
    Sigma_hat0 <- low_rank_part+Sigma_u_hat0
    Sigma_hat0 <- (Sigma_hat0 + t(Sigma_hat0)) / 2
    Sigma_u_hat0 <- (Sigma_u_hat0 + t(Sigma_u_hat0)) / 2
    Exi2 <- estimate_Exi2(Y, mu_hat)$Exi2
    Sigma_hat <- Exi2*Sigma_hat0
    Lambda_hat <- Exi2*Lambda_hat0
    Sigma_u_hat <- Exi2*Sigma_u_hat0
  }
  
  return(list(
    Lambda_hat = Lambda_hat,
    Gamma_hat = Gamma_hat,
    Sigma_u_hat = Sigma_u_hat,
    Sigma_hat = Sigma_hat,
    Sigma_hat0 = Sigma_hat0
  ))
}


evaluate_covariance_estimators <- function(Y, Sigma_true, Sigma_u_true,
                                           method = c("SAMPLE","TylerSPCA","FLW","IPSN"),
                                           estimate_2_part = NULL) {
  method <- match.arg(method)
  p <- ncol(Y)
  m <- 3
  
  eigen_Sigma <- eigen(Sigma_true, symmetric = TRUE)
  Sigma_inv_sqrt <- eigen_Sigma$vectors %*% 
    diag(1 / sqrt(eigen_Sigma$values)) %*% 
    t(eigen_Sigma$vectors)
  
  if (is.null(estimate_2_part)) {
    estimate_2_part <- compute_2_part_nonsparse(Y, m, method = method)
  }
  
  Sigma_hat <- estimate_2_part$Sigma_hat
  
  # ||Sigma_hat - Sigma||_max
  Sigma_max_err <- max(abs(Sigma_hat - Sigma_true))
  
  # ||Sigma_hat - Sigma||_Sigma
  Sigma_rF_err <- norm(
    Sigma_inv_sqrt %*% (Sigma_hat - Sigma_true) %*% Sigma_inv_sqrt,
    type = "F"
  ) / sqrt(p)
  
  # ||Sigma_hat - Sigma||_2
  Sigma_spectral_err <- norm(Sigma_hat - Sigma_true, type = "2")
  
  return(list(
    Sigma_max_err = Sigma_max_err,
    Sigma_rF_err = Sigma_rF_err,
    Sigma_spectral_err = Sigma_spectral_err
  ))
}
run_covariance_simulation_onecase <- function(p = 100,
                                              sims = 10,
                                              mi = 1,
                                              n = 250,
                                              m = 3) {
  methods <- c("SAMPLE", "TylerSPCA", "FLW", "IPSN")
  simulation_results <- list()
  mean_results <- list()
  sd_results <- list()
  
  for (s in 1:sims) {
    cat("mi =", mi, "- Simulation", s, "...\n")
    data_generated <- generate_data_four_cases(p = p, n = n, m = m, mi = mi)
    Sigma_u <- data_generated$Sigma_u0
    Sigma_true <- data_generated$Sigma0
    Y <- data_generated$Y
    
    simulation_results[[s]] <- list()
    
    for (method in methods) {
      loss <- evaluate_covariance_estimators(Y, Sigma_true, Sigma_u,
                                             method = method)
      simulation_results[[s]][[method]] <- list(
        Sigma_max_err = loss$Sigma_max_err,
        Sigma_rF_err = loss$Sigma_rF_err,
        Sigma_spectral_err = loss$Sigma_spectral_err
      )
    }
  }
  for (method in methods) {
    Sigma_max_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$Sigma_max_err)
    Sigma_rF_errs  <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$Sigma_rF_err)
    Sigma_spectral_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$Sigma_spectral_err)
    mean_results[[method]] <- list(
      mean_Sigma_max_err = mean(Sigma_max_errs),
      mean_Sigma_rF_err  = mean(Sigma_rF_errs),
      mean_Sigma_spectral_err = mean(Sigma_spectral_errs)
    )
    
    sd_results[[method]] <- list(
      sd_Sigma_max_err = sd(Sigma_max_errs),
      sd_Sigma_rF_err  = sd(Sigma_rF_errs),
      sd_Sigma_spectral_err = sd(Sigma_spectral_errs)
    )
  }
  
  return(list(
    mean_results = mean_results,
    sd_results = sd_results
  ))
}



set.seed(988)
methods <- c("SAMPLE", "TylerSPCA", "FLW", "IPSN")
metrics <- c("Sigma_max_err", "Sigma_rF_err", "Sigma_spectral_err")
mi_vec <- 1:4
M <- length(methods)
K <- length(metrics)
latex_mats <- vector("list", length(mi_vec))
names(latex_mats) <- paste0("mi", mi_vec)
for (ell in seq_along(mi_vec)) {
  mi <- mi_vec[ell]
  cov_results <- run_covariance_simulation_onecase(p = 500, sims = 100, mi = mi)
  mean_mat <- matrix(NA_real_, nrow = M, ncol = K,
                     dimnames = list(methods, metrics))
  sd_mat   <- matrix(NA_real_, nrow = M, ncol = K,
                     dimnames = list(methods, metrics))
  for (meth in methods) {
    mean_vals <- cov_results$mean_results[[meth]]
    sd_vals   <- cov_results$sd_results[[meth]]
    for (met in metrics) {
      mean_mat[meth, met] <- mean_vals[[paste0("mean_", met)]]
      sd_mat[meth, met]   <- sd_vals[[paste0("sd_", met)]]
    }
  }
  res_mi <- matrix("", nrow = M, ncol = K,
                   dimnames = list(methods, metrics))
  for (i in 1:M) {
    for (j in 1:K) {
      res_mi[i, j] <- sprintf("%.4f(%.4f)",
                              mean_mat[i, j],
                              sd_mat[i, j])
    }
  }
  
  latex_mats[[ell]] <- res_mi
}
latex_mats$mi1
latex_mats$mi2
latex_mats$mi3
latex_mats$mi4