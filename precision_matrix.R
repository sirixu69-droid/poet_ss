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


generate_precision_factor_data_joint <- function(B, 
                                                 n = 250, 
                                                 mu = rep(0, nrow(B)),
                                                 mi = 1) {
  p <- nrow(B)
  m <- ncol(B)
  
  # idiosyncratic precision matrix Vu = (0.4^|i-j|)
  indices <- matrix(1:p, nrow = p, ncol = p)
  V_u_raw <- 0.4^abs(indices - t(indices))
  Sigma_u_raw <- solve(V_u_raw)
  
  # covariance of y before normalization
  Sigma_raw <- B %*% t(B) + Sigma_u_raw
  c0 <- p / sum(diag(Sigma_raw))
  
  Sigma0 <- c0 * Sigma_raw
  Sigma_u0 <- c0 * Sigma_u_raw
  
  # covariance of (f, u)
  Sigma_fu <- matrix(0, nrow = m + p, ncol = m + p)
  Sigma_fu[1:m, 1:m] <- c0 * diag(m)
  Sigma_fu[(m + 1):(m + p), (m + 1):(m + p)] <- c0 * Sigma_u_raw
  
  if (mi == 1) {
    # Scenario I: N(0, Sigma_fu)
    joint_sample <- mvtnorm::rmvnorm(
      n, mean = rep(0, m + p), sigma = Sigma_fu
    )
    
  } else if (mi == 2) {
    # Scenario II: t_4
    df <- 4
    joint_sample <- mvtnorm::rmvt(
      n, sigma = Sigma_fu * (df - 2) / df, df = df
    )
    
  } else if (mi == 3) {
    # Scenario III: t_2.2
    df <- 2.2
    joint_sample <- mvtnorm::rmvt(
      n, sigma = Sigma_fu * (df - 2) / df, df = df
    )
    
  } else if (mi == 4) {
    # Scenario IV: 0.8N(0,Sigma_fu)+0.2N(0,10Sigma_fu)
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
    
    # standardize covariance back to Sigma_fu
    joint_sample <- joint_sample / sqrt(2.8)
    
  } else {
    stop("mi must be 1, 2, 3, or 4.")
  }
  
  f_t <- joint_sample[, 1:m, drop = FALSE]
  u_t <- joint_sample[, (m + 1):(m + p), drop = FALSE]
  Y <- f_t %*% t(B) + u_t
  
  if (length(mu) == p) {
    Y <- Y + matrix(mu, nrow = n, ncol = p, byrow = TRUE)
  }
  
  list(
    data = Y,
    Sigma0 = Sigma0,
    Sigma_u0 = Sigma_u0,
    V0 = solve(Sigma0),
    V_u0 = solve(Sigma_u0)
  )
}
poet_ipsn_hac_test <- function(ri, rn) {
  T_len <- length(ri)
  stopifnot(length(rn) == T_len)
  mu_i <- mean(ri); mu_n <- mean(rn)
  gamma_i <- mean(ri^2); gamma_n <- mean(rn^2)
  Delta_hat <- log(gamma_i - mu_i^2) - log(gamma_n - mu_n^2)
  y_t <- cbind(
    ri - mu_i,
    rn - mu_n,
    ri^2 - gamma_i,
    rn^2 - gamma_n
  )
  
  rho_vec <- apply(y_t, 2, function(yt) {
    yt_lag <- yt[-1]; yt_now <- yt[-length(yt)]
    sum(yt_now * yt_lag) / sum(yt_now^2)
  })
  sigma_sq_vec <- apply(y_t, 2, var)
  
  numerator <- sum((4 * rho_vec^2 * sigma_sq_vec^2) / (1 - rho_vec)^8)
  denominator <- sum(sigma_sq_vec^2 / (1 - rho_vec)^4)
  alpha2_hat <- numerator / denominator
  bw <- 1.3221 * alpha2_hat^(1/5) * T_len^(1/5)
  
  kern_QS <- function(x) {
    if (abs(x) < 1e-8) return(1)
    pi_x <- pi * x
    (25 / (12 * pi_x^2)) * (sin(6 * pi_x / 5) / (6 * pi_x / 5) - cos(6 * pi_x / 5))
  }
  Psi_hat <- matrix(0, 4, 4)
  for (j in -floor(T_len - 1):(T_len - 1)) {
    weight <- kern_QS(j / bw)
    if (abs(weight) < 1e-10) next
    idx_t <- (1 + abs(j)):T_len
    idx_tj <- idx_t - abs(j)
    y_t_j <- y_t[idx_t, , drop = FALSE]
    y_tj_j <- y_t[idx_tj, , drop = FALSE]
    cov_j <- crossprod(y_t_j, y_tj_j) / T_len
    Psi_hat <- Psi_hat + weight * (cov_j + t(cov_j)) / 2
  }
  Psi_hat <- T_len*Psi_hat/(T_len-4)
  
  denom_i <- gamma_i - mu_i^2
  denom_n <- gamma_n - mu_n^2
  dfdv <- c(-2 * mu_i / denom_i, 2 * mu_n / denom_n, 1 / denom_i, -1 / denom_n)
  
  var_Delta_hat <- as.numeric(t(dfdv) %*% Psi_hat %*% dfdv) / T_len
  sd_Delta_hat <- sqrt(var_Delta_hat)
  test_statistic <- Delta_hat / sd_Delta_hat
  p_value <- pnorm(test_statistic, lower.tail = TRUE)
  
  list(
    Delta_hat = Delta_hat,
    test_statistic = test_statistic,
    p_value = p_value,
    sd_Delta_hat = sd_Delta_hat,
    bw = bw,
    rho = rho_vec,
    sigma2 = sigma_sq_vec,
    alpha2_hat = alpha2_hat
  )
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
    scale_full <- mad(y, constant = 1) / 0.6745
    fold_id <- sample(rep(1:nfold, length.out = n))
    cv_loss <- numeric(length(c_candidates))
    
    for (j in seq_along(c_candidates)) {
      c_val <- c_candidates[j]
      H <- c_val * scale_full * sqrt(n / logN) 
      
      fold_errors <- numeric(nfold)
      for (k in 1:nfold) {
        y_train <- y[fold_id != k]
        y_test  <- y[fold_id == k]
        scale_train <- mad(y_train, constant = 1) / 0.6745
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

estimate_Exi2 <- function(Y, mu_hat, epsilon = 0.345, c_candidates = c(0.5, 1.0, 1.5), nfold = 3) {
  n <- nrow(Y)
  p <- ncol(Y)
  xi_sq <- apply(Y, 1, function(y_t) sum((y_t - mu_hat)^2) / p)
  scale_full <- mad(xi_sq, constant = 1) / 0.6745
  fold_id <- sample(rep(1:nfold, length.out = n))
  cv_loss <- numeric(length(c_candidates))
  
  for (j in seq_along(c_candidates)) {
    c_val <- c_candidates[j]
    H <- c_val * n^min(1 / (1 + epsilon / 2), 1/2)
    
    fold_errors <- numeric(nfold)
    for (k in 1:nfold) {
      train_xi <- xi_sq[fold_id != k]
      test_xi  <- xi_sq[fold_id == k]
      scale_train <- mad(train_xi, constant = 1) / 0.6745
      huber_fit <- MASS::huber(train_xi, k = H / scale_train)
      mu_train <- huber_fit$mu
      
      fold_errors[k] <- mean((test_xi - mu_train)^2)
    }
    cv_loss[j] <- mean(fold_errors)
  }
  
  j_best <- which.min(cv_loss)
  c_opt <- c_candidates[j_best]
  H_best <- c_opt * n^min(1 / (1 + epsilon / 2), 1/2)
  huber_fit_final <- MASS::huber(xi_sq, k = H_best / scale_full)
  
  return(list(Exi2 = huber_fit_final$mu, c_opt = c_opt))
}
symmetrize_minabs <- function(mat) {
  mat_sym <- mat
  idx <- upper.tri(mat)
  a <- abs(mat[idx])
  b <- abs(t(mat)[idx])
  use_upper <- a <= b
  mat_sym[idx] <- mat[idx] * use_upper + t(mat)[idx] * (!use_upper)
  mat_sym <- mat_sym + t(mat_sym) - diag(diag(mat_sym))  # 保证对称
  return(mat_sym)
}
symmetrize_Y <- function(Y) {
  n <- nrow(Y)
  m <- floor(n / 2)
  Y_sym <- Y[seq(1, 2 * m, 2), ] - Y[seq(2, 2 * m, 2), ]
  return(Y_sym)
}
TME_svd <- function(Y, tol = 1e-4, maxit = 100) {
  n <- nrow(Y)
  p <- ncol(Y)
  
  if (n < p) {
    svd_Y <- svd(Y)
    U <- svd_Y$u
    D <- diag(svd_Y$d)
    V <- svd_Y$v
    Z <- U %*% D      # n × n
    proj_mat <- V     # p × n
    dim_Z <- n
  } else {
    # 低维情况，直接用原数据
    Z <- Y
    proj_mat <- diag(p)
    dim_Z <- p
  }
  
  # Tyler M-estimator
  Sigma <- diag(dim_Z)
  for (iter in 1:maxit) {
    Sigma_inv <- solve(Sigma)
    updated <- matrix(0, dim_Z, dim_Z)
    for (i in 1:n) {
      xi <- matrix(Z[i, ], ncol = 1)
      denom <- as.numeric(t(xi) %*% Sigma_inv %*% xi)
      updated <- updated + (xi %*% t(xi)) / denom
    }
    Sigma_new <- dim_Z * updated / sum(diag(updated))
    if (norm(Sigma_new - Sigma, type = "F") < tol) break
    Sigma <- Sigma_new
  }
  
  # 投影回 p 维
  Sigma_full <- proj_mat %*% Sigma %*% t(proj_mat)
  return(Sigma_full)
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


tune_factor_number <- function(X, kmax, method = c("SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME"), B = 50) {
  method <- match.arg(method)
  n <- nrow(X)
  p <- ncol(X)
  m <- min(n, p)  
  Y_sym <- symmetrize_Y(X)
  Sigma_full <- switch(method,
                       SAMPLE = cov(X),
                       FLW = RCov(X),
                       TylerSPCA = SCov(X),
                       SPCA = SCov(X),
                       IPSN = RCov(X),
                       RegTME = TME_svd(Y_sym))
  eig_vals <- eigen(Sigma_full, symmetric = TRUE, only.values = TRUE)$values[1:m]
  V <- rev(cumsum(rev(eig_vals)))
  ratio_seq <- log(1 + eig_vals[1:kmax] / V[1:kmax]) /
    log(1 + eig_vals[2:(kmax + 1)] / V[2:(kmax + 1)])
  
  m_hat <- which.max(ratio_seq)
  cat("固定 m̂ =", m_hat, "\n")
  
  return(list(
    best_m_mktcr = m_hat,
    method = method,
    ratio_seq = ratio_seq
  ))
}
compute_2_part_nonsparse <- function(Y, m, method = c("SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME")) {
  method <- match.arg(method)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SAMPLE") {
    Sigma_hat <- cov(Y)
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    Sigma_u_hat <- eigen_projection(R_tilde,epsilon = 1e-3)
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <-Lambda_hat*c 
    #cv_grid <- exp(seq(log(0.01), log(0.5), length.out = 10))
    mu <- sqrt(log(p) / n) +1/sqrt(p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  
  if (method == "FLW") {
    D_hat <- robust_variance_diag(Y)
    estimate <- estimate_Lambda(Y, D_hat, m)
    Lambda_hat <- estimate$Lambda_hat
    Sigma_hat <- estimate$Sigma1_hat
    HK<-RCov(Y)
    eig <- eigen(HK, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    R_tilde1 <- eigen_projection(R_tilde,epsilon = 1e-5)
    max_diff <- max(abs(R_tilde - R_tilde1))
    mu <- sqrt(log(p) / n)
    if (max_diff < 2*mu) {
      R_tilde2 <- R_tilde1
    } else {
      R_tilde2 <- accelerated_proximal_gradient(R_tilde, mu, max_iter = 300, tol = 5e-1)
    }
    Sigma_u_hat <- R_tilde2
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <- Lambda_hat*c
    #cv_grid <- exp(seq(log(0.01), log(0.5), length.out = 10))
    mu <- sqrt(log(p) / n) +1/sqrt(p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  if (method == "SPCA") {
    Sigma_hat <- SCov(Y)*p
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_tilde <- Sigma_hat - low_rank_part
    Sigma_u_hat <- eigen_projection(R_tilde,epsilon = 1e-3)
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <-Lambda_hat*c
    #cv_grid <- exp(seq(log(0.01), log(0.5), length.out = 10))
    mu <- sqrt(log(p) / n) + sqrt(log(n) /n)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  if (method == "TylerSPCA") {
    estimate  <- compute_2_part_nonsparse(Y, m, method = "SPCA")
    Sigma_u_hat0 <- estimate$Sigma_u_hat0
    Gamma_hat    <- estimate$Gamma_hat
    Lambda_hat0  <- estimate$Lambda_hat0
    mu <- sqrt(log(p) / n) + sqrt(log(n) /n)
    #cv_grid <- exp(seq(log(0.15), log(0.5), length.out = 10))
    #cv_result0 <- cv_glasso(Y, cv_grid, K = 2, m = m, method = "SPCA")
    lambda_opt0 <- 0.05 * mu
    #lambda_opt0 <- cv_result0$lambda_opt
    Omega_u_hat0 <- glassoFast(Sigma_u_hat0, lambda_opt0)$wi
    Omega_u_hat0 <- (Omega_u_hat0 + t(Omega_u_hat0)) / 2 
    tmp0 <- solve(solve(Lambda_hat0) + t(Gamma_hat) %*% Omega_u_hat0 %*% Gamma_hat)
    Omega_hat0 <- Omega_u_hat0 - Omega_u_hat0 %*% Gamma_hat %*% tmp0 %*% t(Gamma_hat) %*% Omega_u_hat0
    mu_hat <- spatial.median(Y)
    Y_centered <- sweep(Y, 2, mu_hat, FUN = "-")
    Sigma_S_inv <- Omega_hat0
    quad_forms <- rowSums((Y_centered %*% Sigma_S_inv) * Y_centered) 
    weights <- 1 / quad_forms  
    Sigma_T <- (p / n) * t(Y_centered) %*% diag(weights) %*% Y_centered
    Sigma_hat <- (Sigma_T + t(Sigma_T)) / 2  
    eig_Sigma_T <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig_Sigma_T$vectors[, 1:m, drop = FALSE]
    if(length(eig_Sigma_T$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig_Sigma_T$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig_Sigma_T$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    R_hat_T <- Sigma_hat - low_rank_part
    Sigma_u_hat <- eigen_projection(R_hat_T,epsilon = 1e-3)
    Sigma_u_hat  <- (Sigma_u_hat  + t(Sigma_u_hat )) / 2
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat0 <- Sigma_hat*c
    Sigma_u_hat0 <- Sigma_u_hat*c
    Lambda_hat0 <-Lambda_hat*c
    #cv_grid <- exp(seq(log(0.01), log(0.5), length.out = 10))
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
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
    R_tilde <- Sigma_0_hat - low_rank_part
    Sigma_u_hat0 <- eigen_projection(R_tilde, epsilon = 1e-3)
    Sigma_u_hat0 <- (Sigma_u_hat0 + t(Sigma_u_hat0)) / 2
    Exi2 <- estimate_Exi2(Y, mu_hat)$Exi2
    Lambda_hat <- Exi2*Lambda_hat0
    Sigma_u_hat <- Exi2*Sigma_u_hat0
    #cv_grid <- exp(seq(log(0.01), log(0.5), length.out = 10))
    mu <- sqrt(log(p) / n) + sqrt(log(p) /p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  if (method == "RegTME") {
    Y_sym <- symmetrize_Y(Y)
    Sigma_tme <- regTME_svd(Y_sym, alpha = NULL)
    eig <- eigen(Sigma_tme, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE]
    eigvals <- eig$values[1:m]
    eigvals <- pmax(eigvals, 0)
    if (m == 1) {
      Lambda_hat0 <- matrix(eigvals, nrow = 1, ncol = 1)
    } else {
      Lambda_hat0 <- diag(eigvals)
    }
    n_sym <- nrow(Y_sym)
    Sigma_u_hat0 <- Sigma_tme - Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)
    Sigma_u_hat0 <- eigen_projection(Sigma_u_hat0,epsilon = 1e-3)
    Sigma_u_hat0 <- (Sigma_u_hat0 + t(Sigma_u_hat0)) / 2
    #cv_grid <- exp(seq(log(0.01), log(0.5), length.out = 10))
    mu <- sqrt(log(p) / n_sym)+1/sqrt(p)
    cv_grid <- mu * seq(0.05, 0.1, length.out = 5)
  }
  return(list(
    cv_grid = cv_grid,
    Lambda_hat0 = Lambda_hat0,
    Gamma_hat = Gamma_hat,
    Sigma_u_hat0 = Sigma_u_hat0
  ))
}
cv_glasso <- function(Y, lambda_grid, K = 5, m = 3, method = c("SAMPLE", "SPCA", "TylerSPCA", "FLW","IPSN","RegTME")) {
  n <- nrow(Y)
  folds <- sample(rep(1:K, length.out = n))
  n_lambda <- length(lambda_grid)
  cv_scores <- matrix(NA_real_, nrow = K, ncol = n_lambda)
  
  for (k in 1:K) {
    Y_train <- Y[folds != k, , drop = FALSE]
    Y_valid <- Y[folds == k, , drop = FALSE]
    est_train <- compute_2_part_nonsparse(Y_train, m, method = method)
    est_valid <- compute_2_part_nonsparse(Y_valid, m, method = method)
    Sigma_train_u <- est_train$Sigma_u_hat0
    Sigma_valid_u <- est_valid$Sigma_u_hat0
    for (j in seq_along(lambda_grid)) {
      ll <- lambda_grid[j]
      wi <- glassoFast(Sigma_train_u, ll)$wi
      cv_scores[k, j] <- log(det(wi)) - sum(diag(wi %*% Sigma_valid_u))
    }
  }
  mean_scores <- colMeans(cv_scores)
  best_idx <- which.max(mean_scores)
  return(list(lambda_opt = lambda_grid[best_idx],
              scores = mean_scores))
}

evaluate_precision_estimators <- function(Y, Sigma_true, Sigma_u_true,
                                          method = c("SAMPLE", "SPCA", "TylerSPCA", "FLW","IPSN","RegTME"),
                                          estimate_2_part = NULL) {
  method <- match.arg(method)
  p <- ncol(Y)
  m <- 3
  
  V_true <- solve(Sigma_true)
  V_u_true <- solve(Sigma_u_true)
  
  if (is.null(estimate_2_part)) {
    estimate_2_part <- compute_2_part_nonsparse(Y, m, method = method)
  }
  
  Sigma_u_hat0 <- estimate_2_part$Sigma_u_hat0
  cv_grid <- estimate_2_part$cv_grid
  cv_result0 <- cv_glasso(Y, cv_grid, K = 2, m = m, method = method)
  lambda_opt0 <- cv_result0$lambda_opt
  
  V_u_hat <- glassoFast(Sigma_u_hat0, lambda_opt0)$wi
  V_u_hat <- (V_u_hat + t(V_u_hat)) / 2 
  
  Gamma_hat <- estimate_2_part$Gamma_hat
  Lambda_hat0 <- estimate_2_part$Lambda_hat0
  
  tmp0 <- solve(solve(Lambda_hat0) + t(Gamma_hat) %*% V_u_hat %*% Gamma_hat)
  V_hat <- V_u_hat - V_u_hat %*% Gamma_hat %*% tmp0 %*% t(Gamma_hat) %*% V_u_hat
  
  # Frobenius norm errors
  V_F_err <- norm(V_hat - V_true, type = "F") / sqrt(p)
  V_u_F_err <- norm(V_u_hat - V_u_true, type = "F") / sqrt(p)
  
  # max norm errors
  V_max_err <- max(abs(V_hat - V_true))
  V_u_max_err <- max(abs(V_u_hat - V_u_true))
  
  # spectral norm errors
  V_spectral_err <- norm(V_hat - V_true, type = "2")
  V_u_spectral_err <- norm(V_u_hat - V_u_true, type = "2")
  
  return(list(
    V_F_err = V_F_err,
    V_u_F_err = V_u_F_err,
    V_max_err = V_max_err,
    V_u_max_err = V_u_max_err,
    V_spectral_err = V_spectral_err,
    V_u_spectral_err = V_u_spectral_err
  ))
}

run_precision_simulation_onecase <- function(p = 100,
                                             sims = 10,
                                             mi = 1) {
  methods <- c("SAMPLE", "SPCA", "TylerSPCA", "FLW","IPSN","RegTME")
  simulation_results <- list()
  mean_results <- list()
  sd_results <- list()
  
  for (s in 1:sims) {
    cat("Simulation", s, "...\n")
    
    sm <- c(1, 0.75^2, 0.5^2)
    B <- matrix(0, nrow = p, ncol = 3)
    for (k in 1:3) {
      B[, k] <- rnorm(p, mean = 0, sd = sqrt(sm[k]))
    }
    
    data_generated <- generate_precision_factor_data_joint(
      B, n = 250, mu = rep(0, p), mi = mi
    )
    
    Sigma_u <- data_generated$Sigma_u0
    Sigma_true <- data_generated$Sigma0
    Y <- data_generated$data
    simulation_results[[s]] <- list()
    
    for (method in methods) {
      estimate_2_part <- compute_2_part_nonsparse(Y, 3, method = method)
      loss <- evaluate_precision_estimators(
        Y, Sigma_true, Sigma_u,
        method = method,
        estimate_2_part = estimate_2_part
      )
      
      simulation_results[[s]][[method]] <- list(
        V_F_err = loss$V_F_err,
        V_u_F_err = loss$V_u_F_err,
        V_max_err = loss$V_max_err,
        V_u_max_err = loss$V_u_max_err,
        V_spectral_err = loss$V_spectral_err,
        V_u_spectral_err = loss$V_u_spectral_err
      )
    }
  }
  
  for (method in methods) {
    V_F_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$V_F_err)
    V_u_F_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$V_u_F_err)
    V_max_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$V_max_err)
    V_u_max_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$V_u_max_err)
    V_spectral_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$V_spectral_err)
    V_u_spectral_errs <- sapply(1:sims, function(s) simulation_results[[s]][[method]]$V_u_spectral_err)
    
    mean_results[[method]] <- list(
      mean_V_F_err = mean(V_F_errs),
      mean_V_u_F_err = mean(V_u_F_errs),
      mean_V_max_err = mean(V_max_errs),
      mean_V_u_max_err = mean(V_u_max_errs),
      mean_V_spectral_err = mean(V_spectral_errs),
      mean_V_u_spectral_err = mean(V_u_spectral_errs)
    )
    
    sd_results[[method]] <- list(
      sd_V_F_err = sd(V_F_errs),
      sd_V_u_F_err = sd(V_u_F_errs),
      sd_V_max_err = sd(V_max_errs),
      sd_V_u_max_err = sd(V_u_max_errs),
      sd_V_spectral_err = sd(V_spectral_errs),
      sd_V_u_spectral_err = sd(V_u_spectral_errs)
    )
  }
  
  return(list(
    simulation_results = simulation_results,
    mean_results = mean_results,
    sd_results = sd_results
  ))
}

set.seed(998)
precision_results <- run_precision_simulation_onecase(p = 500, sims = 100, mi = 3)
methods <- c("SAMPLE", "SPCA", "TylerSPCA", "FLW", "IPSN","RegTME")
for (method in methods) {
  cat("=== 方法:", method, "===\n")
  mean_vals <- precision_results$mean_results[[method]]
  sd_vals <- precision_results$sd_results[[method]]
  metrics <- names(mean_vals)
  for (metric in metrics) {
    metric_name <- sub("mean_", "", metric)
    sd_metric <- sub("mean_", "sd_", metric)
    
    cat(sprintf(
      "%-25s: %.4f ± %.4f\n",
      metric_name,
      mean_vals[[metric]],
      sd_vals[[sd_metric]]
    ))
  }
  cat("\n")
}

