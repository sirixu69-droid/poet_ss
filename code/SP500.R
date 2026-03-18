

eigen_projection <- function(A, epsilon = 1e-6) {
  eig <- eigen(A, symmetric = TRUE)
  values <- pmax(eig$values, epsilon)
  return(eig$vectors %*% diag(values) %*% t(eig$vectors))
}

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
      mean(z)  # fallback
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
library(parallel)
library(MASS)
library("SPCAvRP")
library("mvtnorm")
#library("mnormt")
library(ICSNP)
#library(ICS)
library("SpatialNP")
#library("MNM")
#library(flare)
#library(glasso)
library(glassoFast)
library(stats)
#library(tidyverse)
library(lubridate)  
library(dplyr)
library(readr)
# ---- 读取两个表 ----
train_data_raw <- readr::read_csv("D:/R/somework/POET/data/sp500_excess_return_wide.csv") |>
  dplyr::mutate(
    DATE = lubridate::ymd(date),
    year = lubridate::year(DATE),
    month = lubridate::month(DATE)
  ) |>
  dplyr::filter(DATE >= as.Date("1995-01-01") & DATE <= as.Date("2023-12-31")) |>
  dplyr::select(DATE, year, month, dplyr::everything())

test_data_raw <- readr::read_csv("D:/R/somework/POET/data/sp500_anymember_excess_return_wide.csv") |>
  dplyr::mutate(
    DATE = lubridate::ymd(date),
    year = lubridate::year(DATE),
    month = lubridate::month(DATE)
  ) |>
  dplyr::filter(DATE >= as.Date("1995-01-01") & DATE <= as.Date("2023-12-31")) |>
  dplyr::select(DATE, year, month, dplyr::everything())

tune_factor_number <- function(X, kmax, method = c("SAMPLE", "FLW", "SPCA","TylerSPCA","IPSN","RegTME")) {
  method <- match.arg(method)
  n <- nrow(X)
  p <- ncol(X)
  m <- min(n, p) - 1  
  Sigma_full <- Sigma_hat_estimate(X, method = method)
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


Sigma_hat_estimate <- function(Y, method = c("SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME")) {
  method <- match.arg(method)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SAMPLE") {
    Sigma_hat <- cov(Y)
  }
  if (method == "FLW") {
    Sigma_hat <- RCov(Y)
  }
  if (method == "SPCA") {
    Sigma_hat <- SCov(Y)*p
  }
  if (method == "TylerSPCA") {
    Sigma_hat <- SCov(Y)*p
  }
  if (method == "IPSN") {
    Sigma_hat <- RCov(Y)
  }
  if (method == "RegTME") {
    Y_sym <- symmetrize_Y(Y)
    Sigma_hat <- regTME_svd(Y_sym, alpha = NULL)
  }
  return(Sigma_hat)
}

compute_2_part_nonsparse <- function(Y, m, method = c("SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME")) {
  method <- match.arg(method)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "SAMPLE") {
    Sigma_hat <- cov(Y)
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat <- Sigma_hat*c
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    if(length(eig$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    Sigma_u_hat <- Sigma_hat - low_rank_part
    tau <- 0.1*(sqrt(log(p) / n) +1/sqrt(p))
    Sigma_u_hat <- adaptive_threshold(Sigma_u_hat, tau = tau)
    Sigma_hat <- Sigma_u_hat + low_rank_part
    Sigma_hat_inv <- solve(Sigma_hat)
  }
  
  if (method == "FLW") {
    D_hat <- robust_variance_diag(Y)
    estimate <- estimate_Lambda(Y, D_hat, m)
    Lambda_hat <- estimate$Lambda_hat
    Sigma_hat <- estimate$Sigma1_hat
    c<-p/sum(diag(Sigma_hat))
    Sigma_hat <- Sigma_hat*c
    Lambda_hat <- Lambda_hat*c
    HK<-RCov(Y)
    eig <- eigen(HK, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m] 
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    Sigma_u_hat <- Sigma_hat - low_rank_part
    Sigma_u_hat1 <- eigen_projection(Sigma_u_hat,epsilon = 1e-5)
    max_diff <- max(abs(Sigma_u_hat - Sigma_u_hat1))
    mu <- sqrt(log(p) / n)
    if (max_diff < 2*mu) {
      Sigma_u_hat <- Sigma_u_hat1
    } else {
      Sigma_u_hat <- accelerated_proximal_gradient(Sigma_u_hat, mu, max_iter = 300, tol = 5e-3)
    }
    tau <- 0.1*(sqrt(log(p) / n) +1/sqrt(p))
    Sigma_u_hat <- adaptive_threshold(Sigma_u_hat, tau = tau)
    Sigma_hat <- low_rank_part + Sigma_u_hat
    Sigma_hat_inv <- solve(Sigma_hat)
    
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
    mu <- sqrt(log(p) / n)
    tau <- 0.1*(mu +sqrt(log(n)/n))
    Sigma_u_hat <- adaptive_threshold(R_tilde, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat
    Sigma_hat_inv <- solve(Sigma_hat)
  }
  if (method == "TylerSPCA") {
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
    mu <- sqrt(log(p) / n)
    tau <- 0.1*(mu +sqrt(log(n)/n))
    Sigma_u_hat <- adaptive_threshold(R_tilde, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat
    Sigma_hat_inv <- solve(Sigma_hat)
    mu_hat <- spatial.median(Y)
    Y_centered <- sweep(Y, 2, mu_hat, FUN = "-")
    quad_forms <- rowSums((Y_centered %*% Sigma_hat_inv) * Y_centered) 
    weights <- 1 / quad_forms  
    Sigma_T <- (p / n) * t(Y_centered) %*% diag(weights) %*% Y_centered
    c<-p/sum(diag(Sigma_T))
    Sigma_hat <- Sigma_T*c
    eig_Sigma_T <- eigen(Sigma_hat, symmetric = TRUE)
    Gamma_hat <- eig_Sigma_T$vectors[, 1:m, drop = FALSE]
    if(length(eig_Sigma_T$values[1:m]) == 1) {
      Lambda_hat <- matrix(eig_Sigma_T$values[1:m], nrow = 1, ncol = 1)
    } else {
      Lambda_hat <- diag(eig_Sigma_T$values[1:m])
    }
    low_rank_part <- Gamma_hat %*% Lambda_hat %*% t(Gamma_hat)
    Sigma_u_hat <- Sigma_hat - low_rank_part
    Sigma_u_hat <- adaptive_threshold(R_tilde, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat
    Sigma_hat_inv <- solve(Sigma_hat)
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
    Sigma_u_hat0 <- Sigma_0_hat - low_rank_part
    tau <- 0.1*(sqrt(log(p) / n) + sqrt(log(p) /p))
    Sigma_u_hat0 <- adaptive_threshold(Sigma_u_hat0, tau = tau)
    Sigma_hat <- low_rank_part+Sigma_u_hat0
    Sigma_hat_inv <- solve(Sigma_hat)
    
  }
  if (method == "RegTME") {
    Y_sym <- symmetrize_Y(Y)
    Sigma_tme <- regTME_svd(Y_sym, alpha = NULL)
    eig <- eigen(Sigma_tme, symmetric = TRUE)
    Gamma_hat <- eig$vectors[, 1:m, drop = FALSE]
    eigvals <- eig$values[1:m]
    if (m == 1) {
      Lambda_hat0 <- matrix(eigvals, nrow = 1, ncol = 1)
    } else {
      Lambda_hat0 <- diag(eigvals)
    }
    n_sym <- nrow(Y_sym)
    Sigma_u_hat0 <- Sigma_tme - Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)
    tau <- 0.1*(sqrt(log(p) / n_sym)+sqrt(log(p) /p))
    Sigma_u_hat0 <- adaptive_threshold(Sigma_u_hat0, tau = tau)
    Sigma_hat <- Gamma_hat %*% Lambda_hat0 %*% t(Gamma_hat)+Sigma_u_hat0
    Sigma_hat_inv <- solve(Sigma_hat)
  }
  return(list(
    Sigma_hat_inv = Sigma_hat_inv,
    Sigma_hat = Sigma_hat
  ))
}

compute_weights_nonparse <- function(Y, kmax, method = c("SAMPLE","SPCA", "FLW","TylerSPCA","IPSN","RegTME","EW")) {
  method <- match.arg(method)
  Y <- as.matrix(Y)
  p <- ncol(Y)
  n <- nrow(Y)
  if (method == "EW") {
    weights <- rep(1/p, p)  
    weights_matrix <- matrix(weights, ncol = 1)
    rownames(weights_matrix) <- colnames(Y)
    colnames(weights_matrix) <- "weight"
    return(weights_matrix)
  }
  parameter_tuned <- tune_factor_number(Y, kmax, method = method)
  factor_number <- parameter_tuned$best_m_mktcr
  # Step 2: estimate R_hat
  estimate_2_part <- compute_2_part_nonsparse(Y,factor_number, method = method)
  Omega_hat0 <- estimate_2_part$Sigma_hat_inv
  ones <- matrix(1, nrow = ncol(Omega_hat0), ncol = 1)
  numerator <- Omega_hat0 %*% ones
  denominator <- t(ones) %*% Omega_hat0 %*% ones
  denominator <- as.numeric(denominator)  
  if (denominator == 0) {
    stop("分母为0，无法计算权重")
  }
  mvp_weights <- numerator / denominator
  rownames(mvp_weights) <- colnames(Y)
  colnames(mvp_weights) <- "weight"
  return(mvp_weights)
}

# ---- 定义风险计算函数 ----
compute_yearly_risk <- function(train_data_raw, test_data_raw,
                                methods = c("EW","SAMPLE", "FLW","TylerSPCA", "SPCA","IPSN","RegTME")) {
  
  all_months <- train_data_raw %>%
    filter(year >= 2005 & year <= 2023) %>%
    distinct(year, month) %>%
    arrange(year, month) %>%
    mutate(ym = paste(year, month, sep = "-"))
  
  monthly_returns <- list()
  
  for (i in 1:nrow(all_months)) {
    current_year <- all_months$year[i]
    current_month <- all_months$month[i]
    current_ym <- all_months$ym[i]
    cat("处理月份：", current_ym, "\n")
    
    # 训练窗口（前10年）
    train_end_date <- as.Date(paste(current_year, current_month, "01", sep = "-")) - days(1)  # 当月第一天的前一天
    if (month(train_end_date) == 2 && day(train_end_date) == 29) {
      train_end_date1 <- ymd(paste(year(train_end_date), 2, 28, sep = "-"))  # 转为2月28日
    }else{
      train_end_date1 <- train_end_date
    }
    train_start_date <- train_end_date1 - years(10) + days(1)  # 前10年的同一天
    
    train_data <- train_data_raw %>%
      filter(DATE >= train_start_date & DATE <= train_end_date)
    stock_cols <- setdiff(colnames(train_data), c("date","DATE","year","month"))
    
    common_stocks <- stock_cols[
      sapply(train_data[, stock_cols], function(col) all(!is.na(col)))
    ]
    train_data <- train_data %>% dplyr::select(all_of(common_stocks))
    # 测试数据
    test_start_date <- as.Date(paste(current_year, current_month, "01", sep = "-"))
    test_end_date <- ceiling_date(test_start_date, "month") - days(1)
    test_data <- test_data_raw %>%
      filter(DATE >= test_start_date & DATE <= test_end_date)
    
    test_data <- test_data %>% dplyr::select(all_of(common_stocks))
    # ---- 计算权重并组合收益 ----
    kmax <- 10
    for (method in methods) {
      weights_matrix <- compute_weights_nonparse(Y = train_data, kmax, method = method)
      weights <- as.numeric(weights_matrix[,1])
      names(weights) <- rownames(weights_matrix)
      test_matrix <- as.matrix(sapply(test_data, as.numeric))
      weight_mat <- matrix(weights, nrow = nrow(test_matrix), ncol = length(weights), byrow = TRUE)
      avail_mask <- !is.na(test_matrix)
      test_returns <- rowSums(weight_mat * replace(test_matrix, is.na(test_matrix), 0))
      monthly_returns[[paste(current_ym, method, sep = "_")]] <- tibble(
        date = test_data_raw %>% filter(DATE >= test_start_date & DATE <= test_end_date) %>% pull(DATE),
        method = method,
        portfolio_return = as.vector(test_returns)
      )
    }
  }
  
  all_returns <- bind_rows(monthly_returns) %>%
    mutate(year = year(date)) %>%
    filter(year >= 2005 & year <= 2023)
  
  # ---- 年化标准差（按年，长表）----
  yearly_sd <- all_returns %>%
    group_by(year, method) %>%
    summarise(annualized_sd = sd(portfolio_return, na.rm = TRUE) * sqrt(n()),
              .groups = "drop")
  
  # ---- 年度 p-value 检验（长表）----
  yearly_pval <- all_returns %>%
    group_by(year) %>%
    group_split() %>%
    lapply(function(year_data) {
      y <- unique(year_data$year)
      r_sseca <- year_data %>% filter(method == "TylerSPCA") %>% pull(portfolio_return)
      data.frame(
        year = y,
        method = c("SAMPLE", "FLW","SPCA","IPSN","RegTME","EW"),
        pval = c(
          poet_ipsn_hac_test(r_sseca, year_data %>% filter(method=="SAMPLE") %>% pull(portfolio_return))$p_value,
          poet_ipsn_hac_test(r_sseca, year_data %>% filter(method=="FLW") %>% pull(portfolio_return))$p_value,
          poet_ipsn_hac_test(r_sseca, year_data %>% filter(method=="SPCA") %>% pull(portfolio_return))$p_value,
          poet_ipsn_hac_test(r_sseca, year_data %>% filter(method=="IPSN") %>% pull(portfolio_return))$p_value,
          poet_ipsn_hac_test(r_sseca, year_data %>% filter(method=="RegTME") %>% pull(portfolio_return))$p_value,
          poet_ipsn_hac_test(r_sseca, year_data %>% filter(method=="EW") %>% pull(portfolio_return))$p_value
        )
      )
    }) %>% bind_rows()
  
  n_days <- all_returns %>%
    distinct(date) %>%
    mutate(year = year(date)) %>%
    group_by(year) %>%
    summarise(n = n(), .groups = "drop") %>%
    summarise(mean_n = mean(n)) %>%
    pull() %>%
    round()
  
  overall_sd <- all_returns %>%
    group_by(method) %>%
    summarise(annualized_sd = sd(portfolio_return, na.rm = TRUE) * sqrt(n_days),
              .groups = "drop")
  r_sseca_all <- all_returns %>% filter(method == "TylerSPCA") %>% pull(portfolio_return)
  overall_pval <- data.frame(
    method = c("SAMPLE", "FLW","SPCA","IPSN","RegTME","EW"),
    pval = c(
      poet_ipsn_hac_test(r_sseca_all, all_returns %>% filter(method=="SAMPLE") %>% pull(portfolio_return))$p_value,
      poet_ipsn_hac_test(r_sseca_all, all_returns %>% filter(method=="FLW") %>% pull(portfolio_return))$p_value,
      poet_ipsn_hac_test(r_sseca_all, all_returns %>% filter(method=="SPCA")%>% pull(portfolio_return))$p_value,
      poet_ipsn_hac_test(r_sseca_all, all_returns %>% filter(method=="IPSN") %>% pull(portfolio_return))$p_value,
      poet_ipsn_hac_test(r_sseca_all, all_returns %>% filter(method=="RegTME") %>% pull(portfolio_return))$p_value,
      poet_ipsn_hac_test(r_sseca_all, all_returns %>% filter(method=="EW") %>% pull(portfolio_return))$p_value
    )
  )
  
  # ---- 存到文件 ----
  write_csv(yearly_sd, "D:/R/somework/POET/data/yearly_sd.csv")
  write_csv(yearly_pval, "D:/R/somework/POET/data/yearly_pval.csv")
  write_csv(overall_sd, "D:/R/somework/POET/data/overall_sd.csv")
  write_csv(overall_pval, "D:/R/somework/POET/data/overall_pval.csv")
  
  list(
    yearly_sd = yearly_sd,
    yearly_pval = yearly_pval,
    overall_sd = overall_sd,
    overall_pval = overall_pval,
    all_returns = all_returns
  )
}

result <- compute_yearly_risk(train_data_raw, test_data_raw)

print(result$yearly_sd)
print(result$yearly_pval)
print(result$overall_sd)
print(result$overall_pval)
