#' run main microbiome bicluster algorithm.
#'
#'@param W_count The microbiome count matrix that you want to analyze.
#'@param G The number of component
#'@param Q_g The number of biclusters for each component, a vector.
#'@param pi_g A vector of initial guesses of component proportion
#'@param mu_g A list of initial guess of mean vector
#'@param sig_g A list of initial guess of covariance matrix for each component
#'@param m A list of initial guess of variational mean
#'@param V A list of initial guess of variational varaince
#'@param B_g A list of initial guess of bicluster membership
#'@param T_g A list of initial guess of covariance of latent variable: u
#'@param D_g A list of initial guess of error matrix
#'@param cov_str The covaraince structure you choose, there are 16 different models belongs to
#'this family:UUU, UUG, UUD, UUC, UGU, UGG, UGD, UGC, GUU, GUG, GUD, GUC, GGU, GGG, GGD, GGC.
#'@param iter Max iterations, default is 150.
#'@param const the permutation constant in multinomial distribution. Calculated before the main
#'algorithm in order to save computation time.
#'@param beta_g initial guess of covariates coefficients.
#'@param X The regression covariates matrix, which generates by model.matrix.
#'
#'
#'@return z_ig Estimated latent variable z
#'@return cluster Component labels
#'@return mu_g Estimated component mean
#'@return pi_g Estimated component proportion
#'@return B_g Estimated bicluster membership
#'@return T_g Estimated covariance of latent variable u
#'@return D_g Estimated error covariance
#'@return COV Estimated sparsity component covariance
#'@return beta_g Estimated covariates coefficients.
#'@return sigma Estimated original component covariance
#'@return overall_loglik Complete log likelihood value for each iteration
#'@return ICL ICL value
#'@return BIC BIC value
#'@return AIC AIC value



#main algorithm
Mico_bi_jensens <-
  function(W_count,
           G,
           Q_g,
           pi_g,
           mu_g,
           sig_g,
           V,
           m,
           B_g,
           T_g,
           D_g,
           cov_str,
           iter,
           const,
           beta_g,
           X) {

    W_count <- as.matrix(W_count)
    n <- dim(W_count)[1]
    K <- dim(W_count)[2] - 1


    code_split <- stringr::str_split(cov_str, "")[[1]]
    code_B <- code_split[1]
    code_T <- code_split[2]
    code_D <- code_split[3]




    overall_loglik <- rep(-Inf, 2)

    h <- 2

    repeat {
      COV <- list()
      for (g in 1:G) {
        COV[[g]] <- B_g[[g]] %*% T_g[[g]] %*% t(B_g[[g]]) + D_g[[g]]
      }

      #update m and V
      for (g in 1:G) {
        for (i in 1:n) {
          mv <- GD_mvxi(V[[g]][[i]], m[[g]][,i], COV[[g]], W_count[i,], mu_g[[g]][,i], K,
                        B_g[[g]], T_g[[g]], D_g[[g]])
          m[[g]][, i] <- c(mv$m)
          V[[g]][[i]] <- mv$V
        }
      }

      #print(c("mv",p_comp_logw(W_count,xi,alpha,m,V,K,beta_g,B_g,COV,mu_g,pi_g,G,n)))




      log_w_approx <- matrix(0, nrow = n, ncol = G)
      for (i in 1:n) {
      for (g in 1:G) {
          log_w_approx[i, g] <-
            lb_loglik(V[[g]][[i]], K, COV[[g]], m[[g]][, i], mu_g[[g]][,i], W_count[i, ], const[i])
        }
      }
      #print(elbo_fun(pi_g,G,n,log_w_approx))


      #update z_ig
      z_ig <- z_ig_update(pi_g, n, G, log_w_approx)


      #update pi_g
      pi_g <- colMeans(z_ig)
      label <- apply(z_ig, 1, which.max)
      if (length(pi_g[pi_g == 0]) > 0 |
          any(table(label) < 5) | length(unique(label)) < G) {
        break
      }


      #update mu_g, beta_g
      for (g in 1:G) {
        beta_g[[g]]<-ginv(t(X)%*%diag(z_ig[,g])%*%X)%*%t(X)%*%diag(z_ig[,g])%*%t(m[[g]])
        mu_g[[g]]<-t(X%*%beta_g[[g]])
      }
      p_mu=G*K*dim(X)[2]






      #second cycle
      #U_iK
      U_iK <- list()
      for (g in 1:G) {
        U_iK[[g]] <-
          solvecov(t(B_g[[g]]), diag(diag(D_g[[g]]) ^ (-1)), diag(diag(T_g[[g]]) ^
                                                                    (-1), Q_g[g], Q_g[g])) %*%
          t(B_g[[g]]) %*%
          diag(diag(D_g[[g]]) ^ (-1)) %*% (m[[g]] - mu_g[[g]])

      }

      #S_g
      z_ig <- as.matrix(z_ig, K, n)
      S_g <- list()
      for (g in 1:G) {
        S_g[[g]] <-
          (m[[g]] - mu_g[[g]]) %*% diag(z_ig[, g] / sum(z_ig[, g])) %*%
          t(m[[g]] - mu_g[[g]])
      }

      #sig_g
      sig_g <- list()
      for (g in 1:G) {
        beg <- 0
        for (i in 1:n) {
          beg <- beg + z_ig[i, g] * (V[[g]][[i]] + (m[[g]][, i] - mu_g[[g]][,i]) %*%
                                       t(m[[g]][, i] - mu_g[[g]][,i])) / sum(z_ig[, g])
        }
        sig_g[[g]] <- beg
      }

      #e2
      e2 <- list()
      for (g in 1:G) {
        e2[[g]] <-
          solvecov(t(B_g[[g]]), diag(diag(D_g[[g]]) ^ (-1)), diag(diag(T_g[[g]]) ^
                                                                    (-1), Q_g[g], Q_g[g])) +
          solvecov(t(B_g[[g]]), diag(diag(D_g[[g]]) ^ (-1)), diag(diag(T_g[[g]]) ^
                                                                    (-1), Q_g[g], Q_g[g])) %*%
          t(B_g[[g]]) %*%
          diag(diag(D_g[[g]]) ^ (-1)) %*% S_g[[g]] %*% t(diag(diag(D_g[[g]]) ^
                                                                (-1))) %*% B_g[[g]] %*%
          t(solvecov(t(B_g[[g]]), diag(diag(D_g[[g]]) ^ (-1)), diag(diag(T_g[[g]]) ^
                                                                      (-1), Q_g[g], Q_g[g])))
      }








      #T_g
      #E(uu^T|y)
      pre_T_g <- T_g
      for (g in 1:G) {
        T_g[[g]] <- diag(diag(e2[[g]]), Q_g[g], Q_g[g])

      }

      p_T <- sum(Q_g)
      int_T_g <- T_g

      if (code_T == "G" & length(unique(Q_g)) == 1) {
        p_T <- Q_g[1]
        beg_T <- 0
        for (g in 1:G) {
          beg_T <- int_T_g[[g]] * pi_g[g] + beg_T
        }
        for (g in 1:G) {
          T_g[[g]] <- beg_T
        }
      }
      else if (code_T == "D") {
        p_T <- G
        for (g in 1:G) {
          T_g[[g]] <- diag(mean(diag(int_T_g[[g]])), Q_g[g])
        }
      }
      else if (code_T == "C" & length(unique(Q_g)) == 1) {
        p_T <- 1
        beg_T <- 0
        for (g in 1:G) {
          beg_T <- int_T_g[[g]] * pi_g[g] + beg_T
        }
        for (g in 1:G) {
          T_g[[g]] <- diag(mean(diag(beg_T)), Q_g[g])
        }
      }




      #D_g
      pre_D_g <- D_g
      for (g in 1:G) {
        D_g[[g]] <- abs(diag(
          diag(
            sig_g[[g]] - 2 * B_g[[g]] %*%
              solvecov(t(B_g[[g]]), diag(diag(D_g[[g]]) ^
                                           (-1)),
                       diag(diag(T_g[[g]]) ^ (-1), Q_g[g], Q_g[g])) %*%
              t(B_g[[g]]) %*%
              diag(diag(D_g[[g]]) ^ (-1)) %*% S_g[[g]] +
              B_g[[g]] %*% e2[[g]] %*% t(B_g[[g]])
          )
        ))


      }




      p_D <- G * K


      int_D_g <- D_g
      if (code_D == "G") {
        p_D <- K
        beg_D <- 0
        for (g in 1:G) {
          beg_D <- int_D_g[[g]] * pi_g[g] + beg_D
        }
        for (g in 1:G) {
          D_g[[g]] <- beg_D
        }
      }
      else if (code_D == "D") {
        p_D <- K
        for (g in 1:G) {
          D_g[[g]] <- diag(mean(diag(int_D_g[[g]])), K)
        }
      }
      else if (code_D == "C") {
        p_D <- 1
        beg_D <- 0
        for (g in 1:G) {
          beg_D <- int_D_g[[g]] * pi_g[g] + beg_D
        }
        for (g in 1:G) {
          D_g[[g]] <- diag(mean(diag(beg_D)), K)
        }
      }




      pre_B_g <- B_g

      #change B_g
      #H2, loglikelyhood
      if (code_B == "U") {
        p_B <- G * K
        for (g in 1:G) {
          for (r in 1:K) {
            x <- b(r, B_g[[g]])

            #use lapply to calculate
            track_H2 <- lapply(x, function(x, i) {
              invd <- diag(diag(D_g[[i]]) ^ (-1))
              invt <- diag(diag(T_g[[i]]) ^ (-1), Q_g[g], Q_g[g])
              invcov <- solvecov(t(x), invd, invt)
              H2 <- tr(invd %*% x %*% invcov %*% t(x) %*%
                         invd %*% S_g[[i]]) -
                tr(invd %*% x %*% e2[[i]] %*% t(x)) / 2
              return(ifelse(is.na(H2), -Inf, H2))
            }, i = g)
            B_g[[g]] <- x[[which.max(track_H2)]]

          }

        }
      }
      else if (code_B == "G" & length(unique(Q_g)) == 1) {
        p_B <- K
        for (r in 1:K) {
          x <- b(r, B_g[[1]])
          track_H2 <- lapply(x, function(x) {
            H2 <- 0
            for (g in 1:G) {
              invd <- diag(diag(D_g[[g]]) ^ (-1))
              invt <- diag(diag(T_g[[g]]) ^ (-1), Q_g[g], Q_g[g])
              invcov <- solvecov(t(x), invd, invt)
              H2 <- H2 + sum(z_ig[, g]) * (tr(invd %*% x %*% invcov %*% t(x) %*%
                                                invd %*% S_g[[g]]) - tr(invd %*% x %*%
                                                                          e2[[g]] %*% t(x)) / 2)
            }
            return(ifelse(is.na(H2), -Inf, H2))
          })
          B_g[[1]] <- x[[which.max(track_H2)]]
          for (i in 1:G) {
            B_g[[i]] <- B_g[[1]]
          }
          #print(c(r,unlist(track_H2)))
        }
      }


      #instead of ordering, if the overall likelihood decrease, then keep the B_g as last one
      #if (code_B == "G") {
        if (p_comp_logw(W_count, m, V, B_g, T_g, D_g, mu_g, pi_g, G, n, K, const) >= overall_loglik[h -
                                                                                            1]) {
          B_g <- B_g
        }
        else{
          B_g <- pre_B_g
        }


      o_g <- order(pi_g, decreasing = T)
      mu_g <- mu_g[o_g]
      pi_g <- pi_g[o_g]
      z_ig <- z_ig[, o_g]
      B_g <- B_g[o_g]
      D_g <- D_g[o_g]
      T_g <- T_g[o_g]
      Q_g <- Q_g[o_g]
      m <- m[o_g]
      V <- V[o_g]
      U_iK <- U_iK[o_g]
      beta_g<-beta_g[o_g]

      COV <- list()
      for (g in 1:G) {
        COV[[g]] <- B_g[[g]] %*% T_g[[g]] %*% t(B_g[[g]]) + D_g[[g]]
      }



      overall_loglik[h] <-
        p_comp_logw(W_count, m, V, B_g, T_g, D_g, mu_g, pi_g, G, n, K, const)
      if (is.na(overall_loglik[h])) {
        overall_loglik[h] <- overall_loglik[h - 1]
      }

      if (h > 3) {
        a <-
          (overall_loglik[h] - overall_loglik[h - 1]) / (overall_loglik[h - 1] - overall_loglik[h -
                                                                                                  2])
        L_inf <-
          overall_loglik[h - 1] + (overall_loglik[h] - overall_loglik[h - 1]) / (1 -
                                                                                   a)
      } else{
        L_inf <- Inf
      }

      diff <- L_inf - overall_loglik[h]
      if (is.na(diff)) {
        diff <- Inf
      }
      if (abs(diff) < 10 ^ (-2) | h > iter |
          length(pi_g[pi_g == 0]) > 0) {
        break
      }
      h <- h + 1



    }


    if (length(pi_g[pi_g == 0]) > 0 |
        any(table(label) < 5) |
        length(unique(label)) < G) {
      return("Choose fewer components to fit")
    }
    else{
      ICL <-
        overall_loglik[h] * 2 - log(n) * (p_mu + G - 1 + (p_B + p_T + p_D)) + 2 *
        sum(as.vector(z_ig) * log(as.vector(z_ig)))
      BIC <- overall_loglik[h] * 2 - log(n) * (p_mu + G - 1 + (p_B + p_T +
                                                                  p_D))
      AIC <- overall_loglik[h] * 2 - 2 * (p_mu + G - 1 + (p_B + p_T + p_D))

      if (is.vector(z_ig)) {
        lab_rec <- rep(1, n)
      }
      else{
        lab_rec <- apply(z_ig, 1, which.max)
      }


      return(
        list(
          z_ig = z_ig,
          cluster = lab_rec,
          mu_g = mu_g,
          pi_g = pi_g,
          B_g = B_g,
          T_g = T_g,
          D_g = D_g,
          COV = COV,
          beta_g=beta_g,
          sigma = sig_g,
          overall_loglik = overall_loglik,
          ICL = ICL,
          BIC = BIC,
          AIC = AIC
        )
      )
    }


  }


#constant permutation for each i
constant_fun <- function(W, K) {
  funforlb <-
    function(k) {
      if (W[k] == 0) {
        return(0)
      } else{
        return(sum(log(1:W[k])))
      }
    }
  const <- sum(log(1:sum(W))) - sum(sapply(1:(K + 1), funforlb))
  const
}



#calculate ELBO for i and g observation
lb_loglik <- function(V, K, sig, m, mu, W, const) {
  wml <-
    log(det(V)) / 2 + K / 2 - log(det(sig)) / 2 - t(m - mu) %*% MASS::ginv(sig) %*%
    (m - mu) / 2 - tr(MASS::ginv(sig) %*% V) / 2 + W[1:K] %*% m -
    sum(W) * (log(sum(c(exp(
      m + diag(V) / 2
    ), 1)))) + const
  if (is.na(wml)) {
    wml <- -Inf
  }
  return(wml)

}


#complete log likelyhood for lower bound
elbo_fun <- function(pi_g, G, n, log_w_approx) {
 # .Call('elbo_fun', PACKAGE = 'lnmbiclust', pi_g, G, n, log_w_approx)
  forelbo <- matrix(nrow = n, ncol = G)
  for (i in 1:n) {
    for (g in 1:G) {
      if (pi_g[g] * exp(log_w_approx[i, g]) == 0) {
        forelbo[i, g] <- pi_g[g] * .Machine$double.xmin
      }
      else{
        forelbo[i, g] <- pi_g[g] * exp(log_w_approx[i, g])
      }

    }
  }

  mix_den <- rowSums(forelbo)
  return(sum(log(mix_den)))
}

#update z_ig
z_ig_update <- function(pi_g, n, G, log_w_approx) {
  z_ig <- matrix(0,nrow = n, ncol = G)
  for (i in 1:n) {
    for (g in 1:G) {
      if (pi_g[g] * exp(log_w_approx[i, g]) == 0) {
        z_ig[i, g] <- pi_g[g] * .Machine$double.xmin
      }
      else{
        z_ig[i, g] <- pi_g[g] * exp(log_w_approx[i, g])
      }
    }
    if (sum(z_ig[i, ]) == 0) {
      z_ig[i, ] <- pi_g * .Machine$double.xmin
    }
  }

  mix_den <- rowSums(z_ig)
  z_ig <- z_ig / mix_den
  return(z_ig)
}

#update mu_g
mu_update <- function(z_ig, m, g) {
  mu_g <- m[[g]] %*% z_ig[, g] / sum(z_ig[, g])
  return(mu_g)
}


####trace function####
tr <- function(x) {
  sum(diag(x))
}


####update B_g function####
b <- function(r, B_g) {
  #define the number of column for B_g
  q <- length(B_g[r, ])
  B_test <- vector('list', q)
  #create q identical B_g
  for (e in 1:q) {
    B_test[[e]] <- B_g
    B_test[[e]][r, ] <- 0
    B_test[[e]][r, e] <- 1
  }
  B_test
}

####invert factor structure covariance matrix####
solvecov <- function(b, t, d) {
  q <- dim(t)[1]
  p <- dim(d)[1]
  invd <- diag(diag(d) ^ (-1), p, p)
  invt <- diag(diag(t) ^ (-1), q, q)
  invd - invd %*% b %*% MASS::ginv(invt + t(b) %*% invd %*% b) %*% t(b) %*%
    invd
}



#one step newtown raphson to update m and V
GD_mvxi <- function(V,
                    m,
                    sig_g,
                    W_count,
                    mu_g,
                    K,
                    B_g,
                    T_g,
                    D_g) {
  V_old <- sqrt(V)
  m_old <- m


  f_V <- diag(V_old) ^ (-1) - diag(MASS::ginv(sig_g)) * diag(V_old) -
    c(sum(W_count) * diag(V_old) * c(exp(m_old + diag(V_old ^ 2) /
                                                2) /
                                            sum(c(
                                              exp(m_old + diag(V_old ^ 2) / 2), 1
                                            ))))
  s_V <-
    -1 / diag(V_old) ^ 2 - diag(MASS::ginv(sig_g)) - (diag(V_old) ^ 2 +
                                                             1) *
    sum(W_count) * exp(m_old + diag(V_old) ^ 2 / 2) / sum(c(exp(m_old +
                                                                       diag(V_old ^ 2) / 2), 1))


  if (any(is.na(f_V / s_V))) {
    V_new <- V_old
  }
  else{
    V_new <- V_old - diag(f_V / s_V)
  }


  f_m <- -solvecov(B_g, T_g, D_g) %*% (m_old - mu_g) +
    c(W_count[1:K]) - c(sum(W_count) * c(exp(m_old + diag(V_old ^
                                                                    2) / 2) /
                                                   sum(c(
                                                     exp(m_old + diag(V_old ^ 2) / 2), 1
                                                   ))))

  s_m <- -solvecov(B_g, T_g, D_g) -
    sum(W_count) * (diag(exp(m_old + diag(V_old ^ 2) / 2)) / sum(c(exp(
      m_old + diag(V_old ^ 2) / 2
    ), 1)))

  if (det(s_m) == 0 | is.na(det(s_m)) | det(s_m) == Inf) {
    m_new <- m_old
  } else{
    m_new <- m_old - MASS::ginv(s_m) %*% f_m
  }

  list(m = m_new, V = V_new ^ 2)
}




p_comp_logw <- function(W_count,
                        m,
                        V,
                        B_g,
                        T_g,
                        D_g,
                        mu_g,
                        pi_g,
                        G,
                        n,
                        K,
                        const) {
  COV <- list()
  for (g in 1:G) {
    COV[[g]] <- B_g[[g]] %*% T_g[[g]] %*% t(B_g[[g]]) + D_g[[g]]
  }
  #calculate overall log likelyhood.
  log_w_approx <- matrix(0, nrow = n, ncol = G)
  for (i in 1:n) {
  for (g in 1:G) {
      log_w_approx[i, g] <-
        lb_loglik(V[[g]][[i]], K, COV[[g]], m[[g]][, i], mu_g[[g]][,i], W_count[i, ], const[i])
    }
  }
  return(elbo_fun(pi_g, G, n, log_w_approx))
}
