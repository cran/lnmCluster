#' Penalized Logistic Normal Multinomial factor analyzer main estimation process
#'
#' Main function will perform PLNM factor analyzer and return parameters
#' @param W_count The microbiome count matrix
#'
#' @param G All possible number of components. A vector.
#' @param Q_g A specific number of latent dimension.
#' @param pi_g A vector of initial guesses of component proportion
#' @param mu_g A list of initial guess of mean vector
#' @param sig_g A list of initial guess of covariance matrix for each component
#' @param m A list of initial guess of variational mean
#' @param V A list of initial guess of variational varaince
#' @param B_K A list of initial guess of loading matrix.
#' @param T_K A list of identity matrix with dimension q.
#' @param D_K A list of initial guess of error matrix
#' @param cov_str The covaraince structure you choose, there are 2 different models belongs to
#'this family:UUU and GUU. You can choose more than 1 covarance structure to do model selection.
#' @param tuning length G vector with range 0-1, define the tuning parameter for each component
#' @param iter Max iterations, default is 150.
#' @param const the permutation constant in multinomial distribution. Calculated before the main
#'algorithm in order to save computation time.
#'@param beta_g initial guess of covariates coefficients.
#'@param X The regression covariates matrix, which generates by model.matrix.
#'
#'
#'
#'@return z_ig Estimated latent variable z
#'@return cluster Component labels
#'@return mu_g Estimated component mean
#'@return pi_g Estimated component proportion
#'@return B_g Estimated sparsity loading matrix
#'@return D_g Estimated error covariance
#'@return COV Estimated component covariance
#'@return beta_g Estimated covariates coefficients.
#'@return overall_loglik Complete log likelihood value for each iteration
#'@return ICL ICL value
#'@return BIC BIC value
#'@return AIC AIC value
#'@return tuning display the tuning parameter you specified.



#main algorithm
Mico_bi_lasso <-
  function(W_count,
           G,
           Q_g,
           pi_g,
           mu_g,
           sig_g,
           V,
           m,
           B_K,
           T_K,
           D_K,
           cov_str,
           tuning,iter,const,beta_g,X) {
    n <- dim(W_count)[1]
    K <- dim(W_count)[2] - 1

    code_split <- str_split(cov_str, "")[[1]]
    code_B <- code_split[1]
    code_T <- code_split[2]
    code_D <- code_split[3]

    tuning <- tuning / (1 - tuning)



    overall_loglik <- rep(-Inf, 2)

    h <- 2

    repeat {
      COV <- list()
      for (g in 1:G) {
        COV[[g]] <- B_K[[g]] %*% T_K[[g]] %*% t(B_K[[g]]) + D_K[[g]]
      }

      #update m and V
      for (g in 1:G) {
        for (i in 1:n) {
          mv <- GD_mvxi(V[[g]][[i]], m[[g]][, i], COV[[g]], W_count[i, ], mu_g[[g]][,i], K,
                        B_K[[g]], T_K[[g]], D_K[[g]])
          m[[g]][, i] <- c(mv$m)
          V[[g]][[i]] <- mv$V
        }
      }

      #print(c("mv",p_incomp_logw(W_count,alpha,m,V,K,beta_g,B_g,COV,mu_g,pi_g,G,n)))

      log_w_approx <- matrix(0, nrow = n, ncol = G)
      for (g in 1:G) {
        for (i in 1:n) {
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
        # U_iK[[g]]<-T_K[[g]]%*%t(B_K[[g]])%*%
        #   solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%(m[[g]]-c(mu_g[[g]]))


        U_iK[[g]] <-
          solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                    (-1), Q_g, Q_g)) %*%
          t(B_K[[g]]) %*%
          diag(diag(D_K[[g]]) ^ (-1)) %*% (m[[g]] - mu_g[[g]])

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
          beg <- beg + z_ig[i, g] * (V[[g]][[i]] + (m[[g]][, i] - c(mu_g[[g]][,i])) %*%
                                       t(m[[g]][, i] - c(mu_g[[g]][,i]))) / sum(z_ig[, g])
        }
        sig_g[[g]] <- beg
      }

      #e2
      e2 <- list()
      for (g in 1:G) {
        # e2[[g]]<-T_K[[g]]%*%t(B_K[[g]])%*%solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%
        #   S_g[[g]]%*%t(solvecov(B_K[[g]],T_K[[g]],D_K[[g]]))%*%B_K[[g]]%*%T_K[[g]]+
        #   T_K[[g]]-T_K[[g]]%*%t(B_K[[g]])%*%solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%
        #   B_K[[g]]%*%T_K[[g]]


        e2[[g]] <-
          solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                    (-1), Q_g, Q_g)) +
          solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                    (-1), Q_g, Q_g)) %*%
          t(B_K[[g]]) %*%
          diag(diag(D_K[[g]]) ^ (-1)) %*% S_g[[g]] %*% t(diag(diag(D_K[[g]]) ^
                                                                (-1))) %*% B_K[[g]] %*%
          t(solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                      (-1), Q_g, Q_g)))
      }




      #D_K
      pre_D_K <- D_K
      for (g in 1:G) {
        # D_K[[g]]<-diag(diag(sig_g[[g]]-B_K[[g]]%*%T_K[[g]]%*%t(B_K[[g]])%*%
        #                      solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%S_g[[g]]))


        D_K[[g]] <- abs(diag(
          diag(
            sig_g[[g]] - 2 * B_K[[g]] %*%
              solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^
                                           (-1)),
                       diag(diag(T_K[[g]]) ^ (-1), Q_g, Q_g)) %*%
              t(B_K[[g]]) %*%
              diag(diag(D_K[[g]]) ^ (-1)) %*% S_g[[g]] +
              B_K[[g]] %*% e2[[g]] %*% t(B_K[[g]])
          )
        ))


      }




      p_D <- G * K


      int_D_g <- D_K
      if (code_D == "G") {
        p_D <- K
        beg_D <- 0
        for (g in 1:G) {
          beg_D <- int_D_g[[g]] * pi_g[g] + beg_D
        }
        for (g in 1:G) {
          D_K[[g]] <- beg_D
        }
      }
      else if (code_D == "D") {
        p_D <- K
        for (g in 1:G) {
          D_K[[g]] <- diag(mean(diag(int_D_g[[g]])), K)
        }
      }
      else if (code_D == "C") {
        p_D <- 1
        beg_D <- 0
        for (g in 1:G) {
          beg_D <- int_D_g[[g]] * pi_g[g] + beg_D
        }
        for (g in 1:G) {
          D_K[[g]] <- diag(mean(diag(beg_D)), K)
        }
      }




      pre_B_K <- B_K
      beta <- list()
      theta <- list()
      #lasso<-list()
      if (code_B == "U") {
        ng <- pi_g * n
        for (g in 1:G) {
          beta[[g]] <- t(B_K[[g]]) %*% ginv(B_K[[g]] %*% t(B_K[[g]]) + D_K[[g]])
          theta[[g]] <-
            diag(1, Q_g) - beta[[g]] %*% B_K[[g]] + beta[[g]] %*% S_g[[g]] %*% t(beta[[g]])

          for (u in 1:K) {
            for (v in 1:Q_g) {
              x <- ((S_g[[g]] %*% t(beta[[g]]))[u, v] - tuning[g] * D_K[[g]][u, u] / ng[g] -
                      B_K[[g]][u, -v] %*% theta[[g]][-v, v]) / theta[[g]][v, v]
              y <-
                ((S_g[[g]] %*% t(beta[[g]]))[u, v] + tuning[g] * D_K[[g]][u, u] / ng[g] -
                   B_K[[g]][u, -v] %*% theta[[g]][-v, v]) / theta[[g]][v, v]
              if (x > 0) {
                B_K[[g]][u, v] <- x
              }
              else if (y < 0) {
                B_K[[g]][u, v] <- y
              }
              else{
                B_K[[g]][u, v] <- 0
              }
            }
          }
        }
        p_B <- sum(unlist(lapply(B_K, function(x) {
          o_star <- sum(x != 0)
          q_star <- sum(colSums(x) > 0)
          o_star - q_star * (q_star - 1) / 2
        })))#number of parameter that account for orthonormal matrix
      }
      else if (code_B == "G" & length(unique(Q_g)) == 1 &
               length(unique(tuning)) == 1 & (code_D == "G" |
                                              code_D == "C")) {
        beg_S <- 0
        int_S_g <- S_g
        for (g in 1:G) {
          beg_S <- int_S_g[[g]] * pi_g[g] + beg_S
        }

        beta <- t(B_K[[1]]) %*% ginv(B_K[[1]] %*% t(B_K[[1]]) + D_K[[1]])
        theta <- diag(1, Q_g) - beta %*% B_K[[1]] + beta %*% beg_S %*% t(beta)

        for (u in 1:K) {
          for (v in 1:Q_g) {
            x <- ((beg_S %*% t(beta))[u, v] - tuning[1] * D_K[[1]][u, u] / n -
                    B_K[[1]][u, -v] %*% theta[-v, v]) / theta[v, v]
            y <- ((beg_S %*% t(beta))[u, v] + tuning[1] * D_K[[1]][u, u] /
                    n -
                    B_K[[1]][u, -v] %*% theta[-v, v]) / theta[v, v]
            if (x > 0) {
              B_K[[1]][u, v] <- x
            }
            else if (y < 0) {
              B_K[[1]][u, v] <- y
            }
            else{
              B_K[[1]][u, v] <- 0
            }
          }
        }
        for (g in 1:G) {
          B_K[[g]] <- B_K[[1]]
        }
        p_B <- lapply(B_K, function(x) {
          o_star <- sum(x != 0)
          q_star <- sum(colSums(x) > 0)
          o_star - q_star * (q_star - 1) / 2
        })[[1]]#number of parameter that account for orthonormal matrix
      }
      else if (code_B == "G" & length(unique(Q_g)) == 1 &
               length(unique(tuning)) == 1 & (code_D == "D" |
                                              code_D == "U")) {
        ng <- pi_g * n
        Ci <- 0
        for (g in 1:G) {
          beta[[g]] <- t(B_K[[g]]) %*% ginv(B_K[[g]] %*% t(B_K[[g]]) + D_K[[g]])
          theta[[g]] <-
            diag(1, Q_g) - beta[[g]] %*% B_K[[g]] + beta[[g]] %*% S_g[[g]] %*% t(beta[[g]])
          Ci <- Ci + ng[g] * ginv(D_K[[g]]) %*% S_g[[g]] %*% t(beta[[g]])
        }
        for (u in 1:K) {
          for (v in 1:Q_g) {
            ri <- 0
            de <- 0
            for (g in 1:G) {
              ri <- ri + ng[g] * (1 / D_K[[g]][u, u]) * (B_K[[g]][u, -v] %*% theta[[g]][-v, v])
              de <- de + ng[g] * (1 / D_K[[g]][u, u]) * theta[[g]][v, v]
            }
            x <- (Ci[u, v] - ri - tuning[1]) / de
            y <- (Ci[u, v] - ri + tuning[1]) / de
            if (x > 0) {
              B_K[[1]][u, v] <- x
            }
            else if (y < 0) {
              B_K[[1]][u, v] <- y
            }
            else{
              B_K[[1]][u, v] <- 0
            }

          }
        }
        for (g in 1:G) {
          B_K[[g]] <- B_K[[1]]
        }
        p_B <- lapply(B_K,function(x){
          no0rc <- x[rowSums(x)>0,colSums(x)>0]
          zero <- sum(no0rc==0)
          p_star <- dim(no0rc)[1]
          q_star <- dim(no0rc)[2]
          p_star*q_star-q_star*(q_star-1)/2-zero
        })[[1]]#number of parameter that account for orthonormal matrix
      }

      if (code_B == "U") {
        pnt <-
          sum(tuning * unlist(lapply(B_K, function(x) {
            sum(abs(as.vector(x)))
          })))
      }
      else if (code_B == "G") {
        pnt <-
          sum(tuning * unlist(lapply(B_K, function(x) {
            sum(abs(as.vector(x)))
          }))) / G
      }
      overall_loglik[h] <-
        p_comp_logw(W_count, m, V, B_K, T_K, D_K, mu_g, pi_g, G, n, K, const) - pnt
      if (overall_loglik[h] > overall_loglik[h - 1]) {
        B_K <- B_K
      }
      else{
        B_K <- pre_B_K
      }



      o_g <- order(pi_g, decreasing = T)
      mu_g <- mu_g[o_g]
      pi_g <- pi_g[o_g]
      z_ig <- z_ig[, o_g]
      B_K <- B_K[o_g]
      D_K <- D_K[o_g]
      T_K <- T_K[o_g]
      m <- m[o_g]
      V <- V[o_g]
      U_iK <- U_iK[o_g]
      tuning <- tuning[o_g]
      beta_g=beta_g[o_g]


      if (code_B == "U") {
        pnt <-
          sum(tuning * unlist(lapply(B_K, function(x) {
            sum(abs(as.vector(x)))
          })))
      }
      else if (code_B == "G") {
        pnt <-
          sum(tuning * unlist(lapply(B_K, function(x) {
            sum(abs(as.vector(x)))
          }))) / G
      }
      overall_loglik[h] <-
        p_comp_logw(W_count, m, V, B_K, T_K, D_K, mu_g, pi_g, G, n, K, const) - pnt
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
        (overall_loglik[h] + pnt) * 2 - log(n) * (p_mu + G - 1 + (p_B + p_D)) +
        2 * sum(as.vector(z_ig) * log(as.vector(z_ig)))
      BIC <- (overall_loglik[h] + pnt) * 2 - log(n) * (p_mu + G - 1 + (p_B +
                                                                          p_D))
      AIC <- (overall_loglik[h] + pnt) * 2 - 2 * (p_mu + G - 1 + (p_B + p_D))

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
          B_K = B_K,
          D_K = D_K,
          COV = COV,
          beta_g=beta_g,
          #sigma = sig_g,
          overall_loglik = overall_loglik,
          ICL = ICL,
          BIC = BIC,
          AIC = AIC,
          tuning = tuning/(tuning+1)
        )
      )
    }


  }
