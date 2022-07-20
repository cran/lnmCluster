#' run main microbiome Factor Analyzer algorithm.
#'
#'@param W_count The microbiome count matrix that you want to analyze.
#'@param G The number of component
#'@param Q_g The number of latent dimensions for each component, a vector.
#'@param pi_g A vector of initial guesses of component proportion
#'@param mu_g A list of initial guess of mean vector
#'@param sig_g A list of initial guess of covariance matrix for each component
#'@param m A list of initial guess of variational mean
#'@param V A list of initial guess of variational varaince
#'@param B_K A list of initial guess of loading matrix.
#'@param T_K A list of identity matrix with dimension q.
#'@param D_K A list of initial guess of error matrix
#'@param cov_str The covaraince structure you choose, there are 8 different models belongs to
#'this family:UUU, UUG, UUD, UUC, GUU, GUG, GUD, GUC.
#'@param iter Max iterations, default is 150.
#'@param const the permutation constant in multinomial distribution. Calculated before the main
#'algorithm in order to save computation time.
#'@param beta_g initial guess of covariates coefficients.
#'@param X The regression covariates matrix, which generates by model.matrix.
#'
#'@return z_ig Estimated latent variable z
#'@return cluster Component labels
#'@return mu_g Estimated component mean
#'@return pi_g Estimated component proportion
#'@return B_g Estimated loading matix.
#'@return D_g Estimated error covariance
#'@return COV Estimated component covariance
#'@return beta_g Estimated covariates coefficients.
#'@return overall_loglik Complete log likelihood value for each iteration
#'@return ICL ICL value
#'@return BIC BIC value
#'@return AIC AIC value




#main algorithm
Mico_bi_PGMM <-
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
           cov_str, iter, const,beta_g,X) {
    n <- dim(W_count)[1]
    K <- dim(W_count)[2] - 1

    code_split <- str_split(cov_str, "")[[1]]
    code_B <- code_split[1]
    code_T <- code_split[2]
    code_D <- code_split[3]




    overall_loglik <- rep(-Inf, 2)

    h = 2

    repeat {
      COV = list()
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

      #print(c("mv",p_comp_logw(W_count,alpha,m,V,K,beta_g,B_g,COV,mu_g,pi_g,G,n)))

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
      U_iK = list()
      for (g in 1:G) {
        # U_iK[[g]]=T_K[[g]]%*%t(B_K[[g]])%*%
        #   solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%(m[[g]]-c(mu_g[[g]]))


        U_iK[[g]] <-
          solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                    (-1), Q_g[g], Q_g[g])) %*%
          t(B_K[[g]]) %*%
          diag(diag(D_K[[g]]) ^ (-1)) %*% (m[[g]] - mu_g[[g]])

      }

      #S_g
      z_ig <- as.matrix(z_ig, K, n)
      S_g = list()
      for (g in 1:G) {
        S_g[[g]] <-
          (m[[g]] - mu_g[[g]])%*% diag(z_ig[, g] / sum(z_ig[, g])) %*%
          t(m[[g]] - mu_g[[g]])
      }

      #sig_g
      sig_g <- list()
      for (g in 1:G) {
        beg = 0
        for (i in 1:n) {
          beg <- beg + z_ig[i, g] * (V[[g]][[i]] + (m[[g]][, i] - mu_g[[g]][,i]) %*%
                                       t(m[[g]][, i] - mu_g[[g]][,i])) / sum(z_ig[, g])
        }
        sig_g[[g]] <- beg
      }

      #e2
      e2 <- list()
      for (g in 1:G) {
        # e2[[g]]=T_K[[g]]%*%t(B_K[[g]])%*%solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%
        #   S_g[[g]]%*%t(solvecov(B_K[[g]],T_K[[g]],D_K[[g]]))%*%B_K[[g]]%*%T_K[[g]]+
        #   T_K[[g]]-T_K[[g]]%*%t(B_K[[g]])%*%solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%
        #   B_K[[g]]%*%T_K[[g]]


        e2[[g]] <-
          solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                    (-1), Q_g[g], Q_g[g])) +
          solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                    (-1), Q_g[g], Q_g[g])) %*%
          t(B_K[[g]]) %*%
          diag(diag(D_K[[g]]) ^ (-1)) %*% S_g[[g]] %*% t(diag(diag(D_K[[g]]) ^
                                                                (-1))) %*% B_K[[g]] %*%
          t(solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^ (-1)), diag(diag(T_K[[g]]) ^
                                                                      (-1), Q_g[g], Q_g[g])))
      }




      #D_K
      pre_D_K <- D_K
      for (g in 1:G) {
        # D_K[[g]]=diag(diag(sig_g[[g]]-B_K[[g]]%*%T_K[[g]]%*%t(B_K[[g]])%*%
        #                      solvecov(B_K[[g]],T_K[[g]],D_K[[g]])%*%S_g[[g]]))


        D_K[[g]] <- abs(diag(
          diag(
            sig_g[[g]] - 2 * B_K[[g]] %*%
              solvecov(t(B_K[[g]]), diag(diag(D_K[[g]]) ^
                                           (-1)),
                       diag(diag(T_K[[g]]) ^ (-1), Q_g[g], Q_g[g])) %*%
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

      #change B_K
      #H2, loglikelyhood
      beta <- list()
      if (code_B == "U") {
        p_B <- G * (K * Q_g[1] - Q_g[1] * (Q_g[1] - 1) / 2)
        for (g in 1:G) {
          beta[[g]] <- t(B_K[[g]]) %*% solve(B_K[[g]] %*% t(B_K[[g]]) + D_K[[g]])
          B_K[[g]] <- S_g[[g]] %*% t(beta[[g]]) %*%
            solve(diag(1, Q_g[g]) - beta[[g]] %*% B_K[[g]] + beta[[g]] %*% S_g[[g]] %*%
                    t(beta[[g]]))

        }
      }
      else if (code_B == "G" &
               length(unique(Q_g)) == 1 & (code_D == "G" | code_D == "C")) {
        p_B <- K * Q_g[1] - Q_g[1] * (Q_g[1] - 1) / 2
        for (g in 1:G) {
          beta[[g]] <- t(B_K[[g]]) %*% solve(B_K[[g]] %*% t(B_K[[g]]) + D_K[[g]])
          B_K[[g]] <- S_g[[g]] %*% t(beta[[g]]) %*%
            solve(diag(1, Q_g[g]) - beta[[g]] %*% B_K[[g]] + beta[[g]] %*% S_g[[g]] %*%
                    t(beta[[g]]))

        }
        beg_B <- 0
        int_B_g <- B_K
        for (g in 1:G) {
          beg_B <- int_B_g[[g]] * pi_g[g] + beg_B
        }
        for (g in 1:G) {
          B_K[[g]] <- beg_B
        }
      }
      else if (code_B == "G" &
               length(unique(Q_g)) == 1 & (code_D == "D" | code_D == "U")) {
        p_B <- K * Q_g[1] - Q_g[1] * (Q_g[1] - 1) / 2
        ng <- colSums(z_ig)
        ri <- 0
        for (g in 1:G) {
          beta[[g]] <- t(B_K[[g]]) %*% solve(B_K[[g]] %*% t(B_K[[g]]) + D_K[[g]])
          ri <- ri + ng[g] * solve(D_K[[g]]) %*% S_g[[g]] %*% t(beta[[g]])
        }

        lamb <- matrix(0, ncol = Q_g[1], nrow = K)
        for (i in 1:K) {
          dei <- 0
          for (g in 1:G) {
            dei <- dei + (ng[g] / D_K[[g]][i, i]) *
              (diag(1, Q_g[g]) - beta[[g]] %*% B_K[[g]] + beta[[g]] %*% S_g[[g]] %*%
                 t(beta[[g]]))
          }

          lamb[i, ] <- ri[i, ] %*% solve(dei)
        }
        for (g in 1:G) {
          B_K[[g]] <- lamb
        }

      }




      o_g <- order(pi_g, decreasing = T)
      mu_g <- mu_g[o_g]
      pi_g <- pi_g[o_g]
      z_ig <- z_ig[, o_g]
      B_K <- B_K[o_g]
      D_K <- D_K[o_g]
      Q_g <- Q_g[o_g]
      m <- m[o_g]
      V <- V[o_g]
      beta_g<-beta_g[o_g]




      overall_loglik[h] <-
        p_comp_logw(W_count, m, V, B_K, T_K, D_K, mu_g, pi_g, G, n, K, const)
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
        overall_loglik[h] * 2 - log(n) * (p_mu + G - 1 + (p_B + p_D)) + 2 * sum(as.vector(z_ig) *
                                                                                   log(as.vector(z_ig)))
      BIC <- overall_loglik[h] * 2 - log(n) * (p_mu + G - 1 + (p_B + p_D))
      AIC <- overall_loglik[h] * 2 - 2 * (p_mu + G - 1 + (p_B + p_D))

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
          T_K = T_K,
          D_K = D_K,
          COV = COV,
          beta_g=beta_g,
          #sigma=sig_g,
          overall_loglik = overall_loglik,
          ICL = ICL,
          BIC = BIC,
          AIC = AIC
        )
      )
    }


  }

#all functions inside the main function is exactly the same as Micro_bi_jensens
#lb_loglik, elbo_fun, z_ig_update, mu_update, tr, solvecov, GD_mvxi, and p_comp_logw, constant_fun



