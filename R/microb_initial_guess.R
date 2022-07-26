#' Gives default initial guesses for logistic-normal multinomial biclustering algorithm.
#'
#'@param W_count The microbiome count matrix that you want to analyze.
#'@param G The number of component
#'@param Q_g The number of biclusters for each component, a vector.
#'@param cov_str The covaraince structure you choose, there are 16 different models belongs to
#'this family:UUU, UUG, UUD, UUC, UGU, UGG, UGD, UGC, GUU, GUG, GUD, GUC, GGU, GGG, GGD, GGC.
#'@param X The regression covariates matrix, which generated by model.matrix.
#'@return \code{new_pi_g} Initial guess of proportion
#'@return \code{new_mu_g} Initial guess of mean vector
#'@return \code{new_sig_g} Initial guess of covariance matrix for each component
#'@return \code{new_T_g} Initial guess of covariance of latent variable: u
#'@return \code{new_B_g} Initial guess of bicluster membership
#'@return \code{new_D_g} Initial guess of error matrix
#'@return \code{new_m} Initial guess of variational mean
#'@return \code{new_V} Initial guess of variational varaince
#'@return \code{new_beta_g} Initial guess of covariates coefficients.

#'
#@importFrom mvtnorm rmvnorm


initial_variational_gaussian <- function(W_count, G, Q_g, cov_str,X) {
  W_count <- as.matrix(W_count)
  n <- dim(W_count)[1]
  K <- dim(W_count)[2] - 1

  V <- vector("list", G)
  m <- vector("list", G)
  T_g <- list()
  z_ig <- matrix(0, nrow = n, ncol = G)
  pi_g <- rep(0, G)
  B_g <- list()
  D_g <- list()
  mu_g <- list()
  sig_g <- list()
  beta_g=list()

  code <- stringr::str_split(cov_str, "")[[1]]
  code_B <- code[1]
  code_T <- code[2]
  code_D <- code[3]


  if ((code_T == "G" | code_T == "C") & length(unique(Q_g)) != 1)
  {
    return("Q has to be the same for all components when fitting _G_ models")
  }
  else if (code_B == "G" & length(unique(Q_g)) != 1)
  {
    return("Q has to be the same for all components when fitting G__ models")
  }
  else{
    Z <- W_count / rowSums(W_count)
    Z[Z == 0] <- 0.001
    Y <- log(Z[, 1:K] / c(Z[, K + 1]))


      invisible(capture.output(res_class <-pgmm::pgmmEM(
        Y,
        rG = G,
        rq = max(Q_g),
        zstart = 2,
        modelSubset = c("UUU"),
        relax = TRUE
      )))
    lab <- res_class$map

    # res_class <- try(mclust::Mclust(Y, G, modelNames = "VVV"))
    # lab <- res_class$classification

    if(is.null(res_class)){
      res_class <- kmeans(Y,G)
      lab <- res_class$cluster
    }


    #initial for pi_g
    pi_g <- table(lab)/n


    #initial for m
    for (g in 1:G) {
      m[[g]] <- t(Y)
    }

    #transform lab into z_ig
    for(g in 1:G){
      zzz<-lab
      zzz[zzz!=g]<-0
      zzz[zzz!=0]<-1
      z_ig[,g]<-zzz
    }

    #initial for mu_g, beta_g
    for (g in 1:G) {
      beta_g[[g]]<-ginv(t(X)%*%diag(z_ig[,g])%*%X)%*%t(X)%*%diag(z_ig[,g])%*%t(m[[g]])
      mu_g[[g]]<-t(X%*%beta_g[[g]])
    }

    #initial for sig_g
    for (g in 1:G) {
      sig_g[[g]] <- cov(Y[lab == g, ])
    }


    #initial for B_g
    for (g in 1:G) {
      fit <- eigen(cov(scale(Y[lab == g, ])))
      bk <- fit$vectors[, 1:Q_g[g]]
      list_result <- lapply(split(bk, seq(NROW(bk))), function(x, i) {
        y <- rep(0, Q_g[i])
        y[which.max(x)] <- 1
        y
      }, i = g)
      B_g[[g]] <- do.call(rbind, list_result)
    }
    if (code_B == "G" & length(unique(Q_g)) == 1) {
      for (i in 1:G) {
        B_g[[i]] <- B_g[[1]]
      }
    }



    #initial for T_u_g
    for (g in 1:G) {
      T_g[[g]] <- diag(eigen(sig_g[[g]])$values[1:Q_g[g]],
                       Q_g[g], Q_g[g])
    }
    int_T_g <- T_g
    if (code_T == "G" & length(unique(Q_g)) == 1) {
      beg_T <- 0
      for (i in 1:G) {
        beg_T <- int_T_g[[i]] * pi_g[i] + beg_T
      }
      for (i in 1:G) {
        T_g[[i]] <- beg_T
      }
    }
    else if (code_T == "D") {
      for (i in 1:G) {
        T_g[[i]] <- diag(mean(diag(int_T_g[[i]])), Q_g[i])
      }
    }
    else if (code_T == "C" & length(unique(Q_g)) == 1) {
      beg_T <- 0
      for (i in 1:G) {
        beg_T <- int_T_g[[i]] * pi_g[i] + beg_T
      }
      for (i in 1:G) {
        T_g[[i]] <- diag(mean(diag(beg_T)), Q_g[i])
      }
    }




    #initial for D_g
    for (g in 1:G) {
      fit <- eigen(sig_g[[g]])
      s <- matrix(fit$vectors[, 1:Q_g[g]], ncol = Q_g[g], nrow = K)
      l <- diag(fit$values[1:Q_g[g]], Q_g[g])
      D_g[[g]] <- diag(diag(sig_g[[g]] - s %*% l %*% t(s)))
    }

    int_D_g <- D_g
    if (code_D == "G") {
      beg_D <- 0
      for (i in 1:G) {
        beg_D <- int_D_g[[i]] * pi_g[i] + beg_D
      }
      for (i in 1:G) {
        D_g[[i]] <- beg_D
      }
    }
    else if (code_D == "D") {
      for (i in 1:G) {
        D_g[[i]] <- diag(mean(diag(int_D_g[[i]])), K)
      }
    }
    else if (code_D == "C") {
      beg_D <- 0
      for (i in 1:G) {
        beg_D <- int_D_g[[i]] * pi_g[i] + beg_D
      }
      for (i in 1:G) {
        D_g[[i]] <- diag(mean(diag(beg_D)), K)
      }
    }






    #inital for V
    for (g in 1:G) {
      vg <- sig_g[[g]]
      for (i in 1:n) {
        V[[g]][[i]] <- diag(0.1, K)
      }
    }



    return(
      list(
        new_pi_g = pi_g,
        new_T_g = T_g,
        new_B_g = B_g,
        new_D_g = D_g,
        new_m = m,
        new_V = V,
        new_mu_g = mu_g,
        new_sig_g = sig_g,
        new_beta_g=beta_g
      )
    )

    }


}
