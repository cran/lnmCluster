#' Model selections for \code{lnmbicluster}
#'
#' fit several models for lnmbicluster along with 3 criteria values: AIC BIC and ICL
#'
#'@param W_count The microbiome count matrix that you want to analyze.
#'@param range_G All possible number of component groups, a vector.
#'@param range_Q All possible number of bicluster groups Q, a vector.
#'@param model A vector of string that contain cov_str you want to select. Default is all 16 models.
#'@param permutation Only has effect when model contains UUU, UUG, UUD or UUC. If TRUE,
#'it assume the number of biclusters could be different for different components. If FALSE,
#'it assume the number of biclusters are the same cross all components.
#'@param iter Max iterations, defaul is 150.
#'@param const Constant permutation term in multinomial distribution.
#'@param X The regression covariates matrix, which generates from model.matrix.
#'
#'
#'
#'@return A dataframe that contain the cov_str, K, Q, AIC, BIC, ICL values for model. There may
#'be a lot rows if large K and Q, because of lots of combinations: it is a sum of a geometric
#'series with multiplier max(Q) from 1 to max(K).
#'
#'
#'
model_selection <- function(W_count, range_G, range_Q, model,permutation,iter, const,X) {

  df <- W_count
  res <- NULL
  j <- 0

  #registerDoParallel(2)

  ##UGU##
  cov_str <- "UGU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }

  ##UGD##
  cov_str <- "UGD"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }

  ##UGG##
  cov_str <- "UGG"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##UGC##
  cov_str <- "UGC"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GUU##
  cov_str <- "GUU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)
      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  #GUD##
  cov_str <- "GUD"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##GUG##
  cov_str <- "GUG"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GUC##
  cov_str <- "GUC"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GGU##
  cov_str <- "GGU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GGD##
  cov_str <- "GGD"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##GGG##
  cov_str <- "GGG"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##GGC##
  cov_str <- "GGC"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_G) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter, const,X)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }




  ##UUU##
  cov_str <- "UUU"
  if (cov_str %in% model) {
    if(permutation==TRUE){
    model_select <- list()
    h <- 1
    for (K in range_G) {
      comb <-
        gtools::permutations(
          n = length(range_Q),
          r = K,
          v = range_Q,
          repeats.allowed = T
        )
      model_select[[h]] <-
        foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
          Q_K <- comb[j, ]
          paracom(df, K, Q_K, cov_str,iter,const,X)

        }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_G) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter,const,X)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }




  ##UUD##
  cov_str <- "UUD"
  if (cov_str %in% model) {
    if(permutation==TRUE){
      model_select <- list()
      h <- 1
      for (K in range_G) {
        comb <-
          gtools::permutations(
            n = length(range_Q),
            r = K,
            v = range_Q,
            repeats.allowed = T
          )
        model_select[[h]] <-
          foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
            Q_K <- comb[j, ]
            paracom(df, K, Q_K, cov_str,iter,const,X)

          }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_G) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter,const,X)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }



  ##UUG##
  cov_str <- "UUG"
  if (cov_str %in% model) {
    if(permutation==TRUE){
      model_select <- list()
      h <- 1
      for (K in range_G) {
        comb <-
          gtools::permutations(
            n = length(range_Q),
            r = K,
            v = range_Q,
            repeats.allowed = T
          )
        model_select[[h]] <-
          foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
            Q_K <- comb[j, ]
            paracom(df, K, Q_K, cov_str,iter,const,X)

          }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_G) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter, const,X)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }


  ##UUC##
  cov_str <- "UUC"
  if (cov_str %in% model) {
    if(permutation==TRUE){
      model_select <- list()
      h <- 1
      for (K in range_G) {
        comb <-
          gtools::permutations(
            n = length(range_Q),
            r = K,
            v = range_Q,
            repeats.allowed = T
          )
        model_select[[h]] <-
          foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
            Q_K <- comb[j, ]
            paracom(df, K, Q_K, cov_str,iter, const,X)

          }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_G) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter, const,X)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }

  res
}


paracom <- function(df, K, Q_K, cov_str,iter, const,X) {
  initial_guess <- try(initial_variational_gaussian(df, K, Q_K, cov_str,X))
  if (is.character(initial_guess) |
      isTRUE(class(initial_guess) == "try-error")) {
    res <- initial_guess
  }
  else{
    V <- initial_guess$new_V
    m <- initial_guess$new_m
    pi_g <- initial_guess$new_pi_g
    mu_g <- initial_guess$new_mu_g
    sig_g <- initial_guess$new_sig_g
    B_K <- initial_guess$new_B_g
    T_K <- initial_guess$new_T_g
    D_K <- initial_guess$new_D_g
    beta_g<-initial_guess$new_beta_g
    res <-
      try(Mico_bi_jensens(df, K, Q_K, pi_g, mu_g, sig_g, V, m, B_K, T_K, D_K, cov_str,
                          iter, const,beta_g,X))
  }
  Q_hat <- Q_K
  if (isTRUE(class(res) == "try-error") | is.character(res)) {
    f <-
      data.frame(
        name = paste(c(cov_str, " G", K, " Q", Q_hat), collapse = ""),
        AIC = -Inf,
        BIC = -Inf,
        ICL = -Inf
      )
  }
  else {
    f <-
      data.frame(
        name = paste(c(cov_str, " G", K, " Q", Q_hat), collapse = ""),
        AIC = res$AIC,
        BIC = res$BIC,
        ICL = res$ICL
      )
  }
  f
}
