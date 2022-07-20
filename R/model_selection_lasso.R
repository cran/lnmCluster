#' Model selections for \code{plnmfa}
#'
#' fit several models for plnmfa along with 3 criteria values: AIC BIC and ICL
#'
#'@param W_count The microbiome count matrix that you want to analyze.
#'@param K A specific number of component
#'@param Q_K A specific number of latent dimension.
#'@param model A specific model name, UUU or GUU
#'@param range_tuning A range of tuning parameters specified, ranged from 0-1.
#'@param iter Max iterations, defaul is 150.
#'@param const Constant permutation term in multinomial distribution.
#'@param X The regression covariates matrix, which generates from model.matrix.
#'
#'
#'@return A dataframe that contain the cov_str, K, Q, AIC, BIC, ICL values for model. There may
#'be a lot rows if long range of tuning parameters.


#for now, it is not suggested to compare GUU and UUU, so model selection is only for number of G
#in each model.
model_selection_lasso <- function(W_count, K, Q_K, model, range_tuning, iter, const,X){
  df <- W_count
  res <- NULL
  j <- 0

  #UUU
  cov_str<-"UUU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (t in range_tuning) {
      model_select[[h]] <- foreach(j = range_tuning, .combine = rbind) %do%  {
        tuning <- j
        paracom_lasso(df, K, Q_K, cov_str,tuning, iter, const,X)
      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }

  #GUU
  cov_str<-"GUU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (t in range_tuning) {
      model_select[[h]] <- foreach(j = range_tuning, .combine = rbind) %do%  {
        tuning <- j
        paracom_lasso(df, K, Q_K, cov_str,tuning, iter, const,X)
      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }

res
}




paracom_lasso <- function(df, K, Q_K, cov_str,tuning,iter, const,X) {
  initial_guess <- try(initial_variational_lasso(df, K, Q_K, cov_str,X))
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
    tuning=rep(tuning,K)
    res <-
      try(Mico_bi_lasso(df, K, Q_K, pi_g, mu_g, sig_g, V, m, B_K, T_K, D_K, cov_str,tuning, iter,
                        const,beta_g,X))
  }
  Q_hat <- Q_K
  if (isTRUE(class(res) == "try-error") | is.character(res)) {
    f <-
      data.frame(
        name = paste(c(cov_str, " G", K, " Q", Q_hat, " T", tuning[1]), collapse = ""),
        AIC = -Inf,
        BIC = -Inf,
        ICL = -Inf
      )
  }
  else {
    f <-
      data.frame(
        name = paste(c(cov_str, " G", K, " Q", Q_hat, " T", tuning[1]), collapse = ""),
        AIC = res$AIC,
        BIC = res$BIC,
        ICL = res$ICL
      )
  }
  f
}
