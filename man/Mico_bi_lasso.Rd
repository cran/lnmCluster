% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Micro_bi_lasso.R
\name{Mico_bi_lasso}
\alias{Mico_bi_lasso}
\title{Penalized Logistic Normal Multinomial factor analyzer main estimation process}
\usage{
Mico_bi_lasso(
  W_count,
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
  tuning,
  iter,
  const,
  beta_g,
  X
)
}
\arguments{
\item{W_count}{The microbiome count matrix}

\item{G}{All possible number of components. A vector.}

\item{Q_g}{A specific number of latent dimension.}

\item{pi_g}{A vector of initial guesses of component proportion}

\item{mu_g}{A list of initial guess of mean vector}

\item{sig_g}{A list of initial guess of covariance matrix for each component}

\item{V}{A list of initial guess of variational varaince}

\item{m}{A list of initial guess of variational mean}

\item{B_K}{A list of initial guess of loading matrix.}

\item{T_K}{A list of identity matrix with dimension q.}

\item{D_K}{A list of initial guess of error matrix}

\item{cov_str}{The covaraince structure you choose, there are 2 different models belongs to
this family:UUU and GUU. You can choose more than 1 covarance structure to do model selection.}

\item{tuning}{length G vector with range 0-1, define the tuning parameter for each component}

\item{iter}{Max iterations, default is 150.}

\item{const}{the permutation constant in multinomial distribution. Calculated before the main
algorithm in order to save computation time.}

\item{beta_g}{initial guess of covariates coefficients.}

\item{X}{The regression covariates matrix, which generates by model.matrix.}
}
\value{
z_ig Estimated latent variable z

cluster Component labels

mu_g Estimated component mean

pi_g Estimated component proportion

B_g Estimated sparsity loading matrix

D_g Estimated error covariance

COV Estimated component covariance

beta_g Estimated covariates coefficients.

overall_loglik Complete log likelihood value for each iteration

ICL ICL value

BIC BIC value

AIC AIC value

tuning display the tuning parameter you specified.
}
\description{
Main function will perform PLNM factor analyzer and return parameters
}
