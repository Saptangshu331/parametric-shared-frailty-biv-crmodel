#' Observed Log-Likelihood for Bivariate Cure Rate Model
#'
#' Computes the observed log-likelihood of the bivariate cure rate model
#' for given parameter values and observed data.
#'
#' @name Observed_ll_cpp
#' @param gamma1 Parameter vector \eqn{\gamma_1}
#' @param gamma2 Parameter vector \eqn{\gamma_2}
#' @param rho Parameter \eqn{\rho}
#' @param beta1 Parameter vector \eqn{\beta_1}
#' @param beta2 Parameter vector \eqn{\beta_2}
#' @param alpha1 Parameter \eqn{\alpha_1}
#' @param alpha2 Parameter \eqn{\alpha_2}
#' @param q Parameter \eqn{q}
#' @param sigma Parameter \eqn{\sigma}
#' @param t1 Observed lifetimes for event 1
#' @param t2 Observed lifetimes for event 2
#' @param X Covariate matrix
#' @param censor_type_vec Censoring pattern indicator
#'     \describe{
#'       \item{1}{Both events observed}
#'       \item{2}{Event 1 observed, Event 2 censored}
#'       \item{3}{Event 1 censored, Event 2 observed}
#'       \item{4}{Both events censored}
#'     }
#'
#' @return Numeric value representing the observed log-likelihood
#'
#' @details
#' Refer to paper for observed log-likelihood form for the biavraite cure rate model, and further details.
#'
#'
#'
#' @export
NULL
