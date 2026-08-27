#' Simulate Data for Bivariate Cure Rate Model
#'
#' Generates simulated data based on model parameters.
#'
#' @param x_sim Input covariate vector
#' @param gamma1 Parameter vector \eqn{\gamma_1}
#' @param gamma2 Parameter vector \eqn{\gamma_2}
#' @param rho Correlation parameter \eqn{\rho}
#' @param beta1 Parameter vector  \eqn{\beta_1}
#' @param beta2 Parameter vector  \eqn{\beta_2}
#' @param alpha1 Parameter \eqn{\alpha_1}
#' @param alpha2 Parameter \eqn{\alpha_2}
#' @param q Parameter \eqn{q}
#' @param sigma Parameter \eqn{\sigma}
#' @param lambda Parameter \eqn{\lambda}
#' @param lambda_c Parameter \eqn{\kappa}
#'
#' @details
#' Refer to the paper, model assumed for more details on these parameters
#'
#' \eqn{\gamma_1, \gamma_2, \rho} is associated with probabilties of cure combinations
#'
#'  \eqn{\beta_1, \beta_2} are associated with scale parameter   and \eqn{\alpha_1, \alpha_2} are shape parameter in Weibull distribution for conditional lifetime of susceptible individuals
#'
#'  \eqn{q, \sigma, \lambda} in generalized gamma distribution
#'
#'   lambda_c - \eqn{\kappa}. decides censoring proportion, rate parameter in exponential distribution assumed for censoring time distribution
#'
#'
#' @return simulated dataset as a dataframe containing all relevant columns as in observed data
#' \describe{
#'   \item{t1}{Observed lifetime for event 1}
#'   \item{t2}{Observed lifetimes for event 2}
#'   \item{obsvd_status_t1}{Censoring indicators for event 1 (1 = event, 0 = censored)}
#'   \item{obsvd_status_t2}{Censoring indicators for event 2 (1 = event, 0 = censored)}
#'   \item{cure_status_t1}{Cure indicators for event 1 (1 = cured, 0 = susceptible)}
#'   \item{cure_status_t2}{Cure indicators for event 2 (1 = cured, 0 = susceptible)}
#'   \item{x}{Covariate vector}
#'   \item{censor_ind}{Censoring pattern indicators combination of both censoring indicators:
#'     \describe{
#'       \item{1}{Both events observed}
#'       \item{2}{Event 1 observed, Event 2 censored}
#'       \item{3}{Event 1 censored, Event 2 observed}
#'       \item{4}{Both events censored}
#'     }
#'   }
#' }
#'
#' @export
simulate_data <- function(x_sim, gamma1, gamma2, rho,beta1,beta2,alpha1,alpha2,q,sigma,lambda,lambda_c) {


  n_sim=length(x_sim)

  #Data frame for simulated data
  sim_data <- data.frame(
    t1 = numeric(n_sim),
    t2 = numeric(n_sim),
    obsvd_status_t1 = numeric(n_sim),
    obsvd_status_t2 = numeric(n_sim),
    cure_status_t1 = numeric(n_sim),
    cure_status_t2 = numeric(n_sim),
    x = numeric(n_sim),
    stringsAsFactors = FALSE
  )

  sim_data[,]=NA

  sim_data$x=x_sim

  X_sim=cbind(1,x_sim)


  errors=matrix(,n_sim,2)

  for (i in c(1:n_sim)) {
    errors[i,]=mvrnorm(1, mu = c(0, 0), Sigma = matrix(c(1, rho,rho, 1), 2, 2))
  }
  #errors <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, mean(rho), mean(rho), 1), 2, 2))
  eps1 <- errors[, 1]
  eps2 <- errors[, 2]

  # Latent and observed outcomes
  y1_latent <-X_sim %*% gamma1 - eps1
  y2_latent <- X_sim %*% gamma2 - eps2

  y1 <- ifelse(y1_latent <0, 1, 0)
  y2 <- ifelse(y2_latent < 0, 1, 0)

  sim_data$cure_status_t1=as.vector(y1)
  sim_data$cure_status_t2=as.vector(y2)

  t1_given_y=t2_given_y=cens_time=rep(NA,n_sim)

  for (i in c(1:n_sim)) {

    if (sim_data$cure_status_t1[i] == 0 && sim_data$cure_status_t2[i] == 0) {


      if (q != 0) {

        v = rgamma(1, shape = q^-2, rate = q^-2)
        y = (v^(sigma / q)) / lambda
      } else {

        v = rlnorm(1,0,1)
        y = (v^sigma)/lambda
      }

      t1_given_y[i] = rweibull(1, shape = alpha1, scale = y * as.numeric(exp(X_sim[i, ] %*% beta1)))
      t2_given_y[i] = rweibull(1, shape = alpha2, scale = y * as.numeric(exp(X_sim[i, ] %*% beta2)))
      cens_time[i] = rexp(1, lambda_c)

      sim_data[i, "t1"] = ifelse(t1_given_y[i] <= cens_time[i], t1_given_y[i], cens_time[i])
      sim_data[i, "t2"] = ifelse(t2_given_y[i] <= cens_time[i], t2_given_y[i], cens_time[i])
      sim_data[i, "obsvd_status_t1"] = ifelse(t1_given_y[i] <= cens_time[i], 1, 0)
      sim_data[i, "obsvd_status_t2"] = ifelse(t2_given_y[i] <= cens_time[i], 1, 0)

    } else if (sim_data$cure_status_t1[i] == 1 && sim_data$cure_status_t2[i] == 0) {

      if (q != 0){

        v = rgamma(1, shape = q^-2, rate = q^-2)
        y = (v^(sigma / q)) / lambda
      } else {

        v = rlnorm(1,0,1)
        y = (v^sigma)/lambda
      }

      t2_given_y[i] = rweibull(1, shape = alpha2, scale = y * as.numeric(exp(X_sim[i, ] %*% beta2)))
      cens_time[i] = rexp(1, lambda_c)

      sim_data[i, "t1"] = cens_time[i]
      sim_data[i, "t2"] = ifelse(t2_given_y[i] <= cens_time[i], t2_given_y[i], cens_time[i])
      sim_data[i, "obsvd_status_t1"] = 0
      sim_data[i, "obsvd_status_t2"] = ifelse(t2_given_y[i] <= cens_time[i], 1, 0)

    } else if (sim_data$cure_status_t1[i] == 0 && sim_data$cure_status_t2[i] == 1) {

      if (q != 0) {

        v = rgamma(1, shape = q^-2, rate = q^-2)
        y = (v^(sigma / q)) / lambda
      } else {

        v = rlnorm(1,0,1)
        y = (v^sigma)/lambda
      }

      t1_given_y[i] = rweibull(1, shape = alpha1, scale = y * as.numeric(exp(X_sim[i, ] %*% beta1)))
      cens_time[i] = rexp(1, lambda_c)

      sim_data[i, "t1"] = ifelse(t1_given_y[i] <= cens_time[i], t1_given_y[i], cens_time[i])
      sim_data[i, "t2"] = cens_time[i]
      sim_data[i, "obsvd_status_t1"] = ifelse(t1_given_y[i] <= cens_time[i], 1, 0)
      sim_data[i, "obsvd_status_t2"] = 0

    } else {

      cens_time[i] = rexp(1, lambda_c)

      sim_data[i, "t1"] = cens_time[i]
      sim_data[i, "t2"] = cens_time[i]
      sim_data[i, "obsvd_status_t1"] = 0
      sim_data[i, "obsvd_status_t2"] = 0
    }
  }

  #attach(sim_data)

  censor_ind=rep(NA, nrow(sim_data))

  for (i in 1:nrow(sim_data)) {
    if (sim_data$obsvd_status_t1[i] == 1 & sim_data$obsvd_status_t2[i] == 1) {
      censor_ind[i] <- 1
    } else if (sim_data$obsvd_status_t1[i] == 1 & sim_data$obsvd_status_t2[i] == 0) {
      censor_ind[i] <- 2
    } else if (sim_data$obsvd_status_t1[i] == 0 & sim_data$obsvd_status_t2[i] == 1) {
      censor_ind[i] <- 3
    } else {
      censor_ind[i] <- 4
    }
  }

  sim_data$censor_ind=censor_ind

  return(sim_data)

}

#' EM Algorithm for Bivariate Cure Rate Model
#'
#' Estimates model parameters using the Expectation-Maximization (EM) algorithm
#' for the bivariate cure rate model. Refer to the paper for further details
#'
#' @param gamma10 Initial value of parameter vector \eqn{\gamma_1}
#' @param gamma20 Initial value of parameter vector \eqn{\gamma_2}
#' @param rho0 Initial value of parameter \eqn{\rho}
#' @param beta1_init Initial value of parameter vector \eqn{\beta_1}
#' @param beta2_init Initial value of parameter vector \eqn{\beta_2}
#' @param alpha1_init Initial value of parameter \eqn{\alpha_1}
#' @param alpha2_init Initial value of parameter \eqn{\alpha_2}
#' @param q0 Initial value of parameter \eqn{q}
#' @param sigma0 Initial value of parameter \eqn{\sigma}
#' @param t1 Observed lifetimes for event 1
#' @param t2 Observed lifetimes for event 2
#' @param X Covariate matrix
#' @param censor_type_vec Censoring pattern indicator:
#' \describe{
#'   \item{1}{Both events observed}
#'   \item{2}{Event 1 observed, Event 2 censored}
#'   \item{3}{Event 1 censored, Event 2 observed}
#'   \item{4}{Both events censored}
#' }
#' @param print_opt Logical; if TRUE, prints iteration details
#' @param epsilon Convergence tolerance
#' @param conv_obsll Logical; if TRUE, convergence is based on observed log-likelihood
#' @param method1 Optimization method in M-step for \eqn{l_{c_1}} (cure component)
#' @param method2 Optimization method in M-step for \eqn{l_{c_2}} (lifetime component)
#'
#' @details
#' The model assumes:
#'
#' \itemize{
#'   \item \eqn{\gamma_1, \gamma_2, \rho} are associated with probabilities of cure combinations
#'   \item \eqn{\beta_1, \beta_2} are associated with scale parameters of the Weibull distribution
#'   \item \eqn{\alpha_1, \alpha_2} are shape parameters of the Weibull distribution
#'   \item \eqn{q, \sigma} are parameters of the generalized gamma frailty distribution
#' }
#'
#' The EM algorithm iteratively performs:
#' \itemize{
#'   \item E-step: Computes conditional expectations
#'   \item M-step: Updates parameters using by optimizing \eqn{l_{c_1}} and \eqn{l_{c_2}} seperately wrt  to the respective relevant parameters
#' }
#'
#' Convergence is determined either by:
#' \itemize{
#'   \item Change in observed log-likelihood, or
#'   \item Relative change in parameter estimates
#' }
#'
#' @return A named numeric vector containing estimated parameters and its observed log-likelihood:
#' \describe{
#'   \item{gamma10, gamma11}{Estimated \eqn{\gamma_1}}
#'   \item{gamma20, gamma21}{Estimated \eqn{\gamma_2}}
#'   \item{rho}{Estimated \eqn{\rho}}
#'   \item{beta10, beta11}{Estimated \eqn{\beta_1}}
#'   \item{beta20, beta21}{Estimated \eqn{\beta_2}}
#'   \item{alpha1, alpha2}{Estimated  \eqn{\alpha_1}, \eqn{\alpha_1} }
#'   \item{q}{Estimated \eqn{q}}
#'   \item{sigma}{Estimated \eqn{sigma}}
#'   \item{Obs_LL}{Final observed log-likelihood}
#' }
#'
#'
#' @export
EM_func_qunfixed = function(gamma10, gamma20, rho0, beta1_init, beta2_init, alpha1_init, alpha2_init,
                            q0, sigma0, t1, t2, X, censor_type_vec, print_opt = FALSE, epsilon = .0001,conv_obsll=TRUE, method1="Nelder-Mead",method2="Nelder-Mead") {

  gamma1_k = gamma10
  gamma2_k = gamma20
  rho_k = rho0

  beta1_k = beta1_init
  beta2_k = beta2_init
  alpha1_k = alpha1_init
  alpha2_k = alpha2_init
  q_k = q0
  sigma_k = sigma0

  lg = length(gamma1_k)
  lb = length(beta1_k)

  n_iter = 1000
  for (iter in 1:n_iter) {

    phi_k = phi_matrix_cpp(X, gamma1_k, gamma2_k, rho_k)

    phi11_k = phi_k[, 1]
    phi10_k = phi_k[, 2]
    phi01_k = phi_k[, 3]
    phi00_k = phi_k[, 4]

    # E-step
    E_k_c = Expc_k_given_C_cpp(t1, t2, X, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k, censor_type_vec, phi11_k, phi10_k, phi01_k, phi00_k)

    # M-step
    M1_phi_res = optimize_M1_phi(X, E_k_c, censor_type_vec, c(gamma1_k, gamma2_k, rho_k), method1)

    M1_phi_param = M1_phi_res[["parameters"]]

    gamma1_k1 = M1_phi_param[1:lg]
    gamma2_k1 = M1_phi_param[(lg + 1):(2 * lg)]
    rho_k1 = M1_phi_param[2 * lg + 1]

    M2_LC2_res = optimize_M2_lc2(t1, t2, X, E_k_c, censor_type_vec, c(beta1_k, beta2_k, log(alpha1_k), log(alpha2_k), q_k, log(sigma_k)), method2)

    M2_LC2_param = M2_LC2_res[["parameters"]]

    beta1_k1 = M2_LC2_param[1:lb]
    beta2_k1 = M2_LC2_param[(lb + 1):(2 * lb)]
    alpha1_k1 = exp(M2_LC2_param[(2 * lb + 1)])
    alpha2_k1 = exp(M2_LC2_param[(2 * lb + 2)])
    q_k1 = M2_LC2_param[(2 * lb + 3)]
    sigma_k1 = exp(M2_LC2_param[(2 * lb + 4)])


   #  if (max_prd < epsilon) {
   #    print("EM Converges")
   #    break
   #  } else {
   # }
   #
    if (conv_obsll) {
      Obs_LLk = Observed_ll_cpp(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k, t1, t2, X, censor_type_vec)
      Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1, t1, t2, X, censor_type_vec)

      absll_diff = abs(Obs_LLk - Obs_LLk1)

      if (absll_diff < epsilon) {
        #print("EM Converges")
        break
      }
    } else {

      theta_k = c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k)
      theta_k1 = c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1)


      max_prd = max(abs((theta_k - theta_k1) / theta_k))

      if (max_prd < epsilon) {
        # print("EM Converges")
        Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1, t1, t2, X, censor_type_vec)

        break
      }
    }

    if (print_opt) {
      print(paste("E-step No.", iter))
      print(paste("M-step Phi No.", iter))
      print(M1_phi_res)

      print(paste("M-step LC2 No.", iter))
      print(M2_LC2_res)

      cat("gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k -\n",
          paste(c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k), collapse = ", "))
      print("______________________________________________________________________________________")
      cat("gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1 -\n",
          paste(c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1), collapse = ", "))
      cat("\n")


      if (conv_obsll) {


        print(paste("Obs_llk -", Obs_LLk))
        print(paste("Obs_llk1 -", Obs_LLk1))
        print(paste("Observed LL absolute difference", absll_diff))
        print(paste("Observed LL relative difference", absll_diff / abs(Obs_LLk)))

      } else {

        par_dist_sq = (theta_k - theta_k1)^2
        print(paste("Euclidean distance between parameters", sqrt(sum(par_dist_sq))))
        print(paste("Max absolute relative difference in parameters", max_prd))
        print(paste("Index of max absolute relative difference in parameters", which.max(abs(theta_k - theta_k1) / theta_k)))

      }

    }






    gamma1_k = gamma1_k1
    gamma2_k = gamma2_k1
    rho_k = rho_k1
    beta1_k = beta1_k1
    beta2_k = beta2_k1
    alpha1_k = alpha1_k1
    alpha2_k = alpha2_k1
    q_k = q_k1
    sigma_k = sigma_k1


     }


  return(c(gamma10 = gamma1_k1[1], gamma11 = gamma1_k1[2],
           gamma20 = gamma2_k1[1], gamma21 = gamma2_k1[2],
           rho = rho_k1,
           beta10 = beta1_k1[1], beta11 = beta1_k1[2],
           beta20 = beta2_k1[1], beta21 = beta2_k1[2],
           alpha1 = alpha1_k1, alpha2 = alpha2_k1,
           q = q_k1, sigma = sigma_k1,
           Obs_LL = Obs_LLk1
  ))

  # list(
  #   "gamma1" = gamma1_k1,
  #   "gamma2" = gamma2_k1,
  #   "rho" = rho_k1,
  #   "beta1" = beta1_k1,
  #   "beta2" = beta2_k1,
  #   "alpha1" = alpha1_k1,
  #   "alpha2" = alpha2_k1,
  #   "q" = q_k1,
  #   "sigma" = sigma_k1,
  #   "Obs LL" = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1, t1, t2, X, censor_type_vec)
  # )

}


#' EM Algorithm for Bivariate Cure Rate Model (Fixed q)
#'
#' This function is similar to \code{\link{EM_func_qunfixed}}, except that
#' the parameter \eqn{q} is assumed to be fixed and optimization is done
#' w.r.t rest of the parameters using EM Algorithm,
#' this function is mainly used for profile likelihood approach where it is run for grid of
#' values for  \eqn{q} and the MLE is obtained for that value of \eqn{q} at which the final
#'  observed log-likelihood is maximum, see paper for further details
#'
#' @details
#' Refer to \code{\link{EM_func_qunfixed}} for full model specification
#' and parameter interpretation.
#'
#' @seealso \code{\link{EM_func_qunfixed}}
#'
#' @export
EM_qfixed = function(gamma10, gamma20, rho0, beta10, beta20, alpha10, alpha20,
                     q, sigma0, t1, t2, X, censor_type_vec,
                     print_opt = FALSE, epsilon = .0001,conv_obsll=TRUE,method1="Nelder-Mead",method2="Nelder-Mead") {

  gamma1_k = gamma10
  gamma2_k = gamma20
  rho_k = rho0

  beta1_k = beta10
  beta2_k = beta20
  alpha1_k = alpha10
  alpha2_k = alpha20
  sigma_k = sigma0

  lg1=lg2=lb1=lb2=ncol(X)


  n_iter = 1000
  for (iter in 1:n_iter) {

    phi_k = phi_matrix_cpp(X, gamma1_k, gamma2_k, rho_k)

    phi11_k = phi_k[,1]
    phi10_k = phi_k[,2]
    phi01_k = phi_k[,3]
    phi00_k = phi_k[,4]

    # E-step
    E_k_c = Expc_k_given_C_cpp(t1, t2, X, beta1_k, beta2_k, alpha1_k, alpha2_k,
                               q, sigma_k, censor_type_vec, phi11_k, phi10_k, phi01_k, phi00_k)

    # M-step
    M1_phi_res = optimize_M1_phi(X, E_k_c, censor_type_vec, c(gamma1_k, gamma2_k, rho_k) , method1)
    M1_phi_param = M1_phi_res[["parameters"]]

    gamma1_k1 = M1_phi_param[1:lg1]
    gamma2_k1 = M1_phi_param[(lg1 + 1):(lg1 + lg2)]
    rho_k1 = M1_phi_param[(lg1 + lg2 + 1)]

    M2_LC2_res = optimize_M2_lc2_q_fixed(t1, t2, X, E_k_c, censor_type_vec,
                                         q, c(beta1_k, beta2_k, log(alpha1_k), log(alpha2_k), log(sigma_k)), method2)
    M2_LC2_param = M2_LC2_res[["parameters"]]

    beta1_k1 = M2_LC2_param[1:lb1]
    beta2_k1 = M2_LC2_param[(lb1 + 1):(lb1 + lb2)]
    alpha1_k1 = exp(M2_LC2_param[(lb1 + lb2 + 1)])
    alpha2_k1 = exp(M2_LC2_param[(lb1 + lb2 + 2)])
    sigma_k1 = exp(M2_LC2_param[(lb1 + lb2 + 3)])



    # if (max_prd < epsilon) {
    #   print("EM Converges")
    #   break
    # }

    if (conv_obsll) {
      Obs_LLk = Observed_ll_cpp(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q, sigma_k, t1, t2, X, censor_type_vec)
      Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1, t1, t2, X, censor_type_vec)

      absll_diff = abs(Obs_LLk - Obs_LLk1)

      if (absll_diff < epsilon) {
        #print("EM Converges")
        break
      }
    } else {

      theta_k = c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, sigma_k)
      theta_k1 = c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, sigma_k1)

      max_prd = max(abs((theta_k - theta_k1) / theta_k))


      if (max_prd < epsilon) {
       # print("EM Converges")

        Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1, t1, t2, X, censor_type_vec)

        break
      }
    }

    if (print_opt) {
      print(paste("E-step No.", iter))
      print(paste("M-step Phi No.", iter))
      print(M1_phi_res)

      print(paste("M-step LC2 No.", iter))
      print(M2_LC2_res)

      cat("gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q, sigma_k -\n",
          paste(c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q, sigma_k), collapse = ", "))
      print("______________________________________________________________________________________")
      cat("gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1 -\n",
          paste(c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1), collapse = ", "))
      cat("\n")


      if (conv_obsll) {


      print(paste("Obs_llk -", Obs_LLk))
      print(paste("Obs_llk1 -", Obs_LLk1))
      print(paste("Observed LL absolute difference", absll_diff))
      print(paste("Observed LL relative difference", absll_diff / abs(Obs_LLk)))

      } else {

      par_dist_sq = (theta_k - theta_k1)^2
      print(paste("Euclidean distance between parameters", sqrt(sum(par_dist_sq))))
      print(paste("Max absolute relative difference in parameters", max_prd))
      print(paste("Index of max absolute relative difference in parameters", which.max(abs(theta_k - theta_k1) / theta_k)))

      }

    }



    # Update parameters if convergence criteria is not met
    gamma1_k = gamma1_k1
    gamma2_k = gamma2_k1
    rho_k = rho_k1
    beta1_k = beta1_k1
    beta2_k = beta2_k1
    alpha1_k = alpha1_k1
    alpha2_k = alpha2_k1
    sigma_k = sigma_k1






  }

  return(c(gamma10 = gamma1_k1[1], gamma11 = gamma1_k1[2],
    gamma20 = gamma2_k1[1], gamma21 = gamma2_k1[2],
    rho = rho_k1,
    beta10 = beta1_k1[1], beta11 = beta1_k1[2],
    beta20 = beta2_k1[1], beta21 = beta2_k1[2],
    alpha1 = alpha1_k1, alpha2 = alpha2_k1,
    q = q, sigma = sigma_k1,
    Obs_LL = Obs_LLk1
  ))


  # list("gamma1" = gamma1_k1,"gamma2" = gamma2_k1, "rho" = rho_k1, "beta1" = beta1_k1,"beta2" = beta2_k1,
  #   "alpha1" = alpha1_k1,
  #   "alpha2" = alpha2_k1,
  #   "q" = q,
  #   "sigma" = sigma_k1,
  #   "Obs LL" = Obs_LLk1
  # )
}


#' EM Algorithm for Bivariate Cure Rate Model (Gamma distribution as a special case of generalized gamma in the shared frailty)
#'
#' This function executes the EM Algorithm as given in the paper, assuming  gamma  distribution
#' as the special case of generalized gamma distribution, i.e \eqn{q=\sigma}
#' @details
#' Refer to \code{\link{EM_func_qunfixed}} for full model specification  and parameter interpretation.
#' Note that the user should give the initial value of q is not needed and it is assumed to be equal to sigma0 as \eqn{q=\sigma}
#'
#' @seealso \code{\link{EM_func_qunfixed}} and paper for further details
#'
#' @export
EM_func_gamma = function(gamma10, gamma20, rho0, beta10, beta20, alpha10, alpha20,
                          sigma0, t1, t2, X, censor_type_vec,
                         print_opt = FALSE, epsilon = .0001,conv_obsll=TRUE, method1="Nelder-Mead",method2="Nelder-Mead") {

  gamma1_k = gamma10
  gamma2_k = gamma20
  rho_k = rho0

  beta1_k = beta10
  beta2_k = beta20
  alpha1_k = alpha10
  alpha2_k = alpha20

    q_k = sigma0
    sigma_k = sigma0

  n_param = ncol(X)  # Single length assumption

  n_iter = 1000
  for (iter in 1:n_iter) {

    phi_k = phi_matrix_cpp(X, gamma1_k, gamma2_k, rho_k)

    phi11_k = phi_k[, 1]
    phi10_k = phi_k[, 2]
    phi01_k = phi_k[, 3]
    phi00_k = phi_k[, 4]

    # E-step
    E_k_c = Expc_k_given_C_cpp(t1, t2, X, beta1_k, beta2_k, alpha1_k, alpha2_k,
                           q_k, sigma_k, censor_type_vec, phi11_k, phi10_k, phi01_k, phi00_k)

    # M-step
    M1_phi_res = optimize_M1_phi(X, E_k_c, censor_type_vec, c(gamma1_k, gamma2_k, rho_k), method1)
    M1_phi_param = M1_phi_res[["parameters"]]

    gamma1_k1 = M1_phi_param[1:n_param]
    gamma2_k1 = M1_phi_param[(n_param + 1):(2 * n_param)]
    rho_k1 = M1_phi_param[(2 * n_param + 1)]

    M2_LC2_res = optimize_M2_lc2_gamma(t1, t2, X, E_k_c, censor_type_vec,
                                       c(beta1_k, beta2_k, log(alpha1_k), log(alpha2_k), log(sigma_k)), method2)
    M2_LC2_param = M2_LC2_res[["parameters"]]

    beta1_k1 = M2_LC2_param[1:n_param]
    beta2_k1 = M2_LC2_param[(n_param + 1):(2 * n_param)]
    alpha1_k1 = exp(M2_LC2_param[(2 * n_param + 1)])
    alpha2_k1 = exp(M2_LC2_param[(2 * n_param + 2)])
    sigma_k1 = exp(M2_LC2_param[(2 * n_param + 3)])

    q_k1 = sigma_k1  # Ensuring q updates correctly


    if (conv_obsll) {
      Obs_LLk = Observed_ll_cpp(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k, t1, t2, X, censor_type_vec)
      Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1, t1, t2, X, censor_type_vec)

      absll_diff = abs(Obs_LLk - Obs_LLk1)

      if (absll_diff < epsilon) {
        #print("EM Converges")
        break
      }
    } else {

      theta_k = c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k)
      theta_k1 = c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1)

      max_prd = max(abs((theta_k - theta_k1) / theta_k))


      if (max_prd < epsilon) {
        # print("EM Converges")

        Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1, t1, t2, X, censor_type_vec)

        break
      }
    }



    if (print_opt) {
      print(paste("E-step No.", iter))
      print(paste("M-step Phi No.", iter))
      print(M1_phi_res)

      print(paste("M-step LC2 No.", iter))
      print(M2_LC2_res)

      cat("gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k -\n",
          paste(c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k), collapse = ", "))
      print("______________________________________________________________________________________")
      cat("gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1 -\n",
          paste(c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1), collapse = ", "))
      cat("\n")


      if (conv_obsll) {


        print(paste("Obs_llk -", Obs_LLk))
        print(paste("Obs_llk1 -", Obs_LLk1))
        print(paste("Observed LL absolute difference", absll_diff))
        print(paste("Observed LL relative difference", absll_diff / abs(Obs_LLk)))

      } else {

        par_dist_sq = (theta_k - theta_k1)^2
        print(paste("Euclidean distance between parameters", sqrt(sum(par_dist_sq))))
        print(paste("Max absolute relative difference in parameters", max_prd))
        print(paste("Index of max absolute relative difference in parameters", which.max(abs(theta_k - theta_k1) / theta_k)))

      }

    }


    # Update parameters if convergence criteria is not met
    gamma1_k = gamma1_k1
    gamma2_k = gamma2_k1
    rho_k = rho_k1
    beta1_k = beta1_k1
    beta2_k = beta2_k1
    alpha1_k = alpha1_k1
    alpha2_k = alpha2_k1
    sigma_k = sigma_k1
    q_k=q_k1
  }

  return(c(gamma10 = gamma1_k1[1], gamma11 = gamma1_k1[2],
           gamma20 = gamma2_k1[1], gamma21 = gamma2_k1[2],
           rho = rho_k1,
           beta10 = beta1_k1[1], beta11 = beta1_k1[2],
           beta20 = beta2_k1[1], beta21 = beta2_k1[2],
           alpha1 = alpha1_k1, alpha2 = alpha2_k1,
           q = q_k1, sigma = sigma_k1,
           Obs_LL = Obs_LLk1
  ))

  # return(list(
  #   "gamma1" = gamma1_k1,
  #   "gamma2" = gamma2_k1,
  #   "rho" = rho_k1,
  #   "beta1" = beta1_k1,
  #   "beta2" = beta2_k1,
  #   "alpha1" = alpha1_k1,
  #   "alpha2" = alpha2_k1,
  #   "q" = q_k1,
  #   "sigma" = sigma_k1,
  #   "Obs LL" = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
  #                              alpha1_k1, alpha2_k1, q_k1, sigma_k1, t1, t2, X, censor_type_vec)
  # ))
}

#' EM Algorithm for Bivariate Cure Rate Model (Exponential  distribution as a special case of generalized gamma distributionin the shared frailty)
#'
#' This function executes the EM Algorithm as given in the paper, assuming  Exponential  distribution
#' as the special case of generalized gamma distribution, i.e \eqn{q=\sigma=1}
#' @details
#' Refer to \code{\link{EM_func_qunfixed}} for full model specification  and parameter interpretation.
#' Note that the intial values of q and \eqn{\sigma} are not needed and user should give the initial values for only the rest of parameters
#'
#' @seealso \code{\link{EM_func_qunfixed}} and paper for further details
#'
#' @export
EM_func_exp = function(gamma10, gamma20, rho0, beta10, beta20, alpha10, alpha20,
                       t1, t2, X, censor_type_vec,
                       print_opt = FALSE, epsilon = .0001,conv_obsll=TRUE) {

  # Initialize parameters
  gamma1_k = gamma10
  gamma2_k = gamma20
  rho_k = rho0
  beta1_k = beta10
  beta2_k = beta20
  alpha1_k = alpha10
  alpha2_k = alpha20

  # Set q and sigma to 1

  n_param = ncol(X)  # Unified parameter length

  n_iter = 1000
  for (iter in 1:n_iter) {

    # Compute phi matrix
    phi_k = phi_matrix_cpp(X, gamma1_k, gamma2_k, rho_k)

    phi11_k = phi_k[, 1]
    phi10_k = phi_k[, 2]
    phi01_k = phi_k[, 3]
    phi00_k = phi_k[, 4]

    # E-step
    E_k_c = Expc_k_given_C_cpp(t1, t2, X, beta1_k, beta2_k, alpha1_k, alpha2_k,
                               1, 1, censor_type_vec, phi11_k, phi10_k, phi01_k, phi00_k)

    # M-step
    M1_phi_res = optimize_M1_phi(X, E_k_c, censor_type_vec, c(gamma1_k, gamma2_k, rho_k))
    M1_phi_param = M1_phi_res[["parameters"]]

    gamma1_k1 = M1_phi_param[1:n_param]
    gamma2_k1 = M1_phi_param[(n_param + 1):(2 * n_param)]
    rho_k1 = M1_phi_param[(2 * n_param + 1)]

    # Optimize only with respect to beta1, beta2, alpha1, alpha2 (q and sigma are fixed at 1)
    M2_LC2_res = optimize_M2_lc2_exp(t1, t2, X, E_k_c, censor_type_vec,
                                     c(beta1_k, beta2_k, log(alpha1_k), log(alpha2_k)))
    M2_LC2_param = M2_LC2_res[["parameters"]]

    beta1_k1 = M2_LC2_param[1:n_param]
    beta2_k1 = M2_LC2_param[(n_param + 1):(2 * n_param)]
    alpha1_k1 = exp(M2_LC2_param[(2 * n_param + 1)])
    alpha2_k1 = exp(M2_LC2_param[(2 * n_param + 2)])


    # if (max_prd < epsilon) {
    #   print("EM Converges")
    #   break
    # } else {
    # }

    if (conv_obsll) {
      Obs_LLk = Observed_ll_cpp(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, 1,1, t1, t2, X, censor_type_vec)
      Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, 1,1, t1, t2, X, censor_type_vec)

      absll_diff = abs(Obs_LLk - Obs_LLk1)

      if (absll_diff < epsilon) {
        #print("EM Converges")
        break
      }
    } else {

      theta_k = c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k)
      theta_k1 = c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1)

      max_prd = max(abs((theta_k - theta_k1) / theta_k))


      if (max_prd < epsilon) {
        # print("EM Converges")

        Obs_LLk1 = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, 1,1, t1, t2, X, censor_type_vec)

        break
      }
    }



    if (print_opt) {
      print(paste("E-step No.", iter))
      print(paste("M-step Phi No.", iter))
      print(M1_phi_res)

      print(paste("M-step LC2 No.", iter))
      print(M2_LC2_res)

      cat("gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k -\n",
          paste(c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k), collapse = ", "))
      print("______________________________________________________________________________________")
      cat("gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1 -\n",
          paste(c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1), collapse = ", "))
      cat("\n")


       if (conv_obsll) {


        print(paste("Obs_llk -", Obs_LLk))
        print(paste("Obs_llk1 -", Obs_LLk1))
        print(paste("Observed LL absolute difference", absll_diff))
        print(paste("Observed LL relative difference", absll_diff / abs(Obs_LLk)))

      } else {

        par_dist_sq = (theta_k - theta_k1)^2
        print(paste("Euclidean distance between parameters", sqrt(sum(par_dist_sq))))
        print(paste("Max absolute relative difference in parameters", max_prd))
        print(paste("Index of max absolute relative difference in parameters", which.max(abs(theta_k - theta_k1) / theta_k)))

      }


    }



    # Update parameters if convergence criteria is not met
    gamma1_k = gamma1_k1
    gamma2_k = gamma2_k1
    rho_k = rho_k1
    beta1_k = beta1_k1
    beta2_k = beta2_k1
    alpha1_k = alpha1_k1
    alpha2_k = alpha2_k1





  }

  return(c(gamma10 = gamma1_k1[1], gamma11 = gamma1_k1[2],
           gamma20 = gamma2_k1[1], gamma21 = gamma2_k1[2],
           rho = rho_k1,
           beta10 = beta1_k1[1], beta11 = beta1_k1[2],
           beta20 = beta2_k1[1], beta21 = beta2_k1[2],
           alpha1 = alpha1_k1, alpha2 = alpha2_k1,
           q = 1, sigma = 1,
           Obs_LL = Obs_LLk1
  ))

  # return(list(
  #   "gamma1" = gamma1_k1,
  #   "gamma2" = gamma2_k1,
  #   "rho" = rho_k1,
  #   "beta1" = beta1_k1,
  #   "beta2" = beta2_k1,
  #   "alpha1" = alpha1_k1,
  #   "alpha2" = alpha2_k1,
  #   "q" = 1,  # Fixed at 1
  #   "sigma" = 1,  # Fixed at 1
  #   "Obs LL" = Observed_ll_cpp(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
  #                              alpha1_k1, alpha2_k1, 1,1, t1, t2, X, censor_type_vec)
  # ))

  }


#' EM Algorithm for Bivariate Cure Rate Model (Fixed q) using monte carlo summation approximation for integrals
#'
#' This function performs the same thing as \code{\link{EM_qfixed}}, only difference being that all integrals are approximated by monte carlo summation approach, refer to the paper for more details.
#'
#'@param n_mc  sample size of the monte carlo sample used for approximaton of integrals involved using monte carlo summation
#'
#' @details
#' Refer to \code{\link{EM_qfixed}} for full model specification
#' and parameter interpretation.
#'
#' @seealso \code{\link{EM_func_qunfixed}}
#'
#' @export
EM_qfixed_mc = function(gamma10, gamma20, rho0, beta10, beta20, alpha10, alpha20,
                        q, sigma0, t1, t2, X, censor_type_vec,
                        print_opt = FALSE, epsilon = .0001,conv_obsll=TRUE,n_mc=10^3) {

  gamma1_k = gamma10
  gamma2_k = gamma20
  rho_k = rho0

  beta1_k = beta10
  beta2_k = beta20
  alpha1_k = alpha10
  alpha2_k = alpha20
  sigma_k = sigma0

  lg1=lg2=lb1=lb2=ncol(X)



  mcs=if(q!=0){rgamma(n_mc,shape=q^-2,rate=q^-2)
  }else{rlnorm(n_mc)}


  n_iter = 1000
  for (iter in 1:n_iter) {

    phi_k = phi_matrix_cpp(X, gamma1_k, gamma2_k, rho_k)

    phi11_k = phi_k[,1]
    phi10_k = phi_k[,2]
    phi01_k = phi_k[,3]
    phi00_k = phi_k[,4]

    # E-step
    E_k_c = Expc_k_given_C_mc(t1, t2, X, beta1_k, beta2_k, alpha1_k, alpha2_k,
                              q, sigma_k, censor_type_vec, phi11_k, phi10_k, phi01_k, phi00_k,mcs)

    # M-step
    M1_phi_res = optimize_M1_phi(X, E_k_c, censor_type_vec, c(gamma1_k, gamma2_k, rho_k))
    M1_phi_param = M1_phi_res[["parameters"]]

    gamma1_k1 = M1_phi_param[1:lg1]
    gamma2_k1 = M1_phi_param[(lg1 + 1):(lg1 + lg2)]
    rho_k1 = M1_phi_param[(lg1 + lg2 + 1)]

    M2_LC2_res = optimize_M2_lc2_q_fixed_mc(t1, t2, X, E_k_c, censor_type_vec,
                                            q,mcs, c(beta1_k, beta2_k, log(alpha1_k), log(alpha2_k), log(sigma_k)))
    M2_LC2_param = M2_LC2_res[["parameters"]]

    beta1_k1 = M2_LC2_param[1:lb1]
    beta2_k1 = M2_LC2_param[(lb1 + 1):(lb1 + lb2)]
    alpha1_k1 = exp(M2_LC2_param[(lb1 + lb2 + 1)])
    alpha2_k1 = exp(M2_LC2_param[(lb1 + lb2 + 2)])
    sigma_k1 = exp(M2_LC2_param[(lb1 + lb2 + 3)])


    # if (max_prd < epsilon) {
    #   print("EM Converges")
    #   break
    # }

    if (conv_obsll) {
      Obs_LLk = Observed_ll_mc(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q, sigma_k,t1, t2, X, censor_type_vec, mcs)
      Obs_LLk1 = Observed_ll_mc(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1, t1, t2, X, censor_type_vec,mcs)

      absll_diff = abs(Obs_LLk - Obs_LLk1)

      if (absll_diff < epsilon) {
        #print("EM Converges")
        break
      }
    } else {

      theta_k = c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, sigma_k)
      theta_k1 = c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, sigma_k1)

      max_prd = max(abs((theta_k - theta_k1) / theta_k))


      if (max_prd < epsilon) {
        # print("EM Converges")

        Obs_LLk1 = Observed_ll_mc(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1, t1, t2, X, censor_type_vec,mcs)

        break
      }
    }


    if (print_opt) {
      print(paste("E-step No.", iter))
      print(paste("M-step Phi No.", iter))
      print(M1_phi_res)

      print(paste("M-step LC2 No.", iter))
      print(M2_LC2_res)

      cat("gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q, sigma_k -\n",
          paste(c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q, sigma_k), collapse = ", "))
      print("______________________________________________________________________________________")
      cat("gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1 -\n",
          paste(c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q, sigma_k1), collapse = ", "))
      cat("\n")

      if (conv_obsll) {


        print(paste("Obs_llk -", Obs_LLk))
        print(paste("Obs_llk1 -", Obs_LLk1))
        print(paste("Observed LL absolute difference", absll_diff))
        print(paste("Observed LL relative difference", absll_diff / abs(Obs_LLk)))

      } else {

        par_dist_sq = (theta_k - theta_k1)^2
        print(paste("Euclidean distance between parameters", sqrt(sum(par_dist_sq))))
        print(paste("Max absolute relative difference in parameters", max_prd))
        print(paste("Index of max absolute relative difference in parameters", which.max(abs(theta_k - theta_k1) / theta_k)))

      }


    }






    # Update parameters if convergence criteria is not met
    gamma1_k = gamma1_k1
    gamma2_k = gamma2_k1
    rho_k = rho_k1
    beta1_k = beta1_k1
    beta2_k = beta2_k1
    alpha1_k = alpha1_k1
    alpha2_k = alpha2_k1
    sigma_k = sigma_k1






  }

  return(c(gamma10 = gamma1_k1[1], gamma11 = gamma1_k1[2],
           gamma20 = gamma2_k1[1], gamma21 = gamma2_k1[2],
           rho = rho_k1,
           beta10 = beta1_k1[1], beta11 = beta1_k1[2],
           beta20 = beta2_k1[1], beta21 = beta2_k1[2],
           alpha1 = alpha1_k1, alpha2 = alpha2_k1,
           q = q, sigma = sigma_k1,
           Obs_LL = Obs_LLk1
  ))


}



#' EM Algorithm for Bivariate Cure Rate Model (Gamma distribution as a special case of generalized gamma in the shared frailty with Monte Carlo approximation)
#'
#' This function performs the same procedure as \code{\link{EM_func_gamma}}, with the
#' difference that the integrals involved are approximated using
#' the Monte Carlo summation approach.
#'
#' The model assumes a gamma frailty distribution as a special case of the generalized
#' gamma distribution, i.e., \eqn{q = \sigma} , integrals are approximated using
#' the Monte Carlo summation approach.
#'
#' @inheritParams EM_func_gamma
#' @param n_mc Sample size of the Monte Carlo sample used for approximation of integrals
#'
#' @details
#' Refer to \code{\link{EM_func_gamma}} for full model specification and parameter
#' interpretation. Refer to the paper for further details.
#'
#' @seealso \code{\link{EM_func_gamma}}, \code{\link{EM_func_qunfixed}}
#'
#' @export
EM_func_gamma_mc <- function(gamma10, gamma20, rho0, beta10, beta20, alpha10, alpha20,
                             sigma, t1, t2, X, censor_type_vec,
                             print_opt = FALSE, epsilon = .0001, conv_obsll = TRUE, n_mc = 1000) {

  gamma1_k <- gamma10
  gamma2_k <- gamma20
  rho_k <- rho0

  beta1_k <- beta10
  beta2_k <- beta20
  alpha1_k <- alpha10
  alpha2_k <- alpha20
  #q_k <- q0
  #sigma_k <- sigma0

  n_param <- ncol(X)  # Single length assumption

  # Generate MC samples
  mcs <-  rgamma(n_mc, shape =sigma^-2, rate = sigma^-2)

  n_iter <- 1000
  for (iter in 1:n_iter) {




    phi_k <- phi_matrix_cpp(X, gamma1_k, gamma2_k, rho_k)

    phi11_k <- phi_k[, 1]
    phi10_k <- phi_k[, 2]
    phi01_k <- phi_k[, 3]
    phi00_k <- phi_k[, 4]

    # E-step (MC)
    E_k_c <- Expc_k_given_C_mc(t1, t2, X, beta1_k, beta2_k, alpha1_k, alpha2_k,
                               sigma,sigma, censor_type_vec,
                               phi11_k, phi10_k, phi01_k, phi00_k, mcs)

    # M1 step
    M1_phi_res <- optimize_M1_phi(X, E_k_c, censor_type_vec, c(gamma1_k, gamma2_k, rho_k))
    M1_phi_param <- M1_phi_res[["parameters"]]

    gamma1_k1 <- M1_phi_param[1:n_param]
    gamma2_k1 <- M1_phi_param[(n_param + 1):(2 * n_param)]
    rho_k1 <- M1_phi_param[(2 * n_param + 1)]

    # M2 step (MC)
    M2_LC2_res <- optimize_M2_lc2_gamma_mc(t1, t2, X, E_k_c, censor_type_vec, sigma, mcs,
                                           c(beta1_k, beta2_k, log(alpha1_k), log(alpha2_k)))
    M2_LC2_param <- M2_LC2_res[["parameters"]]

    beta1_k1 <- M2_LC2_param[1:n_param]
    beta2_k1 <- M2_LC2_param[(n_param + 1):(2 * n_param)]
    alpha1_k1 <- exp(M2_LC2_param[(2 * n_param + 1)])
    alpha2_k1 <- exp(M2_LC2_param[(2 * n_param + 2)])
    #sigma_k1 <- exp(M2_LC2_param[(2 * n_param + 3)])

    #q_k1 <- sigma_k1  # Ensuring q updates correctly

    if (conv_obsll) {
      Obs_LLk <- Observed_ll_mc(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k,
                                alpha1_k, alpha2_k, sigma,sigma,
                                t1, t2, X, censor_type_vec, mcs)

      Obs_LLk1 <- Observed_ll_mc(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
                                 alpha1_k1, alpha2_k1, sigma,sigma,
                                 t1, t2, X, censor_type_vec, mcs)

      absll_diff <- abs(Obs_LLk - Obs_LLk1)

      if (absll_diff < epsilon) {
        # print("EM Converges")
        break
      }
    } else {

      theta_k <- c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k,
                   alpha1_k, alpha2_k, sigma,sigma)
      theta_k1 <- c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
                    alpha1_k1, alpha2_k1, sigma,sigma)

      max_prd <- max(abs((theta_k - theta_k1) / theta_k))

      if (max_prd < epsilon) {
        # print("EM Converges")

        Obs_LLk1 <- Observed_ll_mc(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
                                   alpha1_k1, alpha2_k1, sigma,sigma,
                                   t1, t2, X, censor_type_vec, mcs)

        break
      }
    }

    if (print_opt) {
      print(paste("E-step No.", iter))
      print(paste("M-step Phi No.", iter))
      print(M1_phi_res)

      print(paste("M-step LC2 No.", iter))
      print(M2_LC2_res)

      cat("gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k, q_k, sigma_k -\n",
          paste(c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k,
                  alpha1_k, alpha2_k, sigma,sigma), collapse = ", "))
      print("______________________________________________________________________________________")
      cat("gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1, q_k1, sigma_k1 -\n",
          paste(c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
                  alpha1_k1, alpha2_k1, sigma,sigma), collapse = ", "))
      cat("\n")

      if (conv_obsll) {
        print(paste("Obs_llk -", Obs_LLk))
        print(paste("Obs_llk1 -", Obs_LLk1))
        print(paste("Observed LL absolute difference", absll_diff))
        print(paste("Observed LL relative difference", absll_diff / abs(Obs_LLk)))
      } else {
        par_dist_sq <- (theta_k - theta_k1)^2
        print(paste("Euclidean distance between parameters", sqrt(sum(par_dist_sq))))
        print(paste("Max absolute relative difference in parameters", max_prd))
        print(paste("Index of max absolute relative difference in parameters", which.max(abs(theta_k - theta_k1) / theta_k)))
      }
    }

    # Update parameters
    gamma1_k <- gamma1_k1
    gamma2_k <- gamma2_k1
    rho_k <- rho_k1
    beta1_k <- beta1_k1
    beta2_k <- beta2_k1
    alpha1_k <- alpha1_k1
    alpha2_k <- alpha2_k1
    # sigma_k <- sigma_k1
    # q_k <- q_k1
  }

  return(c(gamma10 = gamma1_k1[1], gamma11 = gamma1_k1[2],
           gamma20 = gamma2_k1[1], gamma21 = gamma2_k1[2],
           rho = rho_k1,
           beta10 = beta1_k1[1], beta11 = beta1_k1[2],
           beta20 = beta2_k1[1], beta21 = beta2_k1[2],
           alpha1 = alpha1_k1, alpha2 = alpha2_k1,
           q = sigma, sigma = sigma,
           Obs_LL = Obs_LLk1))
}

#' EM Algorithm for Bivariate Cure Rate Model (Exponential  distribution as a special case of generalized gamma distributionin the shared frailty)
#' Monte Carlo summation approach for approximating integrals
#'
#' This function performs the same procedure as \code{\link{EM_func_exp}}, with the
#' difference that the integrals involved in the EM algorithm are approximated using
#' the Monte Carlo summation approach.
#'
#' The model assumes an Exponential distribution as a special case of the generalized
#' gamma distribution, i.e., \eqn{q = \sigma = 1}.
#'
#' @inheritParams EM_func_exp
#' @param n_mc Sample size of the Monte Carlo sample used  in monte carlo summation for approximation of integrals
#'
#' @details
#' Refer to \code{\link{EM_func_exp}} for full model specification and parameter
#' interpretation. Monte Carlo summation is used to approximate integrals involved
#'
#' @seealso \code{\link{EM_func_exp}}, \code{\link{EM_func_qunfixed}}
#'
#' @export
EM_func_exp_mc <- function(gamma10, gamma20, rho0, beta10, beta20, alpha10, alpha20,
                           t1, t2, X, censor_type_vec,
                           print_opt = FALSE, epsilon = .0001, conv_obsll = TRUE, n_mc = 1000) {

  # Initialize parameters
  gamma1_k <- gamma10
  gamma2_k <- gamma20
  rho_k <- rho0
  beta1_k <- beta10
  beta2_k <- beta20
  alpha1_k <- alpha10
  alpha2_k <- alpha20

  n_param <- ncol(X)  # Unified parameter length

  # Set q and sigma to 1
  q_fixed <- 1
  sigma_fixed <- 1

  # Monte Carlo samples for q = 1
  mcs <- rgamma(n_mc, shape = 1, rate = 1)



  n_iter <- 1000
  for (iter in 1:n_iter) {

    # Compute phi matrix
    phi_k <- phi_matrix_cpp(X, gamma1_k, gamma2_k, rho_k)

    phi11_k <- phi_k[, 1]
    phi10_k <- phi_k[, 2]
    phi01_k <- phi_k[, 3]
    phi00_k <- phi_k[, 4]

    # E-step (MC version)
    E_k_c <- Expc_k_given_C_mc(t1, t2, X, beta1_k, beta2_k, alpha1_k, alpha2_k,
                               q_fixed, sigma_fixed, censor_type_vec,
                               phi11_k, phi10_k, phi01_k, phi00_k, mcs)

    # M1-step
    M1_phi_res <- optimize_M1_phi(X, E_k_c, censor_type_vec, c(gamma1_k, gamma2_k, rho_k))
    M1_phi_param <- M1_phi_res[["parameters"]]

    gamma1_k1 <- M1_phi_param[1:n_param]
    gamma2_k1 <- M1_phi_param[(n_param + 1):(2 * n_param)]
    rho_k1 <- M1_phi_param[(2 * n_param + 1)]

    # M2-step (MC version)
    M2_LC2_res <- optimize_M2_lc2_exp_mc(t1, t2, X, E_k_c, censor_type_vec, mcs,
                                         c(beta1_k, beta2_k, log(alpha1_k), log(alpha2_k)))
    M2_LC2_param <- M2_LC2_res[["parameters"]]

    beta1_k1 <- M2_LC2_param[1:n_param]
    beta2_k1 <- M2_LC2_param[(n_param + 1):(2 * n_param)]
    alpha1_k1 <- exp(M2_LC2_param[(2 * n_param + 1)])
    alpha2_k1 <- exp(M2_LC2_param[(2 * n_param + 2)])

    if (conv_obsll) {
      Obs_LLk <- Observed_ll_mc(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k,
                                alpha1_k, alpha2_k, q_fixed, sigma_fixed, t1, t2, X, censor_type_vec, mcs)
      Obs_LLk1 <- Observed_ll_mc(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
                                 alpha1_k1, alpha2_k1, q_fixed, sigma_fixed, t1, t2, X, censor_type_vec, mcs)

      absll_diff <- abs(Obs_LLk - Obs_LLk1)

      if (absll_diff < epsilon) {
        # print("EM Converges")
        break
      }
    } else {

      theta_k <- c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k)
      theta_k1 <- c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1)

      max_prd <- max(abs((theta_k - theta_k1) / theta_k))

      if (max_prd < epsilon) {
        # print("EM Converges")
        Obs_LLk1 <- Observed_ll_mc(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1,
                                   alpha1_k1, alpha2_k1, q_fixed, sigma_fixed, t1, t2, X, censor_type_vec, mcs)
        break
      }
    }

    if (print_opt) {
      print(paste("E-step No.", iter))
      print(paste("M-step Phi No.", iter))
      print(M1_phi_res)

      print(paste("M-step LC2 No.", iter))
      print(M2_LC2_res)

      cat("gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k -\n",
          paste(c(gamma1_k, gamma2_k, rho_k, beta1_k, beta2_k, alpha1_k, alpha2_k), collapse = ", "))
      print("______________________________________________________________________________________")
      cat("gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1 -\n",
          paste(c(gamma1_k1, gamma2_k1, rho_k1, beta1_k1, beta2_k1, alpha1_k1, alpha2_k1), collapse = ", "))
      cat("\n")

      if (conv_obsll) {
        print(paste("Obs_llk -", Obs_LLk))
        print(paste("Obs_llk1 -", Obs_LLk1))
        print(paste("Observed LL absolute difference", absll_diff))
        print(paste("Observed LL relative difference", absll_diff / abs(Obs_LLk)))
      } else {
        par_dist_sq <- (theta_k - theta_k1)^2
        print(paste("Euclidean distance between parameters", sqrt(sum(par_dist_sq))))
        print(paste("Max absolute relative difference in parameters", max_prd))
        print(paste("Index of max absolute relative difference in parameters", which.max(abs(theta_k - theta_k1) / theta_k)))
      }
    }

    # Update parameters
    gamma1_k <- gamma1_k1
    gamma2_k <- gamma2_k1
    rho_k <- rho_k1
    beta1_k <- beta1_k1
    beta2_k <- beta2_k1
    alpha1_k <- alpha1_k1
    alpha2_k <- alpha2_k1
  }

  return(c(gamma10 = gamma1_k1[1], gamma11 = gamma1_k1[2],
           gamma20 = gamma2_k1[1], gamma21 = gamma2_k1[2],
           rho = rho_k1,
           beta10 = beta1_k1[1], beta11 = beta1_k1[2],
           beta20 = beta2_k1[1], beta21 = beta2_k1[2],
           alpha1 = alpha1_k1, alpha2 = alpha2_k1,
           q = 1, sigma = 1,
           Obs_LL = Obs_LLk1))
}


