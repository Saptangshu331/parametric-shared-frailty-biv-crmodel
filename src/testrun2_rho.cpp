#include <RcppArmadillo.h>
#include <pbv.h>
#include <roptim.h>


//#include<Rcpp.h>
using namespace roptim;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(roptim)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(pbv)]]


// [[Rcpp::export]]
double log_md(double x) {
  if (std::isfinite(std::log(x))) {
    return std::log(x);
  } else {
    // Rcpp::Rcerr << "x = " << x << std::endl;

    if (std::isnan(std::log(x))) {
      return NAN;
    } else if (std::log(x) < 0) {
      // Rcpp::Rcerr << "log(x) is -Inf" << std::endl;
      return -1e+200;
    } else {
      // Rcpp::Rcerr << "log(x) is Inf" << std::endl;
      return 1e+200;
    }
  }
}




// [[Rcpp::export]]
arma::mat phi_matrix_cpp(const arma::mat &X, const arma::vec &gamma1, const arma::vec &gamma2, const double &rho) {

  int n = X.n_rows;
  arma::mat phi(n, 4, fill::zeros); // Initialize the phi matrix with zeros

  arma::vec mu1 = X*gamma1; // Compute mu1 = X1 * gamma1
  arma::vec mu2 = X*gamma2; // Compute mu2 = X2 * gamma2


  for (int k = 0; k < n; k++){

    phi(k, 0) = pbv::pbv_rcpp_pbvnorm0(mu1[k], mu2[k], rho);
    phi(k,1)=R::pnorm( mu1[k], 0.0, 1.0, true,  false ) - phi(k, 0);
    phi(k,2)=R::pnorm( mu2[k], 0.0, 1.0, true,  false ) - phi(k, 0);
    phi(k,3)=1-(phi(k,0) + phi(k,1) + phi(k,2));

  }


  return phi;
}



class M1_phi_roptim : public roptim::Functor {
private:
  arma::mat X;
  arma::mat Expc_k_given_C_mat;
  arma::vec censor_type_vec;

public:

  M1_phi_roptim(const arma::mat& X_,const arma::mat& Expc_k_given_C_mat_, const arma::vec& censor_type_vec_)
    : X(X_), Expc_k_given_C_mat(Expc_k_given_C_mat_), censor_type_vec(censor_type_vec_) {}

  // Overloaded operator to compute the objective function
  double operator()(const arma::vec& param) override {
    // Define the number of columns for each gamma vector
    int lg = X.n_cols;


    // Split parameter vector into gamma1, gamma2, gamma3
    arma::vec gamma1 = param.subvec(0, lg - 1);
    arma::vec gamma2 = param.subvec(lg, 2*lg - 1);
    // Changing for better parameter estimation

    double rho = std::isinf(std::exp(param(2 * lg))) ? 1.0 : (std::exp(param(2 * lg)) - 1.0) / (std::exp(param(2 * lg)) + 1.0) ;

    // Compute phi_matrix
    arma::mat phi_mat = phi_matrix_cpp(X, gamma1, gamma2,rho);

    // Objective function (log-likelihood)
    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i)); // Convert censor type to int

      if (censor_type == 1) {
        res += log_md(phi_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(phi_mat(i, 0))* Expc_k_given_C_mat(i, 0) +
          log_md(phi_mat(i, 1))* Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(phi_mat(i, 0))* Expc_k_given_C_mat(i, 0) +
          log_md(phi_mat(i, 2))* Expc_k_given_C_mat(i, 2);
      } else {
        res += log_md(phi_mat(i, 0))* Expc_k_given_C_mat(i, 0) +  log_md(phi_mat(i, 1))* Expc_k_given_C_mat(i, 1) + log_md(phi_mat(i, 2))* Expc_k_given_C_mat(i, 2) + log_md(phi_mat(i, 3))* Expc_k_given_C_mat(i, 3) ;

        // res += arma::dot(Expc_k_given_C_mat.row(i), arma::log(phi_mat.row(i)));
      }
      //Rcout<<" data_point_no."<<i<<"log likelihood contribution"<<res<<std::endl;

    }


    // for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
    //   int censor_type = static_cast<int>(censor_type_vec(i)); // Convert censor type to int
    //
    //   if (censor_type == 1) {
    //     res += std::log(phi_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
    //   } else if (censor_type == 2) {
    //     res += std::log(phi_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
    //       std::log(phi_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
    //   } else if (censor_type == 3) {
    //     res += std::log(phi_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
    //       std::log(phi_mat(i, 2)) * Expc_k_given_C_mat(i, 2);
    //   } else {
    //     res += arma::dot(Expc_k_given_C_mat.row(i), arma::log(phi_mat.row(i)));
    //   }
    //   //Rcout<<" data_point_no."<<i<<"log likelihood contribution"<<res<<std::endl;
    //
    // }

    return -res; // Return negative for minimization
  }
};

//Function to maximize l_c1 in EM

// [[Rcpp::export]]
Rcpp::List optimize_M1_phi(const arma::mat& X,const arma::mat& Expc_k_given_C_mat, const arma::vec& censor_type_vec,arma::vec init_params, std::string method = "Nelder-Mead") {

  M1_phi_roptim obj(X, Expc_k_given_C_mat, censor_type_vec);

  init_params(init_params.n_elem - 1) = std::log((1.0 + init_params(init_params.n_elem - 1)) / (1.0 - init_params(init_params.n_elem - 1)));



  Roptim<M1_phi_roptim> opt(method);
  opt.minimize(obj, init_params);

  arma::vec opt_params=opt.par();
    opt_params(opt_params.n_elem - 1) = (std::exp(opt_params(opt_params.n_elem - 1)) - 1.0) /(std::exp(opt_params(opt_params.n_elem - 1)) + 1.0);

  return Rcpp::List::create(
    Rcpp::Named("parameters") =opt_params ,
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}


// All functions for Montecarlo approximation in lifetime components


// [[Rcpp::export]]
double fmarginal_mc(double ti, arma::vec xi, arma::vec beta,
                    double alpha, double q, double sigma, arma::vec mcs) {

  double mui = exp(arma::accu(beta%xi));

  int n_mc = mcs.n_elem;//50 100000

  double result;


  if(q!=0){
    double lambda = exp(R::lgammafn(pow(q, -2) +( sigma / q)) - R::lgammafn(pow(q, -2)) )/pow(pow(q, -2), (sigma / q)) ;

    // arma::vec mcs = Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma( n_mc, pow(q,-2),  1/ pow(q,-2))));


    arma::vec y= arma::pow(mcs,sigma/q)/lambda;

    arma::vec scale_vec = (y*mui) ;

    NumericVector sample = Rcpp::as<NumericVector>(  Rcpp::wrap(ti/scale_vec)  ) ;

    //NumericVector cond_sample=  Rcpp::dgamma(sample , alpha, 1,  false )/Rcpp::as<NumericVector>( Rcpp::wrap(scale_vec)) ;
    NumericVector cond_sample=  Rcpp::dweibull(sample , alpha, 1,  false )/Rcpp::as<NumericVector>( Rcpp::wrap(scale_vec)) ;

    result = arma::accu(Rcpp::as<arma::vec>(Rcpp::wrap(cond_sample))) / n_mc;
  }else {

    // arma::vec mcs = Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));  // Mean = 0, SD = 1
    double lambda =exp(pow(sigma,2)/2);

    arma::vec y = arma::pow(mcs, sigma)/lambda;

    arma::vec scale_vec = y * mui;

    NumericVector sample = Rcpp::as<NumericVector>(Rcpp::wrap(ti / scale_vec));

    //NumericVector cond_sample=  Rcpp::dgamma(sample , alpha, 1,  false )/Rcpp::as<NumericVector>( Rcpp::wrap(scale_vec)) ;
    NumericVector cond_sample=  Rcpp::dweibull(sample , alpha, 1,  false )/Rcpp::as<NumericVector>( Rcpp::wrap(scale_vec)) ;

    // Compute Monte Carlo estimate
    result = arma::accu(Rcpp::as<arma::vec>(Rcpp::wrap(cond_sample))) / n_mc;
  }


  return result;

}

// [[Rcpp::export]]
double Smarginal_mc(double ti, arma::vec xi, arma::vec beta,
                    double alpha, double q, double sigma,arma::vec mcs) {

  double mui = exp(arma::accu(beta % xi));
  int n_mc = mcs.n_elem ;  //100000
  arma::vec y;

  double lambda;


  if(q!=0){
    lambda = exp(R::lgammafn(pow(q, -2) + (sigma / q)) - R::lgammafn(pow(q, -2))) / pow(pow(q, -2), (sigma / q));

    // arma::vec mcs = Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2))));

    y = arma::pow(mcs, sigma / q) / lambda;

  }else{
    lambda =exp(pow(sigma,2)/2);

     y = arma::pow(mcs, sigma)/lambda;

  }

  arma::vec scale_vec = y * mui;
  NumericVector sample = Rcpp::as<NumericVector>(Rcpp::wrap(ti / scale_vec));
  //  NumericVector cond_values = Rcpp::pgamma(sample, alpha, 1, false, false) ;

  NumericVector cond_values = Rcpp::pweibull(sample, alpha, 1, false, false) ;

  // Approximate the integral
  double result = arma::accu(Rcpp::as<arma::vec>(Rcpp::wrap(cond_values))) / n_mc;



  return result;
}

// [[Rcpp::export]]
double minus_S_t1_mc(double t1i, double t2i, arma::vec xi, arma::vec beta1, arma::vec beta2,
                     double alpha1, double alpha2, double q, double sigma,arma::vec mcs) {

  double mu1i = exp(dot(beta1, xi));
  double mu2i = exp(dot(beta2, xi));

  int n_mc = mcs.n_elem ;//100000;1000

  arma::vec y;

  double lambda;

  if(q!=0){
    lambda = exp(R::lgammafn(pow(q, -2) + (sigma / q)) - R::lgammafn(pow(q, -2))) / pow(pow(q, -2), (sigma / q));

    // arma::vec mcs = Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2))));

    y = arma::pow(mcs, sigma / q) / lambda;

  }else{
    lambda =exp(pow(sigma,2)/2);

     y = arma::pow(mcs, sigma)/lambda;

  }
  // Calculate the integrand
  arma::vec scale_vec1 = y * mu1i;
  arma::vec scale_vec2 = y * mu2i;

  NumericVector sample1 = Rcpp::as<NumericVector>(Rcpp::wrap(t1i / scale_vec1));
  NumericVector sample2 = Rcpp::as<NumericVector>(Rcpp::wrap(t2i / scale_vec2));

  NumericVector term1 = Rcpp::dweibull(sample1, alpha1, 1, false) / Rcpp::as<NumericVector>(Rcpp::wrap(scale_vec1));
  NumericVector term2 = Rcpp::pweibull(sample2, alpha2, 1, false, false);

  double result = arma::accu(Rcpp::as<arma::vec>(Rcpp::wrap(term1 * term2))) / n_mc;

  return result;
}

// [[Rcpp::export]]
double S_mc(double t1i, double t2i, arma::vec xi, arma::vec beta1, arma::vec beta2,
            double alpha1, double alpha2, double q, double sigma, arma::vec mcs) {

  double mu1i = exp(dot(beta1, xi));
  double mu2i = exp(dot(beta2, xi));

  int n_mc = mcs.n_elem ; //100000; 1000

  arma::vec y;

  double lambda;

  if(q != 0){
    lambda = exp(R::lgammafn(pow(q, -2) + (sigma / q)) - R::lgammafn(pow(q, -2))) / pow(pow(q, -2), (sigma / q));
    // arma::vec mcs = Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2))));
    y = arma::pow(mcs, sigma / q) / lambda;

  } else {

    lambda =exp(pow(sigma,2)/2);
    y = arma::pow(mcs, sigma)/lambda;

  }

  arma::vec scale_vec1 = y * mu1i;
  arma::vec scale_vec2 = y * mu2i;

  NumericVector sample1 = Rcpp::as<NumericVector>(Rcpp::wrap(t1i / scale_vec1));
  NumericVector sample2 = Rcpp::as<NumericVector>(Rcpp::wrap(t2i / scale_vec2));

  NumericVector term1 = Rcpp::pweibull(sample1, alpha1, 1, false, false);
  NumericVector term2 = Rcpp::pweibull(sample2, alpha2, 1, false, false);

  return arma::accu(Rcpp::as<arma::vec>(Rcpp::wrap(term1 * term2))) / n_mc;
}

// [[Rcpp::export]]
double f_mc(double t1i, double t2i, arma::vec xi, arma::vec beta1, arma::vec beta2,
            double alpha1, double alpha2, double q, double sigma,arma::vec mcs) {

  double mu1i = exp(dot(beta1, xi));
  double mu2i = exp(dot(beta2, xi));

  int n_mc = mcs.n_elem;//100000;  1000

  arma::vec y;

  double lambda;

  if (q != 0) {
     lambda= exp(R::lgammafn(pow(q, -2) + (sigma / q)) - R::lgammafn(pow(q, -2))) / pow(pow(q, -2), (sigma / q));
    //arma::vec mcs = Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2))));
    y = arma::pow(mcs, sigma / q) / lambda;
  } else {

    lambda =exp(pow(sigma,2)/2);

    y = arma::pow(mcs, sigma)/lambda;

  }

  arma::vec scale_vec1 = y * mu1i;
  arma::vec scale_vec2 = y * mu2i;

  NumericVector sample1 = Rcpp::as<NumericVector>(Rcpp::wrap(t1i / scale_vec1));
  NumericVector sample2 = Rcpp::as<NumericVector>(Rcpp::wrap(t2i / scale_vec2));

  NumericVector term1 = Rcpp::dweibull(sample1, alpha1, 1, false) / Rcpp::as<NumericVector>(Rcpp::wrap(scale_vec1));
  NumericVector term2 = Rcpp::dweibull(sample2, alpha2, 1, false) / Rcpp::as<NumericVector>(Rcpp::wrap(scale_vec2));

  return arma::accu(Rcpp::as<arma::vec>(Rcpp::wrap(term1 * term2))) / n_mc;
}



// Function to calculate E(k|C_i) for one observation
double Expc_k_given_ci_mc(double t1i, double t2i, arma::vec xi,
                          arma::vec beta1, arma::vec beta2, double alpha1, double alpha2,
                          double q, double sigma, int cure_type, double censor_type,
                          double phi11i, double phi10i, double phi01i, double phi00i,arma::vec mcs) {

  if (censor_type == 1) {
    if (cure_type == 1) {
      return 1.0;
    } else {
      return 0.0;
    }
  } else if (censor_type == 2) {
    if (cure_type == 3 || cure_type == 4) {
      return 0.0;
    } else if (cure_type == 1) {
      double temp1 = minus_S_t1_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
      double temp2 = fmarginal_mc(t1i, xi, beta1, alpha1, q, sigma,mcs);
      return phi11i * temp1 / ((phi11i * temp1) + (phi10i * temp2));
    } else {
      double temp1 = minus_S_t1_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
      double temp2 = fmarginal_mc(t1i, xi, beta1, alpha1, q, sigma,mcs);
      return phi10i * temp2 / ((phi11i * temp1) + (phi10i * temp2));
    }
  } else if (censor_type == 3) {
    if (cure_type == 2 || cure_type == 4) {
      return 0;
    } else if (cure_type == 1) {
      double temp1 = minus_S_t1_mc(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma,mcs);
      double temp2 = fmarginal_mc(t2i, xi, beta2, alpha2, q, sigma,mcs);
      return phi11i * temp1 / ((phi11i * temp1) + (phi01i * temp2));
    } else {
      double temp1 = minus_S_t1_mc(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma,mcs);
      double temp2 = fmarginal_mc(t2i, xi, beta2, alpha2, q, sigma,mcs);
      return phi01i * temp2 / ((phi11i * temp1) + (phi01i * temp2));
    }
  } else {

    double temp1 = S_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
    double temp2 = Smarginal_mc(t1i, xi, beta1, alpha1, q, sigma,mcs);
    double temp3 = Smarginal_mc(t2i, xi, beta2, alpha2, q, sigma,mcs);

    double denom = (phi11i * temp1) + (phi10i * temp2) + (phi01i * temp3) + phi00i;

    if (cure_type == 1) {
      return (phi11i * temp1) / denom;
    } else if (cure_type == 2) {
      return (phi10i * temp2) / denom;
    } else if (cure_type == 3) {
      return (phi01i * temp3) / denom;
    } else {
      return phi00i / denom;
    }
  }
}


// Function to calculate E(k|C_i) for all observation and store in a matrix
// [[Rcpp::export]]
arma::mat Expc_k_given_C_mc(arma::vec t1, arma::vec t2, arma::mat X,
                            arma::vec beta1, arma::vec beta2, double alpha1, double alpha2,
                            double q, double sigma, arma::vec censor_type_vec,
                            arma::vec phi11, arma::vec phi10, arma::vec phi01, arma::vec phi00,arma::vec mcs) {

  arma::mat Expc_k_given_C_mat(t1.n_elem, 4, fill::zeros);

  for (size_t i = 0; i < t1.n_elem; ++i) {
    for (int cure = 1; cure <= 4; ++cure) {
      Expc_k_given_C_mat(i, cure - 1) = Expc_k_given_ci_mc(
        t1(i), t2(i), X.row(i).t(), beta1, beta2, alpha1, alpha2, q, sigma, cure,
        censor_type_vec(i), phi11(i), phi10(i), phi01(i), phi00(i), mcs);
    }
  }

  return Expc_k_given_C_mat;
}

// [[Rcpp::export]]
arma::mat L_C2_matrix_mc(arma::vec t1, arma::vec t2,
                         arma::mat X,
                         arma::vec beta1, arma::vec beta2,
                         double alpha1, double alpha2, double q, double sigma,
                         arma::vec censor_type_vec,arma::vec mcs) {

  int n = t1.n_elem;
  arma::mat LC2_dist = arma::mat(n, 4, fill::zeros);

  for (int i = 0; i < n; ++i) {
    double t1i = t1(i);
    double t2i = t2(i);
    arma::vec xi = X.row(i).t();
    double censor_type = censor_type_vec(i);

    if (censor_type == 1) {
      LC2_dist(i, 0) = f_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
    } else if (censor_type == 4) {
      LC2_dist(i, 0) = S_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
      LC2_dist(i, 1) = Smarginal_mc(t1i, xi, beta1, alpha1, q, sigma,mcs);
      LC2_dist(i, 2) = Smarginal_mc(t2i, xi, beta2, alpha2, q, sigma,mcs);
      LC2_dist(i, 3) = 1;
    } else if (censor_type == 2) {
      LC2_dist(i, 0) = minus_S_t1_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
      LC2_dist(i, 1) = fmarginal_mc(t1i, xi, beta1, alpha1, q, sigma,mcs);
    } else {
      LC2_dist(i, 0) = minus_S_t1_mc(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma,mcs);
      LC2_dist(i, 2) = fmarginal_mc(t2i, xi, beta2, alpha2, q, sigma,mcs);
    }
  }

  return LC2_dist;
}


//Observed log likelihood based on montecarlo approx

// [[Rcpp::export]]
double Observed_ll_mc(arma::vec gamma1, arma::vec gamma2, double rho,
                      arma::vec beta1, arma::vec beta2, double alpha1, double alpha2,
                      double q, double sigma, arma::vec t1, arma::vec t2,
                      arma::mat X, arma::vec censor_type_vec,arma::vec mcs) {

  arma::mat phi = phi_matrix_cpp(X, gamma1, gamma2, rho);

  arma::mat Surv_time_dist(t1.n_elem, 4, arma::fill::zeros);
  double res = 0.0;

  for (int i = 0; i < t1.n_elem; ++i) {
    double t1i = t1(i);
    double t2i = t2(i);
    arma::vec xi = X.row(i).t();
    int censor_type = static_cast<int>(censor_type_vec(i));

    if (censor_type == 1) {
      Surv_time_dist(i, 0) = f_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
    } else if (censor_type == 4) {
      Surv_time_dist(i, 0) = S_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
      Surv_time_dist(i, 1) = Smarginal_mc(t1i, xi, beta1, alpha1, q, sigma,mcs);
      Surv_time_dist(i, 2) = Smarginal_mc(t2i, xi, beta2, alpha2, q, sigma,mcs);
      Surv_time_dist(i, 3) = 1.0;
    } else if (censor_type == 2) {
      Surv_time_dist(i, 0) = minus_S_t1_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
      Surv_time_dist(i, 1) = fmarginal_mc(t1i, xi, beta1, alpha1, q, sigma,mcs);
    } else {
      Surv_time_dist(i, 0) = minus_S_t1_mc(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma,mcs);
      Surv_time_dist(i, 2) = fmarginal_mc(t2i, xi, beta2, alpha2, q, sigma,mcs);
    }

    res += std::log(arma::dot(phi.row(i), Surv_time_dist.row(i)));
  }

  return res;
}


// Objective function L_c2 class (q fixed and based on MC approx)

class M2_LC2_qfixed_mc : public Functor {
private:
  arma::mat X;
  arma::vec t1, t2, censor_type_vec, mcs;
  arma::mat Expc_k_given_C_mat;
  double q;

public:
  // Constructor
  M2_LC2_qfixed_mc(const arma::mat& X_, const arma::vec& t1_,
                   const arma::vec& t2_, const arma::mat& Expc_k_given_C_mat_,
                   const arma::vec& censor_type_vec_, const double& q_,const arma::vec& mcs_)
    : X(X_), t1(t1_), t2(t2_), Expc_k_given_C_mat(Expc_k_given_C_mat_),
      censor_type_vec(censor_type_vec_), q(q_),mcs(mcs_) {}

  double operator()(const arma::vec& param) override {
    size_t lb = X.n_cols;  // Single variable instead of lb1 and lb2

    arma::vec beta1 = param.subvec(0, lb - 1);
    arma::vec beta2 = param.subvec(lb, 2 * lb - 1);
    double alpha1 = std::exp(param(2 * lb));   // Ensures positivity
    double alpha2 = std::exp(param(2 * lb + 1));
    double sigma = std::exp(param(2 * lb + 2)); // Ensures positivity

    // Check invalid parameter constraints
    if (q < 0 && sigma >= (-1 / q)) return INFINITY;

    // Compute L_C2_matrix
    arma::mat L_C2_mat = L_C2_matrix_mc(t1, t2, X, beta1, beta2, alpha1, alpha2, q, sigma, censor_type_vec,mcs);

    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i));

      if (censor_type == 1) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2);
      } else {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2) +
          log_md(L_C2_mat(i, 3)) * Expc_k_given_C_mat(i, 3);
      }
    }

    return -res;  // Return negative log-likelihood
  }
};


// [[Rcpp::export]]
Rcpp::List optimize_M2_lc2_q_fixed_mc(arma::vec t1, arma::vec t2, arma::mat X,
                                      arma::mat Expc_k_given_C_mat, arma::vec censor_type_vec,
                                      double q, arma::vec mcs, arma::vec init_params,
                                      std::string method = "Nelder-Mead") {

  M2_LC2_qfixed_mc obj(X, t1, t2, Expc_k_given_C_mat, censor_type_vec, q, mcs);

  Roptim<M2_LC2_qfixed_mc> opt(method);
  opt.control.maxit = 10000;

  opt.minimize(obj, init_params);

  return Rcpp::List::create(
    Rcpp::Named("parameters") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}
// Define the objective function class for gamma mc
class M2_LC2_gamma_mc : public Functor {
private:
  arma::mat X;  // Single covariate matrix
  arma::vec t1, t2, censor_type_vec,mcs;
  arma::mat Expc_k_given_C_mat;
  double sigma;
public:
  // Constructor
  M2_LC2_gamma_mc(const arma::mat& X_, const arma::vec& t1_,
                  const arma::vec& t2_, const arma::mat& Expc_k_given_C_mat_,
                  const arma::vec& censor_type_vec_,const double& sigma_,const arma::vec& mcs_)
    : X(X_), t1(t1_), t2(t2_), Expc_k_given_C_mat(Expc_k_given_C_mat_),
      censor_type_vec(censor_type_vec_),sigma(sigma_),mcs(mcs_) {}

  double operator()(const arma::vec& param) override {
    size_t n_param = X.n_cols; // Unified parameter length

    arma::vec beta1 = param.subvec(0, n_param - 1);
    arma::vec beta2 = param.subvec(n_param, 2 * n_param - 1);
    double alpha1 = std::exp(param(2 * n_param));   // Ensures positivity
    double alpha2 = std::exp(param(2 * n_param + 1));
    //double sigma = std::exp(param(2 * n_param + 2)); // Ensures positivity
    //double q = sigma;  // q is set equal to sigma

    // Compute L_C2_matrix
    arma::mat L_C2_mat = L_C2_matrix_mc(t1, t2, X, beta1, beta2, alpha1, alpha2, sigma, sigma, censor_type_vec,mcs);

    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i));

      if (censor_type == 1) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2);
      } else {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2) +
          log_md(L_C2_mat(i, 3)) * Expc_k_given_C_mat(i, 3);
      }
    }

    return -res;  // Return negative log-likelihood
  }
};


// [[Rcpp::export]]
Rcpp::List optimize_M2_lc2_gamma_mc(arma::vec t1, arma::vec t2, arma::mat X,
                                    arma::mat Expc_k_given_C_mat, arma::vec censor_type_vec,
                                    double sigma, arma::vec mcs,
                                    arma::vec init_params,
                                    std::string method = "Nelder-Mead") {

  M2_LC2_gamma_mc obj(X, t1, t2, Expc_k_given_C_mat, censor_type_vec, sigma, mcs);

  Roptim<M2_LC2_gamma_mc> opt(method);
  opt.control.maxit = 10000;

  opt.minimize(obj, init_params);

  return Rcpp::List::create(
    Rcpp::Named("parameters") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}

// Define the objective function class for the exponential case with Monte Carlo (mc) adjustment
class M2_LC2_exp_mc : public Functor {
private:
  arma::mat X;  // Covariate matrix
  arma::vec t1, t2, censor_type_vec, mcs;
  arma::mat Expc_k_given_C_mat;

public:
  // Constructor
  M2_LC2_exp_mc(const arma::mat& X_, const arma::vec& t1_,
                const arma::vec& t2_, const arma::mat& Expc_k_given_C_mat_,
                const arma::vec& censor_type_vec_, const arma::vec& mcs_)
    : X(X_), t1(t1_), t2(t2_), Expc_k_given_C_mat(Expc_k_given_C_mat_),
      censor_type_vec(censor_type_vec_), mcs(mcs_) {}

  double operator()(const arma::vec& param) override {
    size_t n_param = X.n_cols;

    arma::vec beta1 = param.subvec(0, n_param - 1);
    arma::vec beta2 = param.subvec(n_param, 2 * n_param - 1);
    double alpha1 = std::exp(param(2 * n_param));      // Positive constraint
    double alpha2 = std::exp(param(2 * n_param + 1));

    // Fix q and sigma to 1
    double q = 1.0;
    double sigma = 1.0;

    // Use Monte Carlo based likelihood matrix
    arma::mat L_C2_mat = L_C2_matrix_mc(t1, t2, X, beta1, beta2, alpha1, alpha2, q, sigma, censor_type_vec, mcs);

    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i));

      if (censor_type == 1) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2);
      } else {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2) +
          log_md(L_C2_mat(i, 3)) * Expc_k_given_C_mat(i, 3);
      }
    }

    return -res;
  }
};


// [[Rcpp::export]]
Rcpp::List optimize_M2_lc2_exp_mc(arma::vec t1, arma::vec t2, arma::mat X,
                                  arma::mat Expc_k_given_C_mat, arma::vec censor_type_vec,
                                  arma::vec mcs, arma::vec init_params,
                                  std::string method = "Nelder-Mead") {

  M2_LC2_exp_mc obj(X, t1, t2, Expc_k_given_C_mat, censor_type_vec, mcs);

  Roptim<M2_LC2_exp_mc> opt(method);
  opt.control.maxit = 10000;

  opt.minimize(obj, init_params);

  return Rcpp::List::create(
    Rcpp::Named("parameters") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}



// ================================================
// All functions for Numerical Integration
// ================================================


// Define the integrand function for numerical integration
class fmarginal_Integrand : public Func {
private:
  double ti, mui, alpha, q, sigma, lambda;
  bool use_gamma;

public:
  fmarginal_Integrand(double ti_, double mui_, double alpha_, double q_, double sigma_, bool use_gamma_)
    : ti(ti_), mui(mui_), alpha(alpha_), q(q_), sigma(sigma_), use_gamma(use_gamma_) {


    if (q != 0) {
      lambda = exp( lgamma(pow(q, -2) + (sigma / q)) - lgamma(pow(q, -2)) ) / pow(pow(q, -2), sigma / q);
    }else{

      lambda=exp(pow(sigma,2)/2);

    }
  }

  double operator()(const double& v) const {
    double y = use_gamma ? (pow(v, sigma / q) / lambda) : pow(v, sigma)/lambda;
    double scale = y * mui;

    // Compute density functions using Rcpp's built-in functions
    double density = R::dweibull(ti, alpha, scale, false);
    double mixing = use_gamma ? R::dgamma(v, pow(q, -2), 1.0 / pow(q, -2), false) : R::dlnorm(v, 0, 1, false);

    return density * mixing;
  }
};

// [[Rcpp::export]]
double fmarginal_cpp(double ti, arma::vec xi, arma::vec beta, double alpha, double q, double sigma) {

  double mui = exp(arma::accu(beta%xi)); // Compute mui

  if (q < 0 && sigma >= (-1 / q)) return NA_REAL; // Invalid parameter condition

  bool use_gamma = (q != 0);
  fmarginal_Integrand fmarginal(ti, mui, alpha, q, sigma, use_gamma);

  // Perform numerical integration
  // double err_est;
  // int err_code;
  // double integral_value = integrate(f, 0, R_PosInf, err_est, err_code);
  double integral_value;


  try {
    double err_est;
    int err_code;
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod201;

    integral_value = integrate(fmarginal, 0, R_PosInf, err_est, err_code,rule);
    // if(Rcpp::traits::is_na<REALSXP>(integral_value)||std::isnan(integral_value)){
    //   Rcpp::Rcout<< "numerical integration produces NA or NaN value in f_marginal, integral=" << integral_value<< std::endl;
    //   integral_value =fmarginal_mc(ti, (xi), (beta), alpha, q, sigma);}
  }
  catch (std::exception& e) {  // Catch standard exceptions
    Rcpp::Rcerr << "Integration failed: " << e.what() << "\nSwitching to Monte Carlo method.\n";
    int n_mc=100000;
    arma::vec mcs = use_gamma? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma( n_mc, pow(q,-2),  1/ pow(q,-2)))) :  Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    integral_value = fmarginal_mc(ti, (xi), (beta), alpha, q, sigma,mcs);
  }
  catch (...) {  // Catch any other unknown exception
    Rcpp::Rcerr << "Unknown error occurred. Using Monte Carlo approximation.\n";
    int n_mc=100000;
    arma::vec mcs = use_gamma? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma( n_mc, pow(q,-2),  1/ pow(q,-2)))) :  Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    integral_value = fmarginal_mc(ti, (xi), (beta), alpha, q, sigma,mcs);

  }


  return integral_value;
}


class SmarginalIntegrand : public Func {
private:
  double ti, alpha, q, sigma, lambda;
  arma::vec xi, beta;  //

public:
  SmarginalIntegrand(double ti_, arma::vec xi_, arma::vec beta_, double alpha_, double q_, double sigma_)
    : ti(ti_), xi(xi_), beta(beta_), alpha(alpha_), q(q_), sigma(sigma_) {

    // Compute lambda only if q != 0
    if (q != 0) {
      lambda = exp(lgamma(pow(q, -2) + sigma / q) - lgamma(pow(q, -2))) / pow(pow(q, -2), sigma / q);
    }else{
      lambda=exp(pow(sigma,2)/2);

    }


  }

  double operator()(const double& v) const {
    double mui = exp(dot(beta, xi));  //
    double y = (q != 0) ? (pow(v, sigma / q) / lambda) : (pow(v, sigma)/lambda);
    double scale = y * mui;

    // Compute survival function using Rcpp's built-in functions
    double survival = R::pweibull(ti, alpha, scale, false, false);  // P(X > ti)
    double mixing = (q != 0) ? R::dgamma(v, pow(q, -2), 1.0 / pow(q, -2), false) : R::dlnorm(v, 0, 1, false);

    return survival * mixing;
  }
};

// [[Rcpp::export]]
double Smarginal_cpp(double ti, arma::vec xi, arma::vec beta, double alpha, double q, double sigma) {
  if (q < 0 && sigma >= (-1 / q)) return NA_REAL;  // Invalid parameter condition

  SmarginalIntegrand Smarginal(ti, xi, beta, alpha, q, sigma);

  double integral_value;

  try {
    double err_est;
    int err_code;
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod201;
    integral_value= integrate(Smarginal, 0, R_PosInf, err_est, err_code,rule);

    // if(Rcpp::traits::is_na<REALSXP>(integral_value)||std::isnan(integral_value)){
    //   Rcpp::Rcout<< "numerical integration produces NA or NaN value in Smarginal, integral=" << integral_value<< std::endl;
    //   integral_value =Smarginal_mc(ti, (xi), (beta), alpha, q, sigma);}


  } catch (std::exception& e) {
    Rcpp::Rcerr << "Integration failed: " << e.what() << "\nSwitching to Monte Carlo method.\n";
    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    integral_value =Smarginal_mc(ti, xi, beta, alpha, q, sigma,mcs);  //
  } catch (...) {
    Rcpp::Rcerr << "Unknown error occurred. Using Monte Carlo approximation.\n";
    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    integral_value =Smarginal_mc(ti, xi, beta, alpha, q, sigma,mcs);
  }

  return integral_value;

}


// ========================
// minus_S_t1 Function
// ========================

// Define the integrand function for minus_S_t1

class MinusSIntegrand : public Func {
private:
  double t1i, t2i, alpha1, alpha2, q, sigma, lambda;
  arma::vec xi, beta1, beta2;

public:
  MinusSIntegrand(double t1i_, double t2i_, arma::vec xi_,
                  arma::vec beta1_, arma::vec beta2_, double alpha1_, double alpha2_,
                  double q_, double sigma_)
    : t1i(t1i_), t2i(t2i_), xi(xi_), beta1(beta1_), beta2(beta2_),
      alpha1(alpha1_), alpha2(alpha2_), q(q_), sigma(sigma_) {
    if (q != 0) {
      lambda = exp(lgamma(pow(q, -2) + sigma / q) - lgamma(pow(q, -2))) / pow(pow(q, -2), sigma / q);
    }else{
      lambda=exp(pow(sigma,2)/2);

    }
  }

  double operator()(const double& v) const {
    double mu1i = exp(dot(beta1, xi));
    double mu2i = exp(dot(beta2, xi));

    double y = (q != 0) ? (pow(v, sigma / q) / lambda) : (pow(v, sigma)/lambda) ;
    double scale1 = y * mu1i;
    double scale2 = y * mu2i;

    return R::dweibull(t1i, alpha1, scale1, false) *
      R::pweibull(t2i, alpha2, scale2, false, false) *
      ((q != 0) ? R::dgamma(v, pow(q, -2), 1.0 / pow(q, -2), false) : R::dlnorm(v, 0, 1, false));
  }
};

// [[Rcpp::export]]
double minus_S_t1_cpp(double t1i, double t2i, arma::vec xi,
                      arma::vec beta1, arma::vec beta2, double alpha1, double alpha2,
                      double q, double sigma) {
  if (q < 0 && sigma >= (-1 / q)) return NA_REAL;

  MinusSIntegrand f(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);

  double integral;

  try {
    double err_est;
    int err_code;
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod201;
    integral= integrate(f, 0, R_PosInf, err_est, err_code,rule);
  } catch (std::exception& e) {
    Rcpp::Rcerr << "Integration failed: " << e.what() << "\nSwitching to Monte Carlo method.\n";
    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));

    integral=minus_S_t1_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
  } catch (...) {
    Rcpp::Rcerr << "Unknown error occurred. Using Monte Carlo approximation.\n";
    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    integral= minus_S_t1_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
  }

  return integral;
}

// ========================
// S Function
// ========================

class SIntegrand : public Func {
private:
  double t1i, t2i, alpha1, alpha2, q, sigma, lambda;
  arma::vec xi, beta1, beta2;

public:
  SIntegrand(double t1i_, double t2i_, arma::vec xi_, arma::vec beta1_, arma::vec beta2_,
             double alpha1_, double alpha2_, double q_, double sigma_)
    : t1i(t1i_), t2i(t2i_), xi(xi_), beta1(beta1_), beta2(beta2_),
      alpha1(alpha1_), alpha2(alpha2_), q(q_), sigma(sigma_) {
    if (q != 0) {
      lambda = exp(lgamma(pow(q, -2) + sigma / q) - lgamma(pow(q, -2))) / pow(pow(q, -2), sigma / q);
    }else{
      lambda=exp(pow(sigma,2)/2);

    }
  }

  double operator()(const double& v) const {
    double mu1i = exp(dot(beta1, xi));  //
    double mu2i = exp(dot(beta2, xi));

    double y = (q != 0) ? (pow(v, sigma / q) / lambda) : (pow(v, sigma)/lambda);
    double scale1 = y * mu1i;
    double scale2 = y * mu2i;

    return R::pweibull(t1i, alpha1, scale1, false, false) *
      R::pweibull(t2i, alpha2, scale2, false, false) *
      ((q != 0) ? R::dgamma(v, pow(q, -2), 1.0 / pow(q, -2), false) : R::dlnorm(v, 0, 1, false));
  }
};

// [[Rcpp::export]]
double S_cpp(double t1i, double t2i, arma::vec xi, arma::vec beta1, arma::vec beta2,
             double alpha1, double alpha2, double q, double sigma) {
  if (q < 0 && sigma >= (-1 / q)) return NA_REAL;

  SIntegrand f(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);

  double integral;

  try {
    double err_est;
    int err_code;
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod201;

    integral = integrate(f, 0, R_PosInf, err_est, err_code, rule);

  } catch (std::exception& e) {
    Rcpp::Rcerr << "Integration failed: " << e.what() << "\nSwitching to Monte Carlo method.\n";
    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    integral = S_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
  } catch (...) {
    Rcpp::Rcerr << "Unknown error occurred. Using Monte Carlo approximation.\n";
    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    integral = S_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
  }

  return integral;
}

// ========================
// f Function

class FIntegrand : public Func {
private:
  double t1i, t2i, alpha1, alpha2, q, sigma, lambda;
  arma::vec xi, beta1, beta2;  //

public:
  FIntegrand(double t1i_, double t2i_, arma::vec xi_, arma::vec beta1_, arma::vec beta2_,
             double alpha1_, double alpha2_, double q_, double sigma_)
    : t1i(t1i_), t2i(t2i_), xi(xi_), beta1(beta1_), beta2(beta2_),
      alpha1(alpha1_), alpha2(alpha2_), q(q_), sigma(sigma_) {
    if (q != 0) {
      lambda = exp(lgamma(pow(q, -2) + sigma / q) - lgamma(pow(q, -2))) / pow(pow(q, -2), sigma / q);
    }else{
      lambda=exp(pow(sigma,2)/2);

    }
  }

  double operator()(const double& v) const {
    double mu1i = exp(dot(beta1, xi));  //
    double mu2i = exp(dot(beta2, xi));

    double y = (q != 0) ? (pow(v, sigma / q) / lambda) : pow(v, sigma)/lambda;
    double scale1 = y * mu1i;
    double scale2 = y * mu2i;

    return R::dweibull(t1i, alpha1, scale1, false) *
      R::dweibull(t2i, alpha2, scale2, false) *
      ((q != 0) ? R::dgamma(v, pow(q, -2), 1.0 / pow(q, -2), false) : R::dlnorm(v, 0, 1, false));
  }
};

// [[Rcpp::export]]
double f_cpp(double t1i, double t2i, arma::vec xi, arma::vec beta1, arma::vec beta2,
             double alpha1, double alpha2, double q, double sigma) {
  if (q < 0 && sigma >= (-1 / q)) return NA_REAL;

  FIntegrand f(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);

  double integral;

  try {
    double err_est;
    int err_code;

    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod201;
    integral = integrate(f, 0, R_PosInf, err_est, err_code, rule);

  } catch (std::exception& e) {
    Rcpp::Rcerr << "Integration failed: " << e.what() << "\nSwitching to Monte Carlo method.\n";

    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    return f_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
  } catch (...) {
    Rcpp::Rcerr << "Unknown error occurred. Using Monte Carlo approximation.\n";

    int n_mc=100000;
    arma::vec mcs = (q!=0) ? Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rgamma(n_mc, pow(q, -2), 1 / pow(q, -2)))) : Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::rlnorm(n_mc, 0, 1)));
    return f_mc(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma,mcs);
  }

  return integral;
}



// Function to calculate E(k|C_i)
double Expc_k_given_ci_cpp(double t1i, double t2i, arma::vec xi,
                           arma::vec beta1, arma::vec beta2, double alpha1, double alpha2,
                           double q, double sigma, int cure_type, double censor_type,
                           double phi11i, double phi10i, double phi01i, double phi00i) {

  if (censor_type == 1) {
    if (cure_type == 1) {
      return 1.0;
    } else {
      return 0.0;
    }
  } else if (censor_type == 2) {
    if (cure_type == 3 || cure_type == 4) {
      return 0.0;
    } else if (cure_type == 1) {
      double temp1 = minus_S_t1_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
      double temp2 = fmarginal_cpp(t1i, xi, beta1, alpha1, q, sigma);
      return phi11i * temp1 / ((phi11i * temp1) + (phi10i * temp2));
    } else {
      double temp1 = minus_S_t1_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
      double temp2 = fmarginal_cpp(t1i, xi, beta1, alpha1, q, sigma);
      return phi10i * temp2 / ((phi11i * temp1) + (phi10i * temp2));
    }
  } else if (censor_type == 3) {
    if (cure_type == 2 || cure_type == 4) {
      return 0;
    } else if (cure_type == 1) {
      double temp1 = minus_S_t1_cpp(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma);
      double temp2 = fmarginal_cpp(t2i, xi, beta2, alpha2, q, sigma);
      return phi11i * temp1 / ((phi11i * temp1) + (phi01i * temp2));
    } else {
      double temp1 = minus_S_t1_cpp(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma);
      double temp2 = fmarginal_cpp(t2i, xi, beta2, alpha2, q, sigma);
      return phi01i * temp2 / ((phi11i * temp1) + (phi01i * temp2));
    }
  } else {

    double temp1 = S_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
    double temp2 = Smarginal_cpp(t1i, xi, beta1, alpha1, q, sigma);
    double temp3 = Smarginal_cpp(t2i, xi, beta2, alpha2, q, sigma);

    double denom = (phi11i * temp1) + (phi10i * temp2) + (phi01i * temp3) + phi00i;

    if (cure_type == 1) {
      return (phi11i * temp1) / denom;
    } else if (cure_type == 2) {
      return (phi10i * temp2) / denom;
    } else if (cure_type == 3) {
      return (phi01i * temp3) / denom;
    } else {
      return phi00i / denom;
    }
  }
}

// [[Rcpp::export]]
arma::mat Expc_k_given_C_cpp(arma::vec t1, arma::vec t2, arma::mat X,
                             arma::vec beta1, arma::vec beta2, double alpha1, double alpha2,
                             double q, double sigma, arma::vec censor_type_vec,
                             arma::vec phi11, arma::vec phi10, arma::vec phi01, arma::vec phi00) {

  // Initialize the matrix to store the results (length of t1, and 4 for the 4 cure types)
  arma::mat Expc_k_given_C_mat(t1.n_elem, 4, fill::zeros);

  // Loop over all data points
  for (size_t i = 0; i < t1.n_elem; ++i) {
    for (int cure = 1; cure <= 4; ++cure) {
      Expc_k_given_C_mat(i, cure - 1) = Expc_k_given_ci_cpp(
        t1(i), t2(i), X.row(i).t(), beta1, beta2, alpha1, alpha2, q, sigma, cure,
        censor_type_vec(i), phi11(i), phi10(i), phi01(i), phi00(i));
    }
  }

  return Expc_k_given_C_mat;
}

// [[Rcpp::export]]
arma::mat L_C2_matrix_cpp(arma::vec t1, arma::vec t2,
                          arma::mat X,
                          arma::vec beta1, arma::vec beta2,
                          double alpha1, double alpha2, double q, double sigma,
                          arma::vec censor_type_vec) {

  int n = t1.n_elem;
  arma::mat LC2_dist = arma::mat(n, 4, fill::zeros); // Initialized to 0

  for (int i = 0; i < n; ++i) {
    double t1i = t1(i);
    double t2i = t2(i);
    arma::vec xi = X.row(i).t();
    double censor_type = censor_type_vec(i);

    if (censor_type == 1) {
      LC2_dist(i, 0) = f_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
    } else if (censor_type == 4) {
      LC2_dist(i, 0) = S_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
      LC2_dist(i, 1) = Smarginal_cpp(t1i, xi, beta1, alpha1, q, sigma);
      LC2_dist(i, 2) = Smarginal_cpp(t2i, xi, beta2, alpha2, q, sigma);
      LC2_dist(i, 3) = 1;
    } else if (censor_type == 2) {
      LC2_dist(i, 0) = minus_S_t1_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
      LC2_dist(i, 1) = fmarginal_cpp(t1i, xi, beta1, alpha1, q, sigma);
    } else {
      LC2_dist(i, 0) = minus_S_t1_cpp(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma);
      LC2_dist(i, 2) = fmarginal_cpp(t2i, xi, beta2, alpha2, q, sigma);
    }
  }

  return LC2_dist;
}

// [[Rcpp::export]]
double Observed_ll_cpp(arma::vec gamma1, arma::vec gamma2, double rho,
                       arma::vec beta1, arma::vec beta2, double alpha1, double alpha2,
                       double q, double sigma, arma::vec t1, arma::vec t2,
                       arma::mat X, arma::vec censor_type_vec) {

  // Calculate phi matrix
  arma::mat phi = phi_matrix_cpp(X, gamma1, gamma2, rho);

  arma::mat Surv_time_dist(t1.n_elem, 4, arma::fill::zeros);
  double res = 0.0;

  for (int i = 0; i < t1.n_elem; ++i) {
    double t1i = t1(i);
    double t2i = t2(i);
    arma::vec xi = X.row(i).t();
    int censor_type = static_cast<int>(censor_type_vec(i)); // Convert to integer

    if (censor_type == 1) {
      Surv_time_dist(i, 0) = f_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
    } else if (censor_type == 4) {
      Surv_time_dist(i, 0) = S_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
      Surv_time_dist(i, 1) = Smarginal_cpp(t1i, xi, beta1, alpha1, q, sigma);
      Surv_time_dist(i, 2) = Smarginal_cpp(t2i, xi, beta2, alpha2, q, sigma);
      Surv_time_dist(i, 3) = 1.0;
    } else if (censor_type == 2) {
      Surv_time_dist(i, 0) = minus_S_t1_cpp(t1i, t2i, xi, beta1, beta2, alpha1, alpha2, q, sigma);
      Surv_time_dist(i, 1) = fmarginal_cpp(t1i, xi, beta1, alpha1, q, sigma);
    } else {
      Surv_time_dist(i, 0) = minus_S_t1_cpp(t2i, t1i, xi, beta2, beta1, alpha2, alpha1, q, sigma);
      Surv_time_dist(i, 2) = fmarginal_cpp(t2i, xi, beta2, alpha2, q, sigma);
    }

    // Compute log likelihood
    res += std::log(arma::dot(phi.row(i), Surv_time_dist.row(i)));
  }

  return res;
}



// Define the objective function class
class M2_LC2_nb_qunfixed : public Functor {
private:
  arma::mat X;
  arma::vec t1, t2, censor_type_vec;
  arma::mat Expc_k_given_C_mat;

public:
  // Constructor
  M2_LC2_nb_qunfixed(const arma::mat& X_, const arma::vec& t1_,
                     const arma::vec& t2_, const arma::mat& Expc_k_given_C_mat_,
                     const arma::vec& censor_type_vec_)
    : X(X_), t1(t1_), t2(t2_), Expc_k_given_C_mat(Expc_k_given_C_mat_),
      censor_type_vec(censor_type_vec_) {}

  double operator()(const arma::vec& param) override {
    size_t lb = X.n_cols;  // Single variable instead of lb1 and lb2

    arma::vec beta1 = param.subvec(0, lb - 1);
    arma::vec beta2 = param.subvec(lb, 2 * lb - 1);
    double alpha1 = std::exp(param(2 * lb));   // Ensures positivity
    double alpha2 = std::exp(param(2 * lb + 1));
    double q = param(2 * lb + 2);
    double sigma = std::exp(param(2 * lb + 3)); // Ensures positivity

    // Check invalid parameter constraints
    if (q < 0 && sigma >= (-1 / q)) return INFINITY;

    // Compute L_C2_matrix
    arma::mat L_C2_mat = L_C2_matrix_cpp(t1, t2, X, beta1, beta2, alpha1, alpha2, q, sigma, censor_type_vec);

    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i));

      if (censor_type == 1) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i,2);
      } else {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2) +
          log_md(L_C2_mat(i, 3)) * Expc_k_given_C_mat(i, 3);
      }
    }

    return -res;  // Return negative log-likelihood
  }
};



// [[Rcpp::export]]
Rcpp::List optimize_M2_lc2(arma::vec t1, arma::vec t2, arma::mat X,
                           arma::mat Expc_k_given_C_mat, arma::vec censor_type_vec,
                           arma::vec init_params,
                           std::string method = "Nelder-Mead") {

  M2_LC2_nb_qunfixed obj(X, t1, t2, Expc_k_given_C_mat, censor_type_vec);

  Roptim<M2_LC2_nb_qunfixed> opt(method);
  opt.control.maxit = 10000;

  opt.minimize(obj, init_params);

  return Rcpp::List::create(
    Rcpp::Named("parameters") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}



// Define the objective function class (q fixed case)
class M2_LC2_qfixed : public Functor {
private:
  arma::mat X;
  arma::vec t1, t2, censor_type_vec;
  arma::mat Expc_k_given_C_mat;
  double q;

public:
  // Constructor
  M2_LC2_qfixed(const arma::mat& X_, const arma::vec& t1_,
                const arma::vec& t2_, const arma::mat& Expc_k_given_C_mat_,
                const arma::vec& censor_type_vec_, const double& q_)
    : X(X_), t1(t1_), t2(t2_), Expc_k_given_C_mat(Expc_k_given_C_mat_),
      censor_type_vec(censor_type_vec_), q(q_) {}

  double operator()(const arma::vec& param) override {
    size_t lb = X.n_cols;  // Single variable instead of lb1 and lb2

    arma::vec beta1 = param.subvec(0, lb - 1);
    arma::vec beta2 = param.subvec(lb, 2 * lb - 1);
    double alpha1 = std::exp(param(2 * lb));   // Ensures positivity
    double alpha2 = std::exp(param(2 * lb + 1));
    double sigma = std::exp(param(2 * lb + 2)); // Ensures positivity

    // Check invalid parameter constraints
    if (q < 0 && sigma >= (-1 / q)) return INFINITY;

    // Compute L_C2_matrix
    arma::mat L_C2_mat = L_C2_matrix_cpp(t1, t2, X, beta1, beta2, alpha1, alpha2, q, sigma, censor_type_vec);

    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i));

      if (censor_type == 1) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2);
      } else {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2) +
          log_md(L_C2_mat(i, 3)) * Expc_k_given_C_mat(i, 3);
      }
    }

    return -res;  // Return negative log-likelihood
  }
};


// [[Rcpp::export]]
Rcpp::List optimize_M2_lc2_q_fixed(arma::vec t1, arma::vec t2, arma::mat X,
                                   arma::mat Expc_k_given_C_mat, arma::vec censor_type_vec,
                                   double q, arma::vec init_params,
                                   std::string method = "Nelder-Mead") {

  M2_LC2_qfixed obj(X, t1, t2, Expc_k_given_C_mat, censor_type_vec, q);

  Roptim<M2_LC2_qfixed> opt(method);
  opt.control.maxit = 10000;

  opt.minimize(obj, init_params);

  return Rcpp::List::create(
    Rcpp::Named("parameters") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}

// Define the objective function class
class M2_LC2_gamma : public Functor {
private:
  arma::mat X;  // Single covariate matrix
  arma::vec t1, t2, censor_type_vec;
  arma::mat Expc_k_given_C_mat;

public:
  // Constructor
  M2_LC2_gamma(const arma::mat& X_, const arma::vec& t1_,
               const arma::vec& t2_, const arma::mat& Expc_k_given_C_mat_,
               const arma::vec& censor_type_vec_)
    : X(X_), t1(t1_), t2(t2_), Expc_k_given_C_mat(Expc_k_given_C_mat_),
      censor_type_vec(censor_type_vec_) {}

  double operator()(const arma::vec& param) override {
    size_t n_param = X.n_cols; // Unified parameter length

    arma::vec beta1 = param.subvec(0, n_param - 1);
    arma::vec beta2 = param.subvec(n_param, 2 * n_param - 1);
    double alpha1 = std::exp(param(2 * n_param));   // Ensures positivity
    double alpha2 = std::exp(param(2 * n_param + 1));
    double sigma = std::exp(param(2 * n_param + 2)); // Ensures positivity
    double q = sigma;  // q is set equal to sigma

    // Compute L_C2_matrix
    arma::mat L_C2_mat = L_C2_matrix_cpp(t1, t2, X, beta1, beta2, alpha1, alpha2, q, sigma, censor_type_vec);

    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i));

      if (censor_type == 1) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2);
      } else {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2) +
          log_md(L_C2_mat(i, 3)) * Expc_k_given_C_mat(i, 3);
      }
    }

    return -res;  // Return negative log-likelihood
  }
};



// [[Rcpp::export]]
Rcpp::List optimize_M2_lc2_gamma(arma::vec t1, arma::vec t2, arma::mat X,
                                 arma::mat Expc_k_given_C_mat, arma::vec censor_type_vec,
                                 arma::vec init_params,
                                 std::string method = "Nelder-Mead") {

  M2_LC2_gamma obj(X, t1, t2, Expc_k_given_C_mat, censor_type_vec);

  Roptim<M2_LC2_gamma> opt(method);
  opt.control.maxit = 10000;

  opt.minimize(obj, init_params);

  return Rcpp::List::create(
    Rcpp::Named("parameters") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}

// Define the objective function class
class M2_LC2_exp : public Functor {
private:
  arma::mat X;  // Single covariate matrix
  arma::vec t1, t2, censor_type_vec;
  arma::mat Expc_k_given_C_mat;

public:
  // Constructor
  M2_LC2_exp(const arma::mat& X_, const arma::vec& t1_,
             const arma::vec& t2_, const arma::mat& Expc_k_given_C_mat_,
             const arma::vec& censor_type_vec_)
    : X(X_), t1(t1_), t2(t2_), Expc_k_given_C_mat(Expc_k_given_C_mat_),
      censor_type_vec(censor_type_vec_) {}

  double operator()(const arma::vec& param) override {
    size_t n_param = X.n_cols; // Unified parameter length

    arma::vec beta1 = param.subvec(0, n_param - 1);
    arma::vec beta2 = param.subvec(n_param, 2 * n_param - 1);
    double alpha1 = std::exp(param(2 * n_param));   // Ensures positivity
    double alpha2 = std::exp(param(2 * n_param + 1));

    // Fix q and sigma to 1
    double q = 1.0;
    double sigma = 1.0;

    // Compute L_C2_matrix
    arma::mat L_C2_mat = L_C2_matrix_cpp(t1, t2, X, beta1, beta2, alpha1, alpha2, q, sigma, censor_type_vec);

    double res = 0.0;

    for (size_t i = 0; i < censor_type_vec.n_elem; ++i) {
      int censor_type = static_cast<int>(censor_type_vec(i));

      if (censor_type == 1) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0);
      } else if (censor_type == 2) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1);
      } else if (censor_type == 3) {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2);
      } else {
        res += log_md(L_C2_mat(i, 0)) * Expc_k_given_C_mat(i, 0) +
          log_md(L_C2_mat(i, 1)) * Expc_k_given_C_mat(i, 1) +
          log_md(L_C2_mat(i, 2)) * Expc_k_given_C_mat(i, 2) +
          log_md(L_C2_mat(i, 3)) * Expc_k_given_C_mat(i, 3);
      }
    }

    return -res;  // Return negative log-likelihood
  }
};




// [[Rcpp::export]]
Rcpp::List optimize_M2_lc2_exp(arma::vec t1, arma::vec t2, arma::mat X,
                               arma::mat Expc_k_given_C_mat, arma::vec censor_type_vec,
                               arma::vec init_params,
                               std::string method = "Nelder-Mead") {

  M2_LC2_exp obj(X, t1, t2, Expc_k_given_C_mat, censor_type_vec);

  Roptim<M2_LC2_exp> opt(method);
  opt.control.maxit = 10000;

  opt.minimize(obj, init_params);

  return Rcpp::List::create(
    Rcpp::Named("parameters") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}


double SF_KM_cpp(double x, double y, arma::vec t1, arma::vec t2, arma::vec status_t1, arma::vec status_t2) {
  int n = t1.n_elem;
  arma::vec delta_i_hat = status_t1 % status_t2;
  arma::vec n_i_hat(n, arma::fill::zeros);

  // Compute n_i_hat using a loop
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (std::max(t1[j], t2[j]) >= std::max(t1[i], t2[i])) {
        n_i_hat[i] += 1;
      }
    }
  }

  // Compute num
  double num = arma::mean(arma::conv_to<arma::vec>::from((t1 >= x) && (t2 >= y)));

  // Compute denominator
  arma::uvec idx = arma::find(arma::max(t1, t2) < std::max(x, y));
  arma::vec term = 1.0 - ((1.0 - delta_i_hat.elem(idx)) / n_i_hat.elem(idx));
  double denom = arma::prod(term);

  return num / denom;
}









