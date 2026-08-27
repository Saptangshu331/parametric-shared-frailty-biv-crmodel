#' Real Data used in real data analysis for Bivariate Cure Rate Model
#'
#' This dataset contains observed lifetimes(time to blindness), censoring indicators, and covariates
#' for diabetic retinopathy patients for treated and untreated eye.
#'
#' @format A data frame with 10 variables:
#' \describe{
#'   \item{id}{Unique identifier for each patient}
#'   \item{obstime_trt}{Observed lifetime for the treated eye}
#'   \item{obstime_untrt}{Observed lifetime for the untreated eye}
#'   \item{risk_trt}{Risk indicator for the treated eye}
#'   \item{risk_untrt}{Risk indicator for the untreated eye}
#'   \item{status_trt}{Event(blindness) indicator for the treated eye (1 = event, 0 = censored)}
#'   \item{status_untrt}{Event(blindness) indicator for the untreated eye (1 = event, 0 = censored)}
#'   \item{laser}{Indicator variable for laser treatment}
#'   \item{type}{Categorical variable indicating type/classification}
#'   \item{age}{Age of the patient}
#' }
#'
#' @references
#' Blair, A. L., Hadden, D. R., Weaver, J. A., Archer, D. B., Johnston, P. B.,
#' and Maguire, C. J. (1976).
#' \emph{The 5-year prognosis for vision in diabetes}.
#' \emph{American Journal of Ophthalmology}, 81(3), 383--396.
#'
#'
"diabetic_retinopathy"
