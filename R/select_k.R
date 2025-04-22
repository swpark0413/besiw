#' Select the optimal number of spiked eigenvalues via WAIC
#'
#' Selects the number of spiked eigenvalues \eqn{k} in a covariance matrix by computing the Watanabe-Akaike Information Criterion (WAIC) 
#' under the generalized shrinkage inverse-Wishart (gSIW) prior.
#' 
#' The function evaluates a set of candidate values for \eqn{k} by fitting spiked covariance models 
#' using the gSIW prior through \code{eigenGSIW()}, and computing their WAIC scores.
#' The value of \eqn{k} that yields the lowest WAIC is selected as the optimal number of spiked eigenvalues.
#'
#'
#' @param X An \eqn{n \times p} data matrix.
#' @param kmax  An integer specifying the maximum number of components to be considered. It is automatically adjusted to be no larger than \code{min(n, p) - 1}.
#' @param ncores An integer specifying the number of CPU cores for parallel computation. Default is 2.
#' @param ... Additional arguments passed to \code{eigenGSIW()}, such as prior settings or sampling options.
#'
#' @return An integer giving the selected number of spiked eigenvalues \eqn{k} that minimizes WAIC.
#'
#' @author Sewon Park
#' @seealso \code{eigenGSIW}
#' @keywords spiked covariance
#'
#' @references Kim, S., Lee, K., Park, S., and Lee, J. (2025+), "Eigenstructure inference for high-dimensional covariance with shrinkage inverse-Wishart prior",
#' 	\emph{Arxiv}, URL: \url{https://arxiv.org/abs/2412.10753}.
#'
#' 	Watanabe, S. and Opper, M. (2010). "Asymptotic equivalence of bayes cross validation and
#' 	widely applicable information criterion in singular learning theory", \emph{Journal of machine
#' 	learning research}, 11(12).
#' 	
#' 	Gelman, A., Hwang, J., and Vehtari, A. (2014). "Understanding predictive information criteria for Bayesian models." \emph{Statistics and computing},
#' 	24(6), 997-1016.
#'
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession
#' @importFrom magrittr `%>%`
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # generate a spiked covariance matrix:
#' n <- 10
#' p <- 20
#' k <- 3
#' Sigma0 <-diag(c(100,80,60,rep(1,p-k)))
#' 
#' # generate data
#' set.seed(413)
#' X <- MASS::mvrnorm(n, mu=rep(0,p),Sigma = Sigma0)
#' 
#' # Determine the number of spiked eigenvalues k in the covariance matrix 
#' # using the Watanabe-Akaike Information Criterion (WAIC).
#' H <- 4 *diag(rep(1,p))
#' iter <- 100000; burnin <- floor(iter/2); thin <- 100
#' res <- select_k(X = X, kmax = 10, ncores = 10, H = H, nmcmc = iter, nburn = burnin, nthin = thin)
#' }
#' 

select_k = function(X, kmax = NULL, ncores = 2, ...){
  future::plan(future::multisession, workers = ncores)
  n = nrow(X)
  p = ncol(X)
  kmax = min(min(n,p)-1, kmax)
  kvec = 1:kmax
  res <- kvec %>% furrr::future_map(function(k){eigenGSIW(X = X, k = k, ..., progress = FALSE) %>% compute_waic()},
                                    .options = furrr::furrr_options(seed = 413)) %>% do.call('rbind', .)
  res <- apply(res, 2, function(x){which.min(x)})
  return(min(res))
}


compute_waic = function(res){
  p = ncol(X)
  cov_list = transform_cov(p, res$Lambda, res$Gamma)
  bar_D <- - 2 * sum(sapply(cov_list , function(t){arma_dmvnorm(X,rep(0,p),t,logd=TRUE)}))/length(cov_list);
  LPPD <- sum( log( rowMeans( sapply( cov_list , function(t) arma_dmvnorm(x,rep(0,p),t,log=FALSE) ) ) ) )
  WAIC <- - 2 * LPPD + 2 * sum( apply( log( sapply( cov_list , function(t) arma_dmvnorm(x,rep(0,p),t,log=FALSE) ) ) , 1 , var ) )
  return(WAIC)
}
