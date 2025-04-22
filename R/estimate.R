#' Point and interval estimation of posterior distributions
#'
#' Compute the point estimate (mean) and equal-tailed credible interval to describe posterior distribution.
#'
#' @param object An object from \strong{eigenGSIW}.
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the intervals. The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob.
#' @param orthogonal A logical value indicating whether to ensure the orthogonality of the posterior mean of the eigenvectors.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An array (or matrix) containing the lower and upper bounds of the credible interval, as well as the posterior mean.
#'
#' @author Sewon Park
#' @seealso \code{eigenGSIW}
#' @importFrom magrittr `%>%`
#' @importFrom abind abind
#' @importFrom pracma gramSchmidt
#' @importFrom stats quantile
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' # generate a spiked covariance matrix:
#' n <- 50
#' p <- 100
#' K <- 3
#' leading <- c(50,20,10)
#' remaining <- rep(1, p - K)
#' Sigma0 <- diag(c(leading, remaining), p)
#'
#' # generate data
#' set.seed(413)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' 
#' # Compute both posterior means and credible intervals
#' # for eigenvalues and eigenvectors of a covariance matrix:
#' H <- 4 * diag(rep(1, p))
#' iter <- 100000
#' burnin <- floor(iter / 2)
#' thin <- 100
#' res <- besiw::eigenGSIW(X = X, k = K, H = H, nmcmc = iter, nburn = burnin, nthin = thin)
#' est <- besiw::estimate(res)
#' }
#'

estimate <- function(object, ...) {
  UseMethod("estimate")
}

#' @rdname estimate
#' @export
#'

estimate.besiw <- function(object, prob = 0.95, orthogonal = FALSE, ...) {
  left_tail <- (1 - prob) / 2
  right_tail <- 1 - left_tail
  k = object$k
  
  # posterior samples
  post_evecs <- object$Gamma
  post_evals <- object$Lambda %>% do.call('rbind',.) %>% .[, 1:k]
  
  # posterior estimates for eigenvalues
  post_mean_evals <- colMeans(post_evals)
  credintv_evals <- apply(post_evals, 2, function(x){quantile(x,c(left_tail, right_tail))})
  post_est_evals <- rbind(credintv_evals, post_mean_evals)
  rownames(post_est_evals) <- c('lower','upper','mean')
  
  # posterior estimates for eigenvectors
  temp <- lapply(splitList(post_evecs)[1:k], FUN = function(x){apply(x, 2, function(y){sign(sum(x[,1] * y)) * y})})
  post_mean_evecs = lapply(temp, rowMeans) %>% do.call(cbind,.) %>% apply(., 2, function(x){x/norm(x, type = '2')})
  if(orthogonal){
    post_mean_evecs <- post_mean_evecs %>% pracma::gramSchmidt() %>% .$Q
  }
  credintv_evecs <- apply(post_evecs %>% abind::abind(.,along = 0) %>% .[,,1:k], c(2,3), function(x){quantile(x,c(left_tail, right_tail))})
  dim(post_mean_evecs) <- c(1, dim(post_mean_evecs))
  post_est_evecs <- abind::abind(credintv_evecs, post_mean_evecs, along=1)
  dimnames(post_est_evecs)[[1]] <- c('lower','upper','mean')
  return(list(evals = post_est_evals, evecs = post_est_evecs))
}

