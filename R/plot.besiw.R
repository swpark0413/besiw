#' Plot diagnostics for posterior samples using \pkg{bayesplot} 
#'
#' Provides visual diagnostics of MCMC samples for eigenvalues and eigenvectors, including trace plots, histograms, posterior densities, and autocorrelations.
#'
#' The function supports plotting the posterior samples of either the top \eqn{k} eigenvalues (\code{par = "lambda"}) or the entries of the corresponding eigenvectors (\code{par = "gamma"}).
#' Users can select specific elements to visualize, and combine multiple diagnostic plots using \code{type = "all"}.
#'
#' @param object An object of class \code{"besiw"}, typically returned from \code{eigenGSIW()}.
#' @param par Character, either \code{"lambda"} (default) to plot eigenvalues or \code{"gamma"} to plot eigenvector elements.
#' @param elem A list of indices to specify which elements to plot.  
#' For \code{par = "lambda"}, supply a list like \code{list(1, 3)} to plot the 1st and 3rd eigenvalues.  
#' For \code{par = "gamma"}, supply indices as pairs like \code{list(c(1,1), c(2,3))} to indicate matrix positions.
#' 
#' If \code{elem} is not specified, then:
#' for \code{par = “lambda”}, the top \eqn{k} eigenvalues are plotted;
#' for \code{par = “gamma”}, the diagonal elements \eqn{(i,i)} for \eqn{i = 1,\ldots,k} of the eigenvector matrix are plotted by default.
#' @param type Type of plot to generate. One of \code{"trace"}, \code{"hist"}, \code{"area"}, \code{"acf"}, or \code{"all"} (default).  
#' If \code{"all"}, all four types of diagnostics will be displayed in a 2×2 layout using \pkg{patchwork}.
#'
#' @return A ggplot object (or a patchwork layout of multiple ggplots) displaying the requested diagnostic plot(s).
#'
#' @author Sewon Park
#'
#' @importFrom bayesplot color_scheme_set mcmc_trace mcmc_hist mcmc_areas mcmc_acf_bar
#' @importFrom patchwork plot_layout
#' @importFrom magrittr `%>%`
#' @export
#'
#' @examples
#'
#' \dontrun{
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
#' # estimate eigenvalues and eigenvectors using a generalized shrinkage inverse-Wishart (gSIW) prior:
#' H <- 4 * diag(rep(1, p))
#' iter <- 100000
#' burnin <- floor(iter / 2)
#' thin <- 100
#' res <- besiw::eigenGSIW(X = X, k = K, H = H, nmcmc = iter, nburn = burnin, nthin = thin)
#'
#' # full diagnostics (trace, histogram, density (area), and autocorrelation) for MCMC samples of the k eigenvalues of the covariance matrix
#' plot(res, par = 'lambda', type = 'all')
#'
#' # trace plot for MCMC samples of the first and third eigenvalues
#' plot(res, par = 'lambda', elem = list(1,3), type = 'trace')
#'
#' # full diagnostics (trace, histogram, density (area), and autocorrelation) for the diagonal elements of the orthogonal matrix of k eigenvectors
#' plot(res, par = 'gamma', type = 'all')
#'
#' # autocorrelation plot for the (1,1) and (2,3) elements of the orthogonal matrix of k eigenvectors
#' plot(res, par = 'gamma', elem = list(c(1,1), c(2,3)), type = 'acf')
#'}

plot.besiw = function(object, par = c('lambda', 'gamma'), elem = list(), type = c('all', 'trace', 'hist', 'area', 'acf')){
  k = object$k
  post = extract_posterior(object = object, par = par)
  par = match.arg(par)
  plot_type = match.arg(type)
  if(par == 'lambda'){
    params = sapply(elem, function(x){paste0('lam', x)})
    if(length(elem) == 0){
      params = paste0('lam', 1:k)
    }
  } else if(par == 'gamma'){
    params = sapply(elem, function(x){paste0('gam(', paste(x, collapse = ','), ')')})
    if(length(elem) == 0){
      params = paste0('gam(', 1:k, ',', 1:k, ')')
    }
  }
  bayesplot::color_scheme_set('brightblue')
  if(plot_type == 'trace'){
    p = bayesplot::mcmc_trace(post, pars = params)
  } else if(plot_type == 'hist'){
    p = bayesplot::mcmc_hist(post, pars = params)
  } else if(plot_type == 'area'){
    p = bayesplot::mcmc_areas(post, pars = params, prob = 0.95)
  } else if(plot_type == 'acf'){
    p = bayesplot::mcmc_acf(post, pars = params)
  } else if(plot_type == 'all'){
    bayesplot::color_scheme_set('darkgray')
    p0 = bayesplot::mcmc_trace(post, pars = params)
    bayesplot::color_scheme_set('green')
    p1 = bayesplot::mcmc_hist(post, pars = params)
    bayesplot::color_scheme_set('blue')
    p2 = bayesplot::mcmc_areas(post, pars = params, prob = 0.95)
    bayesplot::color_scheme_set('red')
    p3 = bayesplot::mcmc_acf_bar(post, pars = params)
    p = p0 + p1 + p2 + p3 + patchwork::plot_layout(ncol = 2, nrow = 2)
  }
  return(p)
}

extract_posterior = function(object, par = c('lambda', 'gamma')){
  par = match.arg(par)
  k = object$k
  p = length(object[[1]][[1]])
  niter = length(object[[1]])
  if(par == 'lambda'){
    parname = 'lam'
    coln = paste0(parname, 1:k)
    post_sample = object$Lambda %>% lapply(., function(x){x[1:k]}) %>% do.call('rbind',.)
    total_elem = k
  } else if(par == 'gamma'){
    parname = 'gam'
    coln = paste0('gam', expand.grid(1:p, 1:k) %>% apply(., 1, function(x){sprintf("(%d,%d)", x[1], x[2] )}))
    post_sample = lapply(splitList(object$Gamma)[1:k], FUN = function(x){apply(x, 2, function(y){sign(sum(x[,1] * y)) * y})}) %>% do.call('rbind',.) %>% t()
    total_elem = p * k
  }
  colnames(post_sample) = coln
  posterior = array(post_sample, dim = c(niter, 1, total_elem))
  dimnames(posterior)[[1]] = NULL
  dimnames(posterior)[[2]] = 'chain:1'
  dimnames(posterior)[[3]] = coln
  names(dimnames(posterior)) = c('iterations', 'chains', 'parameters')
  return(posterior)
}



