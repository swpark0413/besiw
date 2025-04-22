#' Bayesian estimation of eigenstructure in spiked covariance models using the gSIW prior
#'
#' Provides posterior samples of eigenvalues and eigenvectors from a spiked covariance matrix using a generalized shrinkage inverse-Wishart (gSIW) prior. 
#' The method supports Bayesian estimation of the eigenstructure in high-dimensional settings.

#'
#' Kim, Lee, Park, and Lee (2025+) proposed a spiked covariance model using a generalized shrinkage inverse-Wishart (gSIW) prior for estimating the spiked eigenstructure of a covariance matrix,
#' focusing on the top \eqn{k} eigenvalues and eigenvectors. The generalized Shrinkage inverse-Wishart (gSIW) prior whose density is defined as 
#' whose density is defined as
#' \deqn{\pi^{gSIW}(\Lambda, \mathrm{U} |a_1,\ldots,a_p,b,H) (d\Lambda) (d\mathrm{U}) \propto   
#' \dfrac{etr(-\dfrac{1}{2}\mathrm{U}\Lambda^{-1} \mathrm{U}^T H)}{\prod\limits_{i=1}^p\lambda_i^{a_i}\prod\limits_{i<j}|\lambda_i-\lambda_j|^{b-1}}(d\Lambda) (d\mathrm{U})}
#' where \eqn{a_1,\ldots,a_p>0}, \eqn{b\in[0,1]}, \eqn{H \in \mathcal{C}_p :=\{ A \in \mathbb{R}^{p\times p}: A \text{ is positive definite} \}}, \eqn{\lambda_1,\cdots,\lambda_p} are the unordered eigenvalues, and \eqn{\mathrm{U}} is the eigenvector matrix of \eqn{\Sigma}.
#' By allowing each \eqn{a_i} to be defined differently for each \eqn{\lambda_i}, we can consider a prior distribution that is more general than the SIW prior (Berger et al., 2020)  to address some of the limitations of the inverse Wishart prior.
#' Unlike the SIW prior, the gSIW prior explicitly distinguishes between the spiked and non-spiked components, leading to improved posterior performance.
#' For simplicity and computational efficiency, we assume \eqn{b=1} throughout this package. 
#' 
#' 
#' We consider the spectral decomposition of \eqn{nS = QWQ^T}, 
#' where \eqn{W = \mathrm{diag}(n\hat{\lambda}_1, \ldots, n\hat{\lambda}_{n \wedge p}, 0, \ldots, 0)} 
#' and \eqn{Q} is a \eqn{p \times p} orthogonal matrix whose \eqn{i}th column is the eigenvector corresponding to the \eqn{i}th eigenvalue. 
#' We define \eqn{\Gamma = Q^T \mathrm{U}} as the rotated eigenvector matrix under the sample eigenspace.
#' Given a generalized shrinkage inverse-Wishart (gSIW) prior, which is conjugate to the multivariate normal distribution,
#' the posterior density is given by
#' \deqn{ \pi(\Lambda,\Gamma\vert \mathbb{X}_n)\propto \prod\limits_{i=1}^p \lambda_i^{-a_i - n/2} \, \mathrm{etr}\left(-\dfrac{1}{2}\Lambda^{-1}\Gamma^T(hI_p + W)\Gamma\right),}
#' where \eqn{\mathbb{X}_n} denotes the \eqn{n} observations and \eqn{h > 0} is a positive constant.
#' The gibbs sampling algorithm for generating posterior samples from \eqn{[\Lambda\vert \Gamma, \mathbb{X}_n]} and \eqn{[\Gamma\vert \Lambda, \mathbb{X}_n]} is given below:
#' \itemize{
#' \item  Independently generate \eqn{\lambda_i \sim \text{Inverse-Gamma}(a_i+n/2-1,c_i/2)}. Here, \eqn{c_i} is the \eqn{(i,i)} element of \eqn{\Gamma^T(hI_p+W)\Gamma}.
#' 
#' \item simulate \eqn{\Gamma} from \eqn{\pi(\Gamma\vert \Lambda,\mathbb{X}_n)\propto etr(-\dfrac{1}{2}\Lambda^{-1}\Gamma^T(hI_p+W)\Gamma)}.
#'
#' }
#' For more details, see Kim, Lee, Park, and Lee (2025+).
#'
#'
#' @param X A \eqn{n \times p} data matrix.
#' @param k An integer specifying the number of spiked eigenvalues.
#' @param a_vec  A numeric vector of length \eqn{p}, specifying prior hyperparameters for the generalized inverse-Wishart prior, which controls the shrinkage behavior of the eigenvalues.
#' @param H A \eqn{p \times p} positive definite matrix used in the prior specification.
#' @param nmcmc The total number of MCMC iterations.
#' @param nburn The number of burn-in iterations to discard.
#' @param nthin The thinning interval for storing posterior samples.
#' @param progress Logical. Whether to display progress updates during sampling.
#'
#' @return A list of class \code{"besiw"} containing:
#' \describe{
#'   \item{\code{Lambda}}{A list of posterior samples of eigenvalues (length: number of retained iterations).}
#'   \item{\code{Gamma}}{A list of posterior samples of eigenvector matrices.}
#'   \item{\code{k}}{The number of spiked eigenvalues used.}
#' }
#'
#' @author Sewon Park
#' @seealso \code{estimate}
#' @keywords spiked covariance eigenvalue eigenvector
#'
#' @references Kim, S., Lee, K., Park, S., and Lee, J. (2025+), "Eigenstructure inference for high-dimensional covariance with shrinkage inverse-Wishart prior",
#' 	\emph{Arxiv}, URL: \url{https://arxiv.org/abs/2412.10753}.
#' 	
#' Berger, J. O., Sun, D., and Song, C. (2020), "Bayesian analysis of the covariance matrix of a multivariate normal distribution with a new class of priors",
#' 	\emph{The Annals of Statistics}, 48(4), 2381–2403.
#' 	
#' @importFrom GPArotation Random.Start
#' @importFrom CholWishart rInvWishart
#' @importFrom progress progress_bar
#' @importFrom purrr map
#' @importFrom magrittr `%>%`
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
#' # Run a MCMC algorithm using a generalized shrinkage inverse-Wishart (gSIW) prior 
#' # for estimating eigenvalues and eigenvectors of a covariance matrix. 
#' H <- 4 * diag(rep(1, p))
#' iter <- 100000
#' burnin <- floor(iter / 2)
#' thin <- 100
#' res <- besiw::eigenGSIW(X = X, k = K, H = H, nmcmc = iter, nburn = burnin, nthin = thin)
#' }
#'



eigenGSIW <- function(X, k, a_vec = NULL, H = NULL,
                  nmcmc = 100000, nburn = floor(nmcmc / 2), nthin = 100,
                  progress = TRUE) {
  
  start_time = Sys.time()
  
  n = nrow(X)
  p = ncol(X)
  is_high = ifelse(n < p, T, F)
  S <- arma_crossprod(X, X)
  
  if(is.null(a_vec)){
    if(is_high){
      a_vec = rep(0, p)
      ed <- eig_decomp(S/n)
      evals_S <- rev(ed$values)
      mean_non_sp <- mean(evals_S[(k + 1):n])
      a_vec[1:k] <- n * mean_non_sp / (evals_S[1:k] - mean_non_sp) / 2 + 2
      a_vec[(k + 1):n] <- p / 2
      a_vec[(n + 1):p] <- p * 2
    } else {
      a_vec = rep(4, p)
    }
  } else {
    msg = check_a_vec(a_vec, k, n, p)
    if(is_high && !is.null(msg)){
      stop(msg)
    }
  }
  
  if(is.null(H)){
    H = diag(4, p)
  } else {
    if(!is_pd(H)) {
      stop("H should be positive definite.")
    }
    if(!is_scalar_identity(H) && is_high){
      stop("H must be a scalar multiple of identity matrix.")
    }
    h = H[1,1]
    if((h >= n) && is_high){
      stop("A scalar h in H = h * I_p must satisfy h < n.")
    }
  }
  
  r_vec <- a_vec + n / 2
  H_0 <- arma_matsum(H, S)
  eig_H0 <- eig_decomp(H_0)
  eigval_H_0 <- rev(eig_H0$values)
  eigvec_H_0 <- eig_H0$vectors[, p:1]
  Gamma0 <- GPArotation::Random.Start(p)
  mod_Gamma_t <- arma_crossprod(eigvec_H_0, Gamma0)
  
  out = gSIW_sampling(mod_Gamma_t, eigval_H_0 = eigval_H_0, eigvec_H_0 = eigvec_H_0,
                      r_vec = r_vec, k = k, iter = nmcmc, burnin = nburn, thin = nthin,
                      progress = progress)
  out$k = k
  out$a_vec = a_vec
  out$H = H
  class(out) = 'besiw'
  
  end_time = Sys.time()
  elapsed = round(end_time - start_time)
  
  if(progress){
    cat(sprintf("MCMC iteration: %s | Elapsed time: %d seconds\n",
                format(nmcmc, big.mark = ",", scientific = FALSE),
                elapsed))
  }
  
  return(out)
}



