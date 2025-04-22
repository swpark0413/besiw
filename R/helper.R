check_a_vec <- function(a_vec, k, n, p) {
  msg <- NULL  
  if (!all(diff(a_vec[1:k]) >= 0)) {
    msg <- "Condition is violated: a_1 <= ... <= a_k must hold."
  } else if (!all(a_vec[(k+1):n] == a_vec[k+1])) {
    msg <- "Condition is violated: a_{k+1}, ..., a_n must be equal."
  } else if (a_vec[k] > a_vec[k+1]) {
    msg <- "Condition is violated: a_k <= a_{k+1} must hold."
  } else if (!all(a_vec[(n+1):p] == a_vec[n+1])) {
    msg <- "Condition is violated: a_{n+1} = ... = a_p must be equal."
  } else if (a_vec[n] > a_vec[n+1]) {
    msg <- "Condition is violated: a_n <= a_{n+1} must hold."
  }
  return(msg)
}

is_scalar_identity<- function(H, tol = 1e-8) {
  d <- diag(H)
  if (sd(d) > tol) return(FALSE)
  if (!Matrix::isDiagonal(H)) return(FALSE)
  return(TRUE)
}