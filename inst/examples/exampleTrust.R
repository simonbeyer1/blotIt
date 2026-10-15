## Minimize a quadratic with minimum at (1, -2)
objfun <- function(x) {
  d <- x - c(1, -2)
  list(value = sum(d^2), gradient = 2 * d, hessian = diag(2, 2))
}
fit <- trust(objfun, parinit = c(a = 0, b = 0), rinit = 1, rmax = 10)
fit$argument
