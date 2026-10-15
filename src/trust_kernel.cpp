// Trust-region optimizer with the trust::trust interface used by blotIt and a
// More-Sorensen subproblem solver. A non-finite Hessian or failed eigensolve
// gives a steepest-descent step instead of an error.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/BLAS.h>     // defines the FCONE macro when USE_FC_LEN_T is set
#include <R_ext/Lapack.h>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <vector>
#include <cmath>
#include <limits>

// Fallback for old R that lacks the hidden Fortran string-length arguments.
#ifndef FCONE
# define FCONE
#endif

using namespace Rcpp;


// Sets the OpenMP thread count (read by MKL) and returns the previous value,
// NA_INTEGER without OpenMP.
//
// [[Rcpp::export]]
int setBlasThreads(int n) {
#ifdef _OPENMP
  int prev = omp_get_max_threads();
  if (n > 0) omp_set_num_threads(n);
  return prev;
#else
  return NA_INTEGER;
#endif
}

namespace {

const double INF = std::numeric_limits<double>::infinity();

// Symmetric eigendecomposition via LAPACK dsyevd: `vals` ascending, `vecs`
// column-major. Returns false instead of throwing on non-finite input or
// LAPACK failure.
bool eigen_sym(const std::vector<double>& A, int K,
               std::vector<double>& vals, std::vector<double>& vecs) {
  for (std::size_t i = 0; i < A.size(); ++i)
    if (!std::isfinite(A[i])) return false;

  vecs = A;                 // dsyevd overwrites the matrix with eigenvectors
  vals.assign(K, 0.0);
  char jobz = 'V', uplo = 'U';
  int n = K, lda = K, info = 0;

  // Workspace query.
  double wkopt = 0.0;
  int iwkopt = 0, lwork = -1, liwork = -1;
  F77_CALL(dsyevd)(&jobz, &uplo, &n, vecs.data(), &lda, vals.data(),
                   &wkopt, &lwork, &iwkopt, &liwork, &info FCONE FCONE);
  if (info != 0) return false;

  lwork = static_cast<int>(wkopt);
  liwork = iwkopt;
  std::vector<double> work(lwork > 0 ? lwork : 1);
  std::vector<int> iwork(liwork > 0 ? liwork : 1);

  F77_CALL(dsyevd)(&jobz, &uplo, &n, vecs.data(), &lda, vals.data(),
                   work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE);
  if (info != 0) return false;

  for (int i = 0; i < K; ++i) if (!std::isfinite(vals[i])) return false;
  return true;
}

// Minimises g^T p + 1/2 p^T H p subject to ||p|| <= r from the eigensystem of
// H. Writes the step to `p_out`, the predicted reduction to `pred_red` and
// whether the interior Newton step was taken to `is_newton`.
void trust_sub(int K, const std::vector<double>& g,
               const std::vector<double>& vals, const std::vector<double>& vecs,
               double r, std::vector<double>& p_out,
               double& pred_red, bool& is_newton) {
  // q = V^T g
  std::vector<double> q(K, 0.0);
  for (int j = 0; j < K; ++j) {
    double s = 0.0;
    for (int i = 0; i < K; ++i) s += vecs[i + (std::size_t) j * K] * g[i];
    q[j] = s;
  }

  const double lam_min = vals[0];
  is_newton = false;

  // Interior Newton step when H is positive definite.
  if (lam_min > 1e-12) {
    std::vector<double> y(K);
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) { y[j] = -q[j] / vals[j]; pn2 += y[j] * y[j]; }
    if (std::sqrt(pn2) <= r) {
      for (int i = 0; i < K; ++i) {
        double s = 0.0;
        for (int j = 0; j < K; ++j) s += vecs[i + (std::size_t) j * K] * y[j];
        p_out[i] = s;
      }
      double pred = 0.0;
      for (int j = 0; j < K; ++j) pred += 0.5 * q[j] * q[j] / vals[j];
      pred_red = pred;
      is_newton = true;
      return;
    }
  }

  // Boundary solution: find shift sigma so that ||(H + sigma I)^-1 g|| = r.
  // sigma must keep all (vals[j] + sigma) > 0, so sigma > -lam_min.
  double sigma_lo = std::max(-lam_min, 0.0) + 1e-12;
  double sigma_hi = sigma_lo + 1.0;
  for (int it = 0; it < 60; ++it) {
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) {
      double d = vals[j] + sigma_hi;
      pn2 += (q[j] / d) * (q[j] / d);
    }
    if (std::sqrt(pn2) <= r) break;
    sigma_hi *= 2.0;
  }
  double sigma = sigma_hi;
  for (int it = 0; it < 80; ++it) {
    sigma = 0.5 * (sigma_lo + sigma_hi);
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) {
      double d = vals[j] + sigma;
      pn2 += (q[j] / d) * (q[j] / d);
    }
    double pn = std::sqrt(pn2);
    if (std::fabs(pn - r) < 1e-10 * r) break;
    if (pn > r) sigma_lo = sigma; else sigma_hi = sigma;
  }

  std::vector<double> y(K);
  for (int j = 0; j < K; ++j) y[j] = -q[j] / (vals[j] + sigma);
  for (int i = 0; i < K; ++i) {
    double s = 0.0;
    for (int j = 0; j < K; ++j) s += vecs[i + (std::size_t) j * K] * y[j];
    p_out[i] = s;
  }
  double pred = 0.0;
  for (int j = 0; j < K; ++j) {
    double d = vals[j] + sigma;
    pred += q[j] * q[j] / d - 0.5 * vals[j] * q[j] * q[j] / (d * d);
  }
  pred_red = pred;
}

bool all_finite(const std::vector<double>& v) {
  for (std::size_t i = 0; i < v.size(); ++i) if (!std::isfinite(v[i])) return false;
  return true;
}

}  // namespace


// [[Rcpp::export]]
List trust_kernel(Function fn,
                  NumericVector parinit,
                  double rinit,
                  double rmax,
                  int    iterlim,
                  double fterm,
                  double mterm,
                  bool   minimize,
                  Nullable<NumericVector> parscale_) {

  const int K = parinit.size();
  const bool has_names = parinit.hasAttribute("names");
  CharacterVector parnames;
  if (has_names) parnames = parinit.names();

  const double sgn = minimize ? 1.0 : -1.0;   // work in a minimisation frame

  bool rescale = parscale_.isNotNull();
  std::vector<double> parscale(K, 1.0);
  if (rescale) {
    NumericVector ps(parscale_);
    for (int i = 0; i < K; ++i) parscale[i] = ps[i];
  }

  auto eval = [&](const std::vector<double>& x, double& val,
                  std::vector<double>& g, std::vector<double>& H) -> bool {
    NumericVector xv(x.begin(), x.end());
    if (has_names) xv.names() = parnames;
    List out;
    try {
      out = as<List>(fn(xv));
    } catch (...) {
      return false;
    }
    double v = as<double>(out["value"]);
    val = sgn * v;
    if (!std::isfinite(val)) return false;   // non-finite value -> reject point
    NumericVector gv = as<NumericVector>(out["gradient"]);
    NumericMatrix Hv = as<NumericMatrix>(out["hessian"]);
    g.assign(K, 0.0);
    H.assign((std::size_t) K * K, 0.0);
    for (int i = 0; i < K; ++i) g[i] = sgn * gv[i];
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i) H[i + (std::size_t) j * K] = sgn * Hv(i, j);
    return true;
  };

  std::vector<double> theta(parinit.begin(), parinit.end());
  double f = 0.0;
  std::vector<double> g, H;
  if (!eval(theta, f, g, H))
    stop("trust: objective not finite at the initial parameters");

  double r = rinit;
  bool converged = false;
  int iter = 0;

  std::vector<double> eigvals, eigvecs, g_s(K), p_s(K), p(K), theta_try(K);

  for (iter = 1; iter <= iterlim; ++iter) {
    if (!all_finite(g)) break;   // no usable gradient -> stop gracefully

    // Subproblem in parscale coordinates: B_s = B / (s s^T), g_s = g / s.
    std::vector<double> B_s = H;
    if (rescale) {
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          B_s[i + (std::size_t) j * K] /= (parscale[i] * parscale[j]);
    }
    for (int i = 0; i < K; ++i) g_s[i] = rescale ? g[i] / parscale[i] : g[i];

    bool ok = eigen_sym(B_s, K, eigvals, eigvecs);
    if (!ok) {
      // Identity surrogate -> steepest-descent step within the trust region.
      eigvals.assign(K, 1.0);
      eigvecs.assign((std::size_t) K * K, 0.0);
      for (int i = 0; i < K; ++i) eigvecs[i + (std::size_t) i * K] = 1.0;
    }

    bool is_newton = false;
    double pred_red = 0.0;
    trust_sub(K, g_s, eigvals, eigvecs, r, p_s, pred_red, is_newton);

    // Back to original coordinates.
    for (int i = 0; i < K; ++i) {
      p[i] = rescale ? p_s[i] / parscale[i] : p_s[i];
      theta_try[i] = theta[i] + p[i];
    }

    double f_try = 0.0;
    std::vector<double> g_try, H_try;
    bool eval_ok = eval(theta_try, f_try, g_try, H_try);

    double actual_red = f - f_try;            // > 0 means improvement
    double rho = (eval_ok && std::isfinite(pred_red) && pred_red > 0.0)
                   ? actual_red / pred_red
                   : -INF;

    bool is_terminate = eval_ok && std::isfinite(pred_red) &&
                        (std::fabs(actual_red) < fterm ||
                         std::fabs(pred_red)   < mterm);

    bool accept;
    if (is_terminate) {
      accept = (f_try < f);
    } else if (rho < 0.25) {
      accept = false;
      r *= 0.25;
    } else {
      accept = true;
      if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
    }

    if (accept && eval_ok) {
      theta = theta_try;
      f = f_try;
      g = g_try;
      H = H_try;
    }
    if (is_terminate) { converged = true; break; }
  }

  NumericVector arg_out(theta.begin(), theta.end());
  if (has_names) arg_out.names() = parnames;

  return List::create(
      Named("argument")   = arg_out,
      Named("value")      = sgn * f,
      Named("converged")  = converged,
      Named("iterations") = iter > iterlim ? iterlim : iter);
}
