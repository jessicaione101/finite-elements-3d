#include <Eigen/Dense>
#include <EigenTypes.h>
#include <limits>

namespace {
  const double tol = std::numeric_limits<double>::epsilon();
  const double p = 0.5;
  const double c = 1e-8;

  template<typename Objective>
  double line_search(const Objective& f, const Eigen::VectorXd& x, const Eigen::VectorXd& d, const Eigen::VectorXd& g) {
    double alpha = 1;
    while (alpha > tol && f(x + alpha*d) > f(x) + c * d.dot(g))
      alpha *= p;
    return alpha;
  }
}

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
void newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrixd> solver;
  Eigen::VectorXd d;

  for (unsigned int i = 0; i < maxSteps; ++i) {
    g(tmp_g, x0);
    if (tmp_g.norm() < tol)
      return;

    H(tmp_H, x0);
    solver.compute(tmp_H);
    if (solver.info() != Eigen::Success)
      return;
    d = solver.solve(-tmp_g);
    if (solver.info() != Eigen::Success)
      return;

    double alpha = line_search(f, x0, d, tmp_g);
    x0 = x0 + alpha*d;
  }
}
