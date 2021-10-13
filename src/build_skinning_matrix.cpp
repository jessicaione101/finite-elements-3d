#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>

bool outside_tetrahedron(Eigen::Ref<const Eigen::Vector4d> phi) {
  for (int i = 0; i < phi.rows(); ++i)
    if (phi(i) < 0)
      return true;
  return false;
}

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                           Eigen::Ref<const Eigen::MatrixXd> V_skin) {
  Eigen::Vector4d phi;
  Eigen::Vector3d x;
  std::vector<Eigen::Triplet<double>> coefficients;
  coefficients.reserve(V_skin.rows() * V.rows());

  for (int i = 0; i < V_skin.rows(); ++i) {
    x = V_skin.row(i);
    for (int j = 0; j < T.rows(); ++j) {
      phi_linear_tetrahedron(phi, V, T.row(j), x);
      if (outside_tetrahedron(phi))
        continue;
      for (int k = 0; k < phi.rows(); ++k)
        coefficients.emplace_back(i, T(j, k), phi(k));
      break;
    }
  }

  N.resize(V_skin.rows(), V.rows());
  N.setFromTriplets(coefficients.begin(), coefficients.end());
}