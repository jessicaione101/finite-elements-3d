#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) {
  Eigen::Matrix1212d dV;
  std::vector<Eigen::Triplet<double>> coefficients;
  coefficients.reserve(q.size()*q.size());

  for (int i = 0; i < T.rows(); ++i) {
    d2V_linear_tetrahedron_dq2(dV, q, V, T.row(i), v0(i), C, D);

    int row, col;

    for (row = 0; row < 3; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 0) * 3 + col, -dV(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 1) * 3 + col - 3, -dV(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 2) * 3 + col - 6, -dV(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 3) * 3 + col - 9, -dV(row, col));
    }

    for (row = 3; row < 6; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 0) * 3 + col, -dV(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 1) * 3 + col - 3, -dV(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 2) * 3 + col - 6, -dV(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 3) * 3 + col - 9, -dV(row, col));
    }

    for (row = 6; row < 9; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 0) * 3 + col, -dV(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 1) * 3 + col - 3, -dV(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 2) * 3 + col - 6, -dV(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 3) * 3 + col - 9, -dV(row, col));
    }

    for (row = 9; row < 12; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 0) * 3 + col, -dV(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 1) * 3 + col - 3, -dV(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 2) * 3 + col - 6, -dV(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 3) * 3 + col - 9, -dV(row, col));
    }
  }

  K.resize(q.size(), q.size());
  K.setFromTriplets(coefficients.begin(), coefficients.end());
}