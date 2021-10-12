#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {
  Eigen::Matrix1212d M_j;
  std::vector<Eigen::Triplet<double>> coefficients;
  coefficients.reserve(qdot.size()*qdot.size());

  for (int i = 0; i < T.rows(); ++i) {
    mass_matrix_linear_tetrahedron(M_j, qdot, T.row(i), density, v0(i));

    int row, col;

    for (row = 0; row < 3; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 0) * 3 + col, M_j(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 1) * 3 + col - 3, M_j(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 2) * 3 + col - 6, M_j(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 0) * 3 + row, T(i, 3) * 3 + col - 9, M_j(row, col));
    }

    for (row = 3; row < 6; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 0) * 3 + col, M_j(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 1) * 3 + col - 3, M_j(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 2) * 3 + col - 6, M_j(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 1) * 3 + row - 3, T(i, 3) * 3 + col - 9, M_j(row, col));
    }

    for (row = 6; row < 9; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 0) * 3 + col, M_j(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 1) * 3 + col - 3, M_j(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 2) * 3 + col - 6, M_j(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 2) * 3 + row - 6, T(i, 3) * 3 + col - 9, M_j(row, col));
    }

    for (row = 9; row < 12; ++row) {
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 0) * 3 + col, M_j(row, col));
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 1) * 3 + col - 3, M_j(row, col));
      for (col = 6; col < 9; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 2) * 3 + col - 6, M_j(row, col));
      for (col = 9; col < 12; ++col)
        coefficients.emplace_back(T(i, 3) * 3 + row - 9, T(i, 3) * 3 + col - 9, M_j(row, col));
    }
  }

  M.resize(qdot.size(), qdot.size());
  M.setFromTriplets(coefficients.begin(), coefficients.end());
}