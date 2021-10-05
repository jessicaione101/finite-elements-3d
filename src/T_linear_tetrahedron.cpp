#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>

void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
  Eigen::Matrix1212d M;
  mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);
  T = 0.5 * qdot.transpose() * M * qdot;
}