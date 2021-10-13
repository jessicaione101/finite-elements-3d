#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>

void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
  Eigen::Matrix1212d M;
  mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);

  Eigen::Vector12d qdot_tet;
  qdot_tet << qdot.segment(element(0)*3, 3),
              qdot.segment(element(1)*3, 3),
              qdot.segment(element(2)*3, 3),
              qdot.segment(element(3)*3, 3);

  T = 0.5 * qdot_tet.transpose() * M * qdot_tet;
}