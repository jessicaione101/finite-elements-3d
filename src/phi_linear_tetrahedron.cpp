#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {
  Eigen::Matrix3d T;
  T.col(0) = V.row(element(1)) - V.row(element(0));
  T.col(1) = V.row(element(2)) - V.row(element(0));
  T.col(2) = V.row(element(3)) - V.row(element(0));

  phi.tail(3) = T.inverse() * (x - V.row(element(0)).transpose());
  phi(0) = 1 - phi(1) - phi(2) - phi(3);
}