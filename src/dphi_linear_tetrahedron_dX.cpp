#include <dphi_linear_tetrahedron_dX.h>

void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
  Eigen::Matrix3d T;
  T.col(0) = V.row(element(1)) - V.row(element(0));
  T.col(1) = V.row(element(2)) - V.row(element(0));
  T.col(2) = V.row(element(3)) - V.row(element(0));
  Eigen::Matrix3d T_inv = T.inverse();

  dphi.row(0) = -T_inv.colwise().sum();
  dphi.bottomRows(3) = T_inv;
}