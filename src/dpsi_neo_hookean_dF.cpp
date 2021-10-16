#include <dpsi_neo_hookean_dF.h>

void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
  Eigen::Vector9d F_flat;
  F_flat << F.row(0).transpose(), F.row(1).transpose(), F.row(2).transpose();

  Eigen::Vector9d F_terms;
  F_terms(0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
  F_terms(1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
  F_terms(2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
  F_terms(3) = F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2);
  F_terms(4) = F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0);
  F_terms(5) = F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1);
  F_terms(6) = F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1);
  F_terms(7) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
  F_terms(8) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);

  double J = F.determinant();
  double tr = (F.transpose() * F).trace();

  dw = 2*D * (J-1) * F_terms + 2*C / pow(J, 2./3) * F_flat - 2*C*tr / (3 * pow(J, 5./3)) * F_terms;
}