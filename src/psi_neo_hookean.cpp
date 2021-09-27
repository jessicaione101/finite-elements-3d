#include <psi_neo_hookean.h>

void psi_neo_hookean(double &psi, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
  double J = F.determinant();
  double tr = (F.transpose() * F).trace();

  psi = C * (tr / pow(J, 2./3) - 3) + D * pow(J-1, 2);
}