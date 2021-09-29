#include <d2psi_neo_hookean_dq2.h>

void d2psi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
  double J = F.determinant();
  double tr = (F.transpose() * F).trace();

  double yellow_term = 2*D;
  double blue_term = -4*C / (3*pow(J, 5./3));
  double green_term = 10*C*tr / (9*pow(J, 8./3));
  double red_term = 2*C*tr / (3*pow(J, 5./3));
  double orange_term = 2*D * (J-1);
  double brown_term = 2*C / pow(J, 2./3);

  double F_single_term, F_sub_term, F_multi_term;

  F_single_term = F(0, 0);
  F_sub_term = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term*F_sub_term;
  ddw(0, 0) = yellow_term*F_multi_term + brown_term - 2*F_single_term*blue_term*F_sub_term + green_term*F_multi_term;








}