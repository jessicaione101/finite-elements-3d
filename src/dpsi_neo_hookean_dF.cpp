#include <dpsi_neo_hookean_dF.h>

void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
  static Eigen::Vector9d F_terms;
  F_terms(0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
//  F_terms(1) = F(, ) * F(, ) - F(, ) * F(, );
//  F_terms(2) = F(, ) * F(, ) - F(, ) * F(, );
//  F_terms(3) = F(, ) * F(, ) - F(, ) * F(, );
//  F_terms(4) = F(, ) * F(, ) - F(, ) * F(, );
//  F_terms(5) = F(, ) * F(, ) - F(, ) * F(, );
//  F_terms(6) = F(, ) * F(, ) - F(, ) * F(, );
//  F_terms(7) = F(, ) * F(, ) - F(, ) * F(, );
//  F_terms(8) = F(, ) * F(, ) - F(, ) * F(, );
}