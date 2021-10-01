#include <d2psi_neo_hookean_dq2.h>

void d2psi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
  double J = F.determinant();
  double tr = (F.transpose() * F).trace();

  // colors in variable names are arbitrary
  double yellow_term = 2*D;
  double blue_term = -4*C / (3*pow(J, 5./3));
  double green_term = 10*C*tr / (9*pow(J, 8./3));
  double red_term = 2*C*tr / (3*pow(J, 5./3));
  double orange_term = 2*D * (J-1);
  double brown_term = 2*C / pow(J, 2./3);

  double F_single_term1, F_single_term2, F_single_term3;
  double F_sub_term1, F_sub_term2;
  double F_multi_term;

  // row 0

  F_single_term1 = F(0, 0);
  F_sub_term1 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(0, 0) = yellow_term*F_multi_term + brown_term - 2*blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 1) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 2) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 3) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_single_term2 = F(2, 2);
  F_single_term3 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 4) = orange_term*F_single_term2 + yellow_term*F_multi_term - red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2 + blue_term*F_single_term3*F_sub_term1 + green_term*F_multi_term;

  F_single_term2 = F(2, 1);
  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 5) = -orange_term*F_single_term2 + yellow_term*F_multi_term + red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 6) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_single_term2 = F(1, 2);
  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 7) = -orange_term*F_single_term2 + yellow_term*F_multi_term + red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_single_term2 = F(1, 1);
  F_single_term3 = F(2, 2);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(0, 8) = orange_term*F_single_term2 + yellow_term*F_multi_term - red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2 + blue_term*F_single_term3*F_sub_term1 + green_term*F_multi_term;

  // row 1

  F_sub_term1 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);

  F_single_term1 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 0) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(1, 1) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 2) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1;

  F_single_term1 = F(2, 2);
  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 3) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 4) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_single_term1 = F(2, 0);
  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 5) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 2);
  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 6) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 7) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 0);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(1, 8) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  // row 2

  F_sub_term1 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);

  F_single_term1 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 0) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 1) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(2, 2) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 1);
  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 3) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 0);
  F_single_term2 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 4) = -orange_term*F_single_term1 + red_term*F_single_term1 + blue_term*F_single_term2*F_sub_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 5) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 1);
  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 6) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 0);
  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 7) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 2);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(2, 8) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  // row 3

  F_sub_term1 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);

  F_single_term1 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 0) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_single_term1 = F(2, 2);
  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 1) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 1);
  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 2) = orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(3, 3) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 4) = yellow_term*F_multi_term + red_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 5) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 6) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 2);
  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 7) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 1);
  F_single_term2 = F(2, 2);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(3, 8) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + blue_term*F_single_term2*F_sub_term1 + green_term*F_multi_term;

  // row 4

  F_single_term1 = F(1, 1);
  F_sub_term1 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);

  F_single_term2 = F(2, 2);
  F_single_term3 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 0) = orange_term*F_single_term2 + yellow_term*F_multi_term - red_term*F_single_term2 + blue_term*F_single_term3*F_sub_term1 + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 1) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_single_term2 = F(2, 0);
  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 2) = -orange_term*F_single_term2 + yellow_term*F_multi_term + red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 3) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(4, 4) = yellow_term*F_multi_term + brown_term +2*blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 5) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_single_term2 = F(0, 2);
  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 6) = -orange_term*F_single_term2 + yellow_term*F_multi_term + red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 7) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  F_single_term2 = F(0, 0);
  F_single_term3 = F(2, 2);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(4, 8) = orange_term*F_single_term2 + yellow_term*F_multi_term - red_term*F_single_term2 + blue_term*F_single_term3*F_sub_term1 + blue_term*F_single_term1*F_sub_term2 + green_term*F_multi_term;

  // row 5

}