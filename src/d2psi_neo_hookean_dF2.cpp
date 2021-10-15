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

  auto brown_pattern = [&](double multi_term, double sub_term, double single_term) {
    return yellow_term*multi_term + 2*blue_term*single_term*sub_term + green_term*multi_term + brown_term;
  };

  auto yellow_blue_green_pattern = [&](double multi_term, double sub_term1, double sub_term2, double single_term1, double single_term2) {
    return yellow_term*multi_term + blue_term*single_term1*sub_term2 + blue_term*single_term2*sub_term1 + green_term*multi_term;
  };

  auto negative_red_pattern = [&](double multi_term, double sub_term1, double sub_term2, double single_term1, double single_term2, double single_term3) {
    return yellow_term*multi_term + blue_term*single_term1*sub_term2 + blue_term*single_term2*sub_term1 + green_term*multi_term + orange_term*single_term3 - red_term*single_term3;
  };

  auto negative_orange_pattern = [&](double multi_term, double sub_term1, double sub_term2, double single_term1, double single_term2, double single_term3) {
    return yellow_term*multi_term + blue_term*single_term1*sub_term2 + blue_term*single_term2*sub_term1 + green_term*multi_term - orange_term*single_term3 + red_term*single_term3;
  };

  double F_single_term, F_sub_term1, F_sub_term2;

  // row 0

  F_single_term = F(0, 0);
  F_sub_term1 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);

  ddw(0, 0) = brown_pattern(F_sub_term1*F_sub_term1, F_sub_term1, F_single_term);

  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  ddw(0, 1) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 1));

  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  ddw(0, 2) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 2));

  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  ddw(0, 3) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 0));

  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  ddw(0, 4) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 1), F(2, 2));

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  ddw(0, 5) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 2), F(2, 1));

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  ddw(0, 6) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 0));

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  ddw(0, 7) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 1), F(1, 2));

  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  ddw(0, 8) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 2), F(1, 1));

  // row 1

  F_single_term = F(0, 1);
  F_sub_term1 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);

  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  ddw(1, 0) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 0));

  ddw(1, 1) = brown_pattern(F_sub_term1*F_sub_term1, F_sub_term1, F_single_term);

  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  ddw(1, 2) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 2));

  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  ddw(1, 3) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 0), F(2, 2));

  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  ddw(1, 4) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 1));

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  ddw(1, 5) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 2), F(2, 0));

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  ddw(1, 6) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 0), F(1, 2));

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  ddw(1, 7) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 1));

  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  ddw(1, 8) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 2), F(1, 0));

  // row 2

  F_single_term = F(0, 2);
  F_sub_term1 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);

  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  ddw(2, 0) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 0));

  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  ddw(2, 1) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 1));

  ddw(2, 2) = brown_pattern(F_sub_term1*F_sub_term1, F_sub_term1, F_single_term);

  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  ddw(2, 3) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 0), F(2, 1));

  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  ddw(2, 4) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 1), F(2, 0));

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  ddw(2, 5) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 2));

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  ddw(2, 6) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 0), F(1, 1));

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  ddw(2, 7) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 1), F(1, 0));

  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  ddw(2, 8) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 2));

  // row 3

  F_single_term = F(1, 0);
  F_sub_term1 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);

  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  ddw(3, 0) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 0));

  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  ddw(3, 1) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 1), F(2, 2));

  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  ddw(3, 2) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(0, 2), F(2, 1));

  ddw(3, 3) = brown_pattern(F_sub_term1*F_sub_term1, F_sub_term1, F_single_term);

  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  ddw(3, 4) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 1));

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  ddw(3, 5) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(1, 2));

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  ddw(3, 6) = yellow_blue_green_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 0));

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  ddw(3, 7) = negative_red_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 1), F(0, 2));

  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  ddw(3, 8) = negative_orange_pattern(F_sub_term1*F_sub_term2, F_sub_term1, F_sub_term2, F_single_term, F(2, 2), F(0, 1));

  // row 4
/*
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

  F_sub_term1 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);

  F_single_term1 = F(2, 1);
  F_single_term2 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 0) = -orange_term*F_single_term1 + red_term*F_single_term1 + blue_term*F_single_term2*F_sub_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 0);
  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 1) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 2) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 3) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 4) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(5, 5) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 1);
  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 6) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 0);
  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 7) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 2);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(5, 8) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  // row 6

  F_sub_term1 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);

  F_single_term1 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 0) = yellow_term*F_multi_term + blue_term*F_single_term1*F_sub_term1 + green_term*F_multi_term;

  F_single_term1 = F(1, 2);
  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 1) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 1);
  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 2) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 3) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 2);
  F_single_term2 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 4) = -orange_term*F_single_term1 + red_term*F_single_term1 + blue_term*F_single_term2*F_sub_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 1);
  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 5) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(6, 6) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 7) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 2);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(6, 8) = blue_term*F_single_term1*F_sub_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  // row 7

  F_sub_term1 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);

  F_single_term1 = F(1, 2);
  F_single_term2 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 0) = -orange_term*F_single_term1 + red_term*F_single_term1 + blue_term*F_single_term2*F_sub_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 1) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 0);
  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 2) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 2);
  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 3) = orange_term*F_single_term1 - red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 4) = blue_term*F_single_term1*F_sub_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(0, 0);
  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 5) = -orange_term*F_single_term1 + red_term*F_single_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 6) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(7, 7) = yellow_term*F_multi_term + green_term*F_multi_term;

  F_single_term1 = F(2, 2);
  F_sub_term2 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(7, 8) = blue_term*F_single_term1*F_sub_term1 + yellow_term*F_multi_term + green_term*F_multi_term;

  // row 8

  F_single_term1 = F(2, 2);
  F_sub_term1 = F(0, 0)*F(1, 1) - F(0, 1)*F(1, 0);

  F_single_term2 = F(1, 1);
  F_single_term3 = F(0, 0);
  F_sub_term2 = F(1, 1)*F(2, 2) - F(1, 2)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 0) = orange_term*F_single_term2 + yellow_term*F_multi_term + green_term*F_multi_term - red_term*F_single_term2 + blue_term*F_single_term3*F_sub_term1 + blue_term*F_single_term1*F_sub_term2;

  F_single_term2 = F(1, 0);
  F_sub_term2 = F(1, 2)*F(2, 0) - F(1, 0)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 1) = -orange_term*F_single_term2 + yellow_term*F_multi_term + green_term*F_multi_term + red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2;

  F_sub_term2 = F(1, 0)*F(2, 1) - F(1, 1)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 2) = yellow_term*F_multi_term + green_term*F_multi_term + blue_term*F_single_term1*F_sub_term2;

  F_single_term2 = F(0, 1);
  F_sub_term2 = F(0, 2)*F(2, 1) - F(0, 1)*F(2, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 3) = -orange_term*F_single_term2 + yellow_term*F_multi_term + green_term*F_multi_term + red_term*F_single_term2 + blue_term*F_single_term1*F_sub_term2;

  F_single_term2 = F(0, 0);
  F_single_term3 = F(1, 1);
  F_sub_term2 = F(0, 0)*F(2, 2) - F(0, 2)*F(2, 0);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 4) = orange_term*F_single_term2 + yellow_term*F_multi_term + green_term*F_multi_term - red_term*F_single_term2 + blue_term*F_single_term3*F_sub_term1 + blue_term*F_single_term1*F_sub_term2;

  F_sub_term2 = F(0, 1)*F(2, 0) - F(0, 0)*F(2, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 5) = yellow_term*F_multi_term + green_term*F_multi_term + blue_term*F_single_term1*F_sub_term2;

  F_sub_term2 = F(0, 1)*F(1, 2) - F(0, 2)*F(1, 1);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 6) = yellow_term*F_multi_term + green_term*F_multi_term + blue_term*F_single_term1*F_sub_term2;

  F_sub_term2 = F(0, 2)*F(1, 0) - F(0, 0)*F(1, 2);
  F_multi_term = F_sub_term1*F_sub_term2;
  ddw(8, 7) = yellow_term*F_multi_term + green_term*F_multi_term + blue_term*F_single_term1*F_sub_term2;

  F_multi_term = F_sub_term1*F_sub_term1;
  ddw(8, 8) = yellow_term*F_multi_term + green_term*F_multi_term + brown_term + 2*blue_term*F_single_term1*F_sub_term1;
  */
}