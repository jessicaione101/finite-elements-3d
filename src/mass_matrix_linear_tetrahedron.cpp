#include <mass_matrix_linear_tetrahedron.h>

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
  auto I = Eigen::Matrix3d::Identity();
  double term = 6 * density * volume;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      double integral = i == j ? 1./60 : 1./120;
      M.block(i*3, j*3, 3, 3) = term * integral * I;
    }
  }
 }