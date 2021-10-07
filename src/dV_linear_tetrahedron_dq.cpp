#include <dV_linear_tetrahedron_dq.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {
  auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    Eigen::Matrix34d x;
    x.col(0) = q.segment(3*element(0), 3);
    x.col(1) = q.segment(3*element(1), 3);
    x.col(2) = q.segment(3*element(2), 3);
    x.col(3) = q.segment(3*element(3), 3);

    Eigen::Matrix43d dphi;
    dphi_linear_tetrahedron_dX(dphi, V, element, X);

    Eigen::Matrix3d F = x * dphi;

    Eigen::Vector9d dw;
    dpsi_neo_hookean_dF(dw, F, C, D);

    Eigen::Matrix<double, 12, 9> B_transpose = Eigen::Matrix<double, 12, 9>::Zero();
    for (int i = 0; i < dphi.rows(); ++i)
      B_transpose.block(i*3, 0, 1, 3) = B_transpose.block(i*3+1, 3, 1, 3) = B_transpose.block(i*3+2, 6, 1, 3) = dphi.row(i);

    dV = B_transpose * dw;
  };

  quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);
}