#include <assemble_forces.h>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) {
  f = Eigen::VectorXd::Zero(q.size());
  Eigen::Vector12d dV;

  for (int i = 0; i < T.rows(); ++i) {
    dV_linear_tetrahedron_dq(dV, q, V, T.row(i), v0(i), C, D);

    f.segment(T(i, 0)*3, 3) -= dV.segment(0, 3);
    f.segment(T(i, 1)*3, 3) -= dV.segment(3, 3);
    f.segment(T(i, 2)*3, 3) -= dV.segment(6, 3);
    f.segment(T(i, 3)*3, 3) -= dV.segment(9, 3);
  }
};