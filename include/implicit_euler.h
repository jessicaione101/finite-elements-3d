#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  energy(q, qdot) -  a function that computes the energy of the FEM system. This takes q and qdot as parameters, returns the energy value.
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_qdot - scratch space for storing velocities 
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename ENERGY, typename FORCE, typename STIFFNESS> 
inline void implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                           const Eigen::SparseMatrixd &mass, ENERGY &energy, FORCE &force, STIFFNESS &stiffness,
                           Eigen::VectorXd &tmp_qdot, Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
  static Eigen::VectorXd tmp_g;
  static Eigen::SparseMatrixd tmp_H;

  auto gradient = [&](Eigen::VectorXd& g, Eigen::VectorXd& x) {
    force(tmp_force, q + dt*x, qdot);
    g = mass * (x - qdot) + dt * -tmp_force;
  };

  auto hessian = [&](Eigen::SparseMatrixd& H, Eigen::VectorXd& x) {
    stiffness(tmp_stiffness, q + dt*x, qdot);
    H = mass + dt*dt * -tmp_stiffness;
  };

  tmp_qdot = qdot;
  newtons_method(tmp_qdot, energy, gradient, hessian, 5, tmp_g, tmp_H);

  qdot = tmp_qdot;
  q = q + dt*qdot;
}