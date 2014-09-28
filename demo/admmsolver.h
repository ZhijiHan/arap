#ifndef _ARAP_DEMO_ADMMSOLVER_H_
#define _ARAP_DEMO_ADMMSOLVER_H_

#include "solver.h"

namespace arap {
namespace demo {

class AdmmSolver : public Solver {
 public:
  AdmmSolver(double rho);

  // Pre computes cotangent weight and left matrix in ADMM.
  void Precompute();

  // This function gives a chance to do preprocessing before starting the
  // solving cycle.
  void SolvePreprocess(const Eigen::MatrixXd& fixed_vertices);
  // Solves for one iteration. This function helps us analyse the algorithm.
  void SolveOneIteration();

  // Given all the current data members, compute the energy.
  double ComputeEnergy() const;

 private:
  // Computes the cotangent angle in one face indicated by |face_id|.
  Eigen::Vector3d ComputeCotangent(int face_id) const;

  // Defines private data members for AdmmSolver.
  // Control parameters in ADMM.
  double rho_;
  // S_ is the variables introduced in ADMM.
  std::vector<Eigen::Matrix3d> S_;
  // T_ is the corresponding dual variables for rotations_ - S_ = 0.
  std::vector<Eigen::Matrix3d> T_;
  // u_ is the dual variables for constraints p - c = 0. The row of u_ is the
  // # of fixed vertices, and the column is 3.
  Eigen::MatrixXd u_;
  // M is the left matrix in the linear solve.
  Eigen::MatrixXD M_;
};

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ADMMSOLVER_H_
