#ifndef _ARAP_DEMO_ARAPSOLVER_H_
#define _ARAP_DEMO_ARAPSOLVER_H_

#include "solver.h"

#include "Eigen/SparseLU"

namespace arap {
namespace demo {

class ArapSolver : public Solver {
 public:
  ArapSolver(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
      const Eigen::VectorXi& fixed, int max_iteration);

  // Pre computes cotangent weight and left matrix in ARAP.
  void Precompute();

  // This function gives a chance to do preprocessing before starting the
  // solving cycle.
  void SolvePreprocess(const Eigen::MatrixXd& fixed_vertices);
  // Solves for one iteration. This function helps us analyse the algorithm.
  void SolveOneIteration();
  // Do post process after all the iterations are done.
  void SolvePostprocess();

  // Given all the current data members, compute the energy.
  Energy ComputeEnergy() const;

 private:
  // Computes the cotangent angle in one face indicated by |face_id|.
  Eigen::Vector3d ComputeCotangent(int face_id) const;

  // The discrete Laplace-Beltrami operator in the left hand side of the linear
  // system defined in equation (9) in the paper.
  Eigen::SparseMatrix<double> lb_operator_;
  // Sparse linear solver.
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
};

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ARAPSOLVER_H_
