#ifndef _ARAP_DEMO_ADMMSOLVER_H_
#define _ARAP_DEMO_ADMMSOLVER_H_

#include "Eigen/CholmodSupport"

#include "solver.h"

namespace arap {
namespace demo {

class AdmmSolver : public Solver {
 public:
  AdmmSolver(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
      const Eigen::VectorXi& fixed, int max_iteration, double rho);

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
  // A helper function to compute |matrix_id|'s |variable_id|-th variable's
  // position.
  int GetMatrixVariablePos(int vertex_num, int matrix_id, int variable_id);

  // A helper function to verify our linear solve. Returns true if the current
  // vertices_updated_ and rotations_ are the optimal solution.
  bool CheckLinearSolve() const;

  // Compute the linear solve energy. Used in CheckLinearSolve.
  double ComputeLinearSolveEnergy(const Eigen::MatrixXd &vertices,
      const std::vector<Eigen::Matrix3d> &rotations) const;

  // Compute the SVD solve energy. Used in CheckSVDSolve.
  double ComputeSVDSolveEnergy() const;

  // Check whether a matrix is in SO(3).
  bool IsSO3(const Eigen::Matrix3d &S) const;

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
  Eigen::SparseMatrix<double> M_;
  // Sparse linear solver. Use Cholmod from SuiteSparse.
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
};

inline int AdmmSolver::GetMatrixVariablePos(int vertex_num, int matrix_id,
    int variable_id) {
  return vertex_num + 3 * matrix_id + variable_id;
}

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ADMMSOLVER_H_
