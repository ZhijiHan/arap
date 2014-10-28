#ifndef _ARAP_DEMO_ADAPTADMMFREESOLVER_H_
#define _ARAP_DEMO_ADAPTADMMFREESOLVER_H_

#include "solver.h"

#include "Eigen/SparseLU"

namespace arap {
namespace demo {

class AdaptAdmmFreeSolver : public Solver {
 public:
  AdaptAdmmFreeSolver(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
      const Eigen::VectorXi& fixed, int max_iteration, double rho);

  // Pre computes cotangent weight and left matrix in ADMM.
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

  double GetRho() const { return rho_; }

 private:
  // A helper function to compute |matrix_id|'s |variable_id|-th variable's
  // position.
  int GetMatrixVariablePos(int matrix_id, int variable_id);

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

  // Compute the squared norm of the primal residual.
  // The primal residual is defined as:
  // r = Ax + Bz - c.
  // and this function returns r^2.
  // where all the notations come from the problem definition:
  // min f(x) + g(z)
  // s.t. Ax + Bz = c
  // For AdaptAdmmFixedSolver, it is extremely easy to compute the primal
  // residual. The return value is just the sum of all the squared norm of
  // R_i - S_i, plus the squared norm of p' - c.
  double ComputePrimalResidual() const;

  // Compute the squared norm of the dual residual.
  // The dual residual is defined as:
  // s = rho * A'B(z^k+1 - z^K).
  // The order for computing all the variables are as follows:
  // Given Rk, Sk, Tk, pk, rhok+1, now we start the k + 1 iteration:
  // Step 1: use rhok+1 to compute Rk+1 and pk+1.
  // Step 2: use rhok+1 to compute Sk+1.
  // Step 3: Use Rk+1 and Sk+1 to compute Tk+1.
  // Step 4: Use Rk+1 and Sk+1 to compute primal dual rk+1.
  // Step 5: Use rhok+1, Sk+1 and Sk to compute sk+1.
  // Step 6: Use rk+1 and sk+1 to compute rhok+2.
  double ComputeDualResidual() const;

  // Defines private data members for AdaptAdmmFreeSolver.
  // Control parameters in ADMM.
  double rho_;
  // S_ is the variables introduced in ADMM.
  std::vector<Eigen::Matrix3d> S_;
  // S_pre_ is the variable used to cache S_ from this iteration, before we
  // update S_. This helps us compute the dual residual.
  std::vector<Eigen::Matrix3d> S_pre_;
  // T_ is the corresponding dual variables for rotations_ - S_ = 0.
  std::vector<Eigen::Matrix3d> T_;
  // u_ is the dual variables for constraints p - c = 0. The row of u_ is the
  // # of fixed vertices, and the column is 3.
  Eigen::MatrixXd u_;
  // M is the left matrix in the linear solve.
  Eigen::SparseMatrix<double> M_;
};

inline int AdaptAdmmFreeSolver::GetMatrixVariablePos(int matrix_id,
    int variable_id) {
  int vertex_num = vertices_.rows();
  return vertex_num + 3 * matrix_id + variable_id;
}

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ADAPTADMMFREESOLVER_H_
