#ifndef _ARAP_DEMO_ARAPBENCHMARKSOLVER_H_
#define _ARAP_DEMO_ARAPBENCHMARKSOLVER_H_

#include "solver.h"

#include <igl/svd3x3/arap.h>

namespace arap {
namespace demo {

// This class wraps up the arap implementation from libigl.
class ArapBenchmarkSolver : public Solver {
 public:
  ArapBenchmarkSolver(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
      const Eigen::VectorXi& fixed, int max_iteration);

  void Precompute();

  // This function gives a chance to do preprocessing before starting the
  // solving cycle.
  void SolvePreprocess(const Eigen::MatrixXd& fixed_vertices);
  // Solves for one iteration. This function helps us analyse the algorithm.
  void SolveOneIteration();

  Energy ComputeEnergy() const;

 private:
  // Wraps up the solver from libigl.
  igl::ARAPData arap_data_;
};

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ARAPBENCHMARKSOLVER_H_
