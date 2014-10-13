#include "arapbenchmarksolver.h"

// C++ standard library
#include <iostream>
#include <vector>

#include "igl/slice.h"
#include "igl/svd3x3/polar_svd3x3.h"

namespace arap {
namespace demo {

ArapBenchmarkSolver::ArapBenchmarkSolver(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces, const Eigen::VectorXi& fixed,
    int max_iteration)
  : Solver(vertices, faces, fixed, max_iteration) {
  arap_data_.max_iter = max_iteration;
  arap_data_.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
  igl::arap_precomputation(vertices, faces, vertices.cols(),
      fixed, arap_data_);
}

void ArapBenchmarkSolver::Precompute() {
  vertices_updated_ = vertices_;
}

void ArapBenchmarkSolver::SolvePreprocess(const Eigen::MatrixXd& fixed_vertices) {
  igl::arap_solve(fixed_vertices, arap_data_, vertices_updated_);
}

void ArapBenchmarkSolver::SolveOneIteration() {
  // Do nothing here.
  // We don't want to expose all the inner states of libigl's implementation.
}

Energy ArapBenchmarkSolver::ComputeEnergy() const {
  return Energy();
}

}  // namespace demo
}  // namespace arap
