#include "arapsolver.h"

namespace arap {
namespace demo {

ArapSolver::ArapSolver() {
}

// If it is the first time to register data, return true; else return false.
bool ArapSolver::RegisterData(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces, const Eigen::VectorXi& selected,
    int max_iteration) {
  static bool registered = false;
  if (!registered) {
    registered = true;
    vertices_ = vertices;
    faces_ = faces;
    selected_ = selected;
    max_iteration_ = max_iteration;
    return true;
  }
  return false;
}

void ArapSolver::Precompute() {

}

void ArapSolver::Solve(const Eigen::MatrixXd& fixed_vertices) {

}

}  // namespace demo
}  // namespace arap
