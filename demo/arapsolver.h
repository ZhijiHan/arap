#ifndef _ARAP_DEMO_ARAPSOLVER_H_
#define _ARAP_DEMO_ARAPSOLVER_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace arap {
namespace demo {

class ArapSolver {
 public:
  ArapSolver();

  // Register data for ArapSolver. Once vertices and faces are set then cannot
  // be reset from outside. The only way to reset vertices is to call Solve,
  // which will run ARAP algorithm based on the current state, and update
  // vertices_.
  bool RegisterData(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
      const Eigen::VectorXi& selected, int max_iteration);

  // Pre compute cotangent weight and left matrix in ARAP.
  void Precompute();

  // Solve ARAP problem. |fixed_vertices| is a # of fixed vertices by 3 matrix.
  // Each row in |fixed_vertices| represents the coordinates of a fixed vertex.
  void Solve(const Eigen::MatrixXd& fixed_vertices);

  // Get vertices_.
  const Eigen::MatrixXd& GetVertices() const { return vertices_; }

  // Get faces_.
  const Eigen::MatrixXi& GetFaces() const { return faces_; }

  // Get indices of the selected vertices.
  const Eigen::VectorXi& GetSelectedIndices() const { return selected_; }

 private:
  Eigen::MatrixXd vertices_;
  Eigen::MatrixXi faces_;
  Eigen::VectorXi selected_;

  // Max number of iterations used to solve ARAP.
  int max_iteration_;
};

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ARAPSOLVER_H_
