#ifndef _ARAP_DEMO_ARAPSOLVER_H_
#define _ARAP_DEMO_ARAPSOLVER_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace arap {
namespace demo {

class ArapSolver {
 public:
  ArapSolver();

  // Registers data for ArapSolver. Once vertices and faces are set then cannot
  // be reset from outside. The only way to reset vertices is to call Solve,
  // which will run ARAP algorithm based on the current state, and update
  // vertices_.
  bool RegisterData(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
      const Eigen::VectorXi& selected, int max_iteration);

  // Pre computes cotangent weight and left matrix in ARAP.
  void Precompute();

  // Solves ARAP problem. |fixed_vertices| is a # of fixed vertices by 3 matrix.
  // Each row in |fixed_vertices| represents the coordinates of a fixed vertex.
  void Solve(const Eigen::MatrixXd& fixed_vertices);

  // Gets vertices_.
  const Eigen::MatrixXd& GetVertices() const { return vertices_; }

  // Gets faces_.
  const Eigen::MatrixXi& GetFaces() const { return faces_; }

  // Gets indices of the selected vertices.
  const Eigen::VectorXi& GetSelectedIndices() const { return selected_; }

 private:
  // Computes the cotangent angle in one face indicated by |face_id|.
  Eigen::Vector3d ComputeCotangent(int face_id) const;

  // vertices_ is a # of vertices by 3 matrix. Each row in vertices_ represents
  // a vertex's position.
  Eigen::MatrixXd vertices_;
  // faces_ is a # of faces by 3 matrix. Each row contains the indices of this
  // face's three vertices.
  Eigen::MatrixXi faces_;
  // selected_ is a vector no longer than the # of vertices. It contains the
  // indices of the vertices that we want to fix during the deformation.
  Eigen::VectorXi selected_;

  // This is a sparse matrix used to store the cotangent weight. The # of rows
  // and columns equals to the # of vertices, i.e., vertices_.rows(). For any i
  // and j, the definition of cot_weight_(i, j) is as follows:
  // If i != j and (i, j) is not an edge, cot_weight_(i, j) = 0.
  // If i != j and (i, j) is an edge, then cot_weight_(i, j) equals to:
  //   cot_weight_(i, j) = 1/2(cot alpha_ij + cot beta_ij).
  // Note that cot_weight_(i, j) = cot_weight_(j, i), or cot_weight_ is
  // symmetric.
  // If i == j, then cot_weight_(i, i) is defined as the sum of all the w_ij.
  Eigen::SparseMatrix<double> cot_weight_;

  // Max number of iterations used to solve ARAP.
  int max_iteration_;
};

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ARAPSOLVER_H_
