#include "arapsolver.h"

// C++ standard library
#include <iostream>
#include <vector>

#include "igl/slice.h"
#include "igl/svd3x3/polar_svd3x3.h"

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
  int vertex_num = vertices_.rows();
  int face_num = faces_.rows();

  // Compute free_.
  int selected_num = selected_.size();
  int free_num = vertex_num - selected_num;
  free_.resize(free_num);
  int j = 0, k = 0;
  for (int i = 0; i < vertex_num; ++i) {
    if (j < selected_num && i == selected_(j)) {
      ++j;
    } else {
      free_(k) = i;
      ++k;
    }
  }
  // Sanity check the sizes of selected_ and free_ are correct.
  if (j != selected_num || k != free_num) {
    std::cout << "Fail to compute free_ in ArapSolver: dimension mismatch."
              << std::endl;
    return;
  }

  // Compute cot_weight_.
  cot_weight_.resize(vertex_num, vertex_num);
  // An index map to help mapping from one vertex to the corresponding edge.
  int index_map[3][2] = { {1, 2}, {0, 2}, {0, 1} };
  // Loop over all the faces.
  for (int f = 0; f < face_num; ++f) {
    // Get the cotangent value with that face.
    Eigen::Vector3d cotangent = ComputeCotangent(f);
    // Loop over the three vertices within the same triangle.
    // i = 0 => A.
    // i = 1 => B.
    // i = 2 => C.
    for (int i = 0; i < 3; ++i) {
      // Indices of the two vertices in the edge corresponding to vertex i.
      int first = faces_(f, index_map[i][0]);
      int second = faces_(f, index_map[i][1]);
      double half_cot = cotangent(i) / 2.0;
      cot_weight_.coeffRef(first, second) += half_cot;
      cot_weight_.coeffRef(second, first) += half_cot;
      // Note that cot_weight_(i, i) is the sum of all the -cot_weight_(i, j).
      cot_weight_.coeffRef(first, first) -= half_cot;
      cot_weight_.coeffRef(second, second) -= half_cot;
    }
  }

  // Compute lb_operator_. This matrix can be computed by extracting free_ rows
  // and columns from -cot_weight_.
  igl::slice(cot_weight_, free_, free_, lb_operator_);
  lb_operator_ *= -1.0;

  // Cholesky factorization.
  solver_.compute(lb_operator_);
  if (solver_.info() != Eigen::Success) {
    // Failed to decompose lb_operator_.
    std::cout << "Fail to do Cholesky factorization." << std::endl;
    return;
  }
}

void ArapSolver::Solve(const Eigen::MatrixXd& fixed_vertices) {
  // The optimization goes alternatively between solving vertices and
  // rotations.
  // Step 0: replace some vertices in vertices_ with fixed_vertices.
  // Step 1: given vertices_, solve the rotations for all the vertices by
  // polar_svd.
  // Step 2: update the rhs in equation (9) with updated rotations.
  // Step 3: repeat Step 1 and Step 2 for max_iteration times.

  // Step 0: replace vertices with fixed_vertices.
  int fixed_num = selected_.size();
  // Sanity check for the # of vertices in fixed_vertices.
  if (fixed_vertices.rows() != fixed_num) {
    std::cout << "Fail to solve: number of fixed vertices mismatch."
              << std::endl;
    return;
  }
  Eigen::MatrixXd vertices_update = vertices_;
  for (int i = 0; i < fixed_num; ++i) {
    vertices_update.row(selected_(i)) = fixed_vertices.row(i);
  }

  int iter = 0;
  int vertex_num = vertices_.rows();
  int face_num = faces_.rows();
  // Initialize rotations_ with Identity matrices.
  rotations_.clear();
  rotations_.resize(vertex_num, Eigen::Matrix3d::Identity());
  // A temporary vector to hold all the edge products for all the vertices.
  // This is the S matrix in equation (5).
  std::vector<Eigen::Matrix3d> edge_product;
  // edge_map is used for looping over three edges in a triangle.
  int edge_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
  while (iter < max_iteration_) {
    // Step 1: solve rotations by polar_svd.
    // Clear edge_product.
    edge_product.clear();
    edge_product.resize(vertex_num, Eigen::Matrix3d::Zero());
    for (int f = 0; f < face_num; ++f) {
      // Loop over all the edges in the mesh.
      for (int e = 0; e < 3; ++e) {
        int first = faces_(f, edge_map[e][0]);
        int second = faces_(f, edge_map[e][1]);
        // Now we have got an edge from first to second.
        Eigen::Vector3d edge = vertices_.row(first) - vertices_.row(second);
        Eigen::Vector3d edge_update
            = vertices_update.row(first) - vertices_update.row(second);
        double weight = cot_weight_.coeff(first, second);
        edge_product[first] += weight * edge * edge_update.transpose();
      }
    }
    for (int v = 0; v < vertex_num; ++v) {
      Eigen::Matrix3d rotation;
      igl::polar_svd3x3(edge_product[v], rotation);
      rotations_[v] = rotation.transpose();
    }

    // Step 2:

    // Update vertices_.
    vertices_ = vertices_update;

    // Increment.
    ++iter;
  }
}

Eigen::Vector3d ArapSolver::ComputeCotangent(int face_id) const {
  Eigen::Vector3d cotangent(0.0, 0.0, 0.0);
  // The triangle is defined as follows:
  //            A
  //           /  -
  //        c /    - b
  //         /        -
  //        /    a      -
  //       B--------------C
  // where A, B, C corresponds to faces_(face_id, 0), faces_(face_id, 1) and
  // faces_(face_id, 2). The return value is (cotA, cotB, cotC).
  // Compute the triangle area first.
  Eigen::Vector3d A = vertices_.row(faces_(face_id, 0));
  Eigen::Vector3d B = vertices_.row(faces_(face_id, 1));
  Eigen::Vector3d C = vertices_.row(faces_(face_id, 2));
  double a_squared = (B - C).squaredNorm();
  double b_squared = (C - A).squaredNorm();
  double c_squared = (A - B).squaredNorm();
  // Compute the area of the triangle. area = 1/2bcsinA.
  double area = (B - A).cross(C - A).norm() / 2;
  // Compute cotA = cosA / sinA.
  // b^2 + c^2 -2bccosA = a^2, or cosA = (b^2 + c^2 - a^2) / 2bc.
  // 1/2bcsinA = area, or sinA = 2area/bc.
  // cotA = (b^2 + c^2 -a^2) / 4area.
  double four_area = 4 * area;
  cotangent(0) = (b_squared + c_squared - a_squared) / four_area;
  cotangent(1) = (c_squared + a_squared - b_squared) / four_area;
  cotangent(2) = (a_squared + b_squared - c_squared) / four_area;
  return cotangent;
}

}  // namespace demo
}  // namespace arap
