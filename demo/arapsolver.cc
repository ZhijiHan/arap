#include "arapsolver.h"

// C++ standard library
#include <iostream>
#include <vector>

#include "igl/slice.h"
#include "igl/svd3x3/polar_svd3x3.h"

namespace arap {
namespace demo {

ArapSolver::ArapSolver(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces, const Eigen::VectorXi& fixed,
    int max_iteration)
  : Solver(vertices, faces, fixed, max_iteration) {
}

void ArapSolver::Precompute() {
  int vertex_num = vertices_.rows();
  int face_num = faces_.rows();

  // Compute weight_.
  weight_.resize(vertex_num, vertex_num);
  // An index map to help mapping from one vertex to the corresponding edge.
  int index_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
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
      weight_.coeffRef(first, second) += half_cot;
      weight_.coeffRef(second, first) += half_cot;
      // Note that weight_(i, i) is the sum of all the -weight_(i, j).
      weight_.coeffRef(first, first) -= half_cot;
      weight_.coeffRef(second, second) -= half_cot;
    }
  }

  // Compute lb_operator_. This matrix can be computed by extracting free_ rows
  // and columns from -weight_.
  igl::slice(weight_, free_, free_, lb_operator_);
  lb_operator_ *= -1.0;

  // Cholesky factorization.
  solver_.compute(lb_operator_);
  if (solver_.info() != Eigen::Success) {
    // Failed to decompose lb_operator_.
    std::cout << "Fail to do Cholesky factorization." << std::endl;
    return;
  }
}

void ArapSolver::SolvePreprocess(const Eigen::MatrixXd& fixed_vertices) {
  // Initialize fixed_vertices_.
  fixed_vertices_ = fixed_vertices;
  // Initialized with vertices_. Note that this is different from the default
  // setting in the ARAP demo, which uses the value from last from for
  // initialization.
  vertices_updated_ = vertices_;
  int fixed_num = fixed_.size();
  for (int i = 0; i < fixed_num; ++i) {
    vertices_updated_.row(fixed_(i)) = fixed_vertices_.row(i);
  }
  // Initialize rotations_ with Identity matrices.
  rotations_.clear();
  int vertex_num = vertices_.rows();
  rotations_.resize(vertex_num, Eigen::Matrix3d::Identity());
}

void ArapSolver::SolveOneIteration() {
  int fixed_num = fixed_.size();
  int vertex_num = vertices_.rows();
  int face_num = faces_.rows();

  // A temporary vector to hold all the edge products for all the vertices.
  // This is the S matrix in equation (5).
  std::vector<Eigen::Matrix3d> edge_product;
  // edge_map is used for looping over three edges in a triangle.
  int edge_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
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
          = vertices_updated_.row(first) - vertices_updated_.row(second);
      double weight = weight_.coeff(first, second);
      edge_product[first] += weight * edge * edge_update.transpose();
    }
  }
  for (int v = 0; v < vertex_num; ++v) {
    Eigen::Matrix3d rotation;
    igl::polar_svd3x3(edge_product[v], rotation);
    rotations_[v] = rotation.transpose();
  }
  // Step 2: compute the rhs in equation (9).
  int free_num = free_.size();
  // The right hand side of equation (9). The x, y and z coordinates are
  // computed separately.
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(free_num, 3);
  for (int f = 0; f < face_num; ++f) {
    // Loop over all the edges in the mesh.
    for (int e = 0; e < 3; ++e) {
      int first = faces_(f, edge_map[e][0]);
      int second = faces_(f, edge_map[e][1]);
      // If first is fixed, skip it.
      if (vertex_info_[first].type == VertexType::Fixed)
        continue;
      int vertex_id = vertex_info_[first].pos;
      double weight = weight_.coeff(first, second);
      Eigen::Vector3d vec = weight / 2.0 *
          (rotations_[first] + rotations_[second]) *
          (vertices_.row(first) - vertices_.row(second)).transpose();
      rhs.row(vertex_id) += vec.transpose();
      if (vertex_info_[second].type == VertexType::Fixed) {
        // If second is fixed, add another term in the right hand side.
        rhs.row(vertex_id) += weight * vertices_updated_.row(second);
      }
    }
  }
  // Solve for free_.
  Eigen::VectorXd solution;
  for (int i = 0; i < 3; ++i) {
    solution = solver_.solve(rhs.col(i));
    if (solver_.info() != Eigen::Success) {
      std::cout << "Fail to solve the sparse linear system." << std::endl;
      return;
    }
    // Sanity check the dimension of solution.
    if (solution.size() != free_num) {
      std::cout << "Fail to write back solution: dimension mismatch."
                << std::endl;
      return;
    }
    // Write back to vertices_updated_.
    for (int j = 0; j < free_num; ++j) {
      int vertex_id = free_(j);
      vertices_updated_(vertex_id, i) = solution(j);
    }
  }
}

Eigen::Vector3d ArapSolver::ComputeCotangent(int face_id) const {
  Eigen::Vector3d cotangent(0.0, 0.0, 0.0);
  // The triangle is defined as follows:
  //            A
  //           /  -
  //        c /     - b
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

Energy ArapSolver::ComputeEnergy() const {
  // Compute the energy.
  int edge_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
  double total = 0.0;
  int face_num = faces_.rows();
  for (int f = 0; f < face_num; ++f) {
    // Loop over all the edges.
    for (int e = 0; e < 3; ++e) {
      int first = faces_(f, edge_map[e][0]);
      int second = faces_(f, edge_map[e][1]);
      double edge_energy = 0.0;
      double weight = weight_.coeff(first, second);
      Eigen::Vector3d vec = (vertices_updated_.row(first) -
          vertices_updated_.row(second)).transpose() -
          rotations_[first] * (vertices_.row(first) -
          vertices_.row(second)).transpose();
      edge_energy = weight * vec.squaredNorm();
      total += edge_energy;
    }
  }
  Energy energy;
  energy.AddEnergyType("Total", total);
  return energy;
}

}  // namespace demo
}  // namespace arap
