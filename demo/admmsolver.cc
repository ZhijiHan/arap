#include "admmsolver.h"

// C++ standard library
#include <iostream>
#include <limits>
#include <vector>

#include "igl/slice.h"
#include "igl/svd3x3/polar_svd3x3.h"

namespace arap {
namespace demo {

const double kMatrixDiffThreshold = 1e-6;

AdmmSolver::AdmmSolver(double rho)
  : rho_(rho) {
}

void AdmmSolver::Precompute() {
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

  // Compute the left matrix M_.
  // TODO

  // Cholesky factorization.
  solver_.compute(M_);
  if (solver_.info() != Eigen::Success) {
    // Failed to decompose M_.
    std::cout << "Fail to do Cholesky factorization." << std::endl;
    return;
  }
}

void AdmmSolver::SolvePreprocess(const Eigen::MatrixXd& fixed_vertices) {
  // Cache fixed_vertices.
  fixed_vertices_ = fixed_vertices;
  // Check the dimension is correct.
  int fixed_num = fixed_.size();
  if (fixed_num != fixed_vertices_.rows()) {
    std::cout << "Wrong dimension in fixed number." << std::endl;
    return;
  }
  // Initialize with vertices_.
  vertices_updated_ = vertices_;
  // Initialize rotations_ with Identity matrices.
  rotations_.clear();
  int vertex_num = vertices_.rows();
  rotations_.resize(vertex_num, Eigen::Matrix3d::Identity());
  // Initialize S_ with Identity matrices.
  S_.clear();
  S_.resize(vertex_num, Eigen::Matrix3d::Identity());
  // Initialize T_ with zero matrices.
  T_.clear();
  T_.resize(vertex_num, Eigen::Matrix3d::Zero());
  // Initialize u_ with zeros.
  u_ = Eigen::MatrixXd::Zero(fixed_num, 3);
}

void AdmmSolver::SolveOneIteration() {
  int fixed_num = fixed_.size();
  int vertex_num = vertices_.rows();
  int face_num = faces_.rows();

  // The iteration contains four steps:
  // Step 1: linear solve.
  // Step 2: SVD solve.
  // Step 3: update u.
  // Step 4: update T.

  // Step 1: linear solve.
  // TODO

  // Step 2: SVD solve.
  // Input: rotations_, S_, T_.
  // Output: S_.
  for (int i = 0; i < vertex_num; ++i) {
    Eigen::Matrix3d rotation;
    igl::polar_svd3x3((rotations_[i] + T_[i]).transpose(), rotation);
    S_[i] = rotation.transpose();
  }

  // Step 3: update u.
  for (int j = 0; j < fixed_num; ++j) {
    u_.rows(j) += vertices_updated_.rows(fixed_(j)) - fixed_vertices_.row(j);
  }

  // Step 4: update T.
  for (int i = 0; i < vertex_num; ++i) {
    T_[i] += rotations_[i] - S_[i];
  }
}

Eigen::Vector3d AdmmSolver::ComputeCotangent(int face_id) const {
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

double AdmmSolver::ComputeEnergy() const {
  // Compute the energy.
  // In order to do early return, let's first test all the S_ matrices to see
  // whether they belong to SO(3).
  double infty = std:numeric_limits<double>::infinity();
  int vertex_num = vertices_.rows();
  Eigen::Matrix3d iden = Eigen::Matrix3d::Identity();
  for (int i = 0; i < vertex_num; ++i) {
    Eigen::Matrix3d si = S_[i];
    double det = si.determinant();
    if ((si * si.transpose() - iden).squaredNorm() > kMatrixDiffThreshold
      || abs(det - 1) > kMatrixDiffThreshold) {
      std::cout << "Si does not belong to SO(3) -- This should never happen.\n"
        << "i = " << i << " det = " << det << std::endl
        << "Si * Si' = \n" << si * si.transpose() << std::endl;
      return infty;
    }
  }
  // Now it passes the indicator function, the energy should be finite.
  int edge_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
  double energy = 0.0;
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
      energy += edge_energy;
    }
  }
  // Add augmented term.
  double half_rho = rho / 2;
  double rotation_aug_energy = 0.0;
  for (int i = 0; i < vertex_num; ++i) {
    rotation_aug_energy += (rotations_[i] - S_[i] + T_[i]).squaredNorm();
  }
  rotation_aug_energy *= half_rho;
  energy += rotation_aug_energy;

  int fixed_num = fixed_.size();
  double vertex_aug_energy = 0.0;
  for (int i = 0; i < fixed_num; ++i) {
    vertex_aug_energy += (vertices_updated_.row(fixed_(i))
      - fixed_vertices_.rows(i) + u_.rows(i)).squaredNorm();
  }
  vertex_aug_energy *= half_rho;
  energy += vertex_aug_energy;
  return energy;
}

}  // namespace demo
}  // namespace arap
