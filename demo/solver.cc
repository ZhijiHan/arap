#include "solver.h"

// C++ standard library
#include <iostream>
#include <vector>

#include "igl/slice.h"
#include "igl/svd3x3/polar_svd3x3.h"

namespace arap {
namespace demo {

Solver::Solver(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces, const Eigen::VectorXi& fixed,
    int max_iteration)
  : vertices_(vertices), faces_(faces), fixed_(fixed),
  max_iteration_(max_iteration) {
  // Compute free_.
  int vertex_num = vertices_.rows();
  int fixed_num = fixed_.size();
  int free_num = vertex_num - fixed_num;
  free_.resize(free_num);
  int j = 0, k = 0;
  for (int i = 0; i < vertex_num; ++i) {
    if (j < fixed_num && i == fixed_(j)) {
      ++j;
    } else {
      free_(k) = i;
      ++k;
    }
  }
  // Sanity check the sizes of fixed_ and free_ are correct.
  if (j != fixed_num || k != free_num) {
    std::cout << "Fail to compute free_ in Solver: dimension mismatch."
              << std::endl;
    return;
  }
  // Compute information for each vertex.
  vertex_info_.resize(vertex_num);
  // Compute information for each fixed vertex.
  for (int i = 0; i < fixed_num; ++i) {
    int vertex_id = fixed_(i);
    VertexInfo info(VertexType::Fixed, i);
    vertex_info_[vertex_id] = info;
  }
  // Compute information for each free vertex.
  for (int i = 0; i < free_num; ++i) {
    int vertex_id = free_(i);
    VertexInfo info(VertexType::Free, i);
    vertex_info_[vertex_id] = info;
  }
  // Sanity check for all the vertices.
  for (int i = 0; i < vertex_num; ++i) {
    VertexInfo info = vertex_info_[i];
    switch (info.type) {
      case Fixed:
        if (fixed_(info.pos) != i) {
          std::cout << "Fail to test vertex info: wrong fixed position."
                    << std::endl;
          return;
        }
        break;
      case Free:
        if (free_(info.pos) != i) {
          std::cout << "Fail to test vertex info: wrong free position."
                    << std::endl;
          return;
        }
        break;
      default:
        std::cout << "Unknown vertex type." << std::endl;
        return;
    }
  }
}

void Solver::Solve(const Eigen::MatrixXd& fixed_vertices) {
  // The optimization goes alternatively between solving vertices and
  // rotations.
  SolvePreprocess(fixed_vertices);
  int iter = 0;
  while (iter < max_iteration_) {
    SolveOneIteration();
    ++iter;
  }
  SolvePostprocess();
}

// RefineRotations will do SVD projection based on given vertices_,
// vertices_updated_ and weight_. It will overwrite rotations_. After calling
// this, the rotation constraints are enforced.
// The typical use can be:
// step 0: call ComputeEnergy to get current energy;
// step 1: cache the current rotations_;
// step 2: call RefineRotations;
// step 3: call ComputeEnergy to update rotations_;
// step 4: if the updated energy is larger, reset rotations_.
void Solver::RefineRotations() {
  int vertex_num = vertices_.rows();
  std::vector<Eigen::Matrix3d> edge_product(vertex_num, Eigen::Matrix3d::Zero());
  for (int i = 0; i < vertex_num; ++i) {
    for (auto& neighbor : neighbors_[i]) {
      int j = neighbor.first;
      double weight = weight_.coeff(i, j);
      Eigen::Vector3d edge = vertices_.row(i) - vertices_.row(j);
      Eigen::Vector3d edge_update =
        vertices_updated_.row(i) - vertices_updated_.row(j);
      edge_product[i] += weight * edge * edge_update.transpose();
    }
  }
  for (int v = 0; v < vertex_num; ++v) {
    Eigen::Matrix3d rotation;
    igl::polar_svd3x3(edge_product[v], rotation);
    rotations_[v] = rotation.transpose();
  }
}

// RefineVertices will update vertices_updated_ based on given rotations_,
// fixed_vertices_, vertices_ and weight_. It will overwrite
// vertices_updated_. After calling this, the fixed vertex constraints are
// enforced.
void Solver::RefineVertices() {
  int fixed_num = fixed_.size();
  // Write back fixed vertices constraints.
  for (int i = 0; i < fixed_num; ++i) {
    vertices_updated_.row(fixed_(i)) = fixed_vertices_.row(i);
  }
  Eigen::SparseMatrix<double> lb_operator;
  igl::slice(weight_, free_, free_, lb_operator);
  lb_operator *= -1.0;
  lb_operator.makeCompressed();

  // LU factorization.
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(lb_operator);
  if (solver.info() != Eigen::Success) {
    // Failed to decompose lb_operator_.
    std::cout << "Fail to do LU factorization." << std::endl;
    exit(EXIT_FAILURE);
  }

  int free_num = free_.size();
  // The right hand side of equation (9). The x, y and z coordinates are
  // computed separately.
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(free_num, 3);
  for (int i = 0; i < free_num; ++i) {
    int i_pos = free_(i);
    for (auto& neighbor : neighbors_[i_pos]) {
      int j_pos = neighbor.first;
      double weight = weight_.coeff(i_pos, j_pos);
      Eigen::Vector3d vec = weight / 2.0
        * (rotations_[i_pos] + rotations_[j_pos])
        * (vertices_.row(i_pos) - vertices_.row(j_pos)).transpose();
      rhs.row(i) += vec;
      if (vertex_info_[j_pos].type == VertexType::Fixed) {
        rhs.row(i) += weight * vertices_updated_.row(j_pos);
      }
    }
  }
  // Solve for free_.
  Eigen::MatrixXd solution = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    std::cout << "Fail to solve the sparse linear system." << std::endl;
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < free_num; ++i) {
    int pos = free_(i);
    vertices_updated_.row(pos) = solution.row(i);
  }
}

// The definition of ARAP energy:
// E = \sum_i \sum_j\in{neighbors_[i]} weight_(i, j)||p'_i - p'_j - R_i(p_i - p_j)||^2

// The four functions below have not been tested. It is highly recommended NOT
// to use them before full testing has been done!

// Given an index in p' and the dimension we are interested in (x, y, or z),
// returns the gradient w.r.t. it.
double Solver::ComputePositionGradient(int point_index, int dim_index) const {
  if (weight_.rows() == 0|| weight_.cols() == 0 || neighbors_.size() == 0) {
    std::cout << "Uninitialized weight_ or neighbors_." << std::endl;
    exit(EXIT_FAILURE);
  }
  double gradient = 0.0;
  // Assume neighborhood relationship is symmetric.
  int vertex_num = vertices_.rows();
  int i = point_index;
  for (auto& neighbor : neighbors_[point_index]) {
    int j = neighbor.first;
    double weight = weight_.coeff(i, j);
    Eigen::Vector3d e = vertices_.row(i) - vertices_.row(j);
    Eigen::Vector3d e_update = vertices_updated_.row(i) -
        vertices_updated_.row(j); 
    // Consider term weight * ||p'_i - p'_j - R_i(p_i - p_j)||^2
    Eigen::Vector3d v = rotations_[i] * e; 
    gradient += 2 * weight * (e_update(dim_index) - v(dim_index));
    // Then, consider weight * ||p'_j - p'_i - R_j(p_j - p_i)||^2
    weight = weight_.coeff(j, i);
    v = rotations_[j] * (-e);
    gradient += 2 * weight * (e_update(dim_index) + v(dim_index));
  }
  return gradient;
}

// Same above, but returns a vertex_num by 3 matrix for the gradients for all
// the points and all the dimensions.
Eigen::MatrixXd Solver::ComputePositionGradient() const {
  if (weight_.rows() == 0|| weight_.cols() == 0 || neighbors_.size() == 0) {
    std::cout << "Uninitialized weight_ or neighbors_." << std::endl;
    exit(EXIT_FAILURE);
  }
  int vertex_num = vertices_.rows();
  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(vertex_num, 3);
  for (int i = 0; i < vertex_num; ++i) {
    for (auto& neighbor : neighbors_[i]) {
      int j = neighbor.first;
      double weight = weight_.coeff(i, j);
      Eigen::Vector3d e = vertices_.row(i) - vertices_.row(j);
      Eigen::Vector3d e_update = vertices_updated_.row(i) -
          vertices_updated_.row(j);
      Eigen::Vector3d v = 2 * weight * (e_update - rotations_[i] * e);
      gradient.row(i) += v;
      gradient.row(j) -= v;
    }
  }
  return gradient; 
}

// Given an index of the rotation, the row and column we are interested in,
// returns the gradient w.r.t. it.
double Solver::ComputeRotationGradient(int rotation_index,
    int row_index, int col_index) const {
  double gradient = 0.0;
  int i = rotation_index;
  for (auto& neighbor : neighbors_[i]) {
    int j = neighbor.first;
    double weight = weight_.coeff(i, j);
    Eigen::Vector3d e = vertices_.row(i) - vertices_.row(j);
    Eigen::Vector3d e_update = vertices_updated_.row(i) -
        vertices_updated_.row(j);
    Eigen::Vector3d v = 2 * weight * (e_update - rotations_[i] * e);
    gradient += v(row_index) * (-e(col_index)); 
  }
  return gradient; 
}

// Same above, but returns a (vertex_num x 3) by 3 matrix for all the
// gradient.
Eigen::MatrixXd Solver::ComputeRotationGradient() const {
  int vertex_num = vertices_.rows();
  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(vertex_num * 3, 3);
  for (int i = 0; i < vertex_num; ++i) {
    for (auto& neighbor : neighbors_[i]) {
      int j = neighbor.first;
      double weight = weight_.coeff(i, j);
      Eigen::Vector3d e = vertices_.row(i) - vertices_.row(j);
      Eigen::Vector3d e_update = vertices_updated_.row(i) -
          vertices_updated_.row(j);
      gradient.block<3, 3>(3 * i, 0) += 2 * weight *
          (e_update - rotations_[i] * e) * (-e.transpose());
    }
  }
  return gradient;
}

}  // namespace demo
}  // namespace arap
