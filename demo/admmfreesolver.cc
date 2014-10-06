#include "admmfreesolver.h"

// C++ standard library
#include <iostream>
#include <limits>
#include <vector>

#include "igl/slice.h"
#include "igl/svd3x3/polar_svd3x3.h"

#define USE_TEST_FUNCTIONS
#define USE_LINEAR_SOLVE_1
//#define USE_LINEAR_SOLVE_2

namespace arap {
namespace demo {

const double kMatrixDiffThreshold = 1e-6;
const double kEnergyTolerance = 0.02;

AdmmFreeSolver::AdmmFreeSolver(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces, const Eigen::VectorXi& fixed,
    int max_iteration, double rho)
  : Solver(vertices, faces, fixed, max_iteration),
  rho_(rho) {
}

void AdmmFreeSolver::Precompute() {
  int vertex_num = vertices_.rows();
  int face_num = faces_.rows();
  int fixed_num = fixed_.size();

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
  // The dimension of M_ should be # of constraints by # of unknowns. i.e.,
  // (4 * vertex_num) * (4 * vertex_num).
  // The arrangement is as follows:
  // variables(columns in M_):
  // col(i): vertex position for p_i.
  // col(vertex_num + 3 * i : vertex_num + 3 * i + 2): the first to third
  // columns in R_i.
  // constraints(rows in M_):
  // row(0 : vertex_num - 1): constraints for each vertex.
  // row(vertex_num : 4 * vertex_num - 1): constraints for rotations.
  // Note that the problem can be decomposed to solve in three dimensions.

  // Here we have two methods to compute M_. The first one is to compute the
  // gradient directly. Unfortunately this is harder to program and check. So
#ifdef USE_LINEAR_SOLVE_1
  // Start of the first method: compute the gradient directly.
  M_.resize(4 * vertex_num, 4 * vertex_num);
  // Loop over all the edges.
  int edge_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
  for (int j = 0; j < fixed_num; ++j) {
    int pos = fixed_(j);
    M_.coeffRef(pos, pos) += rho_;
  }
  for (int i = vertex_num; i < 4 * vertex_num; ++i) {
    M_.coeffRef(i, i) += rho_;
  }
  for (int f = 0; f < face_num; ++f) {
    // Loop over all the edges.
    for (int e = 0; e < 3; ++e) {
      int first = faces_(f, edge_map[e][0]);
      int second = faces_(f, edge_map[e][1]);
      // Each edge is visited twice: i->j and j->i.
      // We treat first as 'i' in our constraints.
      // The i-th constraints are in the i-th and vertex_num + 3 * i :
      // vertex_num + 3 * i + 2 row in M_.
      double weight = weight_.coeff(first, second);
      M_.coeffRef(first, first) +=  2 * weight;
      M_.coeffRef(first, second) -= 2 * weight;
      M_.coeffRef(second, first) -= 2 * weight;
      M_.coeffRef(second, second) += 2 * weight;
      Eigen::Vector3d v = vertices_.row(first) - vertices_.row(second);
      for (int i = 0; i < 3; ++i) {
        M_.coeffRef(first, vertex_num + 3 * first + i) -= 2 * weight * v(i);
        M_.coeffRef(second, vertex_num + 3 * first + i) += 2 * weight * v(i);
      }
      // Rotation constraints.
      Eigen::Matrix3d m = v * v.transpose() * weight * 2;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          M_.coeffRef(vertex_num + 3 * first + i, vertex_num + 3 * first + j)
            += m(i, j);
        }
      }
      for (int i = 0; i < 3; ++i) {
        double val = weight * 2 * v(i);
        M_.coeffRef(vertex_num + 3 * first + i, first) -= val;
        M_.coeffRef(vertex_num + 3 * first + i, second) += val;
      }
    }
  }
  // End of computing the gradient directly. ***/
#endif
  // The second method looks a lot nicer: First let's consider taking the
  // derivatives of f(x) = w||Ax-b||^2:
  // f(x) = w(x'A'-b')(Ax-b) = w(x'A'Ax-2b'Ax+b'b), so
  // \nabla f(x) = w(2A'Ax-2b'A) = (2wA'A)x-2wb'A
  // So we can sum up all these terms to get the normal equation!
  // (\sum 2wA'A)x = \sum 2wA'b
#ifdef USE_LINEAR_SOLVE_2
  // Start of the second method.
  M_.resize(4 * vertex_num, 4 * vertex_num);
  // Loop over all the edges.
  int edge_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
  for (int f = 0; f < face_num; ++f) {
    // Loop over all the edges.
    for (int e = 0; e < 3; ++e) {
      int first = faces_(f, edge_map[e][0]);
      int second = faces_(f, edge_map[e][1]);
      // What is p_1 - p_2?
      Eigen::Vector3d v = vertices_.row(first) - vertices_.row(second);
      // What is the weight?
      double weight = weight_.coeff(first, second);
      // What is the dimension of A?
      // A is a one line sparse row vector!
      Eigen::SparseMatrix<double> A;
      A.resize(1, 4 * vertex_num);
      // Compute A.
      A.coeffRef(0, first) = 1;
      A.coeffRef(0, second) = -1;
      A.coeffRef(0, GetMatrixVariablePos(vertex_num, first, 0)) = -v(0);
      A.coeffRef(0, GetMatrixVariablePos(vertex_num, first, 1)) = -v(1);
      A.coeffRef(0, GetMatrixVariablePos(vertex_num, first, 2)) = -v(2);
      // Add A to M_.
      M_ = M_ + 2 * weight * A.transpose() * A;
      // What is b? b is zero!
    }
  }
  // Add the rotation constraints.
  for (int v = 0; v < vertex_num; ++v) {
    // What is the weight?
    double w = rho_ / 2;
    // What is A?
    Eigen::SparseMatrix<double> A;
    A.resize(3, 4 * vertex_num);
    A.coeffRef(0, GetMatrixVariablePos(vertex_num, v, 0)) = 1;
    A.coeffRef(1, GetMatrixVariablePos(vertex_num, v, 1)) = 1;
    A.coeffRef(2, GetMatrixVariablePos(vertex_num, v, 2)) = 1;
    // Add A to M_.
    M_ = M_ + 2 * w * A.transpose() * A;
  }
  // Add constraints for fixed vertices.
  for (int i = 0; i < fixed_num; ++i) {
    // What is th position of this fixed vertex?
    int pos = fixed_(i);
    // What is the weight w?
    double w = rho_ / 2;
    // What is A?
    Eigen::SparseMatrix<double> A;
    A.resize(1, 4 * vertex_num);
    A.coeffRef(0, pos) = 1;
    // Add A to M_.
    M_ = M_ + 2 * w * A.transpose() * A;
  }
  // End of the second method.
#endif
  // We have compared M_ from both methods, and they're the same!

  // Post-processing: compress M_, factorize it.
  M_.makeCompressed();
  // Cholesky factorization.
  solver_.compute(M_);
  if (solver_.info() != Eigen::Success) {
    // Failed to decompose M_.
    std::cout << "Fail to do Cholesky factorization." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void AdmmFreeSolver::SolvePreprocess(const Eigen::MatrixXd& fixed_vertices) {
  // Cache fixed_vertices.
  fixed_vertices_ = fixed_vertices;
  // Check the dimension is correct.
  int fixed_num = fixed_.size();
  if (fixed_num != fixed_vertices_.rows()) {
    std::cout << "Wrong dimension in fixed number." << std::endl;
    exit(EXIT_FAILURE);
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

void AdmmFreeSolver::SolveOneIteration() {
  int fixed_num = fixed_.size();
  int vertex_num = vertices_.rows();
  int face_num = faces_.rows();

  // The iteration contains four steps:
  // Step 1: linear solve.
  // Step 2: SVD solve.
  // Step 3: update u.
  // Step 4: update T.

  // Step 1: linear solve.
  // Note that the problem can be decomposed in three dimensions.
  // The number of constraints are 4 * vertex_num.
  // The first vertex_num constraints are for vertex, and the remaining
  // 3 * vertex_num constraints are for matrices.
  // Similarly, we implement two methods and cross check both of them.
#ifdef USE_LINEAR_SOLVE_1
  // Method 1: Compute the derivatives directly.
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(4 * vertex_num, 3);
  // Build rhs.
  // For vertex constraints. Since the rhs value for free vertices are 0, we
  // only consider fixed vertices.
  for (int j = 0; j < fixed_num; ++j) {
    rhs.row(fixed_(j)) = rho_ * (fixed_vertices_.row(j) - u_.row(j));
  }
  // For rotation matrix constraints.
  for (int v = 0; v < vertex_num; ++v) {
    rhs.block<3, 3>(vertex_num + 3 * v, 0)
      = rho_ * (S_[v] - T_[v]).transpose();
  }
  // End of Method 1. ***/
#endif

#ifdef USE_LINEAR_SOLVE_2
  // Method 2:
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(4 * vertex_num, 3);
  // f(x) = w||Ax-b||^2:
  // (\sum 2wA'A)x = \sum 2wA'b
  // Add the rotation constraints.
  for (int v = 0; v < vertex_num; ++v) {
    // What is the weight?
    double w = rho_ / 2;
    // What is A?
    Eigen::SparseMatrix<double> A;
    A.resize(3, 4 * vertex_num);
    A.coeffRef(0, GetMatrixVariablePos(vertex_num, v, 0)) = 1;
    A.coeffRef(1, GetMatrixVariablePos(vertex_num, v, 1)) = 1;
    A.coeffRef(2, GetMatrixVariablePos(vertex_num, v, 2)) = 1;
    // What is b?
    Eigen::Matrix3d B = (S_[v] - T_[v]).transpose();
    rhs += (2 * w * A.transpose() * B);
  }
  // Add constraints for fixed vertices.
  for (int i = 0; i < fixed_num; ++i) {
    // What is th position of this fixed vertex?
    int pos = fixed_(i);
    // What is the weight w?
    double w = rho_ / 2;
    // What is A?
    Eigen::SparseMatrix<double> A;
    A.resize(1, 4 * vertex_num);
    A.coeffRef(0, pos) = 1;
    // What is b?
    Eigen::Vector3d b = fixed_vertices_.row(i) - u_.row(i);
    rhs += (2 * w * A.transpose() * b.transpose());
  }
  // End of Method 2. ***/
#endif
  // The two methods have been double checked with each other, it turns out
  // they both give the same rhs! We are free to use either of them (Probably
  // tend to use Method 1 here because it looks shorter).

  // Solve.
  Eigen::MatrixXd solution = solver_.solve(rhs);
  if (solver_.info() != Eigen::Success) {
    std::cout << "Fail to solve the sparse linear system." << std::endl;
    exit(EXIT_FAILURE);
  }
  // Sanity check the dimension of the solution.
  if (solution.rows() != 4 * vertex_num) {
    std::cout << "Fail to write back solution: dimension mismatch."
      << std::endl;
    exit(EXIT_FAILURE);
  }
  // Sanity check the value of the solution.
  if ((M_ * solution - rhs).squaredNorm() > kMatrixDiffThreshold) {
    std::cout << "Sparse linear solver is wrong!" << std::endl;
    exit(EXIT_FAILURE);
  }
  // Write back the solutions.
  vertices_updated_ = solution.topRows(vertex_num);
  for (int v = 0; v < vertex_num; ++v) {
    rotations_[v] = solution.block<3, 3>(vertex_num + 3 * v, 0).transpose();
  }
#ifdef USE_TEST_FUNCTIONS
  // Sanity check whether it is really the optimal solution!
  if (!CheckLinearSolve()) {
    std::cout << "Iteration terminated due to the linear solve error."
      << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  // Step 2: SVD solve.
#ifdef USE_TEST_FUNCTIONS
  // Compute the energy before optimization.
  const double infty = std::numeric_limits<double>::infinity();
  const double energy_before_svd = ComputeSVDSolveEnergy();
#endif
  // Input: R, S_, T_.
  // Output: S_.
  // The solution can be found in wikipedia:
  // http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
  // R + T = U\Sigma V' S = UV'
  // The problem is (for reference):
  // \min_{S} \|S - (R + T)\|^2 s.t. S\in SO(3)
  // i.e., given R + T, find the closest SO(3) matrix.
  for (int i = 0; i < vertex_num; ++i) {
    Eigen::Matrix3d rotation;
    Eigen::Matrix3d res = rotations_[i] + T_[i];
    igl::polar_svd3x3(res, rotation);
    S_[i] = rotation;
  }
#ifdef USE_TEST_FUNCTIONS
  const double energy_after_svd = ComputeSVDSolveEnergy();
  // Test whether the energy really decreases.
  if (energy_before_svd == infty || energy_after_svd == infty
      || energy_before_svd < energy_after_svd - kEnergyTolerance) {
    std::cout << "Iteration terminated due to the svd solve error."
      << std::endl << "Energy before SVD: " << energy_before_svd
      << std::endl << "Energy after SVD: " << energy_after_svd
      << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  // Step 3: update u.
  for (int j = 0; j < fixed_num; ++j) {
    u_.row(j) += vertices_updated_.row(fixed_(j)) - fixed_vertices_.row(j);
  }

  // Step 4: update T.
  for (int i = 0; i < vertex_num; ++i) {
    T_[i] += rotations_[i] - S_[i];
  }
}

Eigen::Vector3d AdmmFreeSolver::ComputeCotangent(int face_id) const {
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

Energy AdmmFreeSolver::ComputeEnergy() const {
  // Compute the energy.
  Energy energy;
  // In order to do early return, let's first test all the S_ matrices to see
  // whether they belong to SO(3).
  double infty = std::numeric_limits<double>::infinity();
  int vertex_num = vertices_.rows();
  for (int i = 0; i < vertex_num; ++i) {
    if (!IsSO3(S_[i])) {
      std::cout << "This should never happen!" << std::endl;
      energy.AddEnergyType("Total", infty);
      return energy;
    }
  }

  // Now it passes the indicator function, the energy should be finite.
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
  energy.AddEnergyType("ARAP", total);

  // Add augmented term.
  double half_rho = rho_ / 2;
  double rotation_aug_energy = 0.0;
  for (int i = 0; i < vertex_num; ++i) {
    rotation_aug_energy += (rotations_[i] - S_[i]).squaredNorm();
  }
  rotation_aug_energy *= half_rho;
  total += rotation_aug_energy;
  energy.AddEnergyType("Rotation", rotation_aug_energy);

  int fixed_num = fixed_.size();
  double vertex_aug_energy = 0.0;
  for (int i = 0; i < fixed_num; ++i) {
    vertex_aug_energy += (vertices_updated_.row(fixed_(i))
      - fixed_vertices_.row(i)).squaredNorm();
  }
  vertex_aug_energy *= half_rho;
  total += vertex_aug_energy;
  energy.AddEnergyType("Vertex", vertex_aug_energy);
  energy.AddEnergyType("Total", total);
  return energy;
}

bool AdmmFreeSolver::CheckLinearSolve() const {
  // Compute the linear solve energy.
  // Don't pollute the solution! Play with a copy instead.
  Eigen::MatrixXd vertices = vertices_updated_;
  std::vector<Eigen::Matrix3d> R = rotations_;
  double optimal_energy = ComputeLinearSolveEnergy(vertices, R);
  std::cout << "Optimal linear energy: " << optimal_energy << std::endl;
  std::cout.precision(15);
  int rows = vertices.rows();
  int cols = vertices.cols();
  double delta = 0.001;
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      // Perturb solution(i, j) a little bit.
      vertices(i, j) += delta;
      double perturbed_enrgy = ComputeLinearSolveEnergy(vertices, R);
      double product = perturbed_enrgy - optimal_energy;
      // Reset value.
      vertices(i, j) = vertices_updated_(i, j);
      // Perturb in another direction.
      vertices(i, j) -= delta;
      perturbed_enrgy = ComputeLinearSolveEnergy(vertices, R);
      perturbed_enrgy *= (perturbed_enrgy - optimal_energy);
      if (product < 0.0) {
        std::cout << "Linear solve check failed!" << std::endl;
        std::cout << "Error occurs in (" << i << ", " << j << ")" << std::endl;
        return false;
      }
      // Reset value.
      vertices(i, j) = vertices_updated_(i, j);
    }
  }
  // Perturb the rotations.
  int vertex_num = vertices_.rows();
  for (int v = 0; v < vertex_num; ++v) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        // Perturb R[v](i, j).
        R[v](i, j) += delta;
        double perturbed_enrgy = ComputeLinearSolveEnergy(vertices, R);
        double product = perturbed_enrgy - optimal_energy;
        R[v](i, j) = rotations_[v](i, j);
        R[v](i, j) -= delta;
        perturbed_enrgy = ComputeLinearSolveEnergy(vertices, R);
        product *= (perturbed_enrgy - optimal_energy);
        if (product < 0.0) {
          std::cout << "Linear solve check failed!" << std::endl;
          std::cout << "Error occurs in (" << v << ", " << i << ", " << j << ")" << std::endl;
          // Print out the parabola.
          std::cout << "Variable\tEnergy\n";
          for (int step = -10; step < 10; ++step) {
            R[v](i, j) = rotations_[v](i, j) + step * delta;
            std::cout << R[v](i, j) << "\t" << ComputeLinearSolveEnergy(vertices, R)
              << std::endl;
          }
          return false;
        }
        R[v](i, j) = rotations_[v](i, j);
      }
    }
  }
  std::cout << "All linear solve tests passed!" << std::endl;
  return true;
}

double AdmmFreeSolver::ComputeLinearSolveEnergy(const Eigen::MatrixXd &vertices,
    const std::vector<Eigen::Matrix3d> &rotations) const {
  // Compute the linear solve energy.
  double energy = 0.0;
  int edge_map[3][2] = { {1, 2}, {2, 0}, {0, 1} };
  int face_num = faces_.rows();
  int vertex_num = vertices_.rows();
  for (int f = 0; f < face_num; ++f) {
    // Loop over all the edges.
    for (int e = 0; e < 3; ++e) {
      int first = faces_(f, edge_map[e][0]);
      int second = faces_(f, edge_map[e][1]);
      double edge_energy = 0.0;
      double weight = weight_.coeff(first, second);
      Eigen::Vector3d vec = (vertices.row(first) -
          vertices.row(second)).transpose() -
          rotations[first] * (vertices_.row(first) -
          vertices_.row(second)).transpose();
      edge_energy = weight * vec.squaredNorm();
      energy += edge_energy;
    }
  }
  // Add up the augmented rotation term.
  double weight = rho_ / 2;
  for (int v = 0; v < vertex_num; ++v) {
    energy += weight * (rotations[v] - S_[v] + T_[v]).squaredNorm();
  }
  // Add up the fixed vertices term.
  int fixed_num = fixed_.size();
  for (int i = 0; i < fixed_num; ++i) {
    int pos = fixed_(i);
    energy += weight * (vertices.row(pos) - fixed_vertices_.row(i)
        + u_.row(i)).squaredNorm();
  }
  return energy;
}

// Compute the SVD solve energy. Used in CheckSVDSolve.
double AdmmFreeSolver::ComputeSVDSolveEnergy() const {
  double infty = std::numeric_limits<double>::infinity();
  int vertex_num = vertices_.rows();
  for (int v = 0; v < vertex_num; ++v) {
    if (!IsSO3(S_[v])) {
      return infty;
    }
  }
  double energy = 0.0;
  for (int v = 0; v < vertex_num; ++v) {
    energy += (rotations_[v] - S_[v] + T_[v]).squaredNorm();
  }
  energy *= rho_ / 2;
  return energy;
}

// Check whether a matrix is in SO(3).
bool AdmmFreeSolver::IsSO3(const Eigen::Matrix3d &S) const {
  double det = S.determinant();
  if ((S * S.transpose() - Eigen::Matrix3d::Identity()).squaredNorm()
      > kMatrixDiffThreshold || abs(det - 1) > kMatrixDiffThreshold) {
    std::cout << "S does not belong to SO(3)" << std::endl;
    std::cout << "S: \n" << S << std::endl;
    return false;
  }
  return true;
}

}  // namespace demo
}  // namespace arap
