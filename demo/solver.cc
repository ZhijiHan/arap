#include "solver.h"

// C++ standard library
#include <iostream>
#include <vector>

#include "igl/slice.h"
#include "igl/svd3x3/polar_svd3x3.h"

namespace arap {
namespace demo {

Solver::Solver() {
}

// If it is the first time to register data, return true; else return false.
bool Solver::RegisterData(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces, const Eigen::VectorXi& fixed,
    int max_iteration) {
  static bool registered = false;
  if (!registered) {
    registered = true;
    vertices_ = vertices;
    faces_ = faces;
    fixed_ = fixed;
    max_iteration_ = max_iteration;
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
      return false;
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
            return false;
          }
          break;
        case Free:
          if (free_(info.pos) != i) {
            std::cout << "Fail to test vertex info: wrong free position."
                      << std::endl;
            return false;
          }
          break;
        default:
          std::cout << "Unknown vertex type." << std::endl;
          return false;
      }
    }
    return true;
  }
  return false;
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
}

}  // namespace demo
}  // namespace arap
