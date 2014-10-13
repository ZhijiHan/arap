#include "adaptadmmfixedsolver.h"
#include "adaptadmmfreesolver.h"
#include "admmfixedsolver.h"
#include "admmfreesolver.h"
#include "arapbenchmarksolver.h"
#include "arapsolver.h"

// C++ standard library
#include <algorithm>
#include <iostream>
#include <vector>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <igl/colon.h>
#include <igl/deform_skeleton.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/dqs.h>
#include <igl/forward_kinematics.h>
#include <igl/lbs_matrix.h>
#include <igl/PI.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/svd3x3/arap.h>
#include <igl/viewer/Viewer.h>

// Vertex matrix. V is the original vertices from .off file, and U is the
// vertices updated in each frame.
Eigen::MatrixXd V, U;
// Face matrix. F is read from .off file.
Eigen::MatrixXi F;
// Fixed vertices.
Eigen::MatrixXd bc;
// Color matrix used to display selected points from S and b below.
Eigen::MatrixXd C;
// S is a column vector representing vertices with predefined coordinates in
// each frame. S is read from .dmat file, whose file format can be found in
// http://igl.ethz.ch/projects/libigl/file-formats/dmat.html. The dimension of
// S is the same as the # of vertices. b contains the indices of those vertices
// whose S is nonnegative.
Eigen::VectorXi S, b;
double anim_t = 0.0;
double anim_t_dir = 0.03;
// Our own implementation of ARAP.
arap::demo::Solver* solver = nullptr;

// Color used to draw precomputed vertices.
static const Eigen::RowVector3d kPurple(80.0 / 255.0,
                                        64.0 / 255.0,
                                        255.0 / 255.0);
static const Eigen::RowVector3d kGold(255.0 / 255.0,
                                      228.0 / 255.0,
                                      58.0 / 255.0);

bool pre_draw(igl::Viewer& viewer) {
  static int iteration = 0;
  if (!viewer.core.is_animating
    || iteration >= solver->GetMaxIteration())
    return false;
  solver->SolveOneIteration();
  arap::demo::Energy energy = solver->ComputeEnergy();
  std::cout << "Iteration: " << iteration << " Energy: " << energy;
  Eigen::MatrixXd solution = solver->GetVertexSolution();
  viewer.data.set_vertices(solution);
  viewer.data.set_points(bc, C);
  viewer.data.compute_normals();
  ++iteration;
  return false;
}

bool key_down(igl::Viewer& viewer, unsigned char key, int mods) {
  switch (key) {
    case ' ':
      viewer.core.is_animating = !viewer.core.is_animating;
      return true;
    default:
      return false;
  }
}

// Usage: ./demo_bin [.off file name] [.dmat file name] [.dmat file name]
// [algorithm name] [iteration number] [rho]
int main(int argc, char *argv[]) {
  if (argc < 6) {
    std::cout << "Not enough input parameters." << std::endl
              << "Usage: demo_bin [.off file name] [.dmat file name] "
                 "[.dmat file name] [algorithm name] [iteration number] "
                 "[rho]" << std::endl;
    return 0;
  }
  // Read V and F from file.
  igl::readOFF(argv[1], V, F);
  // U is initialized with V and gets updated in every frame. V will be the
  // initial state of vertices during the animation.
  U = V;

  // Read S from file. See comments about S above.
  igl::readDMAT(argv[2], S);
  // This works the same as the Matlab code: b = 0 : V.rows() - 1.
  igl::colon<int>(0, V.rows() - 1, b);
  // stable_partition will partition b such that elements with S(i) >= 0 are in
  // front of others. Moreover, if two elements both satisfy S(i) >= 0, their
  // order are kept after the partition.
  // conservativeresize will resize the vector and keep the old values.
  // Example:
  // b = 0 1 2  3 4  5 6
  // S = 0 1 1 -1 1 -1 1
  // After stable_partion:
  // b = 0 1 2 4 6 3 5
  // After conservativeresize:
  // b = 0 1 2 4 6
  b.conservativeResize(std::stable_partition(b.data(), b.data() + b.size(),
    [](int i)->bool { return S(i) >= 0; }) - b.data());
  // Compute the fixed vertices.
  igl::readDMAT(argv[3], bc);

  // Parse the algorithm name.
  std::string algorithm(argv[4]);
  int iter_num = atoi(argv[5]);
  if (algorithm == "arap") {
    std::cout << "Use ArapSolver." << std::endl;
    solver = new arap::demo::ArapSolver(V, F, b, iter_num);
  } else if (algorithm == "admm-fixed") {
    std::cout << "Use AdmmFixedSolver." << std::endl;
    double rho = atof(argv[6]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdmmFixedSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "admm-free") {
    std::cout << "Use AdmmFreeSolver." << std::endl;
    double rho = atof(argv[6]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdmmFreeSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "adapt-admm-fixed") {
    std::cout << "Use AdaptAdmmFixedSolver." << std::endl;
    double rho = atof(argv[6]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdaptAdmmFixedSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "adapt-admm-free") {
    std::cout << "Use AdaptAdmmFreeSolver." << std::endl;
    double rho = atof(argv[6]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdaptAdmmFreeSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "arap-benchmark") {
    std::cout << "Use ArapBenchmarkSolver." << std::endl;
    solver = new arap::demo::ArapBenchmarkSolver(V, F, b, iter_num);
  }

  solver->Precompute();
  // Prepare to solve the problem.
  solver->SolvePreprocess(bc);

  // Set colors for selected vertices.
  C.resize(b.rows(), 3);
  for (int i = 0; i < b.rows(); ++i) {
    C.row(i) = kPurple;
  }
  // Extract selected vertices from U.
  Eigen::MatrixXd U2;
  // U2 = U(b, :).
  igl::slice(U, b, 1, U2);

  igl::Viewer viewer;
  // Draw the mesh.
  viewer.data.set_mesh(U, F);
  // Draw overlays: U2 are selected vertices whose color are set by C.
  viewer.data.set_points(U2, C);
  // pre_draw is called before every frame. Use this function to update
  // vertices.
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core.is_animating = false;
  viewer.core.animation_max_fps = 30.;
  std::cout << "Press [space] to toggle animation" << std::endl;
  viewer.launch();
  return 0;
}
