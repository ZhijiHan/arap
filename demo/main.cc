#include "admmsolver.h"
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

#define USE_ARAP_SOLVER 0
#define USE_ADMM_SOLVER 1

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
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
// Our own implementation of ARAP.
arap::demo::Solver* solver;

// Color used to draw precomputed vertices.
static const Eigen::RowVector3d kPurple(80.0 / 255.0,
                                        64.0 / 255.0,
                                        255.0 / 255.0);
static const Eigen::RowVector3d kGold(255.0 / 255.0,
                                      228.0 / 255.0,
                                      58.0 / 255.0);
// Error threshold.
static const double kErrorPerVertexCoord = 1e-8;

// Given |frame_number|, compute the fixed vertices and return.
Eigen::MatrixXd ComputeFixedVertices(int frame_number) {
  anim_t = anim_t_dir * frame_number;
  Eigen::MatrixXd bc(b.size(), V.cols());
  for(int i = 0; i < b.size(); ++i) {
    bc.row(i) = V.row(b(i));
    switch (S(b(i))) {
      case 0: {
        const double r = mid(0) * 0.25;
        bc(i, 0) += r * sin(0.5 * anim_t * 2. * igl::PI);
        bc(i, 1) -= r + r * cos(igl::PI + 0.5 * anim_t * 2. * igl::PI);
        break;
      }
      case 1: {
        const double r = mid(1) * 0.15;
        bc(i, 1) += r + r * cos(igl::PI + 0.15 * anim_t * 2. * igl::PI);
        bc(i, 2) -= r * sin(0.15 * anim_t * 2. * igl::PI);
        break;
      }
      case 2: {
        const double r = mid(1) * 0.15;
        bc(i, 2) += r + r * cos(igl::PI + 0.35 * anim_t * 2. * igl::PI);
        bc(i, 0) += r * sin(0.35 * anim_t * 2. * igl::PI);
        break;
      }
      default:
        break;
    }
  }
  return bc;
}

bool pre_draw(igl::Viewer& viewer) {
  static int iteration = 0;
  if (!viewer.core.is_animating
    || iteration >= solver->GetMaxIteration())
    return false;
  solver->SolveOneIteration();
  double energy = solver->ComputeEnergy();
  std::cout << "Iteration: " << iteration << " Energy: " << energy << std::endl;
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

// Usage: ./demo_bin [.off file name] [.dmat file name] [frame number]
int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout << "Not enough input parameters." << std::endl
              << "Usage: demo_bin [.off file name] [.dmat file name] "
                 "[frame number]" << std::endl;
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
  // Centroid of the demo, used to compute the positions of selected vertices
  // during the animation.
  mid = 0.5 * (V.colwise().maxCoeff() + V.colwise().minCoeff());
  // Compute the fixed vertices from the frame number.
  bc = ComputeFixedVertices(atoi(argv[3]));

  // Add our own pre computation implementation here.
#if USE_ARAP_SOLVER
  std::cout << "Use ArapSolver." << std::endl;
  solver = new arap::demo::ArapSolver(V, F, b, 500);
#elif USE_ADMM_SOLVER
  std::cout << "Use AdmmSolver." << std::endl;
  solver = new arap::demo::AdmmSolver(V, F, b, 500, 0.1);
#endif
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
