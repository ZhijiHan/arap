#include "adaptadmmfixedsolver.h"
#include "adaptadmmfreesolver.h"
#include "admmfixedsolver.h"
#include "admmfreesolver.h"
#include "arapbenchmarksolver.h"
#include "arapsolver.h"

// C++ standard library
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
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
// Collects data from the solver.
std::ofstream output_file;

// Color used to draw precomputed vertices.
static const Eigen::RowVector3d kPurple(80.0 / 255.0,
                                        64.0 / 255.0,
                                        255.0 / 255.0);
static const Eigen::RowVector3d kGold(255.0 / 255.0,
                                      228.0 / 255.0,
                                      58.0 / 255.0);

static const std::string kDataFolder = "/home/taodu/research/arap/data/";

bool pre_draw(igl::Viewer& viewer) {
  static int iteration = 0;
  static Eigen::MatrixXd last_solution = V;
  if (!viewer.core.is_animating
    || iteration >= solver->GetMaxIteration()) {
    return false;
  }
  // The first time when iteration == 0, we output the energy based on the
  // initialized vertices_updated_ and rotations_. These values should be
  // computed in SolvePreprocess, and should be the same for all the
  // algorithms.
  arap::demo::Energy energy;
  if (iteration == 0) {
    // Get energy types.
    output_file << "iteration\trho\txnorm\txdiffnorm";
    energy = solver->ComputeEnergy();
    std::vector<std::string> types = energy.GetEnergyTypes();
    for (auto it = types.begin(); it != types.end(); ++it) {
      output_file << '\t' << *it;
    }
    output_file << '\n';
  } else {
    solver->SolveOneIteration();
    energy = solver->ComputeEnergy();
  }
  double rho = solver->GetRho();
  const Eigen::MatrixXd solution = solver->GetVertexSolution();
  double solution_diff_norm = (last_solution - solution).norm();
  last_solution = solution;
  std::cout << "After " << iteration << " Iteration:"
            << " rho: " << rho
            << " x norm: " << solution.norm()
            << " x diff norm: " << solution_diff_norm
            << " Energy: " << energy;
  // Write energy back into output data file.
  output_file << iteration << '\t' << rho << '\t'
              << solution.norm() << '\t'
              << solution_diff_norm;
  std::vector<std::string> types = energy.GetEnergyTypes();
  for (auto it = types.begin(); it != types.end(); ++it) {
    output_file << '\t' << energy.GetEnergyValue(*it);
  }
  output_file << '\n';

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

// Usage: ./demo_bin [model file name] [algorithm name] [iteration number] [rho]
// The program will generate data files under the folder kDataFolder.
int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout << "Not enough input parameters." << std::endl
              << "Usage: demo_bin [model file name] "
                 "[algorithm name] [iteration number] "
                 "[rho]" << std::endl;
    return 0;
  }
  std::string model_path = std::string(argv[1]);
  // Get the model name.
  int found = static_cast<int>(model_path.find_last_of("/"));
  std::string model_name = model_path.substr(found + 1);
  std::cout << "model name: " << model_name << std::endl;
  // Build folders under kDataFolder.
  mkdir((kDataFolder + model_name).c_str(), S_IRWXU | S_IRWXG |
        S_IROTH | S_IXOTH);

  // Read V and F from file.
  igl::readOFF(model_path + ".off", V, F);
  // U is initialized with V and gets updated in every frame. V will be the
  // initial state of vertices during the animation.
  U = V;

  // Read S from file. See comments about S above.
  igl::readDMAT(model_path + "-selection.dmat", S);
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
  igl::readDMAT(model_path + ".dmat", bc);

  // Parse the algorithm name.
  std::string algorithm(argv[2]);
  int iter_num = atoi(argv[3]);
  if (algorithm == "arap") {
    std::cout << "Use ArapSolver." << std::endl;
    solver = new arap::demo::ArapSolver(V, F, b, iter_num);
  } else if (algorithm == "admm-fixed") {
    std::cout << "Use AdmmFixedSolver." << std::endl;
    double rho = atof(argv[4]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdmmFixedSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "admm-free") {
    std::cout << "Use AdmmFreeSolver." << std::endl;
    double rho = atof(argv[4]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdmmFreeSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "adapt-admm-fixed") {
    std::cout << "Use AdaptAdmmFixedSolver." << std::endl;
    double rho = atof(argv[4]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdaptAdmmFixedSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "adapt-admm-free") {
    std::cout << "Use AdaptAdmmFreeSolver." << std::endl;
    double rho = atof(argv[4]);
    std::cout << "rho = " << rho << std::endl;
    solver = new arap::demo::AdaptAdmmFreeSolver(V, F, b, iter_num, rho);
  } else if (algorithm == "arap-benchmark") {
    std::cout << "Use ArapBenchmarkSolver." << std::endl;
    solver = new arap::demo::ArapBenchmarkSolver(V, F, b, iter_num);
  }

  // Build output data file name.
  std::string output_file_path = kDataFolder + model_name + "/"
      + algorithm + "-"  // algorithm name.
      + std::string(argv[3]);  // iteration number.
  // If algorithm is admm, add rho into the file name.
  if (algorithm != "arap") {
    output_file_path += "-" + std::string(argv[4]);  // rho.
  }
  output_file_path += ".txt";
  std::cout << "output file path: " << output_file_path << std::endl;
  output_file.open(output_file_path);

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
