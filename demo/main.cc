#include <igl/readOFF.h>
#include <igl/readDMAT.h>
#include <igl/viewer/Viewer.h>

#include <string>

#include <stdlib.h>

// Color for the final destination of constrained vertices.
static const Eigen::RowVector3d kPurple(80.0 / 255.0,
                                        64.0 / 255.0,
                                        255.0 / 255.0);
std::string off_file_name;
std::string dmat_file_name;

bool option = false;
bool has_dmat = false;

bool PreDraw(igl::Viewer& viewer) {
  // Wait for input:
  if (std::cin.get() == 'q') {
    // Enter q to exit.
    exit(0);
  }
  // Read off file again.
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOFF(off_file_name, V, F);
  viewer.data.set_mesh(V, F);

  // Read dmat file again.
  if (has_dmat) {
    Eigen::MatrixXd bc;
    igl::readDMAT(dmat_file_name, bc);
    // Display the final destination.
    Eigen::MatrixXd C = Eigen::VectorXd::Ones(bc.rows()) * kPurple;
    viewer.data.set_points(bc, C);
  }
}

// Displays vertices, faces and possibly vertex constraints.
// Usage: .OFF file [.DMAT file] [option]
// Right now the only option we support is -i, which is used to interactively
// display the object. So if we run with:
// ./build/demo_bin a.off b.dmat -i
// or
// ./build/demo_bin a.off -i
// The program will keep reading from a.off and b.dmat until you press esc.
// All the commands we need to check:
// ./build/demo_bin a.off
// ./build/demo_bin a.off -i
// ./build/demo_bin a.off b.dmat
// ./build/demo_bin a.off b.dmat -i
int main(int argc, char* argv[]) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // Display the object itself.
  off_file_name = std::string(argv[1]);
  igl::readOFF(off_file_name, V, F);
  igl::Viewer viewer;
  viewer.data.set_mesh(V, F);

  // Check whether the last argument is -i.
  if (strcmp(argv[argc - 1], "-i") == 0) {
    std::cout << "Interactive model." << std::endl;
    option = true;
  }
  // Check whether we have dmat file.
  has_dmat = option ? argc == 4 : argc == 3;
  if (has_dmat) {
    dmat_file_name = std::string(argv[option ? argc - 2 : argc - 1]);
  }
  // Sanity check the arguments.
  std::cout << "off file: " << off_file_name << std::endl
            << "dmat file: " << (has_dmat ? dmat_file_name : "N/A") << std::endl
            << "option: " << (option ? "-i" : "N/A") << std::endl;
  if (has_dmat) {
    Eigen::MatrixXd bc;
    igl::readDMAT(dmat_file_name, bc);
    // Display the final destination.
    Eigen::MatrixXd C = Eigen::VectorXd::Ones(bc.rows()) * kPurple;
    viewer.data.set_points(bc, C);
  }

  // If -i is enabled, add a PreDraw handler.
  if (option) {
    viewer.callback_pre_draw = &PreDraw;
  }
  viewer.launch();
  return 0;
}
