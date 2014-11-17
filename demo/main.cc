#include <igl/readOFF.h>
#include <igl/readDMAT.h>
#include <igl/viewer/Viewer.h>

// Color for the final destination of constrained vertices.
static const Eigen::RowVector3d kPurple(80.0 / 255.0,
                                        64.0 / 255.0,
                                        255.0 / 255.0);

// Displays vertices, faces and possibly vertex constraints.
// Usage: .OFF file [.DMAT file]
int main(int argc, char* argv[]) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // Display the object itself.
  igl::readOFF(argv[1], V, F);
  igl::Viewer viewer;
  viewer.data.set_mesh(V, F);

  if (argc > 2) {
    Eigen::MatrixXd bc;
    igl::readDMAT(argv[2], bc);
    // Display the final destination.
    Eigen::MatrixXd C = Eigen::VectorXd::Ones(bc.rows()) * kPurple;
    viewer.data.set_points(bc, C);
  }
  viewer.launch();
  return 0;
}
