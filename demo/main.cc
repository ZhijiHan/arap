#include <igl/readOFF.h>
#include <igl/readDMAT.h>
#include <igl/viewer/Viewer.h>

static const Eigen::RowVector3d kPurple(80.0 / 255.0,
                                        64.0 / 255.0,
                                        255.0 / 255.0);

int main(int argc, char* argv[]) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd bc;

  igl::readOFF(argv[1], V, F);
  igl::readDMAT(argv[2], bc);

  igl::Viewer viewer;
  viewer.data.set_mesh(V, F);

  // Set colors for selected vertices.
  Eigen::MatrixXd C(bc.rows(), 3);
  for (int i = 0; i < bc.rows(); ++i) {
    C.row(i) = kPurple;
  }
  viewer.data.set_points(bc, C);
  viewer.launch();
  return 0;
}
