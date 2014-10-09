#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char* argv[]) {
  igl::readOFF(argv[1], V, F);

  igl::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.launch();
  return 0;
}
