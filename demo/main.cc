#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/svd3x3/arap.h>
#include <igl/viewer/Viewer.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

typedef
  std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> >
  RotationList;

const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
// Vertex matrix. V is the original vertices from .off file, and U is the
// vertices updated in each frame.
Eigen::MatrixXd V, U;
// Face matrix. F is read from .off file.
Eigen::MatrixXi F;
// S is a column vector representing vertices with predefined coordinates in
// each frame. S is read from .dmat file, whose file format can be found in
// http://igl.ethz.ch/projects/libigl/file-formats/dmat.html. The dimension of
// S is the same as the # of vertices. b contains the indices of those vertices
// whose S is nonnegative.
Eigen::VectorXi S, b;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
igl::ARAPData arap_data;

// Color used to draw precomputed vertices.
static const Eigen::RowVector3d kPurple(80.0 / 255.0,
                                        64.0 / 255.0,
                                        255.0 / 255.0);
static const Eigen::RowVector3d kGold(255.0 / 255.0,
                                      228.0 / 255.0,
                                      58.0 / 255.0);

bool pre_draw(igl::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
    MatrixXd bc(b.size(),V.cols());
    for(int i = 0;i<b.size();i++)
    {
      bc.row(i) = V.row(b(i));
      switch(S(b(i)))
      {
        case 0:
        {
          const double r = mid(0)*0.25;
          bc(i,0) += r*sin(0.5*anim_t*2.*igl::PI);
          bc(i,1) -= r+r*cos(igl::PI+0.5*anim_t*2.*igl::PI);
          break;
        }
        case 1:
        {
          const double r = mid(1)*0.15;
          bc(i,1) += r+r*cos(igl::PI+0.15*anim_t*2.*igl::PI);
          bc(i,2) -= r*sin(0.15*anim_t*2.*igl::PI);
          break;
        }
        case 2:
        {
          const double r = mid(1)*0.15;
          bc(i,2) += r+r*cos(igl::PI+0.35*anim_t*2.*igl::PI);
          bc(i,0) += r*sin(0.35*anim_t*2.*igl::PI);
          break;
        }
        default:
          break;
      }
    }
    igl::arap_solve(bc,arap_data,U);
    viewer.data.set_vertices(U);
    MatrixXd C(b.rows(), 3);
    MatrixXd U2;
    igl::slice(U, b, 1, U2);
    for (int v = 0; v < b.rows(); ++v) {
      C.row(v) = kPurple;
    }
    viewer.data.set_points(U2, C);
    viewer.data.compute_normals();
  if(viewer.core.is_animating)
  {
    anim_t += anim_t_dir;
  }
  return false;
}

bool key_down(igl::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core.is_animating = !viewer.core.is_animating;
      return true;
  }
  return false;
}

int main(int argc, char *argv[]) {
  using namespace Eigen;
  using namespace std;
  // Read V and F from file.
  igl::readOFF("model/decimated-knight.off", V, F);
  // U is initialized with V, so V will be the initial state of vertices during
  // the animation.
  U = V;

  // Read S from file. See comments about S above.
  igl::readDMAT("model/decimated-knight-selection.dmat", S);
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
  b.conservativeResize(stable_partition(b.data(), b.data() + b.size(),
    [](int i)->bool { return S(i) >= 0; }) - b.data());
  // Centroid
  mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
  // Precomputation
  arap_data.max_iter = 100;
  igl::arap_precomputation(V,F,V.cols(),b,arap_data);

  // Set colors for each vertices.
  MatrixXd C(b.rows(), 3);
  MatrixXd U2;
  igl::slice(U, b, 1, U2);
  for (int v = 0; v < b.rows(); ++v) {
    C.row(v) = kPurple;
  }

  igl::Viewer viewer;
  viewer.data.set_mesh(U, F);
  viewer.data.set_points(U2, C);
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core.is_animating = false;
  viewer.core.animation_max_fps = 30.;
  cout<<
    "Press [space] to toggle animation"<<endl;
  viewer.launch();
}
