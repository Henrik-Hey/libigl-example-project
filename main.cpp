#include <igl/read_triangle_mesh.h>
#include <igl/loop.h>
#include <igl/upsample.h>
#include <igl/writeOFF.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include "taubin.h"
#include "lifted_loop_wt.h"

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi F, OF;
  Eigen::MatrixXd V, OV;

  // Inline mesh of a tetrahedron
  // OV = (Eigen::MatrixXd(4,3)<<
  //   0.0,0.0,0.0,
  //   1.0,0.0,0.0,
  //   0.66,0.33,1.0,
  //   1.0,1.0,0.0).finished();
  // OF = (Eigen::MatrixXi(4,3)<<
  //   1,2,3,
  //   1,3,4,
  //   1,4,2,
  //   2,4,3).finished().array()-1; // Test

  // Inline mesh of a regular hexagon mesh
  double x = std::sqrt(3)/2.;
  double y = 0.5;
  double l = 1.;
  OV = (Eigen::MatrixXd(7,3)<<
    y,-x,0,
    1,0,0,
    y,x,0,
    -y,x,0,
    -1,0,0,
    -y,-x,0,
    0,0,0
    ).finished();
  OF = (Eigen::MatrixXi(6,3)<<
    2,7,1,
    3,7,2,
    3,4,7,
    4,5,7,
    7,5,6,
    7,6,1).finished().array()-1; // Test

  // string repo_path = "/home/michelle/Documents/LIBIGL/hackathon/libigl-example-project/";
  // const string mesh_path = repo_path + "knightloop.off";
  // igl::read_triangle_mesh(mesh_path,OV,OF);

  V = OV;
  F = OF;

  // cout<<R"(Usage:
  //   1  Restore Original mesh
  //   2  Apply In-plane upsampled mesh
  //   3  Apply Loop subdivided mesh
  //   4  Apply False barycentric subdivision
  // )";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case '1':
        {
          V = OV;
          F = OF;
          break;
        }
        case '2':
        {
          igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
          break;
        }
        case '3':
        {
          igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
          break;
        }
        case '4':
        {
          igl::false_barycentric_subdivision(
            Eigen::MatrixXd(V),Eigen::MatrixXi(F),V,F);
          break;
        }
      }
      viewer.data().clear();
      viewer.data().set_mesh(V,F);
      viewer.data().set_face_based(true);
      return true;
    };

  igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);

  std::vector<int> v_old; // vids from V_in that make up the coarse mesh
  Eigen::MatrixXi F_coarse; // matrix whos e rows are vids in V_in that make up the coarse mesh
  Eigen::MatrixXi fids_covered_by_F_coarse; // #F_coarse x 4 fids in F_in that that tile covers
  if(is_quadrisection(
    F, 
    V, 
    v_old,
    F_coarse,
    fids_covered_by_F_coarse
  )){

    // Display the coarser mesh that taubin extracts
    // F = F_coarse;
    // viewer.data().clear();
    // viewer.data().set_mesh(V,F);
    // viewer.data().set_face_based(true);
    // viewer.launch();

    Eigen::MatrixXd V_copy = V;

    // Call Lifting 1
    fwt_lifting1(
      F, // Connectivity information of the input mesh, referencing V_copy
      v_old, // vids in V_copy that make up F_coarse 
      V_copy, // to get manipulated via the lifting scheme
      F_coarse, // connectivity information between the vids in v_old
      fids_covered_by_F_coarse // F_coarse.rows() by 4 matrix with each row pointing to the 4 fids it covers
    );
  }
  // Henrik takes nice photos

}
