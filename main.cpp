#include <igl/read_triangle_mesh.h>
#include <igl/loop.h>
#include <igl/upsample.h>
#include <igl/writeOFF.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <vector>
#include "taubin.h"
#include "utils.h"
#include "lifted_loop_fwt.h"

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi F, OF;
  Eigen::MatrixXd V, OV;

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

  // string repo_path = "/Users/mihcelle/Documents/wavelets/libigl-example-project/";
  string repo_path = "/home/michelle/Documents/LIBIGL/hackathon/libigl-example-project/";
  const string mesh_path = repo_path + "spot.obj";
  // const string mesh_path = repo_path + "loop_test_close4_6_18.off";
  // const string mesh_path = repo_path + "torus_9_36.off";
  igl::read_triangle_mesh(mesh_path,OV,OF);
  igl::loop( Eigen::MatrixXd(OV), Eigen::MatrixXi(OF), OV,OF);

  V = OV;
  F = OF;

  // Launch viewer
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);

  cout<<R"(Usage:
    1  Restore Original mesh
    2  Apply In-plane upsampled mesh
    3  Apply Loop subdivided mesh
    4  Apply FWT
  )";
  viewer.callback_key_down 
    =[&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
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
        // Make copy of V positions to mutate
        Eigen::MatrixXd V_fine = Eigen::MatrixXd(V);
        Eigen::MatrixXi F_fine = Eigen::MatrixXi(F);
        Eigen::MatrixXi F_coarse;
        if(
          FWT(
            F_fine,
            F_coarse,
            V_fine)
        ) {
          V = V_fine;
          assert(F_coarse.rows() < F.rows());
          F.resize(F_coarse.rows(), F_coarse.cols());
          assert(F.rows() == F_coarse.rows());
          F = F_coarse;
        }
        break;
      }
    }
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_face_based(true);
    return true;
  };

  viewer.launch();

  // Henrik takes nice photos
}
