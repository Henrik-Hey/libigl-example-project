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

  // string repo_path = "/Users/mihcelle/Documents/wavelets/libigl-example-project/";
  string repo_path = "/home/michelle/Documents/LIBIGL/hackathon/libigl-example-project/";
  const string mesh_path = repo_path + "spot.obj";
  // const string mesh_path = repo_path + "loop_test_close4_6_18.off";
  // const string mesh_path = repo_path + "torus_9_36.off";
  igl::read_triangle_mesh(mesh_path,OV,OF);

  V = OV;
  F = OF;

  // cout<<R"(Usage:
  //   1  Restore Original mesh
  //   2  Apply In-plane upsampled mesh
  //   3  Apply Loop subdivided mesh
  //   4  Apply False barycentric subdivision
  // )";
  igl::opengl::glfw::Viewer viewer;

  // viewer.callback_key_down =
  //   [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
  //   {
  //     switch(key)
  //     {
  //       default:
  //         return false;
  //       case '1':
  //       {
  //         V = OV;
  //         F = OF;
  //         break;
  //       }
  //       case '2':
  //       {
  //         igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  //         break;
  //       }
  //       case '3':
  //       {
  //         igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  //         break;
  //       }
  //       case '4':
  //       {
  //         igl::false_barycentric_subdivision(
  //           Eigen::MatrixXd(V),Eigen::MatrixXi(F),V,F);
  //         break;
  //       }
  //     }
  //     viewer.data().clear();
  //     viewer.data().set_mesh(V,F);
  //     viewer.data().set_face_based(true);
  //     return true;
  //   };

  std::cout << "Unleash the wavelets!" << std::endl;
  igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);

  Eigen::MatrixXi v_is_old; // Map vid to 1 if in Vold, and 0 if Vnew
  Eigen::MatrixXi F_coarse; // Matrix whos e rows are vids in V_in that make up the coarse mesh
  Eigen::MatrixXi fids_covered_by_F_coarse; // #F_coarse x 4 fids in F_in that that tile covers
  if(is_quadrisection(
    F, 
    V, 
    v_is_old,
    F_coarse,
    fids_covered_by_F_coarse
  )){

    // Make copy of V positions to mutate
    Eigen::MatrixXd V_fine = Eigen::MatrixXd(V);
    Eigen::MatrixXi F_fine = Eigen::MatrixXi(F);

    // Helper methods
    std::map<std::pair<int,int>, std::vector<int>> edgemap_coarse;
    std::map<std::pair<int,int>, std::vector<int>> edgemap_fine;
    std::map<int, std::vector<int>> neighbours_coarse;
    std::map<int, std::vector<int>> neighbours_fine;
    get_edgemap_and_neighbourhoods(
      F_coarse,
      edgemap_coarse,
      neighbours_coarse
    );
    get_edgemap_and_neighbourhoods(
      F,
      edgemap_fine,
      neighbours_fine
    );
    Eigen::MatrixXi v_is_boundary;
    std::map<int, std::vector<int>> boundary_vold_to_vnew_map;
    std::map<int, std::vector<int>> boundary_vnew_to_vold_map;
    get_boundary_vert_structures(
      edgemap_fine, 
      v_is_old,
      v_is_boundary,
      boundary_vold_to_vnew_map,
      boundary_vnew_to_vold_map
    );

    // Sanity all the pre-lifting helper functions
    int countNumNewBoundVerts = 0;
    int countNumOldBoundVerts = 0;
    for(int i=0; i<V_fine.rows(); i++)
    {
      if(v_is_old(i,0)==1 && v_is_boundary(i,0)==1)
      {
        countNumOldBoundVerts++;
      } 
      else if(v_is_old(i,0)==0 && v_is_boundary(i,0)==1)
      {
        countNumNewBoundVerts++;
      }
    }
    std::cout << "number of boundaries in fine: " << countNumNewBoundVerts+countNumOldBoundVerts << std::endl;
    for(
      std::map<std::pair<int,int>, std::vector<int>>::iterator it = edgemap_fine.begin();
      it!=edgemap_fine.end();
      it++
    ){
      if(it->second.size()>2)
      {
        std::cout << "found irregular edge" << std::endl;
      }
    }
    assert(boundary_vold_to_vnew_map.size()==countNumOldBoundVerts);
    assert(boundary_vnew_to_vold_map.size()==countNumNewBoundVerts);
    assert(neighbours_fine.size()==V_fine.rows());

    // Performing FWT
    fwt_lifting1(
      boundary_vold_to_vnew_map,
      V_fine
    );
    fwt_lifting2(
      boundary_vnew_to_vold_map,
	    V_fine 
    );
    Eigen::MatrixXd V_after_l2 = Eigen::MatrixXd(V_fine);
    fwt_lifting3(
      neighbours_fine,
      v_is_old,
      v_is_boundary,
      V_fine
    );
    Eigen::MatrixXd V_after_l3 = Eigen::MatrixXd(V_fine);
    if(V_after_l2==V_after_l3) std::cout << "L3 did nothing." << std::endl;
    fwt_scaling(
      v_is_old,
      v_is_boundary,
      neighbours_coarse,
      V_fine
    );
    Eigen::MatrixXd V_after_scale = Eigen::MatrixXd(V_fine);
    fwt_lifting4 (
      v_is_old,
      v_is_boundary,
      neighbours_fine,
      neighbours_coarse,
      V_fine
    );
    fwt_lifting5 (
      boundary_vnew_to_vold_map,
      boundary_vold_to_vnew_map,
      V_fine
    );
    fwt_lifting6 (
      v_is_old,
      v_is_boundary,
      neighbours_fine,
      neighbours_coarse,
      V_fine
    );
    Eigen::MatrixXd V_after_l6 = Eigen::MatrixXd(V_fine);
    if(V_after_scale==V_after_l6) std::cout << "L4-L6 collectively did nothing.";

    // Visualize the output of the lifting schemes
    viewer.data().clear();
    viewer.data().set_mesh(V_fine,F_coarse);
    // viewer.data().set_mesh(OV,OF);
    viewer.data().set_face_based(true);
    viewer.launch();
  }
  // Henrik takes nice photos
}
