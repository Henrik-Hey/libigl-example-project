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

  // HELPER FUNCTIONS FOR FWT

  Eigen::MatrixXi v_is_old; // Map vid to 1 if in Vnew, and 0 if Vold
  Eigen::MatrixXi F_coarse; // Matrix whos e rows are vids in V_in that make up the coarse mesh
  Eigen::MatrixXi fids_covered_by_F_coarse; // #F_coarse x 4 fids in F_in that that tile covers
  if(is_quadrisection(
    F, 
    V, 
    v_is_old,
    F_coarse,
    fids_covered_by_F_coarse
  )){

    // Display the coarser mesh that taubin extracts
    // F = F_coarse;
    // viewer.data().clear();
    // viewer.data().set_mesh(V,F);
    // viewer.data().set_face_based(true);
    // viewer.launch();

    // Lifting 1
    // Generate all the edge and neighbour information in F_coarse.
    // Note that all the vids will be referencing V_copy
    std::map<std::pair<int,int>, std::vector<int>> edgemap_coarse;
    std::map<int, std::vector<int>> neighbours_coarse;
    edge_incident_faces(
      F_coarse,
      edgemap_coarse
    );
    // Figure out which Vold vids are boundary verts
    std::vector<int> boundary_vids_coarse;
    get_boundary_vertices(
      edgemap_coarse, 
      boundary_vids_coarse,
      neighbours_coarse
    );

    // Lifting 2
    // Generate a map from Vnew boundary verts
    // to its neighbouring boundary verts in Vold
    std::map<int, std::vector<int>> bound_vnew_to_bound_volds;
    map_bound_vnew_to_bound_vold(
      F,
      fids_covered_by_F_coarse,
      edgemap_coarse,
      neighbours_coarse,
      boundary_vids_coarse,
      bound_vnew_to_bound_volds
    );

    // Lifting 5
    std::map<std::pair<int,int>, std::vector<int>> edgemap_fine;
    std::map<int, std::vector<int>> neighbours_fine;
    std::vector<int> boundary_vids_fine;
    edge_incident_faces(
      F,
      edgemap_fine
    );
    get_boundary_vertices(
      edgemap_fine, 
      boundary_vids_fine, 
      neighbours_fine
    );
    // Assert we have a consistent number of boundary vert ids in vnew
    assert(boundary_vids_fine.size() - boundary_vids_coarse.size() 
            == bound_vnew_to_bound_volds.size());

    // Lifting 6
    std::map<int, std::vector<int>> fig_216f_map;
    get_fig216f_map(
      v_is_old,
      edgemap_fine,
      neighbours_fine,
      fig_216f_map
    );

    // BEGIN FWT

    // Make copy of V positions to mutate
    Eigen::MatrixXd V_copy = Eigen::MatrixXd(V);

    // Call Lifting 1
    fwt_lifting1(
      F,
      fids_covered_by_F_coarse,
      edgemap_coarse,
      neighbours_coarse,
      boundary_vids_coarse,
      V_copy
    );

    fwt_lifting2(
      bound_vnew_to_bound_volds,
	    V_copy 
    );
  
    fwt_lifting5 (
      bound_vnew_to_bound_volds,
      neighbours_coarse,
      V_copy
    );

    fwt_lifting6 (
      fig_216f_map,
      neighbours_fine,
      V_copy
    );

    // END FWT

    // Visualize the output of the lifting schemes
    viewer.data().clear();
    viewer.data().set_mesh(V_copy,F);
    viewer.data().set_face_based(true);
    viewer.launch();

  }
  // Henrik takes nice photos

}
