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

  // string repo_path = "/home/michelle/Documents/LIBIGL/hackathon/libigl-example-project/";
  string repo_path = "/Users/mihcelle/Documents/wavelets/libigl-example-project/";
  const string mesh_path = repo_path + "bobsubdiv1.obj";
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
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);
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
  // igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);

  // Read in visold for knightloop
  // std::vector<int> visoldvec;
  // ifstream visold_in(repo_path+"v_is_old.txt");
  // if (visold_in) {      
  //   double value;
  //   // read the elements in the file into a vector  
  //   while ( visold_in >> value ) {
  //     visoldvec.push_back(value);
  //   }
  // }
  // visold_in.close();

  // HELPER FUNCTIONS FOR FWT
  Eigen::MatrixXi v_is_old; // Map vid to 1 if in Vold, and 0 if Vnew
  Eigen::MatrixXi F_coarse; // Matrix whos e rows are vids in V_in that make up the coarse mesh
  Eigen::MatrixXi fids_covered_by_F_coarse; // #F_coarse x 4 fids in F_in that that tile covers

  // v_is_old.resize(visoldvec.size(),1);
  // int counter = 0;
  // for(int i=0; i<visoldvec.size(); i++)
  // {
  //   v_is_old(i,0) = visoldvec.at(i);
  // }
  // igl::read_triangle_mesh(repo_path+"knightloopcoarse.off",V,F_coarse);
  
  if(is_quadrisection(
    F, 
    V, 
    v_is_old,
    F_coarse,
    fids_covered_by_F_coarse
  )){

    // ofstream visold;
    // ofstream fcoarse;
    // visold.open (repo_path+"v_is_old.txt");
    // fcoarse.open (repo_path+"f_coarse.txt");
    // for(int f=0; f<F_coarse.rows(); f++)
    // {
    //   fcoarse << F_coarse.row(f) << std::endl;
    // }
    // for(int v=0; v<V.rows(); v++)
    // {
    //   visold << v_is_old(v,0) << std::endl;
    // }
    // visold.close();
    // fcoarse.close();
    // assert( visoldvec.size() == v_is_old.rows() );

    // Display the coarser mesh that taubin extracts
    // F = F_coarse;
    // viewer.data().clear();
    // viewer.data().set_mesh(V,F);
    // viewer.data().set_face_based(true);
    // viewer.launch();

    // Make copy of V positions to mutate
    Eigen::MatrixXd V_fine = Eigen::MatrixXd(V);
    Eigen::MatrixXi F_fine = Eigen::MatrixXi(F);

    // Get edge maps
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

    // std::cout << "neighbours_fine: " << std::endl;
    // for(
    //   std::map<int, std::vector<int>>::iterator t = neighbours_fine.begin();
    //   t != neighbours_fine.end();
    //   t++
    // )
    // {
    //   std::cout << "Print fine neighbourhood of vnew vert: " << t->first << std::endl;
    //   for(
    //     std::vector<int>::iterator t2 = t->second.begin();
    //     t2 != t->second.end();
    //     t2++
    //   ){
    //     std::cout << *t2 << std::endl;
    //   }
    //   std::cout << "------" << std::endl;
    // }
    // return 0;


    fwt_lifting1(
      boundary_vold_to_vnew_map,
      V_fine
    );
    // std::cout << V_fine << std::endl;

    fwt_lifting2(
      boundary_vnew_to_vold_map,
	    V_fine 
    );
    // std::cout << V_fine << std::endl;

    // Lifting 3
    fwt_lifting3(
      neighbours_fine,
      v_is_old,
      v_is_boundary,
      V_fine
    );
    // std::cout << V_fine << std::endl;
    // return 0;

    // // Scaling
    fwt_scaling(
      v_is_old,
      v_is_boundary,
      neighbours_coarse,
      V_fine
    );
    // for(int i=0; i<V_fine.rows(); i++)
    // {
    //   if(
    //     std::abs(V_fine(i,0)) > 5 || 
    //     std::abs(V_fine(i,1)) > 5 || 
    //     std::abs(V_fine(i,2)) > 5
    //   ){
    //     std::cout << "vid " << i << " is an exception:" << std::endl; 
    //     if(v_is_old(i,0)==1)
    //     {
    //       std::cout << "it has " << neighbours_coarse[i].size() << " neighbours:" << std::endl; 
    //     }
    //     std::cout << V_fine.row(i) << std::endl; 
    //   }
    // }

    // std::cout << V_fine << std::endl;
    // Lifting 4
    // fwt_lifting4 (
    //   v_is_old,
    //   v_is_boundary,
    //   neighbours_fine,
    //   neighbours_coarse,
    //   V_fine
    // );
    // // std::cout << V_fine << std::endl;
    // // Lifting 5
    // fwt_lifting5 (
    //   boundary_vnew_to_vold_map,
    //   boundary_vold_to_vnew_map,
    //   V_fine
    // );
    // // std::cout << V_fine << std::endl;
    // // Lifting 6 
    // fwt_lifting6 (
    //   v_is_old,
    //   v_is_boundary,
    //   neighbours_fine,
    //   neighbours_coarse,
    //   V_fine
    // );

    // Visualize the output of the lifting schemes
    viewer.data().clear();
    viewer.data().set_mesh(V_fine,F_coarse);
    viewer.data().set_face_based(true);
    viewer.launch();
  }
  // Henrik takes nice photos
}
