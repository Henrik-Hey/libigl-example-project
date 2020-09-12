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
  string repo_path = "/home/gmin7/Documents/wavelets/libigl-example-project/";
  const string mesh_path = repo_path + "spot.obj";
  // const string mesh_path = repo_path + "loop_test_close4_6_18.off";
  // const string mesh_path = repo_path + "torus_9_36.off";
  igl::read_triangle_mesh(mesh_path,OV,OF);
  igl::loop( Eigen::MatrixXd(OV), Eigen::MatrixXi(OF), OV,OF);
  igl::loop( Eigen::MatrixXd(OV), Eigen::MatrixXi(OF), OV,OF);

  V = OV;
  F = OF;

	Eigen::MatrixXi v_is_boundary;
	std::map<std::pair<int, int>, std::vector<int>> edgemap_coarse;
	std::map<std::pair<int, int>, std::vector<int>> edgemap_fine;
	std::map<int, std::vector<int>> neighbours_coarse;
	std::map<int, std::vector<int>> neighbours_fine;
	std::map<int, std::vector<int>> boundary_vold_to_vnew_map;
	std::map<int, std::vector<int>> boundary_vnew_to_vold_map;

	Eigen::MatrixXd V_fine = Eigen::MatrixXd(V);
	Eigen::MatrixXi F_fine = Eigen::MatrixXi(F);
	Eigen::MatrixXi F_coarse;
	Eigen::MatrixXi v_is_old;									// Map vid to 1 if in Vold, and 0 if Vnew
	Eigen::MatrixXi fids_covered_by_F_coarse; // #F_coarse x 4 fids in F_in that that tile covers

  // Launch viewer
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);
  
  Eigen::MatrixXd C = Eigen::MatrixXd(V);
  for(int r=0; r<C.rows(); r++)
  {
    C.row(r) = Eigen::Vector3d(255./255., 119./255., 82./255.);
  }

  // Add per-vertex colors
  viewer.data().set_colors(C);

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

				if (is_quadrisection(
					F_fine,
					V_fine,
					v_is_old,
					F_coarse,
					fids_covered_by_F_coarse))
				{

					get_aux_data(
						V_fine,
						F_fine,
						F_coarse,
						v_is_old,
						v_is_boundary,
						edgemap_coarse,
						edgemap_fine,
						neighbours_coarse,
						neighbours_fine,
						boundary_vold_to_vnew_map,
						boundary_vnew_to_vold_map);

					// FWT time!
					fwd_lifting_scheme(
						neighbours_coarse,
						neighbours_fine,
						boundary_vold_to_vnew_map,
						boundary_vnew_to_vold_map,
						v_is_boundary,
						v_is_old,
						V_fine);

          int countnew = 0;
          for(int v=0; v<V_fine.rows(); v++)
          {
            if(v_is_old(v)==1) countnew++; 
          }

          std::cout << countnew << " hey there" << std::endl;

          Eigen::MatrixXd tempVerts;
          tempVerts.resize(countnew, 3);

          std::map<int, int> vertMap;

          int newCount = 0;
          for(int v=0; v<V_fine.rows(); v++)
          {
            if(v_is_old(v)==1) {
              tempVerts.row(newCount) = V_fine.row(v);
              vertMap.insert({v,newCount});
              std::cout << v << " maps to " << vertMap[v] << " " << std::endl;
              newCount++;
            }
          }

          for(int f=0; f<F_coarse.rows(); f++)
          {


            assert (F_coarse(f,0) < tempVerts.rows());
            assert (F_coarse(f,1) < tempVerts.rows());
            assert (F_coarse(f,2) < tempVerts.rows());

            F_coarse(f,0) = vertMap[F_coarse(f,0)];
            F_coarse(f,1) = vertMap[F_coarse(f,1)];
            F_coarse(f,2) = vertMap[F_coarse(f,2)];
          }

          edgemap_coarse.clear();
          edgemap_fine.clear();
          neighbours_coarse.clear();
          neighbours_fine.clear();
          boundary_vold_to_vnew_map.clear();
          boundary_vnew_to_vold_map.clear();

          V_fine = Eigen::MatrixXd(tempVerts);
          V = Eigen::MatrixXd(tempVerts);

          F_fine = Eigen::MatrixXi(F_coarse);
          F = Eigen::MatrixXi(F_coarse);
				}


        break;
      }
      case '5':
			{

				inv_lifting_scheme(
					neighbours_coarse,
					neighbours_fine,
					boundary_vold_to_vnew_map,
					boundary_vnew_to_vold_map,
					v_is_boundary,
					v_is_old,
					V_fine
				);

				// V = V_fine;
				F = F_fine;

				break;
      }
    }
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_face_based(true);
    Eigen::MatrixXd C = Eigen::MatrixXd(V);
    for(int r=0; r<C.rows(); r++)
    {
      C.row(r) = Eigen::Vector3d(255./255., 119./255., 82./255.);
    }
    viewer.data().set_colors(C);

    return true;
  };

  viewer.launch();

  // Henrik takes nice photos
}
