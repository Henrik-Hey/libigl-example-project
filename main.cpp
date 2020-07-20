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

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi F, OF;
  Eigen::MatrixXd V, OV;

  // inline mesh of a tetrahedron
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

  string repo_path = "/home/henrikhey/Documents/GitHub/libigl-example-project/";
  const string mesh_path = repo_path + "knightloop.off";
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
          cout<<"Number of vertices is: " << V.rows() << endl;
          cout<<"-------"<<endl;
          // igl::writeOFF("/home/michelle/Documents/LIBIGL/hackathon/libigl-example-project/knightloop.off", V, F);
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

  // igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
  // igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);

	Eigen::MatrixXi F_old;
	Eigen::MatrixXd V_old;
	Eigen::MatrixXi F_new;
	Eigen::MatrixXd V_new;
  if(is_quadrisection(F, V, F_old, V_old, F_new, V_new))
  {
    V = V_old;
    F = F_old;
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_face_based(true);
    viewer.launch();
  }

  /**
   * Lifting I:
   * For each boundary vertex in Vold
   *    Grab the neighbouring boundary vertices in Vnew
   *    Update the boundary vertex in Vold with the given equation
  */

  /**
	 * Lifting II
	 * For each boundary vertex in Vnew
	 * 		Get the neighbouring boundary vertices in Vold
	 * 		Update the position of the Vnew boundary vertex by given equation
	*/

  /**                 
	 * Lifting III
	 * Note: we make the assumption that we are updating
	 * both boundary and interior vertices in Vold by 
	 * vertices in Vnew, not just the interior vertices.
	 * Ie., one of the neightburs in Vnew can also be boundary vertices
	 * 
	 * For each vertex in Vold
	 * 		Grab the n neighbours from Vnew 
	 * 		Update the vertex in Vold with the given equation
	*/

  /**                 
	 * Scaling
	 * For each interior vertex in Vold
	 * 		Get the valence of the current interior vertex from Vold
	 * 		Update the position of the vertex by scaling it by the given scalar
	 * 			which is a function of the scalar
	*/

  /**                 
	 * Lifting IV
	 * For each interior vertex in Vnew
	 * 		Grab the 4 vertices in Vold which enclose the interior vertex of Vnew
	 * 		Update the position of the interior vertex in Vnew by the given equation
	*/

  /**                 
	 * Lifting V
	 * For each boundary vertex in Vnew
	 * 		Grab the 4 vertices in Vold specified by Fig2.16e
	 * 		Update the values of those 4 vertices by the given equation
	*/

  /**                 
	 * Lifting VI
	 * For each interior vertex in Vnew
	 * 		Grab the 4 Vold enclosing vertices
	 * 		For each of the enclosing vertices (i=1,2,3,4)
	 * 				Compute the valence
	 * 				Use valence to compute alpha, beta, gamma, delta, and omega
	 * 				Update the position of the enclosing vertex given the equation
	*/

  // Henrik takes nice photos

}
