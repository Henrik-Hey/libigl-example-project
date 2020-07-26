#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <map>
#include "utils.h"
#include "wei_shangyang_244.h"

void fwt_lifting1 (
	const Eigen::MatrixXi& F_fine,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
  const std::map<std::pair<int,int>, std::vector<int>>& edgemap_coarse,
  const std::map<int, std::vector<int>>& neighbours_coarse,
  const std::vector<int>& boundary_vids_coarse,
	Eigen::MatrixXd& V
){
  std::cout << "Begin FWT Lifting 1" << std::endl;

  // Iterate over the vold boundary vertices and 
  // to find its neighbouring Vnew bound verts and update
  // the coords of the vold boundary verts with the neighbours found
  int v1new, v2new, vold;
  Eigen::Vector3d v_prime;
	for(
    std::vector<int>::const_iterator it = boundary_vids_coarse.begin();
    it != boundary_vids_coarse.end();
    it++
  ){
    // Get neighbouring vnew boundary vertex Number 1
    assert(neighbours_coarse.at(*it).size()==2);

    vold = *it;
    v1new = find_boundary_vnew(
      vold,
      neighbours_coarse.at(*it)[0],
      F_fine,
      fids_covered_by_F_coarse,
      edgemap_coarse
    );
    v2new = find_boundary_vnew(
      vold,
      neighbours_coarse.at(*it)[1],
      F_fine,
      fids_covered_by_F_coarse,
      edgemap_coarse
    );

    WT_Lifting_1(
      Eigen::Vector3d(V.row(vold)),
      Eigen::Vector3d(V.row(v1new)),
      Eigen::Vector3d(V.row(v2new)),
      v_prime
    );
    V.row(vold) = v_prime;
  }

  std::cout << "Completed FWT Lifting 1" << std::endl;
};

void fwt_lifting2 (
  const std::map<int, std::vector<int>>& bound_vnew_to_bound_volds,
	Eigen::MatrixXd& V 
){
  std::cout << "Begin FWT Lifting 2" << std::endl;

  int vold1, vold2, vnew;
  Eigen::Vector3d vnew_prime;
  for(
    std::map<int, std::vector<int>>::const_iterator it = bound_vnew_to_bound_volds.begin();
    it != bound_vnew_to_bound_volds.end();
    it++
  ){
    assert(it->second.size()==2);
    
    vnew = it->first;
    vold1 = it->second[0];
    vold2 = it->second[1];
    
    std::cout << "Vnew id: " << it->first << " cleared" << std::endl;
    std::cout << "vold1: " << vold1 << " vold2: " << vold2 << std::endl;

    
    WT_Lifting_2(
      Eigen::Vector3d(V.row(vnew)),
      Eigen::Vector3d(V.row(vold1)),
      Eigen::Vector3d(V.row(vold2)),
      vnew_prime
    );
    V.row(vnew) = vnew_prime;
  }

  std::cout << "Completed FWT Lifting 2" << std::endl;
};

void fwt_lifting3 (
  const std::map<int, std::vector<int>>& neighbours,
  const Eigen::MatrixXi& v_is_old,
  const Eigen::MatrixXi& v_is_boundary,
        Eigen::MatrixXd& vertices
){
  int vold_int, vnew_int;
  Eigen::Vector3d vold;
  Eigen::Vector3d vold_prime;
    std::vector<int> vold_neighbours;
    std::vector<Eigen::Vector3d> vold_int_neighbours;

  for(
    std::map<int, std::vector<int>>::const_iterator it = neighbours.begin();
    it != neighbours.end();
    it++
  ){
    vold_int = it->first;
    if(v_is_old(vold_int, 0) && !v_is_boundary(vold_int, 0))
    {
      vold_neighbours = it->second;
      vold_int_neighbours.clear();
      
      for(
        std::vector<int>::const_iterator n_it = vold_neighbours.begin();
        n_it != vold_neighbours.end();
        n_it++
      ){
        if(!v_is_old(*n_it, 0) && !v_is_boundary(*n_it, 0))
        {
          vold_int_neighbours.push_back(
            Eigen::Vector3d(vertices.row(*n_it))
          );
        }
      }

      if(vold_int_neighbours.size() > 0)
      {
        vold = vertices.row(vold_int);
        WT_Lifting_3(
          vold_int_neighbours,
          vold,
          vold_prime
        );

        vertices.row(vold_int) = vold_prime;
      }
    }
  }
}

void fwt_scaling (
	const Eigen::MatrixXi& v_is_old,
	const Eigen::MatrixXi& v_is_boundary,
  const std::map<int, std::vector<int>>& neighbours_fine,
	Eigen::MatrixXd& V
){
  std::cout << "Begin FWT Scaling" << std::endl;
  Eigen::Vector3d v_prime;
  for(int v=0; v<V.rows(); v++)
  {
    if(v_is_old(v,0)==1&&v_is_boundary(v,0)==0)
    {
      WT_Scaling(
        Eigen::Vector3d(V.row(v)),
        neighbours_fine.at(v).size(),
        v_prime 
      );
      V.row(v) = v_prime;
    }
  }
  std::cout << "Begin FWT Scaling" << std::endl;
};

void fwt_lifting4 (
	const Eigen::MatrixXi& v_is_old,
	const Eigen::MatrixXi& v_is_boundary,
  const std::map<int, std::vector<int>>& neighbours_fine,
	Eigen::MatrixXd& V
){
  std::cout << "Begin FWT Lifting 4" << std::endl;
  Eigen::Vector3d v_prime;
  int v_0, v_1, v_2, v_3;
  for(int v=0; v<V.rows(); v++)
  {
    if(v_is_old(v,0)==0&&v_is_boundary(v,0)==0)
    {
      get_fig216f_map_with_a_splash_of_henrik(
        v,
        neighbours_fine,
        v_is_old,
        v_0,
        v_1,
        v_2,
        v_3
      );

      WT_Lifting_4(
        Eigen::Vector3d(V.row(v)),
        Eigen::Vector3d(V.row(v_0)),
        Eigen::Vector3d(V.row(v_1)),
        Eigen::Vector3d(V.row(v_2)),
        Eigen::Vector3d(V.row(v_3)),
        v_prime
      );
      V.row(v) = v_prime;
    }
  }
  std::cout << "Completed FWT Lifting 4" << std::endl;
};

void fwt_lifting5 (
	const std::map<int, std::vector<int>>& bound_vnew_to_bound_volds,
	const std::map<int, std::vector<int>>& neighbours_coarse,
	Eigen::MatrixXd& V
){
  std::cout << "Begin FWT Lifting 5" << std::endl;

  int vold1, vold2, vold3, vold4, vnew;
  Eigen::Vector3d vold1_prime, vold2_prime, vold3_prime, vold4_prime;
  for(
    std::map<int, std::vector<int>>::const_iterator it = bound_vnew_to_bound_volds.begin();
    it != bound_vnew_to_bound_volds.end();
    it++
  ){
    assert(it->second.size()==2);
    
    vnew = it->first;
    // Following diagram in Fig2.16e
    vold2 = it->second[0];
    vold3 = it->second[1];
    // Find coarse neighbours of vold 2 and vold 3
    assert(neighbours_coarse.at(vold2).size()==2);
    vold1 = neighbours_coarse.at(vold2)[0]==vold3 ? 
            neighbours_coarse.at(vold2)[1] : neighbours_coarse.at(vold2)[0];
    assert(neighbours_coarse.at(vold3).size()==2);
    vold4 = neighbours_coarse.at(vold3)[0]==vold2 ? 
            neighbours_coarse.at(vold3)[1] : neighbours_coarse.at(vold3)[0];

    WT_Lifting_5(
      Eigen::Vector3d(V.row(vnew)),
      Eigen::Vector3d(V.row(vold1)),
      Eigen::Vector3d(V.row(vold2)),
      Eigen::Vector3d(V.row(vold3)),
      Eigen::Vector3d(V.row(vold4)),
      vold1_prime,
      vold2_prime,
      vold3_prime,
      vold4_prime
    );

    V.row(vold1) = vold1_prime;
    V.row(vold2) = vold2_prime;
    V.row(vold3) = vold3_prime;
    V.row(vold4) = vold4_prime;
  }

  std::cout << "Completed FWT Lifting 5" << std::endl;
};

void fwt_lifting6 (
	const Eigen::MatrixXi& v_is_old,
	const Eigen::MatrixXi& v_is_boundary,
  const std::map<int, std::vector<int>>& neighbours_fine,
	Eigen::MatrixXd& V
){
  std::cout << "Begin FWT Lifting 6" << std::endl;

  Eigen::Vector3d v_i_prime;
  int v_0, v_1, v_2, v_3;
  for(int v=0; v<V.rows(); v++)
  {
    if(v_is_old(v,0)==0&&v_is_boundary(v,0)==0)
    {
      get_fig216f_map_with_a_splash_of_henrik(
        v,
        neighbours_fine,
        v_is_old,
        v_0,
        v_1,
        v_2,
        v_3
      );

      std::cout <<  v_0 << std::endl;
      std::cout <<  v_1 << std::endl;
      std::cout <<  v_2 << std::endl;
      std::cout <<  v_3 << std::endl;

      Eigen::Vector4d W;
      WT_Solve_Weights(
        neighbours_fine.at(v_0).size(),
        neighbours_fine.at(v_1).size(),
        neighbours_fine.at(v_2).size(),
        neighbours_fine.at(v_3).size(),
        W
      );

      WT_Lifting_6(
        Eigen::Vector3d(V.row(v_0)), // vold
        W(0),
        Eigen::Vector3d (V.row(v)), // vnew
        v_i_prime
      );
      V.row(v_0) = v_i_prime;

      WT_Lifting_6(
        Eigen::Vector3d(V.row(v_1)), // vold
        W(1),
        Eigen::Vector3d (V.row(v)), // vnew
        v_i_prime
      );
      V.row(v_1) = v_i_prime;

      WT_Lifting_6(
        Eigen::Vector3d(V.row(v_2)), // vold
        W(2),
        Eigen::Vector3d (V.row(v)), // vnew
        v_i_prime
      );
      V.row(v_2) = v_i_prime;

      WT_Lifting_6(
        Eigen::Vector3d(V.row(v_3)), // vold
        W(3),
        Eigen::Vector3d (V.row(v)), // vnew
        v_i_prime
      );
      V.row(v_3) = v_i_prime;
    }
  }

  std::cout << "Completed FWT Lifting 6" << std::endl;
};
