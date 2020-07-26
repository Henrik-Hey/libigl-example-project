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
	const std::map<int, std::vector<int>>& fig_216f_map,
  const std::map<int, std::vector<int>>& neighbours_fine,
	Eigen::MatrixXd& V
){
  std::cout << "Begin FWT Lifting 6" << std::endl;


  std::vector<double> valence;
  for(
    std::map<int, std::vector<int>>::const_iterator it = fig_216f_map.begin();
    it != fig_216f_map.end();
    it++
  ){

    valence.clear();
    Eigen::Vector3d v_i_prime;

    // Iterate over the vold neighbours of current interior vnew
    for(
      std::vector<int>::const_iterator it_n = it->second.begin();
      it_n != it->second.end();
      it_n++
    ){
      valence.emplace_back(neighbours_fine.at(*it_n).size());
    }
    assert(valence.size()==4);

    Eigen::Vector4d W;
    WT_Solve_Weights(
      valence.at(0),
      valence.at(1),
      valence.at(2),
      valence.at(3),
      W
    );

    for(int n=0; n<it->second.size(); n++)
    {
      WT_Lifting_6(
        Eigen::Vector3d(V.row(it->second.at(n))), // vold
        W(n),
        Eigen::Vector3d (V.row(it->first)), // vnew
        v_i_prime
      );

      V.row(it->second.at(n)) = v_i_prime;
    }

  }
  std::cout << "Completed FWT Lifting 6" << std::endl;
};
