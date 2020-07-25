#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <map>
#include "utils.h"
#include "wei_shangyang_244.h"

void fwt_lifting1 (
	const Eigen::MatrixXi& F_fine,
	const std::vector<int>& v_old,
	const Eigen::MatrixXi& F_coarse,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
	Eigen::MatrixXd& V
){
  std::cout << "Begin FWT Lifting 1" << std::endl;

  // Generate all the edge and neighbour information in F_coarse.
  // Note that all the vids will be referencing V_copy
  std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  std::map<int, std::vector<int>> neighbouring_vertices;
  edge_incident_faces(
	  F_coarse,
    incident_faces
  );

  // Figure out which Vold vids are boundary verts
	std::vector<int> boundary_vertices_vold;
  get_boundary_vertices(
    incident_faces, 
    boundary_vertices_vold,
    neighbouring_vertices
  );

  // Iterate over the vold boundary vertices and 
  // to find its neighbouring Vnew bound verts and update
  // the coords of the vold boundary verts with the neighbours found
  int v1new, v2new, vold;
  Eigen::Vector3d v_prime;
	for(
    std::vector<int>::iterator it = boundary_vertices_vold.begin();
    it != boundary_vertices_vold.end();
    it++
  ){
    // Get neighbouring vnew boundary vertex Number 1
    assert(neighbouring_vertices[*it].size()==2);

    vold = *it;
    v1new = find_boundary_vnew(
      vold,
      neighbouring_vertices[*it][0],
      F_fine,
      fids_covered_by_F_coarse,
      incident_faces
    );
    v2new = find_boundary_vnew(
      vold,
      neighbouring_vertices[*it][1],
      F_fine,
      fids_covered_by_F_coarse,
      incident_faces
    );

    // std::cout <<  "vnew1: " << v1 << std::endl;
    // std::cout <<  "vnew2: " << v2 << std::endl;
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
	const Eigen::MatrixXi& F_fine,
	const std::vector<int>& v_old,
	const Eigen::MatrixXi& F_coarse,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
	Eigen::MatrixXd& V 
){
  std::cout << "Begin FWT Lifting 2" << std::endl;

  // Generate a map from Vnew boundary verts
  // to its neighbouring boundary verts in Vold
  std::map<int, std::vector<int>> bound_vnew_to_bound_volds;
  map_bound_vnew_to_bound_vold(
    F_fine,
    F_coarse,
    V,
    v_old,
    fids_covered_by_F_coarse,
    bound_vnew_to_bound_volds
  );

  int vold1, vold2, vnew;
  Eigen::Vector3d vnew_prime;
  for(
    std::map<int, std::vector<int>>::iterator it = bound_vnew_to_bound_volds.begin();
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
