#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <map>
#include "utils.h"
#include "wei_shangyang_244.h"

void fwt_lifting1 (
	const Eigen::MatrixXi& F_in,
	const std::vector<int>& v_old,
	const Eigen::MatrixXi& F_coarse,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
	Eigen::MatrixXd& V
){
  std::cout << "Entered FWT Lifting 1" << std::endl;

  // Generate all the edge and neighbour information in F_coarse.
  // Note that all the vids will be referencing V_copy
  std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  std::map<int, std::vector<int>> neighbouring_vertices;
  edge_incident_faces(
	  F_coarse,
    incident_faces
  );

	std::vector<int> boundary_vertices_vold;
  get_boundary_vertices(
    incident_faces, 
    boundary_vertices_vold,
    neighbouring_vertices
  );

  int v1, v2, v;
  Eigen::Vector3d v_prime;
	for( // Iterate over the vold boundary vertices
    std::vector<int>::iterator it = boundary_vertices_vold.begin();
    it != boundary_vertices_vold.end();
    it++
  ){

    // Get neighbouring vnew boundary vertex Number 1
    assert(neighbouring_vertices[*it].size()==2);

    v = *it;
    v1 = find_boundary_vnew(
      v,
      neighbouring_vertices[*it][0],
      F_in,
      fids_covered_by_F_coarse,
      incident_faces
    );
    v2 = find_boundary_vnew(
      v,
      neighbouring_vertices[*it][1],
      F_in,
      fids_covered_by_F_coarse,
      incident_faces
    );

    std::cout <<  "vnew1: " << v1 << std::endl;
    std::cout <<  "vnew2: " << v2 << std::endl;

    WT_Lifting_1(
      Eigen::Vector3d(V.row(v)),
      Eigen::Vector3d(V.row(v1)),
      Eigen::Vector3d(V.row(v2)),
      v_prime
    );

    V.row(v) = v_prime;

  }

};

void fwt_lifting2 (
	const Eigen::MatrixXi& F_in,
	const std::vector<int>& v_old,
	const Eigen::MatrixXi& F_coarse,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
	Eigen::MatrixXd& V 
){
  
};
