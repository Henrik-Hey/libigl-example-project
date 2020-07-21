#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <map>
#include "utils.h"

void fwt_lifting1 (
	const Eigen::MatrixXi& F_in,
	const std::vector<int>& v_old,
	Eigen::MatrixXd& V,
	Eigen::MatrixXi& F_coarse,
	Eigen::MatrixXi& fids_covered_by_F_coarse
){
  std::cout << "Entered FWT Lifting 1" << std::endl;
  // Generate all the edge and neighbour information in F_coarse.
  // Note that all the vids will be referencing V_copy
  std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  std::map<int, std::vector<int>> neighbouring_vertices;
  get_edges_and_neighbours(
	  F_coarse,
    incident_faces
  );

	std::vector<int> boundary_vertices_vold;
  get_boundary_vertices(
    incident_faces, 
    boundary_vertices_vold,
    neighbouring_vertices
  );

	for(
    std::vector<int>::iterator it = boundary_vertices_vold.begin();
    it != boundary_vertices_vold.end();
    it++
  ){

    // Get neighbouring vnew boundary vertex Number 1
    assert(neighbouring_vertices[*it].size()==2);
    std::cout <<  "vnew1: " << find_boundary_vnew(
      *it,
      neighbouring_vertices[*it][0],
      F_in,
      fids_covered_by_F_coarse,
      incident_faces
    ) << std::endl;


    std::cout << "vnew2: " <<  find_boundary_vnew(
      *it,
      neighbouring_vertices[*it][1],
      F_in,
      fids_covered_by_F_coarse,
      incident_faces
    ) << std::endl;

  }

};
