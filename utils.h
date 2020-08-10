#ifndef UTILS
#define UTILS
#endif

#include <Eigen/Core>
#include <map>


void get_boundary_vert_structures(
	const std::map<std::pair<int,int>, std::vector<int>>& edgemap_fine,
	const Eigen::MatrixXi& v_is_old,
	Eigen::MatrixXi& v_is_boundary,
	std::map<int, std::vector<int>>& boundary_vold_to_vnew_map,
	std::map<int, std::vector<int>>& boundary_vnew_to_vold_map
);


/**
 * Given input mesh, return a map with keys.
 * Each key is an ordered pair of vertex ids (v1 comes before v2 if v1 > v2)
 * and maps to a list of face ids of all the faces 
 * incident to the edge formed by v1 and v2.
*/
void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& edgemap
);


void get_edgemap_and_neighbourhoods(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& edgemap,
	std::map<int, std::vector<int>>& neighbours
);

/**
 * Given map from vid pair to list of incident fids,
 * returns a vector of all the vids which formed a 
 * boundary edge (ie boundary vertices) as well as 
 * a map from the boundary vid to a list of its neighbours.
 * 
 * Inputs:
 * 	edgemap: Map from pair of vids forming edge to list 
 * 									of fids that are incident to it
 * Outputs:
 * 	boundary_vertices: List of vids used in incident faces 
 * 										 that are boundary vertices
 * 	neighbouring_vertices: Map from boundary vid to list of 
 * 												 vids in its neighbourhood
*/
void get_boundary_vertices(
	const std::map<std::pair<int,int>, std::vector<int>>& edgemap,
	std::vector<int>& boundary_vertices,
	std::map<int, std::vector<int>>& neighbouring_vertices
);

/**
 * Given connectivity of a coarse mesh and two vids in the coarse mesh,
 * return the vid from Vnew which would be inserted between them
*/
int find_boundary_vnew(
	const int& vold1,
	const int& vold2,
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
	std::map<std::pair<int,int>, std::vector<int>> edgemap
);

/**
 * Takes the connectivity of a coarse mesh and fine mesh,
 * and vertices of the fine mesh, and returns a map
 * from id of a boundary vert in Vnew to the two neighbouring
 * boundary verts in Vold
*/
void map_bound_vnew_to_bound_vold(
	const Eigen::MatrixXi& F_fine,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
  const std::map<std::pair<int,int>, std::vector<int>>& edgemap_coarse,
  const std::map<int, std::vector<int>>& neighbours_coarse,
  const std::vector<int>& boundary_vids_coarse,
	std::map<int, std::vector<int>>& bound_vnew_to_bound_volds
);

void get_fig216f_map(
	const Eigen::MatrixXi& v_is_old,
	const Eigen::MatrixXi& v_is_boundary,
	const std::map<int, std::vector<int>>& neighbours_fine,
	std::map<int, std::vector<int>>& fig_216f_map
);

void get_fig216f_map_with_a_splash_of_henrik(
	const int v_new_id,
	const std::map<int, std::vector<int>>& neighbours,
  const Eigen::MatrixXi& v_is_old,
	  		int& v_0,
	  		int& v_1,
	  		int& v_2,
	  		int& v_3
);

void get_neighbours(
	const std::map<std::pair<int,int>, std::vector<int>>& edgemap,
	std::map<int, std::vector<int>>& neighbouring_vertices
);

//Sort an array of three entries
void sort3(int arr[]);
