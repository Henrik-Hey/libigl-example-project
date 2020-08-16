#ifndef UTILS
#define UTILS
#endif

#include <Eigen/Core>
#include <map>
#include <vector>


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

void get_fig216f_map_with_a_splash_of_henrik(
	const int v_new_id,
	const std::map<int, std::vector<int>>& neighbours_fine,
	const std::map<int, std::vector<int>>& neighbours_coarse,
  const Eigen::MatrixXi& v_is_old,
	  		int& v_0,
	  		int& v_1,
	  		int& v_2,
	  		int& v_3
);

//Sort an array of three entries
void sort3(int arr[]);
