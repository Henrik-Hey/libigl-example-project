#ifndef UTILS
#define UTILS
#endif

#include <Eigen/Core>
#include <map>

/**
 * Given input mesh, return a map with keys.
 * Each key is an ordered pair of vertex ids (v1 comes before v2 if |v1| > |v2|)
 * and maps to a list of face ids of all the faces 
 * incident to the edge formed by v1 and v2.
*/
void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
);

/**
 * Given map from vid pair to list of incident fids,
 * returns a vector of all the vids which formed a 
 * boundary edge (ie boundary vertices) as well as 
 * a map from the boundary vid to a list of its neighbours.
*/
void get_boundary_vertices(
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces, // map from pair of vids forming edge to list of fids that are incident to it
	std::vector<int>& boundary_vertices, // List of vids used in incident faces that are boundary vertices
	std::map<int, std::vector<int>>& neighbouring_vertices // Map from boundary vid to list of vids in its neighbourhood
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
	std::map<std::pair<int,int>, std::vector<int>> incident_faces
);

//Sort an array of three entries
void sort3(int arr[]);
