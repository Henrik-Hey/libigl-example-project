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

void get_edges_and_neighbours(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
);

void get_boundary_vertices(
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces,
	std::vector<int>& boundary_vertices,
	std::map<int, std::vector<int>>& neighbouring_vertices
);

int find_boundary_vnew(
	const int& vold1,
	const int& vold2,
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
	std::map<std::pair<int,int>, std::vector<int>> incident_faces
);

//Sort an array of three entries
void sort3(int arr[]);

void get_boundary_vertices(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
);