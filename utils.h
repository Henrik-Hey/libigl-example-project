#ifndef UTILS
#define UTILS
#endif

#include <Eigen/Core>
#include <map>

/**
 * Given input mesh, return a map with keys.
 * Each key is an ordered pair of vertex ids (v1 comes before v2 if |v1| > |v2|)
 * and maps to a list of face ids of all the faces incident to the edge formed by v1 and v2.
*/
void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
);

void get_edges_and_neighbours(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces,
	std::map<std::pair<int,int>, std::vector<int>>& neighbouring_vertices
);

//Sort an array of three entries
void sort3(int arr[]);

void get_boundary_vertices(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
);