#ifndef TAUBIN
#define TAUBIN
#endif

#include <Eigen/Core>
#include <vector>
#include <map>

// Adapted from Gabriel Taubin's Loop Inverse Subdivsion Algorithm
// Paper: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.21.1125&rep=rep1&type=pdf


/**
 * Detect whether the input mesh has subdivision connectivity.
 * If so, partitions the input mesh into the coarse mesh and 
 * the new vertices after subdividing.
 * Returns true if has subdivision connectivity,false if not.
*/
bool is_quadrisection(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	std::vector<int>& v_old,
	Eigen::MatrixXi& F_coarse,
	Eigen::MatrixXi& fids_covered_by_F_coarse
);

/**
 * Creates tiles from the input mesh, and the tile
 * sets covered by each.
*/
void covering_mesh(
	const Eigen::MatrixXi& F_in,
 	Eigen::MatrixXi& tiles,
 	Eigen::MatrixXi& covered_faces
);

/**
 * Generate all possible connected components from the tiles.
*/
void connected_components(
	const Eigen::MatrixXi& tiles,
 	std::vector<std::vector<int>>& sub_meshes
);

/**
 * Determine whether input connected component has a 
 * bijection to the original mesh.
*/
void is_equivalence(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	const std::vector<int>& candidate,
	const Eigen::MatrixXi& tiles,
	const Eigen::MatrixXi& covered_faces,
	std::vector<int>& v_old,
	Eigen::MatrixXi& F_coarse,
	Eigen::MatrixXi& fids_covered_by_F_coarse
);
