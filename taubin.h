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
				Eigen::MatrixXi& v_is_old,
				Eigen::MatrixXi& F_coarse,
				Eigen::MatrixXi& fids_covered_by_F_coarse
);

/**
 * Creates tiles from the input mesh, and the tile
 * sets covered by each.
 * 
 * Outputs:
 * tiles - num_tilesx3 matrix of vert ids in V_in
 * covered_faces - # tiles x 4 matrix of fids in F_in
*/
void covering_mesh(
	const Eigen::MatrixXi& F_in,
				Eigen::MatrixXi& tiles,
				Eigen::MatrixXi& covered_faces
);

/**
 * Generate all possible connected components from the tiles.
 * Output:  Vector of vectors containing tile ids in a single connected component found
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
				Eigen::MatrixXi& v_is_old,
				Eigen::MatrixXi& F_coarse,
				Eigen::MatrixXi& fids_covered_by_F_coarse
);

// Hey kid.... want a tetrahedron?
// Inline mesh of a tetrahedron
// OV = (Eigen::MatrixXd(4,3)<<
//   0.0,0.0,0.0,
//   1.0,0.0,0.0,
//   0.66,0.33,1.0,
//   1.0,1.0,0.0).finished();
// OF = (Eigen::MatrixXi(4,3)<<
//   1,2,3,
//   1,3,4,
//   1,4,2,
//   2,4,3).finished().array()-1; // Test
