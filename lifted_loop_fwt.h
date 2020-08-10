#include <Eigen/Core>
#include <vector>

/**
 * Lifting I: Update Vold boundary verts with Vnew boundary verts
 * For each boundary vertex in Vold
 *    Grab the neighbouring boundary vertices in Vnew
 *    Update the boundary vertex in Vold with the given equation
*/
// void fwt_lifting1 (
// 	const Eigen::MatrixXi& F_fine,
// 	const Eigen::MatrixXi& fids_covered_by_F_coarse,
//   const std::map<std::pair<int,int>, std::vector<int>>& edgemap_coarse,
//   const std::map<int, std::vector<int>>& neighbours_coarse,
//   const std::vector<int>& boundary_vids_coarse,
// 	Eigen::MatrixXd& V
// );
void fwt_lifting1 (
	std::map<int, std::vector<int>>& boundary_vold_to_vnew_map,
	Eigen::MatrixXd& V
);

/**
 * Lifting II: Update Vnew boundary verts with Vold boundary verts
 * For each boundary vertex in Vnew
 * 		Get the neighbouring boundary vertices in Vold
 * 		Update the position of the Vnew boundary vertex by given equation
*/
void fwt_lifting2 (
  const std::map<int, std::vector<int>>& bound_vnew_to_bound_volds,
	Eigen::MatrixXd& V
);

/**                 
 * Lifting III: Update interior vertices in Vold with 1-ring neigbours in Vnew
 * Note: we make the assumption that we are updating
 * both boundary and interior vertices in Vold by 
 * vertices in Vnew, not just the interior vertices.
 * Ie., one of the neightburs in Vnew can also be boundary vertices
 * 
 * For each vertex in Vold
 * 		Grab the n neighbours from Vnew 
 * 		Update the vertex in Vold with the given equation
*/
void fwt_lifting3 (
  const std::map<int, std::vector<int>>& neighbours,
  const Eigen::MatrixXi& v_is_old,
  const Eigen::MatrixXi& v_is_boundary,
        Eigen::MatrixXd& V
);

/**                 
 * Scaling: Scale the interior vertices in Vold
 * For each interior vertex in Vold
 * 		Get the valence of the current interior vertex from Vold
 * 		Update the position of the vertex by scaling it by the given scalar
 * 			which is a function of the scalar
*/
void fwt_scaling (
	const Eigen::MatrixXi& v_is_old,
	const Eigen::MatrixXi& v_is_boundary,
  const std::map<int, std::vector<int>>& neighbours_coarse,
	Eigen::MatrixXd& V
);

/**                 
 * Lifting IV: Update interior verices in Vnew with the specified (fig)
 * 						 neighbouring vertices in Vold (either boundary or interior)
 * For each interior vertex in Vnew
 * 		Grab the 4 vertices in Vold which enclose the interior vertex of Vnew
 * 		Update the position of the interior vertex in Vnew by the given equation
*/
void fwt_lifting4 (
	const Eigen::MatrixXi& v_is_old,
	const Eigen::MatrixXi& v_is_boundary,
  const std::map<int, std::vector<int>>& neighbours_fine,
	Eigen::MatrixXd& V
);

/**                 
 * Lifting V: Update boundary verts in Vold with boundary verts in Vnew
 * For each boundary vertex in Vnew
 * 		Grab the 4 vertices in Vold specified by Fig2.16e
 * 		Update the values of those 4 vertices by the given equation
*/
void fwt_lifting5 (
	const std::map<int, std::vector<int>>& boundary_vnew_to_vold_map,
	const std::map<int, std::vector<int>>& boundary_vold_to_vnew_map,
	Eigen::MatrixXd& V
);

/**                 
 * Lifting VI: Update either boundary or interior verts in Vold
 * 						 with an interior vertex in Vnew
 * For each interior vertex in Vnew
 * 		Grab the 4 Vold enclosing vertices
 * 		For each of the enclosing vertices (i=1,2,3,4)
 * 				Compute the valence
 * 				Use valence to compute alpha, beta, gamma, delta, and omega
 * 				Update the position of the enclosing vertex given the equation
*/
void fwt_lifting6 (
	const Eigen::MatrixXi& v_is_old,
	const Eigen::MatrixXi& v_is_boundary,
  const std::map<int, std::vector<int>>& neighbours_fine,
	Eigen::MatrixXd& V
);
