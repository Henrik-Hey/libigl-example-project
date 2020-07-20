#include <Eigen/Core>

/**
 * Lifting I:
 * For each boundary vertex in Vold
 *    Grab the neighbouring boundary vertices in Vnew
 *    Update the boundary vertex in Vold with the given equation
*/
void fwt_lifting1 (
	const Eigen::MatrixXi F_old,
	const Eigen::MatrixXd V_old,
	const Eigen::MatrixXi F_new,
	const Eigen::MatrixXd V_new,
  const Eigen::MatrixXi covered_faces // which faces in F_new 
);

/**
 * Lifting II
 * For each boundary vertex in Vnew
 * 		Get the neighbouring boundary vertices in Vold
 * 		Update the position of the Vnew boundary vertex by given equation
*/

/**                 
 * Lifting III
 * Note: we make the assumption that we are updating
 * both boundary and interior vertices in Vold by 
 * vertices in Vnew, not just the interior vertices.
 * Ie., one of the neightburs in Vnew can also be boundary vertices
 * 
 * For each vertex in Vold
 * 		Grab the n neighbours from Vnew 
 * 		Update the vertex in Vold with the given equation
*/

/**                 
 * Scaling
 * For each interior vertex in Vold
 * 		Get the valence of the current interior vertex from Vold
 * 		Update the position of the vertex by scaling it by the given scalar
 * 			which is a function of the scalar
*/

/**                 
 * Lifting IV
 * For each interior vertex in Vnew
 * 		Grab the 4 vertices in Vold which enclose the interior vertex of Vnew
 * 		Update the position of the interior vertex in Vnew by the given equation
*/

/**                 
 * Lifting V
 * For each boundary vertex in Vnew
 * 		Grab the 4 vertices in Vold specified by Fig2.16e
 * 		Update the values of those 4 vertices by the given equation
*/

/**                 
 * Lifting VI
 * For each interior vertex in Vnew
 * 		Grab the 4 Vold enclosing vertices
 * 		For each of the enclosing vertices (i=1,2,3,4)
 * 				Compute the valence
 * 				Use valence to compute alpha, beta, gamma, delta, and omega
 * 				Update the position of the enclosing vertex given the equation
*/
