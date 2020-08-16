#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <map>
#include "utils.h"
#include "wei_shangyang_244.h"
#include "taubin.h"

void lifting1(
		const std::map<int, std::vector<int>> &boundary_vold_to_vnew_map,
		Eigen::MatrixXd &V,
		const bool isFWT)
{
	// Iterate over the vold boundary vertices,
	// find its neighbouring Vnew bound verts and update
	// the coords of the vold boundary verts with the neighbours found
	int vnew1, vnew2, vold;
	Eigen::Vector3d vold_prime;
	for (
			std::map<int, std::vector<int>>::const_iterator it = boundary_vold_to_vnew_map.begin();
			it != boundary_vold_to_vnew_map.end();
			it++)
	{
		assert(it->second.size() == 2);
		vold = it->first;
		vnew1 = it->second[0];
		vnew2 = it->second[1];
		if (isFWT)
		{
			WT_Lifting_1(
					Eigen::Vector3d(V.row(vold)),
					Eigen::Vector3d(V.row(vnew1)),
					Eigen::Vector3d(V.row(vnew2)),
					vold_prime);
		}
		else
		{
			WT_Lifting_1_inverse(
					Eigen::Vector3d(V.row(vold)),
					Eigen::Vector3d(V.row(vnew1)),
					Eigen::Vector3d(V.row(vnew2)),
					vold_prime);
		}
		V.row(vold) = vold_prime;
	}
	std::cout << "Completed FWT Lifting 1" << std::endl;
};

void lifting2(
		const std::map<int, std::vector<int>> &bound_vnew_to_bound_volds,
		Eigen::MatrixXd &V,
		const bool isFWT)
{
	int vold1, vold2, vnew;
	Eigen::Vector3d vnew_prime;
	for (
			std::map<int, std::vector<int>>::const_iterator it = bound_vnew_to_bound_volds.begin();
			it != bound_vnew_to_bound_volds.end();
			it++)
	{
		assert(it->second.size() == 2);
		vnew = it->first;
		vold1 = it->second[0];
		vold2 = it->second[1];

		if (isFWT)
		{
			WT_Lifting_2(
					Eigen::Vector3d(V.row(vnew)),
					Eigen::Vector3d(V.row(vold1)),
					Eigen::Vector3d(V.row(vold2)),
					vnew_prime);
		}
		else
		{
			WT_Lifting_2_inverse(
					Eigen::Vector3d(V.row(vnew)),
					Eigen::Vector3d(V.row(vold1)),
					Eigen::Vector3d(V.row(vold2)),
					vnew_prime);
		}

		V.row(vnew) = vnew_prime;
	}
	std::cout << "Completed FWT Lifting 2" << std::endl;
};

void lifting3(
		const std::map<int, std::vector<int>> &neighbours,
		const Eigen::MatrixXi &v_is_old,
		const Eigen::MatrixXi &v_is_boundary,
		Eigen::MatrixXd &V,
		const bool isFWT)
{
	int vold_int, vnew_int;
	Eigen::Vector3d vold;
	Eigen::Vector3d vold_prime;
	std::vector<int> vold_all_neighbours;
	std::vector<Eigen::Vector3d> vold_int_neighbours;

	for (
			std::map<int, std::vector<int>>::const_iterator it = neighbours.begin();
			it != neighbours.end();
			it++)
	{
		vold_int = it->first;
		if (v_is_old(vold_int, 0) == 1 && v_is_boundary(vold_int, 0) == 0)
		{
			vold_all_neighbours = it->second;
			vold_int_neighbours.clear();

			for (
					std::vector<int>::const_iterator n_it = vold_all_neighbours.begin();
					n_it != vold_all_neighbours.end();
					n_it++)
			{
				assert(v_is_old(*n_it, 0) == 0);
				vold_int_neighbours.emplace_back(
						Eigen::Vector3d(V.row(*n_it)));
			}

			if (vold_int_neighbours.size() > 0)
			{
				vold = V.row(vold_int);

				if (isFWT)
				{
					WT_Lifting_3(
							vold_int_neighbours,
							vold,
							vold_prime);
				}
				else
				{
					WT_Lifting_3_inverse(
							vold_int_neighbours,
							vold,
							vold_prime);
				}
				V.row(vold_int) = vold_prime;
			}
		}
	}
	std::cout << "Completed FWT Lifting 3" << std::endl;
}

void fwt_scaling(
		const Eigen::MatrixXi &v_is_old,
		const Eigen::MatrixXi &v_is_boundary,
		const std::map<int, std::vector<int>> &neighbours_coarse,
		Eigen::MatrixXd &V,
		const bool isFWT)
{
	Eigen::Vector3d v_prime;
	for (int v = 0; v < V.rows(); v++)
	{
		if (v_is_old(v, 0) == 1 && v_is_boundary(v, 0) == 0)
		{
			if (isFWT)
			{
				WT_Scaling(
						Eigen::Vector3d(V.row(v)),
						neighbours_coarse.at(v).size(),
						v_prime);
			}
			else
			{
				WT_Scaling_inverse(
						Eigen::Vector3d(V.row(v)),
						neighbours_coarse.at(v).size(),
						v_prime);
			}

			V.row(v) = v_prime.transpose();
		}
	}
	std::cout << "Completed FWT Scaling" << std::endl;
};

void lifting4(
		const Eigen::MatrixXi &v_is_old,
		const Eigen::MatrixXi &v_is_boundary,
		const std::map<int, std::vector<int>> &neighbours_fine,
		const std::map<int, std::vector<int>> &neighbours_coarse,
		Eigen::MatrixXd &V,
		const bool isFWT)
{
	Eigen::Vector3d v_prime;
	int v_0, v_1, v_2, v_3;
	for (int v = 0; v < V.rows(); v++)
	{
		if (v_is_old(v, 0) == 0 && v_is_boundary(v, 0) == 0)
		{
			get_fig216f_map_with_a_splash_of_henrik(
					v,
					neighbours_fine,
					neighbours_coarse,
					v_is_old,
					v_0,
					v_1,
					v_2,
					v_3);

			if (isFWT)
			{
				WT_Lifting_4(
						Eigen::Vector3d(V.row(v)),
						Eigen::Vector3d(V.row(v_0)),
						Eigen::Vector3d(V.row(v_1)),
						Eigen::Vector3d(V.row(v_2)),
						Eigen::Vector3d(V.row(v_3)),
						v_prime);
			}
			else
			{
				WT_Lifting_4_inverse(
						Eigen::Vector3d(V.row(v)),
						Eigen::Vector3d(V.row(v_0)),
						Eigen::Vector3d(V.row(v_1)),
						Eigen::Vector3d(V.row(v_2)),
						Eigen::Vector3d(V.row(v_3)),
						v_prime);
			}
			V.row(v) = v_prime;
		}
	}
	std::cout << "Completed FWT Lifting 4" << std::endl;
};

void lifting5(
		const std::map<int, std::vector<int>> &boundary_vnew_to_vold_map,
		const std::map<int, std::vector<int>> &boundary_vold_to_vnew_map,
		Eigen::MatrixXd &V,
		const bool isFWT)
{
	int vold1, vold2, vold3, vold4, vnew, v2temp, v3temp;
	Eigen::Vector3d vold1_prime, vold2_prime, vold3_prime, vold4_prime;
	for (
			std::map<int, std::vector<int>>::const_iterator it = boundary_vnew_to_vold_map.begin();
			it != boundary_vnew_to_vold_map.end();
			it++)
	{
		assert(it->second.size() == 2);

		if (it->second.size() == 2)
		{
			vnew = it->first;
			// Following diagram in Fig2.16e
			vold2 = it->second[0];
			vold3 = it->second[1];

			// Find coarse neighbours of vold 2 and vold 3
			v2temp = boundary_vold_to_vnew_map.at(vold2)[0] == vnew ? boundary_vold_to_vnew_map.at(vold2)[1] : boundary_vold_to_vnew_map.at(vold2)[0];
			vold1 = boundary_vnew_to_vold_map.at(v2temp)[0] == vold2 ? boundary_vnew_to_vold_map.at(v2temp)[1] : boundary_vnew_to_vold_map.at(v2temp)[0];

			v3temp = boundary_vold_to_vnew_map.at(vold3)[0] == vnew ? boundary_vold_to_vnew_map.at(vold3)[1] : boundary_vold_to_vnew_map.at(vold3)[0];
			vold4 = boundary_vnew_to_vold_map.at(v3temp)[0] == vold3 ? boundary_vnew_to_vold_map.at(v3temp)[1] : boundary_vnew_to_vold_map.at(v3temp)[0];

			if (isFWT)
			{
				WT_Lifting_5(
						Eigen::Vector3d(V.row(vnew)),
						Eigen::Vector3d(V.row(vold1)),
						Eigen::Vector3d(V.row(vold2)),
						Eigen::Vector3d(V.row(vold3)),
						Eigen::Vector3d(V.row(vold4)),
						vold1_prime,
						vold2_prime,
						vold3_prime,
						vold4_prime);
			}
			else
			{
				WT_Lifting_5_inverse(
						Eigen::Vector3d(V.row(vnew)),
						Eigen::Vector3d(V.row(vold1)),
						Eigen::Vector3d(V.row(vold2)),
						Eigen::Vector3d(V.row(vold3)),
						Eigen::Vector3d(V.row(vold4)),
						vold1_prime,
						vold2_prime,
						vold3_prime,
						vold4_prime);
			}

			V.row(vold1) = vold1_prime;
			V.row(vold2) = vold2_prime;
			V.row(vold3) = vold3_prime;
			V.row(vold4) = vold4_prime;
		}
	}
	std::cout << "Completed FWT Lifting 5" << std::endl;
};

void lifting6(
		const Eigen::MatrixXi &v_is_old,
		const Eigen::MatrixXi &v_is_boundary,
		const std::map<int, std::vector<int>> &neighbours_fine,
		const std::map<int, std::vector<int>> &neighbours_coarse,
		Eigen::MatrixXd &V,
		const bool isFWT)
{
	Eigen::Vector3d v_i_prime;
	int v_0, v_1, v_2, v_3;
	for (int v = 0; v < V.rows(); v++)
	{
		if (v_is_old(v, 0) == 0 && v_is_boundary(v, 0) == 0)
		{
			get_fig216f_map_with_a_splash_of_henrik(
					v,
					neighbours_fine,
					neighbours_coarse,
					v_is_old,
					v_0,
					v_1,
					v_2,
					v_3);

			Eigen::Vector4d W;
			WT_Solve_Weights(
					neighbours_fine.at(v_0).size(),
					neighbours_fine.at(v_1).size(),
					neighbours_fine.at(v_2).size(),
					neighbours_fine.at(v_3).size(),
					W);

			if (isFWT)
			{
				WT_Lifting_6(
						Eigen::Vector3d(V.row(v_0)), // vold
						W(0),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
				V.row(v_0) = v_i_prime;

				WT_Lifting_6(
						Eigen::Vector3d(V.row(v_1)), // vold
						W(1),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
				V.row(v_1) = v_i_prime;

				WT_Lifting_6(
						Eigen::Vector3d(V.row(v_2)), // vold
						W(2),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
				V.row(v_2) = v_i_prime;

				WT_Lifting_6(
						Eigen::Vector3d(V.row(v_3)), // vold
						W(3),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
			}
			else
			{
				WT_Lifting_6_inverse(
						Eigen::Vector3d(V.row(v_0)), // vold
						W(0),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
				V.row(v_0) = v_i_prime;

				WT_Lifting_6_inverse(
						Eigen::Vector3d(V.row(v_1)), // vold
						W(1),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
				V.row(v_1) = v_i_prime;

				WT_Lifting_6_inverse(
						Eigen::Vector3d(V.row(v_2)), // vold
						W(2),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
				V.row(v_2) = v_i_prime;

				WT_Lifting_6_inverse(
						Eigen::Vector3d(V.row(v_3)), // vold
						W(3),
						Eigen::Vector3d(V.row(v)), // vnew
						v_i_prime);
			}
			V.row(v_3) = v_i_prime;
		}
	}
	std::cout << "Completed FWT Lifting 6" << std::endl;
};

void get_aux_data(
		const Eigen::MatrixXd &V_fine,
		const Eigen::MatrixXi &F_fine,
		const Eigen::MatrixXi &F_coarse,
		const Eigen::MatrixXi &v_is_old,
		Eigen::MatrixXi &v_is_boundary,
		std::map<std::pair<int, int>, std::vector<int>> &edgemap_coarse,
		std::map<std::pair<int, int>, std::vector<int>> &edgemap_fine,
		std::map<int, std::vector<int>> &neighbours_coarse,
		std::map<int, std::vector<int>> &neighbours_fine,
		std::map<int, std::vector<int>> &boundary_vold_to_vnew_map,
		std::map<int, std::vector<int>> &boundary_vnew_to_vold_map)
{
	get_edgemap_and_neighbourhoods(
			F_coarse,
			edgemap_coarse,
			neighbours_coarse);
	get_edgemap_and_neighbourhoods(
			F_fine,
			edgemap_fine,
			neighbours_fine);
	get_boundary_vert_structures(
			edgemap_fine,
			v_is_old,
			v_is_boundary,
			boundary_vold_to_vnew_map,
			boundary_vnew_to_vold_map);

	// Sanity all the pre-lifting helper functions
	int countNumNewBoundVerts = 0;
	int countNumOldBoundVerts = 0;
	for (int i = 0; i < V_fine.rows(); i++)
	{
		if (v_is_old(i, 0) == 1 && v_is_boundary(i, 0) == 1)
		{
			countNumOldBoundVerts++;
		}
		else if (v_is_old(i, 0) == 0 && v_is_boundary(i, 0) == 1)
		{
			countNumNewBoundVerts++;
		}
	}
	std::cout << "number of boundaries in fine: " << countNumNewBoundVerts + countNumOldBoundVerts << std::endl;
	for (
			std::map<std::pair<int, int>, std::vector<int>>::iterator it = edgemap_fine.begin();
			it != edgemap_fine.end();
			it++)
	{
		if (it->second.size() > 2)
			std::cout << "found irregular edge" << std::endl;
	}
	assert(boundary_vold_to_vnew_map.size() == countNumOldBoundVerts);
	assert(boundary_vnew_to_vold_map.size() == countNumNewBoundVerts);
	assert(neighbours_fine.size() == V_fine.rows());
};

void fwd_lifting_scheme(
		const std::map<int, std::vector<int>> &neighbours_coarse,
		const std::map<int, std::vector<int>> &neighbours_fine,
		const std::map<int, std::vector<int>> &boundary_vold_to_vnew_map,
		const std::map<int, std::vector<int>> &boundary_vnew_to_vold_map,
		const Eigen::MatrixXi &v_is_boundary,
		const Eigen::MatrixXi &v_is_old,
		Eigen::MatrixXd &V_fine)
{
	lifting1(
			boundary_vold_to_vnew_map,
			V_fine,
			true);
	lifting2(
			boundary_vnew_to_vold_map,
			V_fine,
			true);
	Eigen::MatrixXd V_after_l2 = Eigen::MatrixXd(V_fine);
	lifting3(
			neighbours_fine,
			v_is_old,
			v_is_boundary,
			V_fine,
			true);
	Eigen::MatrixXd V_after_l3 = Eigen::MatrixXd(V_fine);
	if (V_after_l2 == V_after_l3)
		std::cout << "L3 did nothing." << std::endl;
	fwt_scaling(
			v_is_old,
			v_is_boundary,
			neighbours_coarse,
			V_fine,
			true);
	Eigen::MatrixXd V_after_scale = Eigen::MatrixXd(V_fine);
	lifting4(
			v_is_old,
			v_is_boundary,
			neighbours_fine,
			neighbours_coarse,
			V_fine,
			true);
	lifting5(
			boundary_vnew_to_vold_map,
			boundary_vold_to_vnew_map,
			V_fine,
			true);
	lifting6(
			v_is_old,
			v_is_boundary,
			neighbours_fine,
			neighbours_coarse,
			V_fine,
			true);
	Eigen::MatrixXd V_after_l6 = Eigen::MatrixXd(V_fine);
	if (V_after_scale == V_after_l6)
		std::cout << "L4-L6 collectively did nothing.";
};

void inv_lifting_scheme(
		const std::map<int, std::vector<int>> &neighbours_coarse,
		const std::map<int, std::vector<int>> &neighbours_fine,
		const std::map<int, std::vector<int>> &boundary_vold_to_vnew_map,
		const std::map<int, std::vector<int>> &boundary_vnew_to_vold_map,
		const Eigen::MatrixXi &v_is_boundary,
		const Eigen::MatrixXi &v_is_old,
		Eigen::MatrixXd &V_fine)
{

	lifting6(
			v_is_old,
			v_is_boundary,
			neighbours_fine,
			neighbours_coarse,
			V_fine,
			false);
	lifting5(
			boundary_vnew_to_vold_map,
			boundary_vold_to_vnew_map,
			V_fine,
			false);
	lifting4(
			v_is_old,
			v_is_boundary,
			neighbours_fine,
			neighbours_coarse,
			V_fine,
			false);
	fwt_scaling(
			v_is_old,
			v_is_boundary,
			neighbours_coarse,
			V_fine,
			false);
	lifting3(
			neighbours_fine,
			v_is_old,
			v_is_boundary,
			V_fine,
			false);
	lifting2(
			boundary_vnew_to_vold_map,
			V_fine,
			false);
	lifting1(
			boundary_vold_to_vnew_map,
			V_fine,
			false);
	std::cout << "L4-L6 collectively did nothing.";
};

bool FWT(
		const Eigen::MatrixXi &F_fine,
		Eigen::MatrixXi &F_coarse,
		Eigen::MatrixXd &V_fine)
{
	Eigen::MatrixXi v_is_old;									// Map vid to 1 if in Vold, and 0 if Vnew
	Eigen::MatrixXi fids_covered_by_F_coarse; // #F_coarse x 4 fids in F_in that that tile covers
	if (is_quadrisection(
					F_fine,
					V_fine,
					v_is_old,
					F_coarse,
					fids_covered_by_F_coarse))
	{
		// Aux data structures needed for computing transform
		Eigen::MatrixXi v_is_boundary;
		std::map<std::pair<int, int>, std::vector<int>> edgemap_coarse;
		std::map<std::pair<int, int>, std::vector<int>> edgemap_fine;
		std::map<int, std::vector<int>> neighbours_coarse;
		std::map<int, std::vector<int>> neighbours_fine;
		std::map<int, std::vector<int>> boundary_vold_to_vnew_map;
		std::map<int, std::vector<int>> boundary_vnew_to_vold_map;
		get_aux_data(
				V_fine,
				F_fine,
				F_coarse,
				v_is_old,
				v_is_boundary,
				edgemap_coarse,
				edgemap_fine,
				neighbours_coarse,
				neighbours_fine,
				boundary_vold_to_vnew_map,
				boundary_vnew_to_vold_map);

		// FWT time!
		fwd_lifting_scheme(
				neighbours_coarse,
				neighbours_fine,
				boundary_vold_to_vnew_map,
				boundary_vnew_to_vold_map,
				v_is_boundary,
				v_is_old,
				V_fine);

		inv_lifting_scheme(
				neighbours_coarse,
				neighbours_fine,
				boundary_vold_to_vnew_map,
				boundary_vnew_to_vold_map,
				v_is_boundary,
				v_is_old,
				V_fine);
		return 1;
	}
	else
		return 0;
};
