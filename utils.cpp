#include "utils.h"
#include <iostream>
#include <map>
#include <vector>

void get_boundary_vert_structures(
	const std::map<std::pair<int,int>, std::vector<int>>& edgemap_fine,
	const Eigen::MatrixXi& v_is_old,
				Eigen::MatrixXi& v_is_boundary,
				std::map<int, std::vector<int>>& boundary_vold_to_vnew_map,
				std::map<int, std::vector<int>>& boundary_vnew_to_vold_map
){
	int v1, v2; // one will be in vnew and one in vold
	v_is_boundary = Eigen::MatrixXi::Zero(v_is_old.rows(),1);
	std::map<std::pair<int,int>, std::vector<int>>::const_iterator it = edgemap_fine.begin();
	while (it != edgemap_fine.end())
	{
		// If the current edge only has one incident face
    if(it->second.size()==1)
    {
			v1 = it->first.first;
			v2 = it->first.second;
			v_is_boundary(v1,0) = 1;
			v_is_boundary(v2,0) = 1;
			if(v_is_old(v1,0)==1)
			{
				assert(v_is_old(v2,0)==0);
				// v1 is in vold and v2 is in vnew
				boundary_vold_to_vnew_map[v1].emplace_back(v2);
				boundary_vnew_to_vold_map[v2].emplace_back(v1);

			}
			else if(v_is_old(v2,0)==1)
			{
				assert(v_is_old(v1,0)==0);
				// v2 is in vold and v1 is in vnew
				boundary_vold_to_vnew_map[v2].emplace_back(v1);
				boundary_vnew_to_vold_map[v1].emplace_back(v2);
			}
    }
		it++;
  }
};

void edge_incident_faces(
	const Eigen::MatrixXi& F,
				std::map<std::pair<int,int>, std::vector<int>>& edgemap
){
	for(int f=0; f<F.rows(); f++)
	{
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		edgemap[std::make_pair(std::min(v1,v2),std::max(v1,v2))].emplace_back(f);
		edgemap[std::make_pair(std::min(v2,v3),std::max(v2,v3))].emplace_back(f);
		edgemap[std::make_pair(std::min(v3,v1),std::max(v3,v1))].emplace_back(f);
	}
};

void get_edgemap_and_neighbourhoods(
	const Eigen::MatrixXi& F,
				std::map<std::pair<int,int>, std::vector<int>>& edgemap,
				std::map<int, std::vector<int>>& neighbours
){
	for(int f=0; f<F.rows(); f++)
	{
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		edgemap[std::make_pair(std::min(v1,v2),std::max(v1,v2))].emplace_back(f);
		edgemap[std::make_pair(std::min(v2,v3),std::max(v2,v3))].emplace_back(f);
		edgemap[std::make_pair(std::min(v3,v1),std::max(v3,v1))].emplace_back(f);

		if( std::find(neighbours[v1].begin(),neighbours[v1].end(),v2)==neighbours[v1].end()) 
			neighbours[v1].emplace_back(v2);
		if( std::find(neighbours[v1].begin(),neighbours[v1].end(),v3)==neighbours[v1].end()) 
			neighbours[v1].emplace_back(v3);

		if( std::find(neighbours[v2].begin(),neighbours[v2].end(),v1)==neighbours[v2].end()) 
			neighbours[v2].emplace_back(v1);
		if( std::find(neighbours[v2].begin(),neighbours[v2].end(),v3)==neighbours[v2].end()) 
			neighbours[v2].emplace_back(v3);

		if( std::find(neighbours[v3].begin(),neighbours[v3].end(),v2)==neighbours[v3].end()) 
			neighbours[v3].emplace_back(v2);
		if( std::find(neighbours[v3].begin(),neighbours[v3].end(),v1)==neighbours[v3].end()) 
			neighbours[v3].emplace_back(v1);
	}
};

int find_boundary_vnew(
	const int& vold1,
	const int& vold2,
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
				std::map<std::pair<int,int>, std::vector<int>> edgemap // should be const
){
	
	assert(edgemap[std::make_pair(std::min(vold1, vold2), 
												std::max(vold1, vold2))].size()==1);

	int v1new;
	int coarse_face_id = edgemap[std::make_pair(
											 std::min(vold1, vold2), std::max(vold1, vold2))][0];
	
	int fid1, fid2, vid11, vid12, vid21, vid22;
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(F_in(fids_covered_by_F_coarse(coarse_face_id,i),j)==vold1)
			{
				fid1 = fids_covered_by_F_coarse(coarse_face_id,i);
			}
			if(F_in(fids_covered_by_F_coarse(coarse_face_id,i),j)==vold2)
			{
				fid2 = fids_covered_by_F_coarse(coarse_face_id,i);
			}
		}
	}

	vid11 = -1;
	vid12 = -1;
	for(int j=0; j<3; j++)
	{
		if(F_in(fid1,j)!=vold1 && vid11==-1) vid11 = F_in(fid1,j);
		else if(F_in(fid1,j)!=vold1&& vid12==-1) vid12 = F_in(fid1,j);
	}

	vid21 = -1;
	vid22 = -1;
	for(int j=0; j<3; j++)
	{
		if(F_in(fid2,j)!=vold2 && vid21==-1) vid21 = F_in(fid2,j);
		else if(F_in(fid2,j)!=vold2 && vid22==-1) vid22 = F_in(fid2,j);
	}

	if(vid11==vid21) v1new = vid11;
	else if(vid11==vid22) v1new = vid11;
	else if(vid12==vid21) v1new = vid12;
	else if(vid12==vid22) v1new = vid12;

	return v1new;
};

void get_fig216f_map_with_a_splash_of_henrik(
	const int v_new_id,
	const std::map<int, std::vector<int>>& neighbours_fine,
	const std::map<int, std::vector<int>>& neighbours_coarse,
  const Eigen::MatrixXi& v_is_old,
	  		int& v_0,
	  		int& v_1,
	  		int& v_2,
	  		int& v_3
){
	bool found_v_0 = false;
	bool found_v_2 = false;
	int v_0_id = -1, 
			v_1_id = -1, 
			v_2_id = -1, 
			v_3_id = -1;

	std::vector<int> vnew_relevant_neighbours = neighbours_fine.find(v_new_id)->second;
	for(
    std::vector<int>::iterator n_it = vnew_relevant_neighbours.begin();
    n_it != vnew_relevant_neighbours.end();
    n_it++
  ){
		if(v_is_old(*n_it, 0))
		{
			if(!found_v_0) 
			{
				v_0_id = *n_it;
				found_v_0 = true;
			}
			else
			{
				v_1_id = *n_it;
				break;
			}
		}
	}

	assert(v_0_id != -1 && v_1_id != -1);

	std::vector<int>v0_relevant_neighbours = neighbours_coarse.find(v_0_id)->second;
	std::vector<int>v23_relevant_neighbours;
	for(
    std::vector<int>::iterator n_it = v0_relevant_neighbours.begin();
    n_it != v0_relevant_neighbours.end();
    n_it++
  ){
		if(*n_it != v_1_id && v_is_old(*n_it, 0))
		{
			// std::cout << "got here"  << std::endl;
			v23_relevant_neighbours = neighbours_coarse.find(*n_it)->second;
			if (
				std::find(
					v23_relevant_neighbours.begin(), 
					v23_relevant_neighbours.end(),
					v_1_id
				) != v23_relevant_neighbours.end()
			){
					if(!found_v_2)
					{
						v_2_id = *n_it;
						found_v_2 = true; 
					}
					else 
					{
						v_3_id = *n_it; 
					}
			}
		}
		if(v_2_id != -1 && v_3_id != -1)
		{
			break;
		}
	}

	assert(v_2_id != -1 && v_3_id != -1);

	v_0 = v_0_id;
	v_1 = v_1_id;
	v_2 = v_2_id;
	v_3 = v_3_id;
}

void sort3(int arr[]) 
{ 
	// Insert arr[1] 
	if (arr[1] < arr[0]) 
		std::swap(arr[0], arr[1]); 

	// Insert arr[2] 
	if (arr[2] < arr[1]) 
	{ 
		std::swap(arr[1], arr[2]); 
		if (arr[1] < arr[0]) 
			std::swap(arr[1], arr[0]); 
	} 
};
