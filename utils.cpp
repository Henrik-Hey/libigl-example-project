#include "utils.h"
#include <iostream>
#include <map>
#include <vector>

void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
){
	for(int f=0; f<F.rows(); f++)
	{
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		incident_faces[std::make_pair(std::min(v1,v2),std::max(v1,v2))].emplace_back(f);
		incident_faces[std::make_pair(std::min(v2,v3),std::max(v2,v3))].emplace_back(f);
		incident_faces[std::make_pair(std::min(v3,v1),std::max(v3,v1))].emplace_back(f);
	}
};

void get_boundary_vertices(
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces,
	std::vector<int>& boundary_vertices,
	std::map<int, std::vector<int>>& neighbouring_vertices
){
	std::map<std::pair<int,int>, std::vector<int>>::iterator it = incident_faces.begin();
	while (it != incident_faces.end())
	{
    if(it->second.size()==1)
    {
      int v1 = it->first.first;
			if( std::find(boundary_vertices.begin(), boundary_vertices.end(), v1) 
					== boundary_vertices.end() 
			){
				boundary_vertices.emplace_back(v1);
			}
      int v2 = it->first.second;
			if( std::find(boundary_vertices.begin(), boundary_vertices.end(), v2) 
					== boundary_vertices.end() 
			){
				boundary_vertices.emplace_back(v2);
			}
			neighbouring_vertices[v1].emplace_back(v2);
			neighbouring_vertices[v2].emplace_back(v1);
    }
		it++;
  }
};

int find_boundary_vnew(
	const int& vold1,
	const int& vold2,
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
	std::map<std::pair<int,int>, std::vector<int>> incident_faces // should be const
){
	
	assert(incident_faces[std::make_pair(std::min(vold1, vold2), 
												std::max(vold1, vold2))].size()==1);

	int v1new;
	int coarse_face_id = incident_faces[std::make_pair(
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

void map_bound_vnew_to_bound_vold(
	const Eigen::MatrixXi& F_fine,
	const Eigen::MatrixXi& fids_covered_by_F_coarse,
  const std::map<std::pair<int,int>, std::vector<int>>& edgemap_coarse,
  const std::map<int, std::vector<int>>& neighbours_coarse,
  const std::vector<int>& boundary_vids_coarse,
	std::map<int, std::vector<int>>& bound_vnew_to_bound_volds
){

	// Iterate over each Vold boundary vert and 
	// find the two Vnew verts that are its neighbours.
	// While doing so, keep track of which boundary verts in Vold 
	// neighbour the boundary Vnew verts

  int v1new, v2new, vold;
	for( // Iterate over the vold boundary vertices
    std::vector<int>::const_iterator it = boundary_vids_coarse.begin();
    it != boundary_vids_coarse.end();
    it++
  ){

    // Get neighbouring vnew boundary vertex Number 1
    assert(neighbours_coarse.at(*it).size()==2);

    vold = *it;
    v1new = find_boundary_vnew(
      vold,
      neighbours_coarse.at(*it)[0],
      F_fine,
      fids_covered_by_F_coarse,
      edgemap_coarse
    );
    v2new = find_boundary_vnew(
      vold,
      neighbours_coarse.at(*it)[1],
      F_fine,
      fids_covered_by_F_coarse,
      edgemap_coarse
    );

		bound_vnew_to_bound_volds[v1new].emplace_back(vold);
		bound_vnew_to_bound_volds[v2new].emplace_back(vold);
  }
};

void get_fig216f_map(
	const Eigen::MatrixXi& v_is_old,
	const std::map<std::pair<int,int>, std::vector<int>>& edgemap_fine,
	const std::map<int, std::vector<int>>& neighbours_fine,
	std::map<int, std::vector<int>>& fig_216f_map
){

	int edge_v1, edge_v2;
	for(
    std::map<std::pair<int,int>, std::vector<int>>::const_iterator it = edgemap_fine.begin();
    it != edgemap_fine.end();
    it++
  ){
		if(it->second.size()==2)
		{
			// Apply to both verts that form the 
			// current regular edge.
			// First in pair
			edge_v1 = v_is_old(it->first.first,0);
			edge_v2 = v_is_old(it->first.second,0);
			if(v_is_old(edge_v1,0)==0)
			{
				// Iterate over its neighbours and 
				// store the ones that are in Vold
				for(
					std::vector<int>::const_iterator it_n = neighbours_fine.at(edge_v1).begin();
					it_n != neighbours_fine.at(edge_v1).end();
					it_n++
				){
					if(v_is_old(*it_n,0)==1)
					{
						fig_216f_map[edge_v1].emplace_back(*it_n);
					}
				}
				assert(fig_216f_map[edge_v1].size()==4);
			}
			// Do the same thing for the other vert in the curr reg edge
			if(v_is_old(edge_v2,0)==0)
			{
				// Iterate over its neighbours and 
				// store the ones that are in Vold
				for(
					std::vector<int>::const_iterator it_n = neighbours_fine.at(edge_v2).begin();
					it_n != neighbours_fine.at(edge_v2).end();
					it_n++
				){
					if(v_is_old(*it_n,0)==1)
					{
						fig_216f_map[edge_v2].emplace_back(*it_n);
					}
				}
				assert(fig_216f_map[edge_v2].size()==4);
			}
		}
  }
};

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

