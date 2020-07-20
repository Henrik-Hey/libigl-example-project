#include "utils.h"

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

void get_edges_and_neighbours(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces,
	std::map<int, std::vector<int>>& neighbouring_vertices
){
	for(int f=0; f<F.rows(); f++)
	{
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		incident_faces[std::make_pair(std::min(v1,v2),std::max(v1,v2))].emplace_back(f);
		incident_faces[std::make_pair(std::min(v2,v3),std::max(v2,v3))].emplace_back(f);
		incident_faces[std::make_pair(std::min(v3,v1),std::max(v3,v1))].emplace_back(f);

    neighbouring_vertices[v1].emplace_back(v2);
    neighbouring_vertices[v1].emplace_back(v3);
    neighbouring_vertices[v2].emplace_back(v1);
    neighbouring_vertices[v2].emplace_back(v3);
    neighbouring_vertices[v3].emplace_back(v1);
    neighbouring_vertices[v3].emplace_back(v2);
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

