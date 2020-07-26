#include "taubin.h"
#include "utils.h"
#include <igl/loop.h>
#include <iostream>
#include <vector>
#include <map>

bool is_quadrisection(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	Eigen::MatrixXi& v_is_old,
	Eigen::MatrixXi& F_coarse,
	Eigen::MatrixXi& fids_covered_by_F_coarse
){
  std::cout << "Num verts in input mesh: " << V_in.rows() << std::endl;
  std::cout << "Num faces in input mesh: " << F_in.rows() << std::endl;
  
	// Construct a bunch of tiles
	std::cout << "Begin covering mesh" << std::endl;
  Eigen::MatrixXi tiles; // A tile soup!
  Eigen::MatrixXi covered_faces; // row i is tile i from tiles to 4 fids from F_in
  covering_mesh(F_in, tiles, covered_faces);
	std::cout << "Completed covering mesh" << std::endl;
	// Generate sets of tiles that connect
	std::cout << "Begin connected components" << std::endl;
	std::vector<std::vector<int>> sub_meshes; // Sets of tile soups!
  connected_components(tiles, sub_meshes);
	std::cout << "Completed connected components" << std::endl;
	std::cout << "Number of candidates: " << sub_meshes.size() << std::endl;

	// Find a candidate connected component
	for(auto it=sub_meshes.begin(); it!=sub_meshes.end(); it++)
	{
		// Try to construct a bijection from current
		// connect component to the input mesh
		is_equivalence( F_in, V_in, *it, tiles, covered_faces, v_is_old, F_coarse, fids_covered_by_F_coarse );
		// Stop early as soon as we find a
		// successful candidate
		if(F_coarse.rows()>0) 
			return true;
	}
	std::cout << "Mesh does not have subdivision connectivity" << std::endl;
	return false;
};

void covering_mesh(
	const Eigen::MatrixXi& F_in,
 	Eigen::MatrixXi& tiles, // num_tilesx3 matrix of vert ids in V_in
 	Eigen::MatrixXi& covered_faces // # tiles x 4 matrix of fids in F_in
){
	std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  edge_incident_faces(F_in, incident_faces);
  std::cout << "Num edges: " << incident_faces.size() << std::endl;

	std::vector<std::tuple<int, int, int>> found_tiles;
	std::vector<std::tuple<int, int, int, int>> neighbouring_faces;
	for(int f=0; f<F_in.rows(); f++)
	{
		int v1 = F_in(f,0);
		int v2 = F_in(f,1);
		int v3 = F_in(f,2);

		int e1v1 = std::min(v1,v2);
		int e1v2 = std::max(v1,v2);
		int e2v1 = std::min(v2,v3);
		int e2v2 = std::max(v2,v3);
		int e3v1 = std::min(v3,v1);
		int e3v2 = std::max(v3,v1);

		std::pair<int, int> e1(e1v1, e1v2);
		std::pair<int, int> e2(e2v1, e2v2);
		std::pair<int, int> e3(e3v1, e3v2);

		bool isRegular = 
			   incident_faces[e1].size()==2
			&& incident_faces[e2].size()==2
			&& incident_faces[e3].size()==2;

		if(isRegular)
		{
			// Neighbouring triangles to current triangle
			int nt1 = incident_faces[e1][0]==f ? incident_faces[e1][1] : incident_faces[e1][0];
			int nt2 = incident_faces[e2][0]==f ? incident_faces[e2][1] : incident_faces[e2][0];
			int nt3 = incident_faces[e3][0]==f ? incident_faces[e3][1] : incident_faces[e3][0];

			// Get tile vertices
			int v1p, v2p, v3p;
			for(int i=0; i<3; i++)
			{
				if(F_in(nt1, i)!=e1v1 && F_in(nt1, i)!=e1v2) v1p = F_in(nt1, i);
				if(F_in(nt2, i)!=e2v1 && F_in(nt2, i)!=e2v2) v2p = F_in(nt2, i);
				if(F_in(nt3, i)!=e3v1 && F_in(nt3, i)!=e3v2) v3p = F_in(nt3, i);
			}
			found_tiles.push_back(std::tuple<int, int, int>(v1p,v2p,v3p));
			neighbouring_faces.push_back(std::tuple<int, int, int, int>(f,nt1,nt2,nt3));
		}

		tiles.setIdentity(found_tiles.size(),3);
		covered_faces.setIdentity(found_tiles.size(),4);
		for(int t=0; t < found_tiles.size(); t++)
		{
			tiles(t, 0) = std::get<0>(found_tiles[t]);
			tiles(t, 1) = std::get<1>(found_tiles[t]);
			tiles(t, 2) = std::get<2>(found_tiles[t]);

			covered_faces(t, 0) = std::get<0>(neighbouring_faces[t]);
			covered_faces(t, 1) = std::get<1>(neighbouring_faces[t]);
			covered_faces(t, 2) = std::get<2>(neighbouring_faces[t]);
			covered_faces(t, 3) = std::get<3>(neighbouring_faces[t]);
		}
	}
	std::cout << "Num tiles: " << tiles.rows() << std::endl;
};

void connected_components(
	const Eigen::MatrixXi& tiles,
 	std::vector<std::vector<int>>& sub_meshes // Vector of vectors containing tile ids in a single connected component found
){
	std::map<int, std::vector<int>*> where_are_you; // Which partition is the face in
	for(int tid=0; tid<tiles.rows(); tid++)
	{
		// At first, the partition of the mesh is 
		// a collection of singletons of the fids
		where_are_you[tid] = new std::vector<int>;
		where_are_you[tid]->emplace_back(tid);
	}
	assert(where_are_you.size()==tiles.rows());

	// Get all the edges in the mesh
	std::map<std::pair<int,int>, std::vector<int>> incident_tiles;
  edge_incident_faces(tiles, incident_tiles);

	std::map<std::pair<int,int>, std::vector<int>>::iterator it = incident_tiles.begin();
	while (it != incident_tiles.end())
	{
		// assert((it->second).size()<=2);
		if((it->second).size()==2)
		{
			int tid1 = it->second[0];
			int tid2 = it->second[1];
			if(where_are_you[tid1]!=where_are_you[tid2])
			{
				assert(where_are_you[tid1]!=NULL && where_are_you[2]!=NULL);

				// where_are_you[tid1]->resize(where_are_you[tid1]->size()+where_are_you[tid2]->size());
				where_are_you[tid1]->insert(
					where_are_you[tid1]->end(),
					where_are_you[tid2]->begin(),
					where_are_you[tid2]->end()
				);

				for(int i=0; i<where_are_you[tid1]->size(); i++)
				{
					if(where_are_you[tid1]->at(i)!=tid2 && where_are_you[tid1]->at(i)!=tid1)
					{
						where_are_you[where_are_you[tid1]->at(i)] = where_are_you[tid1];
					}
				}

				delete where_are_you[tid2];
				where_are_you[tid2] = where_are_you[tid1];
			}
			else
			{
				assert(where_are_you[tid1]->size() == where_are_you[tid2]->size());
			}
		}
		it++;
	}

	std::cout << "Find unique submeshes" << std::endl;
	for(int tid=0; tid<tiles.rows(); tid++)
	{
		if(where_are_you[tid] != NULL)
		{
			std::vector<int>* temp = where_are_you[tid];
			sub_meshes.emplace_back(*temp);
			for(size_t i=0; i<temp->size(); i++)
			{
				where_are_you[ temp->at(i) ] = NULL;
			}
			delete temp;
		}
	}
};

void is_equivalence(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	// Tile indices that make up the candidate mesh.
	const std::vector<int>& candidate,
	// #tiles by 3 matrix of indices into V_in that
	// make up each tile.
	const Eigen::MatrixXi& tiles,
	const Eigen::MatrixXi& covered_faces,
	Eigen::MatrixXi& v_is_old,
	Eigen::MatrixXi& F_coarse,
	Eigen::MatrixXi& fids_covered_by_F_coarse
){
	// First test
	if(candidate.size()*4==F_in.rows())
	{
		std::cout << "--Begin analyzing the candidate connected component--" << std::endl;

		// Turn candidate into the proper F, V 
		// matrix format so that we can quadrisect it
		Eigen::MatrixXi submesh; // "Faces" from canadidate
		Eigen::MatrixXi submesh_covered_faces; 
		std::vector<int> V_i; // All the corners of the candidates
		submesh.setIdentity(candidate.size(),3);
		submesh_covered_faces.setIdentity(candidate.size(),4);
		int t=0;
		v_is_old = Eigen::MatrixXi::Zero(V_in.rows(),1);
		for(auto it2=candidate.begin(); it2!=candidate.end(); it2++)
		{
			for(int i=0; i<3; i++)
			{
				// Populate submesh as a matrix with
				// rows being vids of V that form the tiles
				// of the candidate.
				submesh(t, i) = tiles(*it2,i);

				// Populate V_i as a set of unique vids 
				// encountered while iterating over
				// the tiles which form the candidate.
				if( std::find(V_i.begin(), V_i.end(), tiles(*it2,i)) == V_i.end() )
				{
					v_is_old(tiles(*it2,i),0) = 1;
					V_i.emplace_back(tiles(*it2,i));
				}
			}

			submesh_covered_faces(t, 0) = covered_faces(*it2,0);
			submesh_covered_faces(t, 1) = covered_faces(*it2,1);
			submesh_covered_faces(t, 2) = covered_faces(*it2,2);
			submesh_covered_faces(t, 3) = covered_faces(*it2,3);

			t++;
		}

		// At this point:
		// V_i has the vids of V covered by candidate.
		// submesh has #tilesincandidate rows
		// where each row has the 3 vids from 
		// V which make up that tile
		// assert(submesh.rows()==6);

		// We need num verts and num edges in candidate
		// connected component for the second test
		std::map<std::pair<int,int>, std::vector<int>> incident_tiles;
		edge_incident_faces(submesh, incident_tiles);

		// Second test
		if(V_i.size()+incident_tiles.size()==V_in.rows())
		{
			// Make sure that every vert in the subdivided
			// candidate can be found in the original V
			int found = true;
			int temp = false;
			int visited_cout = 0;
			std::map<int,bool> f_visited;
			for(int f=0;f<F_in.rows();f++)
			{
				f_visited[f] = false;
			}

			for(auto it2=candidate.begin(); it2!=candidate.end(); it2++)
			{
				
				for(int i=0; i<4; i++)
				{
					if(!f_visited[covered_faces(*it2,i)])
					{
						f_visited[covered_faces(*it2,i)] = true;
						visited_cout++;
						temp = true;
					}
				}
				if(temp==true)
				{
					temp = false;
				} else {
					std::cout << "Candidate failed" << std::endl;
					found = false;
					break;
				}
			}

			// Third (final) test
			if(found && visited_cout==F_in.rows()) 
			{ 
				std::cout << "Gagnant!" << std::endl; 

				F_coarse = submesh;
				fids_covered_by_F_coarse = submesh_covered_faces;
			}
		} 
		else
		{ 
			std::cout << "Second test in is_equivalence failed." << std::endl; 
		}
	}
	else 
	{
		std::cout << "First test in is_equivalence failed." << std::endl; 
	}
};
