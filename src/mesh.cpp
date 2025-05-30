#include "mesh.h"
#include <set>
#include "igl/per_face_normals.h"
#include "igl/principal_curvature.h"
#include "igl/readSTL.h"
#include "drogon/drogon.h"
#include "igl/remove_duplicate_vertices.h"

#include "time.h"
// #include "igl/decimate.h"

namespace CrownGen
{

	Mesh::Mesh() {}

	void Mesh::curvatureExtraction()
	{
		int NUM_NEIGHBORS = 4;
		Eigen::MatrixXd decimated_vertices;
		Eigen::MatrixXi decimated_faces;
		// igl::decimate(this->Vertices, this->Faces, 12000, decimated_vertices, decimated_faces);
		igl::principal_curvature(this->Vertices, this->Faces, this->verticesNormal, this->facesNormal_, this->K1, this->K2, this->bad_vertices, NUM_NEIGHBORS, false);
		// std::cout<< "bad vertices:"<<std::endl;
		// for (int i : this->bad_vertices)
		// {
		// 	std::cout << i<< std::endl;
		// }
		// std::cout<<std::endl;
		// igl::per_face_normals(this->Vertices, this->Faces, this->facesNormal);
		std::cout << "Number of bad vertices: " << this->bad_vertices.size() << std::endl;
		std::cout << "Curvatures calculated for mesh" << std::endl;
	}

	void Mesh::setNumVertices(int numVertices)
	{
		this->numVertices = numVertices;
	}

	void Mesh::setNumFaces(int numFaces)
	{
		this->numFaces = numFaces;
	}

	void Mesh::calculateAdjacencyVerticesFaces()
	{
		LOG_INFO << "Calculating Adjacency Vertices and Faces";
		// add adjacency_vertices for each vertex
		this->adjacency_vertices.resize(this->numVertices);
		for (int i = 0; i < this->numFaces; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				this->adjacency_vertices[this->Faces(i, j)].push_back(this->Faces(i, (j + 1) % 3));
				this->adjacency_vertices[this->Faces(i, j)].push_back(this->Faces(i, (j + 2) % 3));
			}
		}
		// remove duplicates
		for (auto &adj : this->adjacency_vertices)
		{
			std::set<int> unique_adj(adj.begin(), adj.end());
			adj.assign(unique_adj.begin(), unique_adj.end());
		}

		// add adjacency_faces for each vertices
		this->adjacency_faces.resize(this->numVertices);
		for (int i = 0; i < this->numFaces; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				this->adjacency_faces[this->Faces(i, j)].push_back(i);
			}
		}
		// remove duplicates
		for (auto &adj : this->adjacency_faces)
		{
			std::set<int> unique_adj(adj.begin(), adj.end());
			adj.assign(unique_adj.begin(), unique_adj.end());
		}
	}

	void Mesh::initSTL(std::string filename)
	{
		FILE *stlfile = fopen("data.stl", "rb");
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXd N;
		igl::readSTL(stlfile, V, F, N);

		igl::remove_duplicate_vertices(V,F,10e-3,this->Vertices,this->Faces);

		for (int i = 0; i < 3; i++)
		{
			std::cout << "Vertex " << i << " " << this->Vertices.row(i) << std::endl;
			std::cout << "Faces " << i << " " << this->Faces.row(i) << std::endl;
		}
		this->setNumVertices(this->Vertices.rows());
		LOG_INFO << " Number of vertices: " << this->numVertices;
		this->setNumFaces(this->Faces.rows());
		LOG_INFO << " Number of faces: " << this->numFaces;
		this->calculateAdjacencyVerticesFaces();
		time_t prevtime = time(NULL);
		this->curvatureExtraction();
		time_t nowtime = time(NULL);
		LOG_INFO << "time elapsed: " << nowtime - prevtime;
		// fclose(stlfile);
	}

}
