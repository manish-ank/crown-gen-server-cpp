#include "mesh.h"
#include "igl/readSTL.h"
#include "igl/readPLY.h"
#include "igl/readOBJ.h"
#include <set>
#include "igl/per_face_normals.h"
#include "igl/principal_curvature.h"
#include "igl/remove_duplicate_vertices.h"
#include "stdlib.h"
#include "igl/decimate.h"
#include "math.h"
#include "igl/octree.h"
#include "igl/knn.h"
#include <iostream>
#include "igl/writeOFF.h"
#include <drogon/drogon.h>

#include "time.h"
#include "igl/decimate.h"

#include "nlohmann/json.hpp"

namespace CrownGen
{
	Mesh::Mesh() {}

	void Mesh::init(std::string filename, std::string out_filename,std::string out_color_filename)
	{
		// opening and loading PLY files
		FILE *file = fopen((char *)filename.c_str(), "rb");
		Eigen::MatrixXd temp_V;
		Eigen::MatrixXi temp_F;
		Eigen::MatrixXd temp_N;
		std::string extension = filename.substr(filename.find_last_of(".") + 1);
		if (extension == "stl")
			igl::readSTL(file, temp_V, temp_F, temp_N);
		else if (extension == "ply")
			igl::readPLY(file, temp_V, temp_F);
		else if (extension == "obj")
			igl::readOBJ(filename, temp_V, temp_F);
		else
			return;
		this->fileName = filename;
		this->out_fileName = out_filename;
		this->out_color_fileName = out_color_filename;
		// removing duplicate vertices.
		Eigen::VectorXi SVI, SVJ;
		igl::remove_duplicate_vertices(temp_V, temp_F, 10e-10, this->Vertices, SVI, SVJ, this->Faces);

		this->numFaces = this->Faces.rows();
		this->numVertices = this->Vertices.rows();

		// calculating perface normal
		igl::per_face_normals(this->Vertices, this->Faces, this->faceNormal);

		// calculating per vertex normal;
		igl::per_vertex_normals(this->Vertices, this->Faces, this->verticesNormal);
	}

	void Mesh::decimate()
	{
		this->decimated = true;
		this->pre_decimatedFaces = this->Faces;
		this->pre_decimatedVertices = this->Vertices;
		Eigen::MatrixXd decimated_vertices;
		Eigen::MatrixXi decimated_faces;
		int max_faces = this->Faces.rows() * 0.7;

		time_t prev_time = time(NULL);
		igl::decimate(this->Vertices, this->Faces, max_faces, decimated_vertices, decimated_faces, this->decimated_faces_map, this->decimated_vertex_map);
		time_t current_time = time(NULL);
		std::cout << "Time elapsed: " << current_time - prev_time << std::endl;

		// removing duplicate vertices.
		Eigen::VectorXi SVI, SVJ;
		igl::remove_duplicate_vertices(decimated_vertices, decimated_faces, 10e-10, this->Vertices, SVI, SVJ, this->Faces);

		this->numFaces = this->Faces.rows();
		this->numVertices = this->Vertices.rows();

		// calculating perface normal
		igl::per_face_normals(this->Vertices, this->Faces, this->faceNormal);

		// calculating per vertex normal;
		igl::per_vertex_normals(this->Vertices, this->Faces, this->verticesNormal);
	}

	void Mesh::curvatureCalculation(float min_concave_threshold = 0.0)
	{
		std::cout << "Running Curvature calculation..." << std::endl;

		int NUM_NEIGHBORS = 4;
		time_t prev_time = time(NULL);
		igl::principal_curvature(this->Vertices, this->Faces, this->maximum_curvature_direction, this->minimum_curvature_direction, this->maximum_curvature_value, this->minimum_curvature_value, this->bad_vertices, NUM_NEIGHBORS, false);
		time_t current_time = time(NULL);
		std::cout << "Time elapsed: " << current_time - prev_time << std::endl;
		// taking minimum curvature and setting threshold such that value represents concave features.

		this->concave_curvature_value = -this->minimum_curvature_value;
		int count = this->minimum_curvature_value.rows();
		for (size_t i = 0; i < concave_curvature_value.size(); i++)
		{
			if (this->concave_curvature_value[i] < min_concave_threshold)
			{
				this->concave_curvature_value[i] = 0;
				count--;
			}
		}
		std::cout << "Number of Values higher than threshold: " << count << std::endl;

		if (this->decimated)
		{
			Eigen::VectorXd decimated_curvature = this->concave_curvature_value;
			this->concave_curvature_value.setConstant(pre_decimatedVertices.rows(), 0.0);

			for (size_t i = 0; i < this->decimated_vertex_map.size(); i++)
			{
				int indx = this->decimated_vertex_map[i];
				this->concave_curvature_value[indx] = 0.5 * (decimated_curvature[i] + this->concave_curvature_value[indx]);
			}
		}
	}

	void Mesh::curvatureCorrection(float curv_thresh, float smooth_thres)
	{
		this->smoothness_threshold = smooth_thres;
		this->curvature_threshold = curv_thresh;
		this->curvatureCorrection();
	}

	void Mesh::curvatureCorrection()
	{
		std::cout << "Running Curvature Correction..." << std::endl;
		Eigen::VectorXd corrected_curvature;
		corrected_curvature.resize(this->numVertices);
		for (int v = 0; v < numVertices; v++)
		{
			const auto &adj_faces = adjacency_faces[v];
			// iterating over all the faces taking a pair at a time
			for (size_t i = 0; i < adj_faces.size(); i++)
			{
				for (size_t j = i + 1; j < adj_faces.size(); j++)
				{
					// gets the index of adjacent faces in this->Faces
					int f1 = adj_faces[i];
					int f2 = adj_faces[j];

					// get the normal vectors from this->faceNormal. value of this->faceNormal at index i gives the normal vector of face at index i of this->Faces
					Eigen::RowVectorXd n1 = this->faceNormal.row(f1).normalized();
					Eigen::RowVectorXd n2 = this->faceNormal.row(f2).normalized();

					// find the dot product between normalized vectors and get angle by cos inverse of that dot product
					float dot = std::clamp(n1.dot(n2), -1.0, 1.0);
					float angle = std::acos(dot);

					// normals are facing opposite direction
					bool facing_away = (n1.dot(n2) < 0);

					// check concave angle threshold
					if (angle > this->smoothness_threshold)
					{
						// deep valley like surface if facing away else normal concave
						if (facing_away)
						{
							// increase own curvature value
							corrected_curvature[v] += this->curvature_threshold + 0.1f;
							// increase curvature value for all the neighbouring vertices in the adjacecy_vertex.
							for (int neighbour : this->adjacency_vertices[v])
							{
								corrected_curvature[neighbour] += curvature_threshold + 0.05;
							}
						}
						else
						{
							corrected_curvature[v] += curvature_threshold;
						}
					}
				}
			}
		}

		this->corrected_concave_curvature_value = corrected_curvature;

		if (this->decimated)
		{
			Eigen::VectorXd decimated_curvature = this->corrected_concave_curvature_value;
			this->corrected_concave_curvature_value.setConstant(this->pre_decimatedVertices.rows(), 0.0);
			for (size_t i = 0; i < this->decimated_vertex_map.size(); i++)
			{
				int org_idx = this->decimated_vertex_map[i];
				std::cout << this->Vertices.row(i) - this->pre_decimatedVertices.row(org_idx) << std::endl;
				this->corrected_concave_curvature_value[org_idx] = decimated_curvature[i];
			}
		}
	}

	double Mesh::angle_between(Eigen::Vector3d n1, Eigen::Vector3d n2)
	{
		double dot = std::max(-1.0, std::min(1.0, n1.dot(n2)));
		return acos(dot) * 180.0 / M_PI;
	}

	void Mesh::regionGrowingSegmentation(int neighbourNumbers, float curvature_threshold_filter, float angle_threshold_filter, int min_cluster_size, int max_cluster_size)
	{
		time_t prev = time(NULL);

		// variable to store the cluster labels for each point and also denotes if the vertex is processed previously or not
		Eigen::VectorXi point_labels = Eigen::VectorXi::Constant(this->numVertices, -1);

		// build an octree for KNN
		std::vector<std::vector<int>> point_indices;
		Eigen::MatrixXi CH;
		Eigen::MatrixXd CN;
		Eigen::VectorXd W;
		igl::octree(this->Vertices, point_indices, CH, CN, W);

		// find nearest Neighbour indices
		Eigen::MatrixXi knn_indices;
		igl::knn(this->Vertices, neighbourNumbers, point_indices, CH, CN, W, knn_indices);

		// Create an indices list and sort based on increasing curvature
		std::vector<int> indices(this->numVertices);
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&](int i, int j)
				  { return this->concave_curvature_value(i) < this->concave_curvature_value(j); });

		int current_cluster = 0;
		// int boundary_cluster = -10;
		// bool isboundary_cluster = false;
		for (int idx : indices)
		{
			// if the points are processed or the curvature value is less than threshold skip the vertex
			if (point_labels(idx) != -1 || this->concave_curvature_value(idx) > curvature_threshold_filter)
				continue;
			std::vector<int> cluster;
			std::queue<int> Q;

			// vertex at idx is seed point and its neighbours will be looked for addtion to cluster
			Q.push(idx);
			point_labels(idx) = current_cluster;

			while (!Q.empty())
			{
				// take the element in front of the queue and add it to cluster
				int seed = Q.front();
				Q.pop();
				cluster.push_back(seed);

				// check all the k nearest neighbours
				for (int n = 0; n < neighbourNumbers; n++)
				{
					int neighbour = knn_indices(seed, n);
					if (point_labels(neighbour) != -1)
						continue;
					double angle = this->angle_between(this->verticesNormal.row(seed), this->verticesNormal.row(neighbour));
					if (angle > angle_threshold_filter)
						continue;
					if (this->concave_curvature_value(neighbour) < curvature_threshold_filter)
					{
						point_labels(neighbour) = current_cluster;
						Q.push(neighbour);
					}
				}
			}
			// Filter by size of the cluster
			// if cluster is too small or too large then change the points to unlabeled/ unprocessed.
			if (cluster.size() < min_cluster_size)
			{
				for (int i : cluster)
					point_labels(i) = -1;
			}
			else if (cluster.size() > max_cluster_size)
			{
				// how to split cluster?

				for (int i : cluster)
					point_labels(i) = -1;
			}
			else
			{
				std::cout << "cluster size: " << cluster.size() << std::endl;
				current_cluster++;
			}
		}
		std::cout << "Detected " << current_cluster << " clusters.\n";

		current_cluster = -10;
		std::sort(indices.begin(), indices.end(), [&](int i, int j)
				  { return this->concave_curvature_value(i) > this->concave_curvature_value(j); });
		for (int idx : indices)
		{
			// if the points are processed or the curvature value is less than threshold skip the vertex
			if (point_labels(idx) != -1 || this->concave_curvature_value(idx) < curvature_threshold_filter)
				continue;

			std::vector<int> cluster;
			std::queue<int> Q;

			// vertex at idx is seed point and its neighbours will be looked for addtion to cluster
			Q.push(idx);
			point_labels(idx) = current_cluster;
			while (!Q.empty())
			{
				// take the element in front of the queue and add it to cluster
				int seed = Q.front();
				Q.pop();
				cluster.push_back(seed);

				// check all the k nearest neighbours
				for (int n = 0; n < neighbourNumbers; n++)
				{
					int neighbour = knn_indices(seed, n);
					if (point_labels(neighbour) != -1)
						continue;
					double angle = this->angle_between(this->verticesNormal.row(seed), this->verticesNormal.row(neighbour));
					if (angle < angle_threshold_filter)
						continue;
					if (this->concave_curvature_value(neighbour) > curvature_threshold_filter)
					{
						point_labels(neighbour) = current_cluster;
						Q.push(neighbour);
					}
				}
			}
			// Filter by size of the cluster
			// if cluster is too small or too large then change the points to unlabeled/ unprocessed.
			if (cluster.size() < min_cluster_size)
			{
				for (int i : cluster)
					point_labels(i) = -1;
			}
			else if (cluster.size() > max_cluster_size)
			{
				// how to split cluster?

				for (int i : cluster)
					point_labels(i) = -1;
			}
			else
			{
				std::cout << "cluster size: " << cluster.size() << std::endl;
				current_cluster--;
			}
		}
		time_t diff = time(NULL) - prev;
		std::cout << "execution_time: " << diff << std::endl;

		Eigen::MatrixXd colors(this->numVertices, 3);
		std::map<int, Eigen::RowVector3d> cluster_colors;
		nlohmann::json j_array = nlohmann::json::array();

		auto get_random_color = []()
		{
			return Eigen::RowVector3d::Random().cwiseAbs(); // RGB in [0,1]
		};

		for (int i = 0; i < point_labels.size(); ++i)
		{
			int cluster_id = point_labels[i];
			if (cluster_id == -1)
			{
				colors.row(i) = Eigen::RowVector3d(0.7, 0.7, 0.7); // Gray for unclustered
			}
			else
			{
				if (cluster_colors.find(cluster_id) == cluster_colors.end())
				{
					cluster_colors[cluster_id] = get_random_color();
				}
				colors.row(i) = cluster_colors[cluster_id];
			}
		}
		Eigen::MatrixXi F(0, 3); // No faces

		for (const auto &[key, color] : cluster_colors)
		{
			nlohmann::json color_obj = nlohmann::json::array({color[0], color[1], color[2]});
			j_array.push_back(color_obj);
		}


		std::ofstream file(this->out_color_fileName);
		if (file.is_open())
		{
			file << j_array.dump(4); // pretty-print with 4-space indent
			file.close();
			std::cout << "JSON file written successfully." << std::endl;
		}
		else
		{
			std::cerr << "Failed to open file for writing." << std::endl;
		}
		igl::writeOFF(this->out_fileName, this->Vertices, this->Faces, colors);
	}

	void Mesh::setCurvatureCorrectionThresholds(double curvature_threshold, double smoothness_threshold)
	{
		this->curvature_threshold = curvature_threshold;
		this->smoothness_threshold = smoothness_threshold;
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
		std::cout << "Calculating Adjacency Vertices and Faces" << std::endl;
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
}
