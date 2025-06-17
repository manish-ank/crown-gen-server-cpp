#ifndef __MESH_H
#define __MESH_H

#include <Eigen/Core>
#include <vector>

namespace CrownGen
{

	class Mesh
	{
	public:
		void init(std::string filename, std::string out_filename, std::string out_color_filename);
		void curvatureCalculation(float min_concave_threshold);
		void calculateAdjacencyVerticesFaces();
		void curvatureCorrection();

		void setNumVertices(int numVertices);
		void setNumFaces(int numFaces);
		void setCurvatureCorrectionThresholds(double, double);
		void curvatureCorrection(float curv_thresh, float smooth_thres);
		void regionGrowingSegmentation(int neighbourNumbers, float curvature_threshold_filter, float angle_threshold_filter, int min_cluster_size, int max_cluster_size);

		void decimate();

		Mesh();

		Eigen::MatrixXd Vertices, pre_decimatedVertices;
		Eigen::MatrixXi Faces, pre_decimatedFaces;
		Eigen::MatrixXd verticesNormal, faceNormal;
		Eigen::VectorXi decimated_vertex_map, decimated_faces_map;

	private:
		Eigen::VectorXd concave_curvature_value;
		Eigen::VectorXd corrected_concave_curvature_value;

		Eigen::MatrixXd maximum_curvature_direction, minimum_curvature_direction;
		Eigen::VectorXd maximum_curvature_value, minimum_curvature_value; // principle curvature values K1 and K2 respectively.
		std::vector<unsigned int> bad_vertices;

		std::vector<std::vector<int>> adjacency_vertices;
		std::vector<std::vector<int>> adjacency_faces;

		int numVertices, numFaces;

		std::string fileName, out_fileName, out_color_fileName;

		double curvature_threshold, smoothness_threshold;

		bool decimated = false;

		double angle_between(Eigen::Vector3d n1, Eigen::Vector3d n2);
	};
}
#endif
