#ifndef __MESH_H
#define __MESH_H

#include <Eigen/Core>
#include <vector>

namespace CrownGen {

class Mesh {
public:
	void curvatureExtraction();
	void calculateAdjacencyVerticesFaces();

	void setNumVertices(int numVertices);
	void setNumFaces(int numFaces);
	void initSTL(std::string filename);

	Mesh();

	Eigen::MatrixXd Vertices;
	Eigen::MatrixXi Faces;
	Eigen::MatrixXd Normal;

private:
	Eigen::MatrixXd verticesNormal, facesNormal, facesNormal_;
	Eigen::VectorXd K1, K2;
	std::vector<unsigned int> bad_vertices;

	std::vector<std::vector<int>> adjacency_vertices;
    std::vector<std::vector<int>> adjacency_faces;

	int numVertices, numFaces;

	std::vector<double> _curvatures;
};

}



#endif
