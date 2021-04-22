#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <memory>
#include <Eigen>

class Mesh
{
public:
    Mesh();
    ~Mesh();
     void preProcess(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F);

    Eigen::SparseMatrix<double> _cotangentWeigths;
    std::vector<std::vector<int>> _incidentFaces;
    std::vector<std::vector<int>> _vertexIndicesIncidents;
    Eigen::MatrixXd _normals;
    std::vector<Eigen::MatrixXi> _vertexList;

    std::vector<Eigen::VectorXd> _weigthsVecList;
    std::vector<Eigen::MatrixXd> _dv;
    Eigen::MatrixXd _areaXVertex;

    int _numVertices;
    int _numFaces;

    double _relaviteDisplacement;
};

#endif // MESH_H
