#include "Mesh.h"

#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/triangle_triangle_adjacency.h>
Mesh::Mesh()
{

}

Mesh::~Mesh()
{

}

void Mesh::preProcess(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces)
{
    _numVertices = vertices.size()/3;
    _numFaces= faces.size()/3;

    igl::cotmatrix(vertices,faces,_cotangentWeigths);
    Eigen::MatrixXd normals;
    igl::per_vertex_normals(vertices,faces,_normals);


    igl::vertex_triangle_adjacency(vertices.rows(),faces,_incidentFaces,_vertexIndicesIncidents);


    _vertexList.resize(vertices.rows());

    _weigthsVecList.resize(faces.rows());
    _dv.resize(vertices.rows());

    //the barycentric area around vertex
    Eigen::SparseMatrix<double> tmp;
    igl::massmatrix(vertices,faces,igl::MASSMATRIX_TYPE_BARYCENTRIC,tmp);
    _areaXVertex = tmp.diagonal();

    Eigen::VectorXi b;

    std::vector<int> adjacentFaces;
    for(int i = 0 ; i <  vertices.rows(); i++)
    {
     adjacentFaces = _incidentFaces[i];
     int ns = adjacentFaces.size() * 3;
     _vertexList[i].resize(ns, 2);
     _weigthsVecList[i].resize(ns);

     for (size_t j=0; j<adjacentFaces.size(); j++)
     {
         int v0 = faces(adjacentFaces[j],0);
         int v1 = faces(adjacentFaces[j],1);
         int v2 = faces(adjacentFaces[j],2);

         _vertexList[i](3*j,0) = v0;
         _vertexList[i](3*j,1) = v1;
         _vertexList[i](3*j+1,0) = v1;
         _vertexList[i](3*j+1,1) = v2;
         _vertexList[i](3*j+2,0) = v2;
         _vertexList[i](3*j+2,1) = v0;


         _weigthsVecList[i](3*j) = _cotangentWeigths.coeff(v0,v1);
         _weigthsVecList[i](3*j+1) = _cotangentWeigths.coeff(v1,v2);
         _weigthsVecList[i](3*j+2) = _cotangentWeigths.coeff(v2,v0);
     }

     _dv[i].resize(3, adjacentFaces.size()*3);
     Eigen::MatrixXd hV0, hV1;

     Eigen::MatrixXi m =_vertexList[i].col(0);
//      std::cout << "Here is the matrix m:" << std::endl   << m << std::endl;


     igl::slice(vertices,m,1,hV0);
     igl::slice(vertices,_vertexList[i].col(1),1,hV1);
     _dv[i] = (hV1 - hV0).transpose();
    }

}
